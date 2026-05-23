nextflow.enable.dsl=2

params.reads  = "fastq/*_{1,2}.fastq.gz"
params.outdir = "results"

process DOWNLOAD_REFERENCE {
    publishDir "${params.outdir}/reference", mode: 'copy'
    output:
    path "reference.fasta", emit: reference
    script:
    """
    wget -q \
      "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_003197.2&rettype=fasta&retmode=text" \
      -O reference.fasta
    """
}

process FASTQC_RAW {
    tag "$sample_id"
    publishDir "${params.outdir}/fastqc_raw", mode: 'copy'
    input:
    tuple val(sample_id), path(reads)
    output:
    path "*_fastqc.{html,zip}", emit: fastqc
    script:
    """
    fastqc ${reads}
    """
}

process FASTP {
    tag "$sample_id"
    publishDir "${params.outdir}/trimmed", mode: 'copy'
    input:
    tuple val(sample_id), path(reads)
    output:
    tuple val(sample_id), path("${sample_id}_1.trimmed.fastq.gz"), path("${sample_id}_2.trimmed.fastq.gz"), emit: trimmed
    script:
    """
    fastp \
        -i ${reads[0]} \
        -I ${reads[1]} \
        -o ${sample_id}_1.trimmed.fastq.gz \
        -O ${sample_id}_2.trimmed.fastq.gz \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --length_required 36 \
        --thread ${task.cpus}
    """
}

process SNIPPY {
    tag "$sample_id"
    publishDir "${params.outdir}/snippy", mode: 'copy'
    input:
    tuple val(sample_id), path(r1), path(r2)
    path reference
    output:
    path "${sample_id}/", emit: snippy_dir
    script:
    """
    snippy \
        --outdir ${sample_id} \
        --ref ${reference} \
        --R1 ${r1} --R2 ${r2} \
        --cpus ${task.cpus}
    """
}

process SNIPPY_CORE {
    publishDir "${params.outdir}/snippy_core", mode: 'copy'
    input:
    path snippy_dirs
    path reference
    output:
    path "core.full.aln", emit: core_aln
    script:
    """
    snippy-core --ref ${reference} --prefix core ${snippy_dirs}
    """
}

process IQTREE {
    publishDir "${params.outdir}/phylogeny", mode: 'copy'
    input:
    path alignment
    output:
    path "mytree.*", emit: tree
    script:
    """
    iqtree -s ${alignment} -m GTR+G -bb 1000 -nt ${task.cpus} -pre mytree
    """
}

workflow {
    reads_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    DOWNLOAD_REFERENCE()
    ref_ch = DOWNLOAD_REFERENCE.out.reference
    FASTQC_RAW(reads_ch)
    FASTP(reads_ch)
    SNIPPY(FASTP.out.trimmed, ref_ch)
    SNIPPY_CORE(SNIPPY.out.snippy_dir.collect(), ref_ch)
    IQTREE(SNIPPY_CORE.out.core_aln)
}
