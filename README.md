PHYLOGENETIC AND GENOMIC ANALYSIS OF MULTI RESISTANCE IN SALMONELLA ENTERICA 
The study investigates how phylogenetic relationships shape the distribution and evolution of multidrug resistance genes in Salmonella enterica from diverse African sources. It combines genomic characterization with evolutionary analysis to identify resistance determinants and their evolutionary relationships.

Methodology Overview
-This repository contains a comprehensive bioinformatics pipeline for analyzing antimicrobial resistance (AMR) in Salmonella enterica isolates from African sources. The analysis integrates whole-genome sequencing, resistance gene profiling, phylogenetic reconstruction, and comparative genomics to elucidate the distribution and evolution of multidrug resistance mechanisms.
-Pipeline Architecture
The entire workflow is implemented using Nextflow, a workflow management system designed for scalable and reproducible bioinformatics analysis. All scripts and configuration files are located in the ./scripts directory.

-1. Data Retrieval
Source

Database: NCBI Pathogen Detection database
Isolates: 20 Salmonella enterica strains from diverse African countries (Tanzania, Nigeria, South Africa, Tunisia, Malawi, Ethiopia, Uganda, Kenya, Mauritius, Democratic Republic of the Congo, Burkina Faso)
Sequencing Platform: Illumina (short-read data)
Metadata Collected:

BioProject accession numbers
BioSample identifiers
SRA (Sequence Read Archive) accessions
Strain identifiers
Assembly accessions
Isolation source (clinical specimens, food animals, food products)
Geographic origin and year of isolation



Data Preparation
Genome sequences and associated metadata are downloaded and organized prior to analysis. A manifest file (Table 1 in manuscript) contains all strain information including size, region, source, and sequencing platform details.
Nextflow Script: scripts/01_data_retrieval.nf

2. Genome Annotation and Resistance Gene Identification
2.1 Genome Annotation
Tool: Rapid Annotation using Subsystem Technology toolkit (RASTtk)
Platform: Bacterial and Viral Bioinformatics Resource Center (BV-BRC, version 3.55.17)
The assembled genomes are submitted to BV-BRC for comprehensive functional annotation, which includes:

Protein-coding gene prediction
Functional role assignment
Subsystem classification
Quality assessment metrics

Nextflow Script: scripts/02_genome_annotation.nf
2.2 Antibiotic Resistance Gene Screening
Tool: BV-BRC Resistance Annotation Pipeline
Database: Comprehensive AMR database curated at BV-BRC
The annotated genomes are screened against the resistance database to identify:

β-lactamase genes (e.g., TEM, OXA, CTX-M, CMY)
Fluoroquinolone resistance determinants (e.g., mutations in gyrA, parC)
Plasmid-mediated quinolone resistance (PMQR) genes (qnr alleles)
Efflux pump genes (RND, MFS, MATE families)
Regulatory genes (marA, marR, soxS, robA, ramA, cpxR, baeR, baeS)
Other aminoglycoside and tetracycline resistance markers

2.3 Biocide and Heavy Metal Resistance Gene Detection
Tool: BacMet 2.0 Database
URL: http://bacmet.biomedicine.gu.se/
Database Type: Manually curated resistance genes to biocides and heavy metals
Annotated genomes are screened against BacMet to identify resistance determinants to:

Biocides: Acriflavine, phenol, triclosan, benzalkonium chloride (BAC), tetrachlorosalicylanilide (TCS), sodium deoxycholate (SDC), sodium dodecyl sulfate (SDS)
Heavy Metals: Gold (Au), zinc (Zn), iron (Fe), silver (Ag), mercury (Hg), cadmium (Cd), tungsten (W)

Integration: Results from both BV-BRC and BacMet analyses are combined to identify genes conferring multi-resistance and cross-resistance to both antibiotics and non-antibiotic antibacterials.
Nextflow Script: scripts/02_resistance_gene_identification.nf

3. Phylogenetic Reconstruction and SNP Analysis
3.1 SNP Calling and Core Genome Extraction
Tool: Snippy v4.6.0
Components:

Read Alignment: Burrows-Wheeler Aligner (BWA)
Variant Calling: FreeBayes
Variant Annotation: SnpEff, SAMtools, BCFtools

Workflow:

FASTQ reads from all 20 Salmonella enterica strains are aligned to the reference genome (NC_011294)
SNPs are called across the entire genome
Core genome SNPs are extracted using snippy-core, excluding variable regions and mobile elements
A whole-genome SNP alignment (FASTA format) is generated for downstream phylogenetic analysis

Output: Multi-sequence alignment of core genome SNPs suitable for phylogenetic reconstruction
Nextflow Script: scripts/03_snp_calling.nf
3.2 Maximum-Likelihood Phylogenetic Inference
Tool: IQ-TREE v2
Parameters:

Substitution Model: GTR+G (General Time Reversible with Gamma rate heterogeneity)
Bootstrap Replicates: 1,000 ultrafast bootstrap replicates
Convergence: Automatic convergence criteria

Methodology:

The core genome SNP alignment is used as input
The phylogenetic tree is rooted using the reference genome NC_011294 (Salmonella enterica subsp. enterica serovar Typhimurium str. LT2)
Branch support is quantified using bootstrap values (0-100%)
Bootstrap values of 100% indicate maximum statistical support for internal nodes

Interpretation:

High bootstrap support (BS ≥ 95%) indicates strong phylogenetic signal
Tree topology reveals evolutionary relationships and clonal structure
Branch lengths are proportional to the number of SNP substitutions per site

Nextflow Script: scripts/03_phylogenetic_reconstruction.nf
3.3 Pairwise SNP Distance Analysis
Tool: Custom Python/R script using SNP alignment from Snippy
Method: Hamming distance calculation between all isolate pairs
Output: Pairwise SNP distance matrix showing:

Minimum distances: 8 SNPs (between closely related pairs)
Maximum distances: >60,000 SNPs (between divergent lineages)
Distance ranges identify clonal clusters and distinct evolutionary lineages

Nextflow Script: scripts/03_snp_distance_matrix.nf

4. SNP-Type Pattern Analysis
Methodology
Tool: Custom R script with hierarchical clustering
SNP Categories:

No SNP (0): No variation detected
Synonymous SNP only (1): Silent mutations not affecting amino acid sequence
Non-synonymous SNP only (2): Missense mutations altering protein sequence
Both synonymous and non-synonymous (3): Mixed mutation types at the locus

Analysis:

SNP-type patterns are computed for 33 efflux-associated resistance genes across all 20 isolates
Hierarchical clustering using Ward's D2 linkage method is applied
Isolates and genes are clustered based on SNP pattern similarity

Biological Interpretation:

Synonymous mutations: May affect gene expression through mRNA secondary structure, transcript stability, or codon usage bias
Non-synonymous mutations: Directly alter protein function, substrate specificity, or efflux efficiency
Clustering patterns: Reveal distinct mutational profiles reflecting different selective pressures and evolutionary trajectories

Visualization: Heatmap with color-coded SNP types showing genome-wide mutation distribution
Nextflow Script: scripts/04_snp_type_analysis.nf

5. Comparative Genomics and Pangenome Analysis
5.1 BRIG Circular Genome Comparison
Tool: BLAST Ring Image Generator (BRIG)
Reference Genome: NC_011294 (4,685,848 bp)
Methodology:

Each of the 20 African isolates is compared against the reference genome using BLAST
Results are visualized as concentric circular rings
Each ring represents one isolate
Sequence identity at each genomic position is color-coded:

Fully colored regions: High sequence similarity (>90% identity)
White/absent regions: Low similarity, deletions, or absent genomic regions

Interpretation:

Identifies conserved core genome
Reveals discrete regions of absence or reduced coverage
Highlights genomic islands, prophage insertions, and horizontally acquired elements
Shows increasing genomic divergence toward outer rings

Output: Circular comparative map with genomic coordinate scale and GC content/skew tracks
Nextflow Script: scripts/05_brig_comparison.nf
5.2 LASTZ Pairwise Whole-Genome Alignment
Tool: Geneious Prime (alignment engine: LASTZ)
Method: Whole-genome pairwise alignment against reference
Workflow:

Each isolate is aligned to the reference genome NC_011294 using LASTZ algorithm
Alignment tracks are generated showing:

Red segments: Forward-strand alignments (same orientation as reference)
Blue segments: Reverse-strand (complementary) alignments
White gaps: Absent or unaligned sequences



Analysis Features:

Synteny Visualization: Reveals conserved gene order and structure
Inversions/Rearrangements: Identified as blue segments interspersed within red regions
Nucleotide Identity Track: Shows overall sequence conservation across genome length
Gene Annotation Track: Displays reference genome features (genes, coding sequences)

Biological Significance:

Predominantly red alignment indicates maintained core genomic architecture
Blue segments suggest mobile element integration or structural rearrangements
Alignment fragmentation in gene-rich regions suggests preferential insertion of accessory elements

Nextflow Script: scripts/05_lastz_alignment.nf
5.3 Average Nucleotide Identity (ANI) Analysis
Tool: JSpecies (ANI calculation via BLAST)
Method: Pairwise genome comparison
Calculation:

ANI is calculated from pairwise BLAST comparisons of assembled genome contigs
Formula: Percentage of nucleotide positions with sequence identity across alignments
Normalized to account for variable genome sizes

Species Boundary:

Threshold: ≥95% ANI = same species
Result: All 20 isolates show ≥97% ANI, confirming species identity

Analysis Depth:

Fine-scale discrimination: 98-100% ANI values reveal subtle strain-level differences
Clustering: Hierarchical clustering groups closely related isolates (≥99.9% ANI)
Strain Relatedness: Low ANI values indicate greater genomic divergence and independent evolutionary trajectories

Visualization: Heatmap with color gradient (dark green = 100% identity, light/white = ~97-98% identity)
Nextflow Script: scripts/05_ani_analysis.nf

6. Data Visualization and Statistical Analysis
6.1 Resistance Gene Association Heatmaps
Tools: R (RStudio v4.4.3)
Packages Used:

readxl: Data import
dplyr: Data manipulation
tidyr: Data reshaping
pheatmap: Heatmap generation
ggplot2: Graphics
RColorBrewer: Color palettes

Visualizations:
Figure 1: Efflux Gene-Antibiotic Associations

X-axis: Antibiotic or drug combinations
Y-axis: Efflux genes and regulatory systems
Bubble Size: Represents detection frequency
Color Intensity: Darker (red) = higher frequency; lighter (yellow/orange) = lower frequency
Purpose: Reveals substrate-specific pump associations across drug classes

Figure 2: Efflux Gene-Biocide Associations

X-axis: Biocide and disinfectant compounds
Y-axis: Resistance genes/efflux systems
Interpretation: Demonstrates cross-resistance between antibiotics and chemical disinfectants
Biological Relevance: Indicates co-selection from routine disinfectant use

Figure 3: Efflux Gene-Metal Associations

X-axis: Metal compounds (Au, Zn, Fe, Ag, Hg, Cd, W)
Y-axis: Resistance genes
Note: Uniform detection frequency (18) indicates lineage-specific distribution
Significance: Reveals environmental co-selection pressures

Figure 4: Gene-Protein Product Associations

Purpose: Systems-level view of efflux organization
Shows: Structural components and regulatory genes functioning as integrated networks
Interpretation: Confirms coordinated resistance mechanisms rather than individual gene expression

Nextflow Script: scripts/06_visualization.nf
6.2 SNP Distance Heatmap
Method: Pairwise distance calculation from SNP alignment
Visualization:

Matrix format with rows and columns representing isolates
Cell color indicates SNP distance (darker = greater divergence)
Diagonal cells show self-comparisons (distance = 0)
Outliers highlighted separately

Nextflow Script: scripts/06_snp_distance_heatmap.nf
6.3 ANI Heatmap
Method: Hierarchical clustering of ANI values
Color Scale:

Red: ≥99.9% identity (highest similarity)
White: ~99% identity (intermediate)
Blue: ~98% identity (lower similarity)

Clustering Algorithm: Complete linkage hierarchical clustering
