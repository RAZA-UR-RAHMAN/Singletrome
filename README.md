# Singletrome

Singletrome, is an enhanced human lncRNA genome annotation for scRNA-seq analysis, by merging protein-coding and lncRNA databases with additional filters for quality control. Using Singletrome, we were able to cluster cell types based solely on lncRNAs expression, providing evidence of the depth and diversity of lncRNA reads contained in scRNA-seq data. Our analysis identified lncRNAs differentially expressed in specific cell types with development of liver fibrosis. Importantly, lncRNAs alone were able to predict cell types and human disease pathology through the application of machine learning. This comprehensive annotation will allow mapping of lncRNA expression across cell types of the human body facilitating the development of an atlas of human lncRNAs in health and disease. 

### The scripts for creating Singletrome and for processing single cell RNA-Seq 10X Genomics Data are uploaded and explained below.

#### Genome_scripts 
The scripts in the Genome_scripts were used to download the genome fasta files, genome annotation from GENCODE and LncExpDB, scripts to build genome indices and create BED file format for comparing protein-coding and lncRNA genes as well as to perform RSEQC analysis. The gene_and_transcript_length_analysis.R has the code to perform correlation analysis of gene and transcript length for both protein-coding genes and lncRNA genes.

#### GTF_Processesing 
GtfProcessor.ipynb contains all the code and functions to read and write a GTF file with Genes, Transcripts and Exons. create_ULGA_and_TLGA_gtfs.ipynb contains all the code to create TLGA and ULGA genome annotations. This notebook deletes lncRNA genes that overlap protein-coding genes on the same strand and trim lncRNA exons that overlap protein-coding genes on the opposite strand.

#### ANALYSIS_PIPELINES 
The scripts in the ANALYSIS_PIPELINES were used to perform analysis and create plots such as CellrangerCompiler.Rmd to compile and compare the mapping stats produced by the cellranger, Expression_Comparison_and_Filtering.Rmd to filter lncRNAs by comparing the expression of lncRNAs in TLGA and ULGA, Read_Distribution_Filtering.Rmd filter lncRNAs by assessing the read distribution across the transcripts, SeuratPipeline.Rmd contains the analysis of protein-coding and lncRNA genes across different genome annotations and Filtered_DownStream_Analysis.R to perform clustering and differential expression of the lncRNAs that passed the quality control filter. The geneBody_coverage_modified.py is the modified script form RSEQC to report the coverage for each transcript before normalizing to 1.  CustomTheme.R is the theme that is used to plot figures.

#### Utility
The Utility directory contains FUNCTIONS.R which has the functions to perform analysis. This script is used by many other scripts and contain FUNCTIONS that are called by other scripts.

#### Machine_Learning
The Machine_Learning directory contains scripts that were used to convert the Seurat objects into matrix form for cell type and disease prediction. The directory also includes the code to split data into training and testing sets as well as all the models for cell type and disease prediction. 
