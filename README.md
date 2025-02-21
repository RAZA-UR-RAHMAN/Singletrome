# Singletrome

Singletrome, is an enhanced human lncRNA genome annotation for scRNA-seq analysis, by merging protein-coding and lncRNA databases with additional filters for quality control. Using Singletrome, we were able to cluster cell types based solely on lncRNAs expression, providing evidence of the depth and diversity of lncRNA reads contained in scRNA-seq data. Our analysis identified lncRNAs differentially expressed in specific cell types with development of liver fibrosis. Importantly, lncRNAs alone were able to predict cell types and human disease pathology through the application of machine learning. This comprehensive annotation will allow mapping of lncRNA expression across cell types of the human body facilitating the development of an atlas of human lncRNAs in health and disease. 

## Overview
This repository provides a collection of scripts designed for single-cell RNA sequencing (scRNA-seq) analysis, with a focus on incorporating long noncoding RNAs (lncRNAs) into transcriptomic studies. The toolkit includes scripts for generating enhanced genome annotations, creating BED files, and running Cell Ranger for scRNA-seq processing and visualization.



## Installation and Dependencies
1. Ensure you have [Apptainer](https://apptainer.org/) installed.
2. Download the prebuilt SIF container or rebuild it using `build.apptainer`.
---
## Prerequisites

The toolkit requires [Apptainer](https://apptainer.org/) to run the provided containerized environment.

A prebuilt SIF container is available for [download](https://www.dropbox.com/scl/fi/oy2b68i3j2vxg9rhi8ghp/singletrome_apptainer_cellranger602.sif?rlkey=qirzit6khced16svn4pdqk1md&st=h8wdx66r&dl=0).

### Download the Prebuilt SIF Container
You can download the prebuilt Apptainer container from the following link: [[Dropbox](https://www.dropbox.com/scl/fi/oy2b68i3j2vxg9rhi8ghp/singletrome_apptainer_cellranger602.sif?rlkey=qirzit6khced16svn4pdqk1md&st=h8wdx66r&dl=0)]

If you need to rebuild the image, use the `build.apptainer` script and refer to the LSF submission example in `Apptainer.lsf`.

---


## Scripts and Usage

### 1. `Singletrome.py`: Generate Singletrome GTF
`Singletrome.py` generates the Singletrome GTF, an enhanced human lncRNA genome annotation for scRNA-seq analysis.

#### Usage:
```bash
python Singletrome.py --help
```

#### Options:
```
-h, --help            Show this help message and exit.
-v, --version         Show program's version number and exit.
-lncbook_path, --pathToLncBook
                        Path to the LncExpDB GTF file containing lncRNAs.
                        Default: ftp://download.big.ac.cn/lncexpdb/0-ReferenceGeneModel/1-GTFFiles/LncExpDB_OnlyLnc.tar.gz
-gencode_path, --pathToGencode
                        Path to the Gencode GTF file containing gene annotations.
                        Default: http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz
-o, --output_dir      Directory to save output files. Default: ./Singletrome_output
-m, --mkref          Flag to indicate whether to create a cellranger reference genome. Use T (True) or F (False). Default: F
-fasta_path, --pathToFasta
                        Path to the reference genome FASTA file, required if `mkref` is set to T.
                        Default: https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
```

---

### 2. `run_cellranger.py`: Run Cell Ranger and generate merged BAM and BigWig files for visualization
This script facilitates running Cell Ranger, merging BAM files for RSeQC, and creating BigWig for visualizations.

#### Usage:
```bash
python run_cellranger.py <input_file> <transcriptome_path> [--localcores=<cores>] [--localmem=<memory>] [--create-bam=<true|false>]
```

#### Options:
- `<input_file>`: Path to input FASTQ files.
- `<transcriptome_path>`: Path to the reference transcriptome.
- `--localcores`: Number of cores to use (optional).
- `--localmem`: Amount of memory to allocate (optional).
- `--create-bam`: Whether to create BAM files (`true` or `false`).

---

### 3. `CreateBedFilesForGTF.py`: Convert GTF to BED
This script converts a GTF file to a BED12 format by first converting to GenePred and then to BED12.

#### Usage:
```bash
python CreateBedFilesForGTF.py --pathToGTF <GTF_FILE> --outDir <OUTPUT_DIR>
```

#### Options:
```
-h, --help            Show this help message and exit.
--pathToGTF           Path to the input GTF file.
--outDir              Directory to store output files.
```

---

## ANALYSIS_PIPELINES
The scripts in the ANALYSIS_PIPELINES directory were used to perform analysis and create plots, including:
- `CellrangerCompiler.Rmd`: Compiles and compares mapping stats produced by Cell Ranger.
- `Expression_Comparison_and_Filtering.Rmd`: Filters lncRNAs by comparing expression in TLGA and ULGA.
- `Read_Distribution_Filtering.Rmd`: Filters lncRNAs by assessing read distribution across transcripts.
- `SeuratPipeline.Rmd`: Analyzes protein-coding and lncRNA genes across different genome annotations.
- `Filtered_DownStream_Analysis.R`: Performs clustering and differential expression analysis of lncRNAs passing QC filters.
- `geneBody_coverage_modified.py`: A modified RSeQC script reporting coverage for each transcript before normalization.
- `CustomTheme.R`: Custom theme used for plotting figures.

---

## Utility
The `Utility` directory contains `FUNCTIONS.R`, which includes reusable functions for analysis. This script is used by multiple other scripts in the repository.

---

## Machine_Learning
The `Machine_Learning` directory contains scripts for converting Seurat objects into matrix form for cell type and disease prediction. It includes code to split data into training/testing sets and models for classification.

### `Create_Matrices_For_ML.R`
- Creates two matrices from Seurat Object:
  1. `Matrix.csv`: Contains expression counts.
  2. `metadata.csv`: Contains metadata, including cell type or disease annotations.

### `PreProcessing.ipynb`
- Takes `Matrix.csv` and `metadata.csv` as input.
- Merges them based on cell barcode.
- Identifies and removes columns with all zero values.
- Creates dummy variables for outcome variables (Cell Type/Disease).
- Saves the resultant dataframe as CSV for XGBoost classifier input.

### `CellType-Prediction-Classifier.ipynb`
- Uses XGBoost classifier to predict cell types.

### `Disease-Prediction-Classifier.ipynb`
- Uses XGBoost classifier to predict diseases.


---

## Citation

If you use this toolkit in your research, please cite:
Rahman, R., Ahmad, I., Sparks, R., Ben Saad, A., & Mullen, A. (2022). [Singletrome: A method to analyze and enhance the transcriptome with long noncoding RNAs for single cell analysis]. bioRxiv. https://doi.org/10.1101/2022.10.31.514182
---

