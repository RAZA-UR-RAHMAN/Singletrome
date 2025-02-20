#!/usr/bin/env python
import os
import re
import subprocess
import pandas as pd
import numpy as np
import sys
import itertools
from copy import deepcopy
import copy
import logging
import argparse
import textwrap
import Fun_lib
import datetime
import time

def Singletrome_gtf_generate(LncBook_file, genCod_file, outputDir):
    logging.info("Starting Singletrome GTF generation process...")
    logging.info("-----------------------------------------")

    # Define BED files
    BED_file1 = outputDir + '/LncBook_Version2.0_OnlyLnc_hg38_Exons.bed'
    BED_file2 = outputDir + '/gencode.v32.primary_assembly.annotation.gtf.filtered_Exons.bed'

    # Create exon BED files
    logging.info("Step 7 - Creating exon BED file for LncBook...")
    Fun_lib.createExonBedFile(LncBook_file, BED_file1)
    logging.info("LncBook exon BED file created.")
    logging.info("")

    logging.info("Step 8 - Creating exon BED file for Gencode...")
    Fun_lib.createExonBedFile(genCod_file, BED_file2)
    logging.info("Gencode exon BED file created.")
    logging.info("")

    # Define overlap BED files
    senseExonOverlapBedFile = outputDir + "/strand_LnRNA_ExonOverlap_ProteinCodingExons.bed"
    antiSenseExonOverlapBedFile = outputDir + "/antiStrand_LnRNA_ExonOverlap_ProteinCodingExons.bed"

    # Run Bedtools intersect for sense strand overlap
    logging.info("Step 9 - Running Bedtools intersect for sense strand overlap...")
    cmd_1 = 'intersect -s -wo -a {BED_file1} -b {BED_file2}'
    cmd_1_R = cmd_1.format(BED_file1=BED_file1, BED_file2=BED_file2)
    f1 = open(senseExonOverlapBedFile, "w")
    cmd_1_Run = subprocess.Popen(['bedtools'] + cmd_1_R.split(), stdout=f1)
    cmd_1_Run.communicate()
    f1.close()
    logging.info("Sense strand overlap BED file created.")
    logging.info("")

    # Run Bedtools intersect for antisense strand overlap
    logging.info("Step 10 - Running Bedtools intersect for antisense strand overlap...")
    cmd_2 = 'intersect -S -wo -a {BED_file1} -b {BED_file2}'
    cmd_2_R = cmd_2.format(BED_file1=BED_file1, BED_file2=BED_file2)
    f2 = open(antiSenseExonOverlapBedFile, "w")
    cmd_2_Run = subprocess.Popen(['bedtools'] + cmd_2_R.split(), stdout=f2)
    cmd_2_Run.communicate()
    f2.close()
    logging.info("Antisense strand overlap BED file created.")
    logging.info("")

    # Generate cleaned_Same_Strand_Exons.gtf
    logging.info("Step 11 - Generating cleaned same-strand exons GTF file...")
    LncBookGenes, LncBookTranscripts, LncBookExons = Fun_lib.read_gtf_file(LncBook_file)
    uniqueGeneIdsOverlappingSenseStrand = Fun_lib.deleteLncRNAGenesOverlappingSameExonsSameStrand(senseExonOverlapBedFile, LncBookGenes, outputDir)
    logging.info("Cleaned same-strand exons GTF file generated.")
    logging.info("")

    # Generate cleanedForSense_and_Antisense.gtf
    logging.info("Step 12 - Generating cleaned GTF file for sense and antisense strands...")
    uniqueGeneIdsOverlappingAntisense, df = Fun_lib.antiSenseExonOverlapCompressedBedFile_antiSenseExonsChoppedDfFilePath(antiSenseExonOverlapBedFile, outputDir)
    logging.info("Antisense overlap compressed BED file and chopped dataframe generated.")
    logging.info("")

    logging.info("Step 13 - Identifying remaining antisense genes for deletion...")
    alreadyDeletedInSenseFile = outputDir + "/alreadyDeletedInSenseFile.txt"
    remainingAntiSenseToBeDeleted = Fun_lib.getRemainingAntiSenseToBeDeleted(uniqueGeneIdsOverlappingSenseStrand, uniqueGeneIdsOverlappingAntisense, alreadyDeletedInSenseFile)
    logging.info("Remaining antisense genes to be deleted identified.")
    logging.info("")

    cleanedSameStrandExonsGTFFile = outputDir + "/cleaned_Same_Strand_Exons.gtf"
    CleanedSameStrandGenes, CleanedSameStrandTranscripts, CleanedSameStrandExons = Fun_lib.read_gtf_file(cleanedSameStrandExonsGTFFile)
    cleanedForSenseAndAntisenseFile = outputDir + "/cleanedForSense_and_Antisense.gtf"
    Fun_lib.cleaneAntiSenseOverlappingExons(remainingAntiSenseToBeDeleted, CleanedSameStrandGenes, cleanedForSenseAndAntisenseFile, df)

    # Generate Singletrome annotation file
    logging.info("Step 14 - Final GTF file generation for Singletrome annotations...")
    cmd_3 = "cat " + genCod_file + " " + outputDir + "/cleanedForSense_and_Antisense.gtf >" + outputDir + "/singletrome.gtf"
    os.system(cmd_3)
    logging.info("Singletrome annotation GTF file generated.")
    logging.info("")

    # Completion log
    TimeNow = str(datetime.datetime.now())
    logging.info("Singletrome GTF generation process completed successfully at %s", TimeNow)
    logging.info("")

def Singletrome_gtf_preprocess(pathToLncBook, pathTogenCod, outputDir):
    # Download LncExpDB GTF
    logging.info("Step 1 - Downloading LncExpDB GTF file...")
    logging.info("Downloading from %s", pathToLncBook)
    current_directory = os.getcwd()
    os.chdir(outputDir)
    cmd_1 = "curl -o LncExpDB_OnlyLnc.tar.gz " + pathToLncBook
    os.system(cmd_1)
    logging.info("Successfully downloaded LncExpDB GTF file.")
    logging.info("")

    # Extract and filter LncExpDB GTF
    logging.info("Step 2 - Extracting LncExpDB GTF and removing specific gene entries...")
    cmd_2 = "tar -xzvf LncExpDB_OnlyLnc.tar.gz && cat LncExpDB_OnlyLnc.gtf | grep -v 'HSALNG0056858\|HSALNG0059740\|HSALNG0078365\|HSALNG0092690\|HSALNG0093062\|HSALNG0089130\|HSALNG0089954\|HSALNG0095105' >LncBook_corrected.gtf && rm -f LncExpDB_OnlyLnc.gtf && rm -f LncExpDB_OnlyLnc.tar.gz"
    os.system(cmd_2)
    logging.info("LncExpDB GTF file processed and filtered successfully.")
    os.chdir(current_directory)
    logging.info("")

    # Download Gencode GTF
    logging.info("Step 3 - Downloading Gencode GTF file...")
    logging.info("Downloading from %s", pathTogenCod)
    cmd_3 = "curl -sS " + pathTogenCod + " | zcat >" + outputDir + "/gencode.primary_assembly.annotation.gtf"
    os.system(cmd_3)
    logging.info("Successfully downloaded Gencode GTF file.")
    logging.info("")

    # Modify Gencode GTF with sed
    logging.info("Step 4 - Modifying Gencode GTF file by adding version information...")
    id_pattern = r'"(ENS(MUS)?[GTE][0-9]+)\.([0-9]+)"'
    cmd_4 = (
        f"cat {outputDir}/gencode.primary_assembly.annotation.gtf | "
        f"sed -E 's/gene_id {id_pattern};/gene_id \"\\1\"; gene_version \"\\3\";/' | "
        f"sed -E 's/transcript_id {id_pattern};/transcript_id \"\\1\"; transcript_version \"\\3\";/' | "
        f"sed -E 's/exon_id {id_pattern};/exon_id \"\\1\"; exon_version \"\\3\";/' "
        f"> {outputDir}/gencode.primary_assembly.annotation.gtf.modified"
    )
    os.system(cmd_4)
    logging.info("Gencode GTF file modified successfully.")
    os.remove(f"{outputDir}/gencode.primary_assembly.annotation.gtf")
    logging.info("Gencode GTF file modified successfully.")
    logging.info("")

    # Construct the gene ID allowlist
    logging.info("Step 5 - Constructing gene ID allowlist for protein-coding genes...")
    BIOTYPE_PATTERN = "(protein_coding)"
    GENE_PATTERN = f'gene_type "{BIOTYPE_PATTERN}"'
    TX_PATTERN = f'transcript_type "{BIOTYPE_PATTERN}"'
    READTHROUGH_PATTERN = 'tag "readthrough_transcript"'
    PAR_PATTERN = 'tag "PAR"'

    cmd_5 = (
        f"cat {outputDir}/gencode.primary_assembly.annotation.gtf.modified | "
        f"awk '$3 == \"transcript\"' | "
        f"grep -E '{GENE_PATTERN}' | "
        f"grep -E '{TX_PATTERN}' | "
        f"grep -Ev '{READTHROUGH_PATTERN}' | "
        f"grep -Ev '{PAR_PATTERN}' | "
        f"sed -E 's/.*(gene_id \"[^\"]+\").*/\\1/' | "
        f"sort | uniq > {outputDir}/gene_allowlist_protein_coding"
    )
    os.system(cmd_5)
    logging.info("Gene ID allowlist for protein-coding genes constructed successfully.")
    logging.info("")

    # Filter the GTF file based on the gene allowlist
    logging.info("Step 6 - Filtering GTF file for protein-coding annotations...")
    cmd_6 = (
        f"grep -E '^#' {outputDir}/gencode.primary_assembly.annotation.gtf.modified > "
        f"{outputDir}/gencode.primary_assembly.annotation.gtf.protein_coding.gtf && "
        f"grep -Ff {outputDir}/gene_allowlist_protein_coding {outputDir}/gencode.primary_assembly.annotation.gtf.modified >> "
        f"{outputDir}/gencode.primary_assembly.annotation.gtf.protein_coding.gtf"
    )
    os.system(cmd_6)
    logging.info("GTF file filtered successfully to retain protein-coding annotations.")
    logging.info("")

    # Final log message
    TimeNow = str(datetime.datetime.now())
    logging.info("GTF preprocessing completed successfully at %s", TimeNow)
    logging.info("")

def Singletrome_mkref(pathToLncBook, pathTogenCod, outputDir, pathToFasta):
    LOG_FILENAME = outputDir + "/logfile.log"
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    logging.basicConfig(filename=LOG_FILENAME,level=logging.INFO)
    
    logging.info("=========================================")
    logging.info("Log Timestamp: " + str(datetime.datetime.now()))
    logging.info("Starting Singletrome GTF Creation Process")
    logging.info("-------------------------------------------------")
    logging.info("Downloading necessary GTF files...")
    logging.info("LncRNA GTF file: " + pathToLncBook)
    logging.info("Protein-coding GTF file: " + pathTogenCod)
    logging.info("Reference genome FASTA file: " + pathToFasta)
    logging.info("Output directory: " + outputDir)
    logging.info("-----------------------------------------")
    logging.info("")

    Singletrome_gtf_preprocess(pathToLncBook, pathTogenCod, outputDir)
    LncBook_file = outputDir + "/LncBook_corrected.gtf"
    genCod_file = outputDir + "/gencode.primary_assembly.annotation.gtf.protein_coding.gtf"
    Singletrome_gtf_generate(LncBook_file, genCod_file, outputDir)
    
    TimeNow = str(datetime.datetime.now())
    logging.info("Final Singletrome GTF file generated successfully.")
    logging.info("Output file is located at: " + outputDir + "/singletrome.gtf")
    logging.info("")

    # Download Human reference
    cmd_7 = "curl -C - " + pathToFasta + " -o " + outputDir + "/refdata-gex-GRCh38-2020-A.tar.gz"
    cmd_8 = "tar -xzvf " + outputDir + "/refdata-gex-GRCh38-2020-A.tar.gz -C " + outputDir
    cmd_9 = "rm -f " + outputDir + "/refdata-gex-GRCh38-2020-A.tar.gz"
    os.system(cmd_7)
    os.system(cmd_8)
    os.system(cmd_9)
   
    # Build genome index for cellranger
    logging.info("Building genome index for Cell Ranger...")
    current_directory = os.getcwd()
    os.chdir(outputDir)
    cmd_10 = 'mkref --genome=Singletrome_Genome_index --fasta=refdata-gex-GRCh38-2020-A/fasta/genome.fa --genes=singletrome.gtf'
    cmd_10_R = cmd_10.format(outputDir=outputDir)
    cmd_10_Run = subprocess.Popen(['/opt/software/cellranger-3.1.0/cellranger'] + cmd_10_R.split())
    cmd_10_Run.communicate()
    os.chdir(current_directory)
    
    TimeNow = str(datetime.datetime.now())
    logging.info("Genome index for Cell Ranger built successfully.")
    logging.info("Output file is located at: " + outputDir + "/Singletrome_Genome_index")
    logging.info("")
    logging.info("=========================================")
    logging.info("Process completed successfully!")

def Singletrome_no_mkref(pathToLncBook, pathTogenCod, outputDir):
    LOG_FILENAME = outputDir + "/logfile.log"
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    logging.basicConfig(filename=LOG_FILENAME,level=logging.INFO)

    logging.info("=========================================")
    logging.info("Log Timestamp: " + str(datetime.datetime.now()))
    logging.info("Starting Singletrome GTF Creation Process")
    logging.info("-------------------------------------------------")
    logging.info("Downloading necessary GTF files...")
    logging.info("LncRNA GTF file: " + pathToLncBook)
    logging.info("Protein-coding GTF file: " + pathTogenCod)
    logging.info("Output directory: " + outputDir)
    logging.info("-----------------------------------------")
    logging.info("")
    Singletrome_gtf_preprocess(pathToLncBook, pathTogenCod, outputDir)
    LncBook_file = outputDir + "/LncBook_corrected.gtf"
    genCod_file = outputDir + "/gencode.primary_assembly.annotation.gtf.protein_coding.gtf"
    Singletrome_gtf_generate(LncBook_file, genCod_file, outputDir)
   
    TimeNow = str(datetime.datetime.now())
    logging.info("Final Singletrome GTF file generated successfully.")
    logging.info("Output file is located at: " + outputDir + "/singletrome.gtf")
    logging.info("")
    logging.info("=========================================")
    logging.info("Process completed successfully!")

# Define shared parser
shared_parser = argparse.ArgumentParser(
    usage="python Singletrome.py [options]",
    formatter_class=argparse.RawTextHelpFormatter,
    description="Singletrome: A method to analyze and enhance the transcriptome with long noncoding RNAs for single cell analysis."
)

# Add version argument
shared_parser.add_argument(
    '-v', '--version',
    action='version',
    version='Singletrome v1.0'
)

# Add lncbook path argument
shared_parser.add_argument(
    '-lncbook_path', '--pathToLncBook',
    type=str,
    default='',
    help=textwrap.dedent('''
        Path to the LncExpDB GTF file. This file contains information on long non-coding RNAs (lncRNAs) from the LncExpDB database.
        Default: ftp://download.big.ac.cn/lncexpdb/0-ReferenceGeneModel/1-GTFFiles/LncExpDB_OnlyLnc.tar.gz
    ''')
)

# Add gencode path argument
shared_parser.add_argument(
    '-gencode_path', '--pathToGencode',
    type=str,
    default='',
    help=textwrap.dedent('''
        Path to the Gencode GTF file. This file contains gene annotations from the Gencode project.
        Default: http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz
    ''')
)

# Add output directory argument
shared_parser.add_argument(
    '-o', '--output_dir',
    type=str,
    default='./Singletrome_output',
    help=textwrap.dedent('''
        Directory where the output files will be saved.
        Default: ./Singletrome_output
    ''')
)

# Add mkref flag argument
shared_parser.add_argument(
    '-m', '--mkref',
    type=str,
    choices=['T', 'F'],
    default='F',
    help='Flag to indicate whether to create a genome reference. Use T for True or F for False. Default: F'
)

# Add fasta path argument
shared_parser.add_argument(
    '-fasta_path', '--pathToFasta',
    type=str,
    default='',
    required=False,
    help=textwrap.dedent('''
        Path to the reference genome FASTA file. This file is required if the mkref flag is set to T.
        Default: https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
    ''')
)

# Parsing arguments
args = shared_parser.parse_args()

# Assigning paths for LncBook and Gencode resources
path_to_lncbook = args.pathToLncBook if args.pathToLncBook else "ftp://download.big.ac.cn/lncexpdb/0-ReferenceGeneModel/1-GTFFiles/LncExpDB_OnlyLnc.tar.gz"
path_to_gencode = args.pathToGencode if args.pathToGencode else "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz"

# Creating output directory if it does not exist
output_dir = args.output_dir if args.output_dir else "./Singletrome_output"
if not os.path.isdir(output_dir):
    os.makedirs(output_dir)

mkref_flag = args.mkref == 'T'
pathToFasta = args.pathToFasta if args.pathToFasta else "https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz"

start_time = time.time()  # Record start time
if mkref_flag:
    Singletrome_mkref(path_to_lncbook, path_to_gencode, output_dir, pathToFasta)
else:
    Singletrome_no_mkref(path_to_lncbook, path_to_gencode, output_dir)


cmd_gtfToGenePred_lncRNA = [
        "gtfToGenePred",
        "-genePredExt",
        "-geneNameAsName2",
        output_dir + "/cleanedForSense_and_Antisense.gtf",
        output_dir + "/cleanedForSense_and_Antisense.gtf-genePredExt"
    ]
subprocess.run(cmd_gtfToGenePred_lncRNA, check=True)
Fun_lib.createBed12FormateFromgenePredExt(output_dir + "/cleanedForSense_and_Antisense.gtf-genePredExt", output_dir + "/cleanedForSense_and_Antisense.gtf-bed12.bed")

cmd_gtfToGenePred_pc = [
        "gtfToGenePred",
        "-genePredExt",
        "-geneNameAsName2",
        output_dir + "/gencode.primary_assembly.annotation.gtf.protein_coding.gtf",
        output_dir + "/gencode.primary_assembly.annotation.gtf-genePredExt"
    ]
subprocess.run(cmd_gtfToGenePred_pc, check=True)
Fun_lib.createBed12FormateFromgenePredExt(output_dir + "/gencode.primary_assembly.annotation.gtf-genePredExt", output_dir + "/gencode.primary_assembly.annotation.gtf-bed12.bed")

cmd_1 = "rm -f " + output_dir + "/alreadyDeletedInSenseFile.txt"
cmd_2 = "rm -f " + output_dir + "/antiSenseExonOverlapBedFileCompressed.csv"
cmd_3 = "rm -f " + output_dir + "/antiSenseExonsChoppedDf.csv"
cmd_4 = "rm -f " + output_dir + "/antiStrand_LnRNA_ExonOverlap_ProteinCodingExons.bed"
cmd_5 = "rm -f " + output_dir + "/cleanedForSense_and_Antisense.gtf"
cmd_6 = "rm -f " + output_dir + "/cleaned_Same_Strand_Exons.gtf"
cmd_7 = "rm -f " + output_dir + "/gencode.primary_assembly.annotation.gtf.protein_coding.gtf"
cmd_8 = "rm -f " + output_dir + "/gencode.v32.primary_assembly.annotation.gtf.filtered_Exons.bed"
cmd_9 = "rm -f " + output_dir + "/gene_allowlist_protein_coding"
cmd_10 = "rm -f " + output_dir + "/LncBook_corrected.gtf"
cmd_11 = "rm -f " + output_dir + "/LncBook_Version2.0_OnlyLnc_hg38_Exons.bed"
cmd_12 = "rm -f " + output_dir + "/same_Strand_Deleted_Genes.gtf"
cmd_13 = "rm -f " + output_dir + "/strand_LnRNA_ExonOverlap_ProteinCodingExons.bed"
cmd_14 = "rm -rf " + output_dir + "/refdata-gex-GRCh38-2020-A"
cmd_15 = "rm -f " + output_dir + "/Log.out"
cmd_16 = "rm -f " + output_dir + "/cleanedForSense_and_Antisense.gtf-genePredExt"
cmd_17 = "rm -f " + output_dir + "/gencode.primary_assembly.annotation.gtf-genePredExt"

os.system(cmd_1)
os.system(cmd_2)
os.system(cmd_3)
os.system(cmd_4)
os.system(cmd_5)
os.system(cmd_6)
os.system(cmd_7)
os.system(cmd_8)
os.system(cmd_9)
os.system(cmd_10)
os.system(cmd_11)
os.system(cmd_12)
os.system(cmd_13)
os.system(cmd_14)
os.system(cmd_15)
os.system(cmd_16)
os.system(cmd_17)

end_time = time.time()  # Record end time
running_time = end_time - start_time
running_time_minutes = running_time / 60  # Convert seconds to minutes
logging.info(f"Running time: {running_time_minutes:.2f} minutes")
