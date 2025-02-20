#!/usr/bin/env python3
"""
CreateBedFilesForGTF.py

Example usage:
    python CreateBedFilesForGTF.py --pathToGTF singletrome.gtf --outDir outputDirectory

This script runs:
  1) gtfToGenePred -genePredExt -geneNameAsName2 {GTF} {outDir}/{basename}.gtf-genePredExt.tmp
  2) python ConvertGenePredToBed12.py {outDir}/{basename}.gtf-genePredExt.tmp {outDir}/{basename}.gtf-bed12.bed
"""

import argparse
import os
import subprocess

def listToStringWithoutBrackets(list1):
    s=str(list1).replace('[','').replace(']','')
    return str(s).replace(' ','')
def getUpdatedLine(line):
    exonStart=line[8].split(",")
    exonEnd=line[9].split(",")
    strand=line[2]
    start=int(line[3])
    end=int(line[4])
    nummberOfExons=int(line[7])
    while("" in exonStart) :
        exonStart.remove("")
        exonEnd.remove("")
    exonStart = [int(i) for i in exonStart]
    exonEnd = [int(i) for i in exonEnd]

    block_size = []
    zip_object = zip(exonStart, exonEnd)
    for exonStart_i, exonEnd_i in zip_object:
        if(strand=='+'):
            block_size.append(exonEnd_i-exonStart_i)
        else:
            block_size.append(exonEnd_i-exonStart_i)
    block_size=listToStringWithoutBrackets(block_size)
    block_start=exonStart
    for i in range(len(exonStart)):
        if(strand=='+'):
            block_start[i] = exonStart[i]-start
        else:
            block_start[i] =exonStart[i]-start
    block_start=listToStringWithoutBrackets(block_start)
    lineToReplace=str(line[1]+" "+str(start)+" "+str(end)+" "+line[0]+" "+str(0)+" "+strand+" "+str(start)+" "+str(end)+" "+'0'+" "+str(nummberOfExons)+" "+block_size+" "+block_start+"\n")
    return lineToReplace

def createBed12FormateFromgenePredExt(inputFilePath,outputFilePath):
    with open(inputFilePath, 'r') as file:
        # read a list of lines into data
        data = file.readlines()

    for x in range(len(data)):
        line=getUpdatedLine(data[x].split())
        data[x]=line

    # and write everything back
    with open(outputFilePath, 'w') as file:
        file.writelines( data )

def main():
    parser = argparse.ArgumentParser(
        description="Convert GTF to BED12 by first converting to GenePred, then to BED12."
    )
    parser.add_argument(
        "--pathToGTF",
        required=True,
        help="Path to the input GTF file."
    )
    parser.add_argument(
        "--outDir",
        required=True,
        help="Directory where output files will be stored."
    )
    args = parser.parse_args()

    # Make sure the output directory exists
    os.makedirs(args.outDir, exist_ok=True)

    # Extract the base name from the GTF file (e.g., 'singletrome' from 'singletrome.gtf')
    gtf_basename = os.path.basename(args.pathToGTF)
    base_name = os.path.splitext(gtf_basename)[0]

    # Construct output paths
    genePred_file = os.path.join(args.outDir, f"{base_name}.gtf-genePredExt.tmp")
    bed12_file = os.path.join(args.outDir, f"{base_name}.gtf-bed12.bed")

    # Command 1: Convert GTF to GenePred
    cmd_gtfToGenePred = [
        "gtfToGenePred",
        "-genePredExt",
        "-geneNameAsName2",
        args.pathToGTF,
        genePred_file
    ]


    # Run command 1
    print("Running:", " ".join(cmd_gtfToGenePred))
    subprocess.run(cmd_gtfToGenePred, check=True)

    # Run command 2
    print("Running: createBed12FormateFromgenePredExt " + genePred_file + " " + bed12_file)
    createBed12FormateFromgenePredExt(genePred_file,bed12_file)

    os.system("rm -f " + genePred_file)
    print("All steps completed successfully.")
    print(f"GenePred file: {genePred_file}")
    print(f"BED12 file: {bed12_file}")

if __name__ == "__main__":
    main()

