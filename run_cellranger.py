import subprocess
import csv
import os
import sys

# Ensure the user provides an input file
def main():
    if len(sys.argv) < 3:
        print("Usage: python run_cellranger.py <input_file> <transcriptome_path> [--localcores=<cores>] [--localmem=<memory>] [--create-bam=<true|false>]")
        sys.exit(1)

    input_file = sys.argv[1]
    transcriptome_path = sys.argv[2]

    # Set default values for optional parameters
    local_cores = 4
    local_mem = 124

    # Parse optional parameters
    for arg in sys.argv[3:]:
        if arg.startswith('--localcores='):
            local_cores = int(arg.split('=')[1])
        elif arg.startswith('--localmem='):
            local_mem = int(arg.split('=')[1])

    # Check if the input file exists
    if not os.path.exists(input_file):
        print(f"Error: File '{input_file}' not found.")
        sys.exit(1)

    # List to store generated BAM file paths
    bam_files = []

    # Read the input file and execute cellranger count for each sample
    with open(input_file, 'r') as file:
        csv_reader = csv.reader(file, delimiter='\t')
        for row in csv_reader:
            if len(row) < 3:
                print("Error: Each line in the input file must have at least 3 fields: <sample_name> <fastq_path> <sample_id>")
                continue

            sample_name = row[0]
            fastq_path = row[1]
            sample_id = row[2]

            # Construct the cellranger command
            command = [
                "/opt/software/cellranger-6.0.2/cellranger",
                "count",
                f"--id={sample_name}",
                f"--fastqs={fastq_path}",
                f"--transcriptome={transcriptome_path}",
                f"--localcores={local_cores}",
                f"--localmem={local_mem}",
                f"--sample={sample_id}"
            ]

            # Run the command
            try:
                print(f"Running cellranger count for sample: {sample_name}")
                subprocess.run(command, check=True)

                bam_file = os.path.join(sample_name, "outs", "possorted_genome_bam.bam")
                if os.path.exists(bam_file):
                    bam_files.append(bam_file)
                else:
                    print(f"Warning: BAM file not found for sample '{sample_name}' at {bam_file}")
            except subprocess.CalledProcessError as e:
                print(f"Error: cellranger count failed for sample '{sample_name}'. Details: {e}")

    # Merge BAM files
    if bam_files:
        merged_bam = "merged.bam"
        merge_bam_files(bam_files, merged_bam)

        # Index the merged BAM file
        index_bam_file(merged_bam)

        # Convert BAM to BigWig
        bigwig_file = "merged.bw"
        convert_bam_to_bigwig(merged_bam, bigwig_file)

def merge_bam_files(bam_files, output_file):
    try:
        print("Merging BAM files...")
        command = ["samtools", "merge", "-@", "100", output_file] + bam_files
        subprocess.run(command, check=True)
        print(f"Merged BAM file created: {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error: Failed to merge BAM files. Details: {e}")

def index_bam_file(bam_file):
    try:
        print("Indexing BAM file...")
        command = ["samtools", "index", "-@", "100", bam_file]
        subprocess.run(command, check=True)
        print(f"BAM file indexed: {bam_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error: Failed to index BAM file. Details: {e}")

def convert_bam_to_bigwig(bam_file, output_file):
    try:
        print("Converting BAM to BigWig...")
        command = [
            "bamCoverage",
            "--verbose",
            "-of", "bigwig",
            "--numberOfProcessors", "max",
            "-b", bam_file,
            "-o", output_file
        ]
        subprocess.run(command, check=True)
        print(f"BigWig file created: {output_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error: Failed to convert BAM to BigWig. Details: {e}")

if __name__ == "__main__":
    main()
