## Genomes to build and their paths
genomesBase=/data/Mullen_1/Raza/scRNA/Genomes/
DefaultFrom10X=$genomesBase/DefaultFrom10X

# STEP1 Download genomes and GTF
# download DefaultFrom10X
sh /data/Mullen_1/Raza/scRNA/Analysis/Scripts/Genome_scripts/download_10xGenomes.sh $DefaultFrom10X

# download fasta and GTF files
source=$genomesBase/reference_sources
# create dirs
mkdir $source
cd $source

# fasta
fasta_url="http://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
fasta_in="${source}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
# Gencode GTF
gtf_url="http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz"
gtf_in="${source}/gencode.v32.primary_assembly.annotation.gtf"
# LncExpDB GTF
LncExpDB_gtf_url="ftp://download.big.ac.cn/lncexpdb/0-ReferenceGeneModel/1-GTFFiles/LncExpDB_OnlyLnc.tar.gz"
LncExpDB_gtf_in="${source}/LncBook_Version2.0_OnlyLnc_hg38.gtf"

if [ ! -f "$fasta_in" ]; then
    curl -sS "$fasta_url" | zcat > "$fasta_in"
fi
if [ ! -f "$gtf_in" ]; then
    curl -sS "$gtf_url" | zcat > "$gtf_in"
fi
if [ ! -f "$LncExpDB_gtf_in" ]; then
    cd $source
    wget -c $LncExpDB_gtf_url -O - | tar -xz && mv $source/LncBook_Version2.0_onlylnc.gtf $LncExpDB_gtf_in
    ## delete genes with Invalid exons in trancript and delete genes where exons of transcript are not stored in ascending order.
    #HSALNG0056858,HSALNG0059740,HSALNG0078365,HSALNG0092690,HSALNG0093062---HSALNG0089130,HSALNG0089954,HSALNG0095105
    correctedLncGtf="${source}/$(basename "$LncExpDB_gtf_in")_corrected.gtf"
    cat $LncExpDB_gtf_in | grep -v 'HSALNG0056858\|HSALNG0059740\|HSALNG0078365\|HSALNG0092690\|HSALNG0093062\|HSALNG0089130\|HSALNG0089954\|HSALNG0095105'> $correctedLncGtf
fi


# Modify sequence headers in the Ensembl FASTA to match the file
# "GRCh38.primary_assembly.genome.fa" from GENCODE. Unplaced and unlocalized
# sequences such as "KI270728.1" have the same names in both versions.
#
# Input FASTA:
#   >1 dna:chromosome chromosome:GRCh38:1:1:248956422:1 REF
#
# Output FASTA:
#   >chr1 1
fasta_modified="$source/$(basename "$fasta_in").modified"
# sed commands:
# 1. Replace metadata after space with original contig name, as in GENCODE
# 2. Add "chr" to names of autosomes and sex chromosomes
# 3. Handle the mitochrondrial chromosome
cat "$fasta_in" \
    | sed -E 's/^>(\S+).*/>\1 \1/' \
    | sed -E 's/^>([0-9]+|[XY]) />chr\1 /' \
    | sed -E 's/^>MT />chrM /' \
    > "$fasta_modified"

# Remove version suffix from transcript, gene, and exon IDs in order to match
# previous Cell Ranger reference packages
#
# Input GTF:
#     ... gene_id "ENSG00000223972.5"; ...
# Output GTF:
#     ... gene_id "ENSG00000223972"; gene_version "5"; ...
gtf_modified="$source/$(basename "$gtf_in").modified"
# Pattern matches Ensembl gene, transcript, and exon IDs for human or mouse:
ID="(ENS(MUS)?[GTE][0-9]+)\.([0-9]+)"
cat "$gtf_in" \
| sed -E 's/gene_id "'"$ID"'";/gene_id "\1"; gene_version "\3";/' \
| sed -E 's/transcript_id "'"$ID"'";/transcript_id "\1"; transcript_version "\3";/' \
| sed -E 's/exon_id "'"$ID"'";/exon_id "\1"; exon_version "\3";/' \
> "$gtf_modified"

# Define string patterns for GTF tags
# NOTES:
# - Since GENCODE release 31/M22 (Ensembl 97), the "lincRNA" and "antisense"
#   biotypes are part of a more generic "lncRNA" biotype.
# - These filters are relevant only to GTF files from GENCODE. The GTFs from
#   Ensembl release 98 have the following differences:
#   - The names "gene_biotype" and "transcript_biotype" are used instead of
#     "gene_type" and "transcript_type".
#   - Readthrough transcripts are present but are not marked with the
#     "readthrough_transcript" tag.
#   - Only the X chromosome versions of genes in the pseudoautosomal regions
#     are present, so there is no "PAR" tag.
##  Removed by Raza
#lncRNA| IG_C_gene|IG_D_gene|IG_J_gene|IG_LV_gene|IG_V_gene|IG_V_pseudogene|IG_J_pseudogene|IG_C_pseudogene|TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene|TR_V_pseudogene|TR_J_pseudogene
BIOTYPE_PATTERN="(protein_coding)"
GENE_PATTERN="gene_type \"${BIOTYPE_PATTERN}\""
TX_PATTERN="transcript_type \"${BIOTYPE_PATTERN}\""
READTHROUGH_PATTERN="tag \"readthrough_transcript\""
PAR_PATTERN="tag \"PAR\""


# Construct the gene ID allowlist. We filter the list of all transcripts
# based on these criteria:
#   - allowable gene_type (biotype)
#   - allowable transcript_type (biotype)
#   - no "PAR" tag (only present for Y chromosome PAR)
#   - no "readthrough_transcript" tag
# We then collect the list of gene IDs that have at least one associated
# transcript passing the filters.
cat "$gtf_modified" \
    | awk '$3 == "transcript"' \
    | grep -E "$GENE_PATTERN" \
    | grep -E "$TX_PATTERN" \
    | grep -Ev "$READTHROUGH_PATTERN" \
    | grep -Ev "$PAR_PATTERN" \
    | sed -E 's/.*(gene_id "[^"]+").*/\1/' \
    | sort \
    | uniq \
    > "${source}/gene_allowlist_protein_coding"


# Filter the GTF file based on the gene allowlist
gtf_protein_coding="${source}/$(basename "$gtf_in").protein_coding.gtf"
# Copy header lines beginning with "#"
grep -E "^#" "$gtf_modified" > "$gtf_protein_coding"
# Filter to the gene allowlist
grep -Ff "${source}/gene_allowlist_protein_coding" "$gtf_modified" >> "$gtf_protein_coding"

### Step 1.1 (option) Also create a GTF for Gencode lncRNAs only for overlapping with Gencode protein coding (just for the numbers)
BIOTYPE_PATTERNLncRNA="(lncRNA)"
GENE_PATTERNLncRNA="gene_type \"${BIOTYPE_PATTERNLncRNA}\""
cat "$gtf_modified" | awk '$3 == "transcript"' | grep -E "$GENE_PATTERNLncRNA" | sed -E 's/.*(gene_id "[^"]+").*/\1/' | sort | uniq > "${source}/gene_allowlist_gencode_lncRNA"
gtf_geneCodeLncRNAPath="${source}/$(basename "$gtf_in").gencode_lncRNA.gtf"
grep -E "^#" "$gtf_modified" > "$gtf_geneCodeLncRNAPath"
grep -Ff "${source}/gene_allowlist_gencode_lncRNA" "$gtf_modified" >> "$gtf_geneCodeLncRNAPath"

# Step 2.
# clean lncRNA GTF in Python from (gtf_protein_coding) gencode.v32.primary_assembly.annotation.gtf.protein_coding.gtf and ($correctedLncGtf) LncBook_Version2.0_OnlyLnc_hg38.gtf_corrected.gtf (created above) and get cleanedLNC
cleanedLncRNAGTF=$genomesBase/OverlapCleansing/cleanedForSense_and_Antisense.gtf
#Merge $cleanedLncRNAGTF and gencode.v32.primary_assembly.annotation.gtf.protein_coding.gtf
completeCleanGTF=${source}/LncBook_Version2.0__gencode.v32.gtf
cat $gtf_protein_coding $cleanedLncRNAGTF > $completeCleanGTF


## for nonchopped lncRNAs but cleaned for sense strand
cleanedForSenseStrand=$genomesBase/OverlapCleansing/cleaned_Same_Strand_Exons.gtf
completeCleanForSenseStrandGTF=${source}/cleanedForSenseStrand_LncBook_Version2.0__gencode.v32.gtf
# $gtf_protein_coding /data/Mullen_1/Raza/scRNA/Genomes//reference_sources/gencode.v32.primary_assembly.annotation.gtf.protein_coding.gtf
cat $gtf_protein_coding $cleanedForSenseStrand > $completeCleanForSenseStrandGTF



### extract the gene names from both GTFs to be used later on for analysis reasons
# protein coding
cat $gtf_protein_coding |  awk 'OFS="\t" {if ($3=="gene") {print $16}}' | tr -d '";' | sort | uniq > ${source}/protein_coding_genenames_unique.txt
cat $gtf_protein_coding |  awk 'OFS="\t" {if ($3=="gene") {print $10}}' | tr -d '";' | sort | uniq > ${source}/protein_coding_geneids_unique.txt
cat ${source}/protein_coding_genenames_unique.txt ${source}/protein_coding_geneids_unique.txt > ${source}/protein_coding_genenames_and_ids_unique.txt
# long non coding RNAs of GENCODE
cat $gtf_geneCodeLncRNAPath  |  awk 'OFS="\t" {if ($3=="gene") {print $16}}' | tr -d '";' | sort | uniq > ${source}/GENECODE_lncRNA_genenames_unique.txt
cat $gtf_geneCodeLncRNAPath  |  awk 'OFS="\t" {if ($3=="gene") {print $10}}' | tr -d '";' | sort | uniq > ${source}/GENECODE_lncRNA_geneids_unique.txt
cat ${source}/GENECODE_lncRNA_genenames_unique.txt ${source}/GENECODE_lncRNA_geneids_unique.txt > ${source}/GENECODE_lncRNA_genenames_and_ids_unique.txt
# long non coding RNAs of Singletrome
cat $cleanedForSenseStrand |  awk 'OFS="\t" {if ($3=="gene") {print $12}}' | tr -d '";' | sort | uniq > ${source}/Singletrome_lncRNA_genenames_unique.txt
cat $cleanedForSenseStrand |  awk 'OFS="\t" {if ($3=="gene") {print $10}}' | tr -d '";' | sort | uniq > ${source}/Singletrome_lncRNA_geneids_unique.txt
cat ${source}/Singletrome_lncRNA_genenames_unique.txt ${source}/Singletrome_lncRNA_geneids_unique.txt  > ${source}/Singletrome_genenames_and_ids_unique.txt






# 5' expression cleaned GTF
inhouseGTF=$genomesBase/OverlapCleansing/cleanedFor5PrimeExpressionFile.gtf



# Step 3.
# create gnome indexes
# scripts and tools paths
cellrangerPath=/data/Mullen_1/Raza/installations/Softwares/cellranger/cellranger-3.1.0/cellranger
cd $genomesBase
bsub -q bigmem $cellrangerPath mkref --ref-version="GRCh38-proteinCoding" --genome="GRCh38-proteinCoding" --fasta="$fasta_modified" --genes="$gtf_protein_coding"
bsub -q bigmem $cellrangerPath mkref --ref-version="GRCh38-lncRNA" --genome="GRCh38-lncRNA" --fasta="$fasta_modified" --genes="$cleanedLncRNAGTF"
bsub -q bigmem $cellrangerPath mkref --ref-version="GRCh38-proteinLncRNA" --genome="GRCh38-proteinLncRNA" --fasta="$fasta_modified" --genes="$completeCleanGTF"
bsub -q bigmem $cellrangerPath mkref --ref-version="GRCh38-proteinLncRNA-inhouse" --genome="GRCh38-proteinLncRNA-inhouse" --fasta="$fasta_modified" --genes="$inhouseGTF"
bsub -q bigmem $cellrangerPath mkref --ref-version="GRCh38-LncExpDB" --genome="GRCh38-LncExpDB" --fasta="$fasta_modified" --genes="$completeCleanForSenseStrandGTF"


echo '-Genome generation done-'

#### Step 4. Check all the four genomes (default and 3 created against the small test data from cellranger) ################   Run test datasets
cellrangerPath=/data/Mullen_1/Raza/installations/Softwares/cellranger/cellranger-6.0.2/cellranger
dataDir=/data/Mullen_1/Raza/installations/Softwares/cellranger/cellranger-6.0.2/external/cellranger_tiny_fastq/
testDir=/data/Mullen_1/Raza/scRNA/Results/Test/Test_genomes/
mkdir $testDir
cd ${testDir}
$cellrangerPath testrun --id=tinygex
$cellrangerPath count --id=refdata-gex-GRCh38-2020-A --transcriptome=/data/Mullen_1/Raza/scRNA/Genomes//DefaultFrom10X/refdata-gex-GRCh38-2020-A --fastqs=${dataDir} --sample=tinygex --chemistry=SC3Pv3
$cellrangerPath count --id=GRCh38-lncRNA --transcriptome=/data/Mullen_1/Raza/scRNA/Genomes/GRCh38-lncRNA/ --fastqs=${dataDir} --sample=tinygex --chemistry=SC3Pv3
$cellrangerPath count --id=GRCh38-proteinCoding --transcriptome=/data/Mullen_1/Raza/scRNA/Genomes/GRCh38-proteinCoding/ --fastqs=${dataDir} --sample=tinygex --chemistry=SC3Pv3
$cellrangerPath count --id=GRCh38-proteinLncRNA --transcriptome=/data/Mullen_1/Raza/scRNA/Genomes/GRCh38-proteinLncRNA/ --fastqs=${dataDir} --sample=tinygex --chemistry=SC3Pv3
