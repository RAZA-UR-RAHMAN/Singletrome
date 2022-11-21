genomeDir=$1
mkdir $genomeDir
cd $genomeDir
# download Human reference  from
GRCh38_url=https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
GRCh38_in="${genomeDir}/refdata-gex-GRCh38-2020-A"
# download mouse reference  from
mm10_url=https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz
mm10_in="${genomeDir}/refdata-gex-mm10-2020-A"

if [ ! -d "$GRCh38_in" ]; then
    wget -P $genomeDir "$GRCh38_url" && tar -xzvf $genomeDir/refdata-gex-GRCh38-2020-A.tar.gz -C $genomeDir
fi
if [ ! -d "$mm10_in" ]; then
    wget -P $genomeDir "$mm10_url" && tar -xzvf $genomeDir/refdata-gex-mm10-2020-A.tar.gz -C $genomeDir
fi
