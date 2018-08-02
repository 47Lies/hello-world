GENECODE_VERSION=28
GRCH_VERSION=38

echo "`date` Rsem_prepare start"
source "/data/users/ltaing/.bashrc"
source activate PROTEOMEGENERATOR

RSEM_RESSOURCES_DIR=/data/tmp/ltaing/SUPPA/RSEM/RESSOURCES/
RSEM_INDEX_DIR=${RSEM_RESSOURCES_DIR}GENCODE/GencodeV${GENECODE_VERSION}/INDEX/
mkdir --parents ${RSEM_INDEX_DIR}

GENOME_URL="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENECODE_VERSION}/GRCh${GRCH_VERSION}.primary_assembly.genome.fa.gz"
GENOME="${RSEM_RESSOURCES_DIR}GENCODE/GencodeV${GENECODE_VERSION}/GRCh${GRCH_VERSION}.primary_assembly.genome.fa"

if [ ! -f ${GENOME} ];
then
 wget ${GENOME_URL} \
 --output-document ${GENOME}.gz
 gunzip ${GENOME}.gz
fi

GTF_URL="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENECODE_VERSION}/gencode.v${GENECODE_VERSION}.annotation.gtf.gz"
GTF="${RSEM_RESSOURCES_DIR}GENCODE/GencodeV${GENECODE_VERSION}/gencode.v${GENECODE_VERSION}.annotation.gtf"

if [ ! -f ${GTF} ];
then
 wget ${GTF_URL} \
 --output-document ${GTF}.gz
 gunzip ${GTF}.gz
fi
echo "`date` DL done"

RSEM_INDEX_NAME=${RSEM_INDEX_DIR}GRCh${GRCH_VERSION}_x_gencode.v${GENECODE_VERSION}
rsem-prepare-reference --gtf ${GTF} --star -p 10 ${GENOME} ${RSEM_INDEX_NAME}

