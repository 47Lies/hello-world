#!/bin/bash
## Torque Configuration
#PBS -l walltime=100:00:00
#PBS -l mem=120gb
#PBS -l nodes=1:ppn=12
#PBS -q batch
#PBS -N GetGen
#PBS -j oe

#Ligne Supposee :            <h1>Release 28 (GRCh38.p12)</h1>
GENECODE_VERSION=`wget --quiet https://www.gencodegenes.org/releases/current.html -O - | grep "<h1>" | sed -e "s/.*Release //" | sed -e "s/ .*//"`
echo "Gencode version from www.gencodegenes.org : ${GENECODE_VERSION}"

GRCH_VERSION=`wget --quiet https://www.gencodegenes.org/releases/current.html -O - | grep "<h1>" | sed -e "s/.*(//" | sed -e "s/).*//"`
echo "GRCH version from www.gencodegenes.org : ${GRCH_VERSION}"

ResDir=RESSOURCES/GENCODE/GencodeV${GENECODE_VERSION}/
mkdir --parents ${ResDir}

GENOME_URL="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENECODE_VERSION}/${GRCH_VERSION}.genome.fa.gz"
GENOME="${ResDir}${GRCH_VERSION}.genome.fa"

if [ ! -f ${GENOME} ];
then
 wget ${GENOME_URL} \
 --output-document ${GENOME}.gz
 gunzip ${GENOME}.gz
else
 echo "${GENOME} file already present"
fi

GTF_URL="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENECODE_VERSION}/gencode.v${GENECODE_VERSION}.annotation.gtf.gz"
GTF="${ResDir}/gencode.v${GENECODE_VERSION}.annotation.gtf"

if [ ! -f ${GTF} ];
then
 wget ${GTF_URL} \
 --output-document ${GTF}.gz
 gunzip ${GTF}.gz
else
 echo "${GTF} file already present"
fi
echo "`date` DL done"

source "/data/users/ltaing/.bashrc"
source activate PROTEOMEGENERATOR
RSEM_INDEX_DIR=RESSOURCES/GENCODE/GencodeV${GENECODE_VERSION}/INDEX/RSEM/
mkdir --parents ${RSEM_INDEX_DIR}
RSEM_INDEX_NAME=${RSEM_INDEX_DIR}${GRCH_VERSION}_x_gencode.v${GENECODE_VERSION}
rsem-prepare-reference --gtf ${GTF} -p 8 --star ${GENOME} ${RSEM_INDEX_NAME}

