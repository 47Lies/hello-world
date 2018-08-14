#!/bin/bash
## Torque Configuration
#PBS -l walltime=100:00:00
#PBS -l mem=60gb
#PBS -l nodes=1:ppn=10
#PBS -q batch
#PBS -N RsemCalc
#PBS -j oe

echo "start at `date` on `hostname`"

source "/data/users/ltaing/.bashrc"
source activate PROTEOMEGENERATOR

if [ -z "${Sample}" ];
then
 Sample=MB25
fi

GENECODE_VERSION="28"
GRCH_VERSION="GRCh38.p12"
FQ_DIR=/bioinfo/users/ltaing/ProteoGenomics/data/2000553/


RSEM_INDEX_DIR=RESSOURCES/GENCODE/GencodeV${GENECODE_VERSION}/INDEX/RSEM/
RSEM_INDEX_NAME=${RSEM_INDEX_DIR}${GRCH_VERSION}_x_gencode.v${GENECODE_VERSION}
rsem-calculate-expression --star \
	--star-gzipped-read-file \
	--paired-end \
	-p 10 \
	${FQ_DIR}${Sample}_R1.fastq.gz \
	${FQ_DIR}${Sample}_R2.fastq.gz \
	${RSEM_INDEX_NAME} \
	${Sample}.RSEM
