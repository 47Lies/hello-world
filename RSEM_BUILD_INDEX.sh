#!/bin/bash
## Torque Configuration
#PBS -l walltime=100:00:00
#PBS -l mem=120gb
#PBS -l nodes=1:ppn=12
#PBS -q batch
#PBS -N GetGen
#PBS -j oe
source ./GetGencode.sh

source "/data/users/ltaing/.bashrc"
source activate PROTEOMEGENERATOR
RSEM_INDEX_DIR=RESSOURCES/GENCODE/GencodeV${GENECODE_VERSION}/INDEX/RSEM/
mkdir --parents ${RSEM_INDEX_DIR}
RSEM_INDEX_NAME=${RSEM_INDEX_DIR}${GRCH_VERSION}_x_gencode.v${GENECODE_VERSION}
rsem-prepare-reference --gtf ${GTF} -p 12 --star ${GENOME} ${RSEM_INDEX_NAME}

