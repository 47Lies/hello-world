#!/bin/bash
## Torque Configuration
#PBS -l walltime=240:00:00
#PBS -l mem=150
#PBS -l nodes=1:ppn=24
#PBS -q q_freak
#PBS -N MaxQuant_SM_QS
#PBS -j oe
source /bioinfo/users/ltaing/.bashrc
source activate PROTEOMEGENERATOR


cd /bioinfo/users/ltaing/DATA_TMP/ltaing/SUPPA/RSEM/SCRIPTS/hello-world


snakemake --cores 24
