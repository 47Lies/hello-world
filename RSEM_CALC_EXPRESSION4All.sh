LogsDir=/data/users/ltaing/DATA_TMP/ltaing/SUPPA/RSEM/SCRIPTS/hello-world/Logs/
for file in MB01 \
MB02 \
MB03 \
MB04 \
MB05 \
MB06 \
MB07 \
MB08 \
MB09 \
MB13 \
MB14 \
MB15 \
MB16 \
MB17 \
MB19 \
MB20 \
MB22 \
MB24 \
MB25 \
MB30 \
MB31 \
MB34 \
MB36 \
MB38 \
MB39 \
MB40 \
MB41 \
MB42 \
MB43 \
MT1179 \
MT1219 \
MT1364 \
MT1377 \
MT1402 \
MT161 \
MT229 \
MT2381 \
MT314 \
MT435
do
 qsub -v Sample=${file} RSEM_CALC_EXPRESSION4QSUB.sh -o ${LogsDir}${file}.RSEM_CALC_EXPR.QSubLog.txt 
done
