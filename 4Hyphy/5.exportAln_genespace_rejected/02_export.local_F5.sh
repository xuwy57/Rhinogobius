#!/bin/bash

#SET THIS TO THE CORRECT NUMBER
TOTALPARTS=16
END=$(( $TOTALPARTS - 1 ))

mkdir -p output_F5
mkdir -p logs_F5

for i in $( seq 0 $END ); do 

php 02_exportAln_multiref_withorigseq_F5.php -N ${TOTALPARTS} -f ${i} > logs_F5/export.${i}.log 2>&1 &

done
wait
exit

