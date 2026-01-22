#IN=../genespace/rundir/orthofinder/Results_Sep19/Orthogroups/Orthogroups.GeneCount.tsv
IN=Orthogroups.GeneCount.tsv
cat $IN | cut -f1-26 | awk '{if ($1=="Orthogroup") {print "Desc\t"$0 } else {print "(null)\t"$0 } }' > orthogroups.genecount.txt
