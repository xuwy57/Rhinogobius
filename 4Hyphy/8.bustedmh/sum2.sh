clade=Clade7
#clade=Callopanchax
#clade=Aphyosemion
#clade=Scriptaphyosemion
#clade=Archiaphyosemion
#clade=Epiplatys
#clade=Fundulopanchax
#clade=PronothoNothos
#clade=AKN
suffix=
clade=${clade}${suffix}

cat ../busted*/clade7/Relax_${clade}/ret*.txt | grep -v "error" | grep -v "GeneName" | sort -t $'\t' -k 5n,5n  > allsum_${clade}.txt


