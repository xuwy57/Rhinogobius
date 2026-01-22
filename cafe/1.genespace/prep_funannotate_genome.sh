#for i in ~/APL/compare/ref_data/ensembl_ref/*; do
#{	sStem=$(basename $i)
#	if [ -d "$i" ]
#	then
#		echo $sStem
#		php prep_funannotate_genome.php $sStem $i/$sStem.longest_isoform.prot.fa
#	else
#		echo "this file is not directory"
#	fi
#} &	
#done
#wait

for i in ~/APL/compare/ref_data/ref/*; do
{       sStem=$(basename $i)
        if [ -d "$i" ]
        then
#               echo $sStem
                php prep_funannotate_genome.php $sStem $i/$sStem.noisoform_annotated.fa
        else
                echo "this file is not directory"
        fi
} &
done
wait
