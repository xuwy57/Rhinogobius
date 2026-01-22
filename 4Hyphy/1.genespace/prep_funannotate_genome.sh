for i in /fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/ref/addition_0812/*; do
{	sStem=$(basename $i)
	if [ -d "$i" ]
	then
		echo $sStem
		php prep_funannotate_genome.php $sStem $i/$sStem.longest_isoform.prot.fa
	else
		echo "this file is not directory"
	fi
} &	
done
wait

