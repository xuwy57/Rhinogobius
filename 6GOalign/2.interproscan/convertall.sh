mkdir -p ref_GOout

for i in ./ref/*gaf*; do
	sSrcDir=`dirname $i`
	sSp=Boleophthalmus_pectinirostris
	sGAF=$i
	sGFF=/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/ref/outgroups/Boleophthalmus_pectinirostris.gff
	echo $sSp
	php gff_rna2geneid.php $sGFF $sGAF ref_GOout/$sSp.tsv
done
