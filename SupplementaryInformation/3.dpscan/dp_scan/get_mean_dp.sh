module load bcftools
vcf=$1
sample=$2
nSampleLines=$3 #921612192

mkdir -p meandp
doneflag=meandp/$sample.done
if [ -f $doneflag ]; then
echo $doneflag
else
bcftools query -f '%DP' $1 | head -n $nSampleLines |  awk '{sitecount+=1; if ($1!=".") totaldp+=$1} END{print "'$sample'\t"sitecount"\t"totaldp/sitecount}' > meandp/$sample.mean.dp.txt && touch $doneflag
fi
