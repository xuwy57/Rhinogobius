module load bcftools
vcf=$1
sample=$2

mkdir -p normdp
doneflag1=meandp/$sample.done
doneflag2=normdp/$sample.done
if [ -f $doneflag1 ]; then
if [ -f $doneflag2 ]; then
echo $doneflag2
else
nMeanDP=`cut -f3 meandp/$sample.mean.dp.txt`
bcftools query -f '%CHROM	%POS	%DP' $vcf | awk 'BEGIN{print "CHROM\tPOS\t'$sample'"} {print $1"\t"$2"\t"$3/'$nMeanDP'}' | gzip -c > normdp/$sample.normdp.gz && touch $doneflag2
fi
else
 echo $doneflag1 not found compute mean cov first
fi
