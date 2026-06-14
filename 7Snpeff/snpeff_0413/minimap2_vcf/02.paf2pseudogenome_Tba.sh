module load bcftools
module load bedtools
module load samtools

genome1=/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/ref/ref/Rform.fasta
genome2=/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/ref/outgroups/Tba_sort.fa
sp1=Rform
sp2=Tridentiger_barbatus

cd ./rform_Tba
bcftools view -H -i 'TYPE!="snp"' out.$sp1.$sp2.vcf.gz | awk '{ len=length($4); start=$2-1; end=start+len;  print $1"\t"start"\t"end }' > indels.mask.$sp1.$sp2.bed
cat covered.regions.$sp1.$sp2.txt | awk '{print $1"\t"($2-1)"\t"$3}' > covered.regions.$sp1.$sp2.bed
cut -f1-2 $genome1.fai | sort -k1,1 -k2,2n > genome.chrlen.$sp1.txt
bedtools complement -i  covered.regions.$sp1.$sp2.bed -g genome.chrlen.$sp1.txt >  uncovered.regions.$sp1.$sp2.bed
cat indels.mask.$sp1.$sp2.bed uncovered.regions.$sp1.$sp2.bed | sort -k1,1 -k2,2n > masked.regions.$sp1.$sp2.bed
bcftools consensus -i 'TYPE=="snp"'  -f $genome1 out.$sp1.$sp2.vcf.gz > _pseudogenome.$sp1.$sp2.fa
bedtools maskfasta -mc N -fi _pseudogenome.$sp1.$sp2.fa  -fo pseudogenome.$sp1.$sp2.fa -bed  masked.regions.$sp1.$sp2.bed
samtools faidx pseudogenome.$sp1.$sp2.fa
rm _pseudogenome.$sp1.$sp2.fa
