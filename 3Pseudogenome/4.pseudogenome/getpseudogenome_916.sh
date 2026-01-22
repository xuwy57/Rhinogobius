#!/usr/bin/bash
#SBATCH -p long,gpu,blade
#SBATCH -c 4
#SBATCH --mem=150G

# Getting alignments using genome reference for concatenated analyses and bucky (allowing missing data):
# -w 1000 break off threshold 1000bp; if more than 1000bp gap exists between two alignment blocks, break them into two genes
# -c list containing references to ref genome and the bcf files
# -R character requirements
# -s 0.5: allow 50% missing per column
# -l 500, final alignment length at least 500bp
# -d 5, coverage < 5x turned to "N"
# -N no,  when applying the final divergence filter, don't count "N" as a difference.
# -o output folder

# -D no = turn off divergence filter. -n include N positions in the ref sequence. -I yes: code indels or multiple non-ref alleles as N, so that coordinate identical with ref.
# -r ratio of reads supporting the alternative allele, for example, 1 means 0 reads supporting ref, 1 supporting alt, or 0 to 2, 0 to 3 etc. if doesn't pass this, masked as N

module load bcftools
module load clusterbasics

mkdir ./testaln_916
php getAlignmentAllowN_GATK.php -n yes -q 1 -D no -w 1000000000 -c list_916.txt -R charRequirements.txt -s 1 -l 10 -d 1 -N no -r 0.5 -o ./testaln_916 -I yes -b no -M yes
[ $? -ne 0 ] &&  exit 1 ;



arrInput=( `cut -f1 list_916.txt` )
#arrInput=( `cut -f1 list1taxa.txt` )
sRefGenome=`cut -f2 list_916.txt | head -n 1`
#sRefGenome=`cut -f2 list1taxa.txt | head -n 1`


#mkdir $outFolder
for i in "${arrInput[@]}" 
do 
echo converting $i ...

( php phy2fasta.php -R $sRefGenome -L ${i} -o ${i}_pseudogenome.fasta ; fastahack -i ${i}_pseudogenome.fasta ) &

done
wait
