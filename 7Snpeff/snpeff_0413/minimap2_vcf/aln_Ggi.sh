#!/usr/bin/bash
#SBATCH -p blade,himem
#SBATCH -c 36
#SBATCH --mem=50G

module load minimap2
module load bcftools
module load samtools

ref=/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/ref/ref/Rform.fasta
sp2=Glossogobius_giuris
sp1=Rform
genome2=/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/ref/outgroups/Ggi_sort.fa

paftools=paftools.js

minimap2 -t36  --secondary=no -x asm20 --cs $ref $genome2 | sort -k6,6 -k8,8n --parallel=4 > pairwise.aln.$sp1.$sp2.paf

 cut -f1,3,4,5,6,8,9,12 pairwise.aln.$sp1.$sp2.paf | awk '$8>=60 {print}' > pairwise.simp.$sp1.$sp2.txt

( $paftools call -s $sp2 -q 60 -f $ref -L200 -l200 pairwise.aln.$sp1.$sp2.paf | bgzip -c | bcftools norm -w 99999999 -O z -f $ref  - > out.$sp1.$sp2.vcf.gz ) 2>out.vcf.info
tabix out.$sp1.$sp2.vcf.gz

