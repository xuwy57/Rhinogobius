#!/usr/bin/bash
#SBATCH -p blade
#SBATCH -c 1
#SBATCH -J snpEff

module load java
module load samtools
module load bcftools

sSNPEFF="/public/software/snpEff/exec/snpeff";
ref=/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/ref/ref/Rform.fasta
gff=/fast3/group_crf/home/b20xuwy57/Rhinogobius_spp/ref/ref/Rform.gff3

mkdir -p snpeff_data/genomes
mkdir -p snpeff_data/ref

ln -sf `realpath $ref` snpeff_data/genomes/ref.fa
ln -sf `realpath $gff` snpeff_data/ref/genes.gff

echo "data.dir = ./snpeff_data/" > snpeff.config
echo "ref.genome : ref" >> snpeff.config

$sSNPEFF build -c snpeff.config -dataDir `pwd`/snpeff_data/  -gff3 -noCheckProtein -noCheckCds -v ref > ./snpeff.build.log 2>&1
$sSNPEFF ann -c snpeff.config ref 'final_merged.vcf.gz' | bgzip -c > 'ann.final_merged.vcf.gz'
tabix ann.final_merged.vcf.gz

