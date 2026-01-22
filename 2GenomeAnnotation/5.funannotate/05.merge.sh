#!/user/bin/bash
#SBATCH -c 32

updated_gff=/public3/group_crf/home/g21mazhy26/R.formosanus_data/10_funannotate/funannotate/new_train/traineid_miniprot/update_results/Rhinogobius_formosanus_DWC2021.gff3
origin_gff=/public3/group_crf/home/g21mazhy26/R.formosanus_data/10_funannotate/funannotate/new_train/traineid_miniprot/predict_results/Rhinogobius_formosanus_DWC2021.gff3
fasta=/public3/group_crf/home/g21mazhy26/R.formosanus_data/10_funannotate/funannotate/new_train/scf.softmasked.fa
sp=orig
module load gffread
 
grep -vP "^#" $origin_gff | \
sed 's/FUN_/'$sp'/g' | \
sed 's/novel_/'$sp'_/g' \
> orig.gff3

cat $updated_gff orig.gff3 > combined.gff

gffread -g $fasta --t-adopt --keep-genes -V -F -M -d dup.info -K -w out.mrna.fa -x out.cds.fa -y out.prot.fa -S -o merged.gff combined.gff
