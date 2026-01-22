genome=/public3/group_crf/home/g21mazhy26/R.formosanus_data/07_nextpolish/01_rundir/genome.nextpolish.fasta

genome=`realpath $genome`

mkdir -p references
ln -sf $genome references/ref.fa
bwa index references/ref.fa
fastahack -i references/ref.fa

mkdir -p restriction_sites
cut -f1,2 references/ref.fa.fai > restriction_sites/ref.chrom.sizes

cd restriction_sites
python /data/software/juicer/misc/generate_site_positions.py DpnII ref ../references/ref.fa

cd ..
mkdir fastq
cd fastq
ln -s /data/projects/zma/R.formosanus_data/02_clean.data/hic/cleaned_chopped_1.fq.gz ./read_R1.fastq.gz
ln -s /data/projects/zma/R.formosanus_data/02_clean.data/hic/cleaned_chopped_2.fq.gz ./read_R2.fastq.gz
