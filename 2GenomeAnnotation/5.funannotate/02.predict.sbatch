#!/usr/bin/bash
#SBATCH -p long
#SBATCH -c 64
#SBATCH --mem=100G

source /public/apps/miniconda3/etc/profile.d/conda.sh
conda activate funannotate1.8.15

source genomedef.sh

CPU=$SLURM_CPUS_PER_TASK


 funannotate predict --force  -i $GENOME \
            -o trained -s "$SPECIES" --strain $STRAIN --cpus $CPU \
--max_intronlen 200000 \
--weights augustus:10 hiq:10 pasa:10 snap:0 genemark:0 glimmerhmm:0 proteins:10 transcripts:1 \
--optimize_augustus \
--busco_db actinopterygii \
--busco_seed_species zebrafish \
--organism other \
--protein_evidence ~/bahaha_assembly/annotations/protein_refs/filtered.prot.evidence.fa \
--stringtie $STRINGTIE \
--keep_no_stops \
--min_protlen 10 \
--repeat_filter none
#--repeats2evm \
#--no-evm-partitions \
# StartWeights = {'augustus': 1, 'hiq': 2, 'genemark': 1, 'pasa': 6,'codingquarry': 0, 'snap': 1, 'glimmerhmm': 1, 'proteins': 1, 'transcripts': 1}
