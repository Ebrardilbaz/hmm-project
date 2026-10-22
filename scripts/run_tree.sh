#!/bin/bash
#SBATCH --job-name=grm_tree
#SBATCH --account=adelab
#SBATCH --partition=genomics
#SBATCH --qos=adelab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=tree_%j.log

module load conda/miniconda_20250420
source $(conda info --base)/etc/profile.d/conda.sh
conda activate grm_tree

mkdir -p tree_validation
rm -f tree_validation/id_map.tsv tree_validation/all_grm_unaligned.fasta

for f in data/raw/*/GRM*_noallgaps.fasta; do
    fam=$(basename "$f" | cut -d_ -f1)
    awk -v fam="$fam" -v map=tree_validation/id_map.tsv '
        /^>/ {n++; id=fam"_"n; print id "\t" substr($0,2) >> map; print ">" id; next}
        {gsub("-",""); print}' "$f" >> tree_validation/all_grm_unaligned.fasta
done

mafft --auto --thread 4 tree_validation/all_grm_unaligned.fasta > tree_validation/all_grm_aligned.fasta
FastTree -lg tree_validation/all_grm_aligned.fasta > tree_validation/grm_tree.nwk
