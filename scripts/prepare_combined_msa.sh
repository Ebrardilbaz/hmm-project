#!/bin/bash
#SBATCH --job-name=combined_msa
#SBATCH --account=adelab
#SBATCH --partition=genomics
#SBATCH --qos=adelab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/cta/users/ebrar.dilbaz/hmm-project/combined_msa_%j.log

module load conda/miniconda_20250420
source $(conda info --base)/etc/profile.d/conda.sh
conda activate grm_tree

cd /cta/users/ebrar.dilbaz/hmm-project

# Gap temizle ve birleştir
for f in data/train/target/*.fasta data/train/closest/*.fasta data/train/rest/*.fasta; do
    fam=$(basename "$f" | cut -d_ -f1)
    awk -v fam="$fam" '/^>/{print ">" fam "_" ++n} !/^>/{gsub("-",""); print}' "$f"
done > results/combined_unaligned.fasta

# MAFFT ile hizala
mafft --auto --thread 4 results/combined_unaligned.fasta > results/combined_aligned.fasta

echo "Done!"
