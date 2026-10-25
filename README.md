GRM1 Subfamily-Specific Profile HMM
Overview

This project constructs a subfamily-specific profile Hidden Markov Model (HMM) for GRM1 (mGluR1, Metabotropic Glutamate Receptor 1) following the methodology of Bircan et al. (2024), applied here to the mGluR family for the first time. Ortholog sequence data is sourced from Selçuk et al. (2026).

Key result: The GRM1-specific HMM achieves perfect discrimination of GRM1 from all other mGluR subfamilies on held-out test sequences (Sensitivity = Specificity = Precision = 1.000, bit score threshold ≥ 1000).

Usage — How to Test a Sequence Against the GRM1 HMM
Requirements
HMMER 3.x installed (conda install -c bioconda hmmer)
Step 1 — Clone the repository
bash
git clone https://github.com/Ebrardilbaz/hmm-project.git
cd hmm-project
Step 2 — Prepare your sequence(s)

Your input must be a plain FASTA file with no gap characters:

>my_protein
MWSMMRMNCVLFLSLLLLCATASNIYERSLVPRAATRSVARMDGDIIIGA...

If your sequence comes from an MSA (contains - characters), remove them first:

bash
awk '/^>/{print} !/^>/{gsub("-",""); print}' your_msa.fasta > your_sequence.fasta
Step 3 — Run hmmsearch
bash
hmmsearch --tblout results.txt \
          results/hmm_profiles/GRM1_UpdatedProfiles.hmm \
          your_sequence.fasta
Step 4 — Interpret results
bash
grep -v "^#" results.txt

The output table contains one row per hit. The key column is bit score (column 6):

sequence_name   -   Updated   -   E-value   bit_score   ...
Bit Score	Interpretation
≥ 1000	✓ Sequence belongs to GRM1
< 1000	✗ Sequence is NOT GRM1
Example — Python script for batch testing
python
import subprocess
import sys

def is_grm1(fasta_file, hmm_file, threshold=1000):
    result = subprocess.run(
        ["hmmsearch", "--tblout", "/tmp/hmm_out.txt", hmm_file, fasta_file],
        capture_output=True
    )
    hits = {}
    with open("/tmp/hmm_out.txt") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            seq_id = parts[0]
            bit_score = float(parts[5])
            hits[seq_id] = bit_score >= threshold
    return hits

# Usage
results = is_grm1("your_sequence.fasta", "results/hmm_profiles/GRM1_UpdatedProfiles.hmm")
for seq, is_hit in results.items():
    print(f"{seq}: {'GRM1 ✓' if is_hit else 'NOT GRM1 ✗'}")
Family Selection
Role	Protein	Group	Coupling	Train	Test
Target	GRM1 (mGluR1)	Group I	Gq	154	39
Closest	GRM5 (mGluR5)	Group I	Gq	223	56
Rest	GRM2 (mGluR2)	Group II	Gi	186	47
Rest	GRM4 (mGluR4)	Group III	Gi	183	46

Rationale: Based on the Class C GPCR phylogenetic tree (Bircan et al. 2024, Figure 2A), GRM1 and GRM5 share the closest evolutionary relationship within Group I. GRM2 (Group II) and GRM4 (Group III) cover the remaining phylogenetic diversity of mGluRs. This selection was validated by reconstructing a phylogenetic tree from the ortholog sequences, which was consistent with Bircan et al. Figure 2A.

Biological Background

GRM1 (mGluR1) is a Class C GPCR that senses glutamate, the primary excitatory neurotransmitter in the brain. It couples to Gq proteins and activates phospholipase C (PLC), leading to intracellular calcium release. It plays critical roles in synaptic plasticity, learning, and memory, and is implicated in cerebellar ataxia and schizophrenia.

GRM5 (mGluR5) is the closest relative of GRM1, also coupling to Gq. Together they form Group I mGluRs and are often co-expressed in the brain.

GRM2 (mGluR2) and GRM4 (mGluR4) belong to Groups II and III respectively, both coupling to Gi proteins and inhibiting adenylyl cyclase.

Data Sources

Ortholog MSA files from Selçuk et al. (2026) Zenodo repository (DOI: 10.5281/zenodo.16505397):

data/raw/target/  → GRM1_26-10-2022_orthologs_orderedMSA_noallgaps.fasta  (193 seq)
data/raw/closest/ → GRM5_26-10-2022_orthologs_orderedMSA_noallgaps.fasta  (279 seq)
data/raw/rest/    → GRM2_26-10-2022_orthologs_orderedMSA_noallgaps.fasta  (233 seq)
                  → GRM4_26-10-2022_orthologs_orderedMSA_noallgaps.fasta  (229 seq)

Files contain pre-aligned ortholog sequences (orderedMSA) with fully-gapped columns removed (noallgaps). Note: partial gaps within sequences remain and are removed separately before HMM search.

Pipeline
Step 1 — Phylogenetic Validation

All raw sequences (934 total) were merged, gaps removed (awk), re-aligned with MAFFT, and a maximum likelihood tree built with FastTree (LG model). The tree was visualized in iTOL and confirmed to be consistent with Bircan et al. (2024) Figure 2A.

Tools: awk, MAFFT v7, FastTree 2.2, iTOL
Output: tree_validation/grm_tree.nwk

Step 2 — Train/Test Split (80% / 20%)

Each family's sequences were randomly shuffled (seed=42) and split into 80% training and 20% test sets. The HMM is built exclusively on training data; test data is held out for unbiased evaluation.

Tools: Python (train_test_split.py)
Output: data/train/, data/test/

Step 3 — Combined MSA + Trim

Training sequences from all 4 families were merged, gaps removed, and jointly aligned with MAFFT (6051 positions total). The combined MSA was then trimmed to retain only positions where GRM1 has amino acids (1172 positions), following Bircan et al. (2024). This ensures all families share the same coordinate system for identity score comparison.

Tools: awk, MAFFT (sbatch, Tosun HPC), Python (trim_msa.py)
Output: results/*_trimmed_msa.fasta (1172 positions each)

Step 4 — Identity Score Calculation

For each alignment position, the frequency of the most common amino acid in each family is computed from the trimmed MSA. Positions conserved only in GRM1 receive high scores; positions conserved across all families receive low scores. These scores drive position weighting in the HMM.

Tools: Python (compute_identity_score_grm.py, fasta_dict.py)
Output: results/identity_scores_grm1.csv (1172 positions)

Step 5 — Standard HMM Construction

A standard profile HMM is built from the GRM1 trimmed MSA using HMMER. --pnone disables pseudocount priors so observed frequencies are used directly (emission probabilities will be modified in the next step). --symfrac 0 ensures all 1172 positions become match states, matching the identity score coordinates exactly.

Tools: HMMER (hmmbuild --pnone --symfrac 0)
Output: results/hmm_profiles/GRM1.hmm (1172 match states)

Step 6 — Subfamily-Specific HMM Construction

The standard HMM's emission probabilities are modified by identity scores following the Bircan et al. (2024) algorithm. Positions are classified into four types (Type I → highest weight, Type II → lowest weight) based on conservation patterns across GRM1, GRM5, GRM2, and GRM4. BLOSUM62 is used to assess amino acid physicochemical similarity when classifying positions.

Tools: R (construct_HMM.R — Bircan et al., run_construct_hmm.R — wrapper), Biostrings (BLOSUM62), Tosun HPC compute node
Output: results/hmm_profiles/GRM1_UpdatedProfiles.hmm

Step 7 — Performance Evaluation

Test sequences (188 total) had gaps removed and were renamed (GRM1_1, GRM5_1... etc.) before HMM search. The GRM1-specific HMM was searched against all test sequences using hmmsearch. A bit score threshold of ≥ 1000 was selected based on the clear separation between GRM1 scores (~2835–3167) and GRM5 scores (~631–711).

Tools: awk (gap removal + rename), HMMER (hmmsearch --tblout), Python (evaluate_hmm.py)
Output: results/hmmsearch_results.txt

Results
Internal Test Set (20% held-out)
	HMM HIT (≥1000)	HMM NO HIT (<1000)
GRM1 (+)	39 (TP)	0 (FN)
Other (−)	0 (FP)	149 (TN)
Metric	Value
Sensitivity	1.000
Specificity	1.000
Precision	1.000
Threshold	Bit score ≥ 1000
External Validation (UniProt sequences, not in training data)
Sequence	UniProt	Bit Score	Result
Human GRM1	Q13255	3019.7	✓ HIT
Mouse GRM1	P97772	2955.2	✓ HIT
Rat GRM1	P23385	2951.3	✓ HIT
Human GRM5	P41594	631.8	✓ NO HIT
Human GRM3	Q14832	—	✓ NO HIT
Human GRM6	O15303	17.8	✓ NO HIT
Human GRM7	Q14831	—	✓ NO HIT
Human GRM8	O00222	6.4	✓ NO HIT
Repository Structure
hmm-project/
├── data/
│   ├── raw/
│   │   ├── target/          # GRM1 MSA from Zenodo (193 seq)
│   │   ├── closest/         # GRM5 MSA from Zenodo (279 seq)
│   │   └── rest/            # GRM2, GRM4 MSA from Zenodo
│   ├── train/               # 80% training sequences
│   │   ├── target/
│   │   ├── closest/
│   │   └── rest/
│   └── test/                # 20% held-out test sequences
│       ├── target/
│       ├── closest/
│       └── rest/
├── results/
│   ├── hmm_profiles/
│   │   ├── GRM1.hmm                  # Standard HMM
│   │   └── GRM1_UpdatedProfiles.hmm  # Subfamily-specific HMM
│   ├── GRM1_trimmed_msa.fasta        # Trimmed MSA (1172 pos)
│   ├── GRM2_trimmed_msa.fasta
│   ├── GRM4_trimmed_msa.fasta
│   ├── GRM5_trimmed_msa.fasta
│   ├── combined_aligned.fasta        # Combined MSA (6051 pos)
│   ├── combined_unaligned.fasta
│   ├── identity_scores_grm1.csv      # Per-position identity scores
│   ├── hmmsearch_results.txt         # HMM search results
│   ├── test_sequences.fasta          # Prepared test sequences
│   └── evaluation/                   # External UniProt test files
├── tree_validation/
│   ├── grm_tree.nwk                  # Phylogenetic tree (Newick)
│   ├── all_grm_aligned.fasta
│   ├── all_grm_unaligned.fasta
│   ├── id_map.tsv
│   └── colors.txt
├── scripts/
│   ├── train_test_split.py           # 80/20 split
│   ├── prepare_combined_msa.sh       # Combined MSA (MAFFT, sbatch)
│   ├── trim_msa.py                   # Trim to GRM1 positions
│   ├── compute_identity_score_grm.py # Identity score calculation
│   ├── fasta_dict.py                 # FASTA reader helper
│   ├── construct_HMM.R               # Bircan et al. original script
│   ├── run_construct_hmm.R           # GRM1 wrapper script
│   ├── evaluate_hmm.py               # TP/FP/FN/TN evaluation
│   ├── run_tree.sh                   # Phylogenetic tree pipeline
│   ├── compute_identity_score.py     # Original CaSR script (reference)
│   └── Test_HMM.R                    # Bircan et al. original script
└── notebook/
References
Bircan A, Kuru N, Dereli O, Selçuk B, Adebali O (2024)
Evolutionary history of calcium-sensing receptors unveils hyper/hypocalcemia-causing mutations
PLoS Computational Biology 20(11): e1012591
https://doi.org/10.1371/journal.pcbi.1012591
Selçuk B et al. (2026)
British Journal of Pharmacology
https://doi.org/10.1111/bph.70538
Selçuk B et al. (2022)
Life Science Alliance
https://doi.org/10.26508/lsa.202201439
Zenodo Dataset
https://doi.org/10.5281/zenodo.16505397