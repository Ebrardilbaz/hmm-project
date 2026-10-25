# GRM1 Subfamily-Specific Profile HMM

A subfamily-specific profile Hidden Markov Model for **GRM1 (mGluR1)**, 
built following the methodology of Bircan et al. (2024) and applied here 
to the mGluR family for the first time.

## Family Selection

| Role    | Protein | Group     | Train | Test |
|---------|---------|-----------|-------|------|
| Target  | GRM1    | Group I   | 154   | 39   |
| Closest | GRM5    | Group I   | 223   | 56   |
| Rest    | GRM2    | Group II  | 186   | 47   |
| Rest    | GRM4    | Group III | 183   | 46   |

## Pipeline Overview

```
Raw MSA (Zenodo) → Train/Test Split (80/20) → Combined MSA (MAFFT)
→ Trim to GRM1 positions → Identity Scores → hmmbuild → construct_HMM.R
→ GRM1_UpdatedProfiles.hmm → hmmsearch → Evaluation
```

## Repository Structure

```
hmm-project/
├── data/
│   ├── raw/          # Original MSA files from Zenodo
│   ├── train/        # 80% training sequences
│   └── test/         # 20% held-out test sequences
├── results/
│   ├── hmm_profiles/
│   │   ├── GRM1.hmm                  # Standard HMM
│   │   └── GRM1_UpdatedProfiles.hmm  # Subfamily-specific HMM ← use this
│   ├── identity_scores_grm1.csv
│   ├── hmmsearch_results.txt
│   └── evaluation/                   # External validation files
├── tree_validation/                  # Phylogenetic tree files
└── scripts/
    ├── train_test_split.py
    ├── prepare_combined_msa.sh
    ├── trim_msa.py
    ├── compute_identity_score_grm.py
    ├── construct_HMM.R               # Bircan et al. original
    ├── run_construct_hmm.R           # GRM1 wrapper
    └── evaluate_hmm.py
```

## Usage — Test a Sequence Against the GRM1 HMM

### Requirements
```bash
conda install -c bioconda hmmer
```

### 1. Clone the repository
```bash
git clone https://github.com/Ebrardilbaz/hmm-project.git
cd hmm-project
```

### 2. Prepare your sequence (remove gaps if any)
```bash
awk '/^>/{print} !/^>/{gsub("-",""); print}' input.fasta > query.fasta
```

### 3. Run hmmsearch
```bash
hmmsearch --tblout results.txt \
          results/hmm_profiles/GRM1_UpdatedProfiles.hmm \
          query.fasta
```

### 4. Interpret results
```bash
grep -v "^#" results.txt
```

| Bit Score | Interpretation         |
|-----------|------------------------|
| ≥ 1000    | ✓ Sequence is GRM1     |
| < 1000    | ✗ Sequence is not GRM1 |

## References

- Bircan et al. (2024) PLoS Comput Biol — https://doi.org/10.1371/journal.pcbi.1012591
- Selçuk et al. (2026) Br J Pharmacol — https://doi.org/10.1111/bph.70538
- Zenodo Dataset — https://doi.org/10.5281/zenodo.16505397
