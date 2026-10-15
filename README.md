# Profile HMM Project

## Purpose
To build a subfamily-specific profile HMM for mGluR1 (GRM1)
using phylogenetic information and pre-processed ortholog 
sequence data from Selçuk et al. (2026), following the 
subfamily-specific HMM construction methodology from 
Bircan et al. (2024). Performance is evaluated on held-out 
test sequences.

## Family Selection
- **Target:** GRM1 (mGluR1 — Metabotropic Glutamate Receptor 1)
- **Closest:** GRM5 (mGluR5 — Group I, closest relative)
- **Rest:** GRM2 (Group II representative) + GRM4 (Group III representative)

## Rationale
Metabotropic glutamate receptors (mGluRs) are divided into 
three functional groups based on sequence similarity, 
pharmacology and signal transduction:

    Group I:   GRM1, GRM5 → couple to Gq, activate PLC
    Group II:  GRM2, GRM3 → couple to Gi, inhibit AC
    Group III: GRM4, GRM6, GRM7, GRM8 → couple to Gi, inhibit AC

Based on the Class C GPCR phylogenetic tree from Bircan et al. 
(2024) Figure 2A, GRM1 and GRM5 share the closest evolutionary 
relationship within Group I, making GRM5 the natural choice 
for the closest group. GRM2 and GRM4 are selected as rest 
group representatives from Group II and Group III respectively, 
ensuring coverage of the full phylogenetic diversity of mGluR 
subfamilies while keeping the analysis tractable.

## Biological Background
**GRM1 (mGluR1)** is a Class C GPCR that senses glutamate,
the primary excitatory neurotransmitter in the brain. It 
couples to Gq proteins and activates phospholipase C (PLC),
leading to intracellular calcium release. It plays critical 
roles in synaptic plasticity, learning and memory, and is 
implicated in neurological disorders including cerebellar 
ataxia and schizophrenia.

**GRM5 (mGluR5)** is the closest relative of GRM1, also 
coupling to Gq and activating PLC. Together they form 
Group I mGluRs and are often co-expressed in the brain.

**GRM2 (mGluR2)** belongs to Group II and couples to Gi 
proteins, inhibiting adenylyl cyclase. It is involved in 
modulating glutamate release and is a target for anxiety 
and schizophrenia treatments.

**GRM4 (mGluR4)** belongs to Group III and also couples 
to Gi proteins. It is primarily located presynaptically, 
regulating neurotransmitter release, and is implicated 
in Parkinson's disease and pain modulation.

## Data Sources
Ortholog MSA files are obtained from Selçuk et al. (2026) 
Zenodo repository (DOI: 10.5281/zenodo.16505397). For each 
GRM subfamily, the following file is used:

    GRM1_26-10-2022_orthologs_orderedMSA_noallgaps.fasta → target
    GRM5_26-10-2022_orthologs_orderedMSA_noallgaps.fasta → closest
    GRM2_26-10-2022_orthologs_orderedMSA_noallgaps.fasta → rest
    GRM4_26-10-2022_orthologs_orderedMSA_noallgaps.fasta → rest

These files contain pre-aligned ortholog sequences with gap 
positions removed based on the human sequence, following the 
same preprocessing approach used in Bircan et al. (2024).

## References
- Bircan et al. (2024) PLoS Computational Biology
  https://doi.org/10.1371/journal.pcbi.1012591
- Selçuk et al. (2022) Life Science Alliance
  https://doi.org/10.26508/lsa.202201439
- Selçuk et al. (2026) British Journal of Pharmacology
  https://doi.org/10.1111/bph.70538

## Pipeline
1. Data collection from Selçuk et al. (2026) Zenodo
2. Train/Test split (80%/20%) per subfamily
3. Identity score calculation (compute_identity_score.py)
4. Subfamily-specific HMM construction (construct_HMM.R)
5. Performance evaluation (Test_HMM.R)
6. Results reported as TP/FP/FN/TN table

## Repository Structure

    hmm-project/
    ├── data/
    │   ├── raw/
    │   │   ├── target/     # GRM1 MSA from Zenodo
    │   │   ├── closest/    # GRM5 MSA from Zenodo
    │   │   └── rest/       # GRM2, GRM4 MSA from Zenodo
    │   ├── train/
    │   │   ├── target/
    │   │   ├── closest/
    │   │   └── rest/
    │   └── test/
    │       ├── target/
    │       ├── closest/
    │       └── rest/
    ├── results/
    │   ├── hmm_profiles/   # Constructed HMM profiles
    │   └── evaluation/     # Performance results (TP/FP/FN/TN)
    ├── scripts/            # Python and R scripts
    └── notebook/           # Analysis and visualization

## Expected Output
A GRM1-specific profile HMM that correctly identifies 
GRM1 sequences while rejecting GRM5 (closest) and 
GRM2, GRM4 (rest), evaluated using True Positive, False 
Positive, False Negative and True Negative metrics, 
following the format of Table 1 in Bircan et al. (2024).