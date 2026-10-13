# Profile HMM-project
Profile HMM pipeline for GPCRs.

## Purpose
To build a subfamily-specific profile HMM for TAS1R1 
using phylogenetic tree information and evaluate its 
performance on held-out test sequences.

## Family Selection
- **Target:** TAS1R1 
- **Closest:** TAS1R2 
- **Rest:** TAS1R3 and GPRC6A
## Rationale
Based on the Class C GPCR phylogenetic tree from 
Bircan et al. (2024), TAS1R1 and TAS1R2 share a common 
ancestor from taxonomy clade 117571-Euteleostomi, making 
TAS1R2 the closest relative to TAS1R1. TAS1R3 and GPRC6A evolved 
earlier (7776 Gnathostomata) and is therefore selected 
as the rest group.

## Biological Background
TAS1R1 is a Class C GPCR involved in umami taste 
perception. It forms a heterodimer with TAS1R3 to detect 
L-glutamate and other L-amino acids. Unlike TAS1R2 (sweet 
taste) and TAS1R3 (shared subunit), TAS1R1 has a distinct 
ligand-binding profile focused on savory/umami stimuli. GPRC6A is a Class C GPCR that responds to a broad 
spectrum of ligands including basic L-amino acids 
(L-Arg, L-Lys, L-Orn), divalent cations (calcium, 
magnesium, zinc), osteocalcin, and testosterone. It plays 
roles in bone metabolism, testosterone production, and 
insulin secretion. Unlike TAS1R family members which are 
highly specialized for taste perception, GPRC6A acts as 
a metabolic sensor responding to multiple ligand types

## References
- Bircan et al. (2024) PLoS Computational Biology
- Selçuk et al. (2022) Life Science Alliance
- Selçuk et al. (2026) British Journal of Pharmacology

## Pipeline
1. Sequence collection (NCBI)
2. Train/Test split (80%/20%)
3. Multiple sequence alignment (MAFFT)
4. HMM construction (hmmbuild)
5. HMM search (hmmsearch)
6. Performance evaluation (TP/FP/FN/TN)

## Repository Structure

    hmm-project/
    ├── data/
    │   ├── raw/
    │   │   ├── target/     # Raw TAS1R1 sequences from NCBI
    │   │   ├── closest/    # Raw TAS1R2 sequences from NCBI
    │   │   └── rest/       # Raw TAS1R3 and GPRC6A sequences from NCBI
    │   ├── train/
    │   │   ├── target/
    │   │   ├── closest/
    │   │   └── rest/
    │   └── test/
    │       ├── target/
    │       ├── closest/
    │       └── rest/
    ├── results/
    │   ├── alignments/     # MAFFT output files
    │   ├── hmm_profiles/   # Constructed HMM profiles
    │   └── evaluation/     # Performance results
    ├── scripts/            # Python and bash scripts
    └── notebook/           # Analysis and visualization

## Expected Output
A TAS1R1-specific profile HMM that correctly identifies 
TAS1R1 sequences while rejecting TAS1R2 and TAS1R3, GPRC6A
evaluated using True Positive, False Positive, 
False Negative and True Negative metrics.