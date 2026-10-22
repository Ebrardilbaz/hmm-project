import sys
import os
import csv
import operator
from fasta_dict import *

def get_aa_dict(msa_dict, position):
    aa_dict = {}
    for k, v in msa_dict.items():
        aa = v[position]
        if not aa in aa_dict.keys():
            aa_dict[aa] = 1
        else:
            aa_dict[aa] += 1
    if not "-" in aa_dict.keys():
        aa_dict["-"] = 0
    return aa_dict

def compute_identity_score(msa_dict, gap_threshold):
    seq_number = len(list(msa_dict.keys()))
    aln_length = len(list(msa_dict.values())[0])
    identity_dict = {}
    for i in range(aln_length):
        aa_dict = get_aa_dict(msa_dict, i)
        gap_prop = aa_dict["-"] / seq_number
        max_aa = max(aa_dict.items(), key=operator.itemgetter(1))[0]
        if max_aa == "-":
            identity_score = 0
            identity_dict[i+1] = {'identity_score': identity_score, 'consensus_aa': max_aa}
        elif gap_prop <= gap_threshold:
            identity_score = aa_dict[max_aa] / (seq_number - aa_dict["-"])
            identity_dict[i+1] = {'identity_score': identity_score, 'consensus_aa': max_aa}
        else:
            identity_score = aa_dict[max_aa] / seq_number
            identity_dict[i+1] = {'identity_score': identity_score, 'consensus_aa': max_aa}
    return identity_dict

def write_csv(score_file, grm1_file, grm5_file, grm2_file, grm4_file):
    grm1_dict = get_fasta_dict(grm1_file)
    grm1_identity = compute_identity_score(grm1_dict, 0.3)

    grm5_dict = get_fasta_dict(grm5_file)
    grm5_identity = compute_identity_score(grm5_dict, 0.3)

    grm2_dict = get_fasta_dict(grm2_file)
    grm2_identity = compute_identity_score(grm2_dict, 0.3)

    grm4_dict = get_fasta_dict(grm4_file)
    grm4_identity = compute_identity_score(grm4_dict, 0.3)

    with open(score_file, 'w') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(["position",
                         "GRM1", "GRM1_aa",
                         "GRM5", "GRM5_aa",
                         "GRM2", "GRM2_aa",
                         "GRM4", "GRM4_aa"])
        for i in range(1, len(grm1_identity) + 1):
            print(i)
            writer.writerow([i,
                grm1_identity[i]['identity_score'], grm1_identity[i]['consensus_aa'],
                grm5_identity[i]['identity_score'], grm5_identity[i]['consensus_aa'],
                grm2_identity[i]['identity_score'], grm2_identity[i]['consensus_aa'],
                grm4_identity[i]['identity_score'], grm4_identity[i]['consensus_aa']])

# Train setlerini kullan
write_csv(
    "results/identity_scores_grm1.csv",
    "data/train/target/GRM1_26-10-2022_orthologs_orderedMSA_noallgaps_train.fasta",
    "data/train/closest/GRM5_26-10-2022_orthologs_orderedMSA_noallgaps_train.fasta",
    "data/train/rest/GRM2_26-10-2022_orthologs_orderedMSA_noallgaps_train.fasta",
    "data/train/rest/GRM4_26-10-2022_orthologs_orderedMSA_noallgaps_train.fasta"
)
