import csv
import operator
import sys
sys.path.insert(0, 'scripts')
from fasta_dict import get_fasta_dict

def get_aa_dict(msa_dict, position):
    aa_dict = {}
    for k, v in msa_dict.items():
        aa = v[position]
        if aa not in aa_dict:
            aa_dict[aa] = 1
        else:
            aa_dict[aa] += 1
    if "-" not in aa_dict:
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

# Trimmed MSA dosyalarını oku
grm1 = get_fasta_dict('results/GRM1_trimmed_msa.fasta')
grm5 = get_fasta_dict('results/GRM5_trimmed_msa.fasta')
grm2 = get_fasta_dict('results/GRM2_trimmed_msa.fasta')
grm4 = get_fasta_dict('results/GRM4_trimmed_msa.fasta')

print(f"GRM1: {len(grm1)} seq")
print(f"GRM5: {len(grm5)} seq")
print(f"GRM2: {len(grm2)} seq")
print(f"GRM4: {len(grm4)} seq")

# Identity score hesapla
grm1_identity = compute_identity_score(grm1, 0.3)
grm5_identity = compute_identity_score(grm5, 0.3)
grm2_identity = compute_identity_score(grm2, 0.3)
grm4_identity = compute_identity_score(grm4, 0.3)

# CSV yaz
with open("results/identity_scores_grm1.csv", "w") as f:
    writer = csv.writer(f, delimiter=',')
    writer.writerow(["position", "GRM1", "GRM1_aa", "GRM5", "GRM5_aa", "GRM2", "GRM2_aa", "GRM4", "GRM4_aa"])
    for i in range(1, len(grm1_identity) + 1):
        writer.writerow([i,
            grm1_identity[i]['identity_score'], grm1_identity[i]['consensus_aa'],
            grm5_identity[i]['identity_score'], grm5_identity[i]['consensus_aa'],
            grm2_identity[i]['identity_score'], grm2_identity[i]['consensus_aa'],
            grm4_identity[i]['identity_score'], grm4_identity[i]['consensus_aa']])

print("Done! identity_scores_grm1.csv created.")
