import os

# hmmsearch sonuçlarını oku
results = {}
with open("results/hmmsearch_results.txt") as f:
    for line in f:
        if line.startswith("#") or line.strip() == "":
            continue
        parts = line.split()
        seq_id = parts[0]
        evalue = float(parts[4])
        score = float(parts[5])
        fam = seq_id.split("_")[0]
        results[seq_id] = {"family": fam, "evalue": evalue, "score": score}

# Bit score threshold
threshold_score = 1000

TP = FP = FN = TN = 0

# Test setindeki tüm sequenceleri say
all_seqs = {}
for fname in os.listdir("data/test/target"):
    n = 0
    for line in open(f"data/test/target/{fname}"):
        if line.startswith(">"):
            n += 1
            all_seqs[f"GRM1_{n}"] = "GRM1"

for fname in os.listdir("data/test/closest"):
    n = 0
    for line in open(f"data/test/closest/{fname}"):
        if line.startswith(">"):
            n += 1
            all_seqs[f"GRM5_{n}"] = "GRM5"

for fname in os.listdir("data/test/rest"):
    fam = fname.split("_")[0]
    n = 0
    for line in open(f"data/test/rest/{fname}"):
        if line.startswith(">"):
            n += 1
            all_seqs[f"{fam}_{n}"] = fam

print(f"Toplam test sequence: {len(all_seqs)}")

# TP/FP/FN/TN hesapla
for seq_id, fam in all_seqs.items():
    is_grm1 = (fam == "GRM1")
    hit = seq_id in results and results[seq_id]["score"] >= threshold_score

    if is_grm1 and hit:
        TP += 1
    elif not is_grm1 and hit:
        FP += 1
    elif is_grm1 and not hit:
        FN += 1
    elif not is_grm1 and not hit:
        TN += 1

print(f"\nBit Score Threshold: >= {threshold_score}")
print(f"TP: {TP}")
print(f"FP: {FP}")
print(f"FN: {FN}")
print(f"TN: {TN}")
if (TP+FN) > 0:
    print(f"\nSensitivity: {TP/(TP+FN):.3f}")
if (TN+FP) > 0:
    print(f"Specificity: {TN/(TN+FP):.3f}")
if (TP+FP) > 0:
    print(f"Precision:   {TP/(TP+FP):.3f}")
