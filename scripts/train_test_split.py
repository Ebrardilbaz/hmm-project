import random
from pathlib import Path

# Tekrarlanabilir sonuçlar için seed
random.seed(42)

# Dosya yolları
files = {
    "target": "data/raw/target/GRM1_26-10-2022_orthologs_orderedMSA_noallgaps.fasta",
    "closest": "data/raw/closest/GRM5_26-10-2022_orthologs_orderedMSA_noallgaps.fasta",
    "rest_grm2": "data/raw/rest/GRM2_26-10-2022_orthologs_orderedMSA_noallgaps.fasta",
    "rest_grm4": "data/raw/rest/GRM4_26-10-2022_orthologs_orderedMSA_noallgaps.fasta"
}

output = {
    "target":    ("data/train/target",  "data/test/target"),
    "closest":   ("data/train/closest", "data/test/closest"),
    "rest_grm2": ("data/train/rest",    "data/test/rest"),
    "rest_grm4": ("data/train/rest",    "data/test/rest")
}

def read_fasta(filepath):
    sequences = []
    current_header = ""
    current_seq = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_header:
                    sequences.append((current_header, "".join(current_seq)))
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)
        if current_header:
            sequences.append((current_header, "".join(current_seq)))
    return sequences

def write_fasta(sequences, filepath):
    Path(filepath).parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, "w") as f:
        for header, seq in sequences:
            f.write(f"{header}\n{seq}\n")

for name, filepath in files.items():
    print(f"\nProcessing {name}...")
    sequences = read_fasta(filepath)
    total = len(sequences)
    random.shuffle(sequences)
    split_idx = int(total * 0.8)
    train_seqs = sequences[:split_idx]
    test_seqs = sequences[split_idx:]
    print(f"  Total: {total} | Train: {len(train_seqs)} | Test: {len(test_seqs)}")
    train_dir, test_dir = output[name]
    filename = Path(filepath).name
    write_fasta(train_seqs, f"{train_dir}/{filename.replace('.fasta', '_train.fasta')}")
    write_fasta(test_seqs, f"{test_dir}/{filename.replace('.fasta', '_test.fasta')}")
    print(f"  Saved to {train_dir}/ and {test_dir}/")

print("\nDone!")
