import sys
sys.path.insert(0, 'scripts')
from fasta_dict import get_fasta_dict

# Combined MSA'yı oku
combined = get_fasta_dict('results/combined_aligned.fasta')

# GRM1_1'i referans al, gap olmayan pozisyonları bul
ref = combined['GRM1_1']
keep_positions = [i for i, c in enumerate(ref) if c != '-']
print(f"Toplam pozisyon: {len(ref)}")
print(f"GRM1 referans sonrası kalan pozisyon: {len(keep_positions)}")

# Tüm sequenceleri bu pozisyonlara göre kısalt
trimmed = {}
for header, seq in combined.items():
    trimmed[header] = ''.join(seq[i] for i in keep_positions)

# 4 aileye böl ve kaydet
families = ['GRM1', 'GRM5', 'GRM2', 'GRM4']
for fam in families:
    fam_seqs = {k: v for k, v in trimmed.items() if k.startswith(fam)}
    outfile = f"results/{fam}_trimmed_msa.fasta"
    with open(outfile, 'w') as f:
        for header, seq in fam_seqs.items():
            f.write(f">{header}\n{seq}\n")
    print(f"{fam}: {len(fam_seqs)} seq, {len(list(fam_seqs.values())[0])} pozisyon → {outfile}")

print("Done!")
