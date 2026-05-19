# =========================================================
# Codon usage and RSCU analysis for CDS FASTA files
# Colab-friendly script
# =========================================================

# If running in Google Colab, keep the next two lines.
# If running locally, comment them and place your FASTA files in the same folder.
!pip install biopython -q

import re
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict, Counter
from Bio import SeqIO
from google.colab import files

# Upload your CDS FASTA file(s) in Colab
uploaded = files.upload()

# =========================================================
# USER CONFIGURATION
# Option 1: analyze one FASTA file
# Option 2: add multiple entries to compare multiple files/species/genomes
# =========================================================

FASTA_FILES = {
    "sample_1": "sample_cds.fasta",
    # "sample_2": "another_sample_cds.fasta",
}

OUTPUT_PREFIX = "codon_usage"
EXCLUDE_INTERNAL_STOPS = False

# =========================================================
# CODON TABLE: standard genetic code
# =========================================================

CODON_TO_AA = {
    'TTT':'Phe', 'TTC':'Phe', 'TTA':'Leu', 'TTG':'Leu',
    'TCT':'Ser', 'TCC':'Ser', 'TCA':'Ser', 'TCG':'Ser',
    'TAT':'Tyr', 'TAC':'Tyr', 'TAA':'End', 'TAG':'End',
    'TGT':'Cys', 'TGC':'Cys', 'TGA':'End', 'TGG':'Trp',
    'CTT':'Leu', 'CTC':'Leu', 'CTA':'Leu', 'CTG':'Leu',
    'CCT':'Pro', 'CCC':'Pro', 'CCA':'Pro', 'CCG':'Pro',
    'CAT':'His', 'CAC':'His', 'CAA':'Gln', 'CAG':'Gln',
    'CGT':'Arg', 'CGC':'Arg', 'CGA':'Arg', 'CGG':'Arg',
    'ATT':'Ile', 'ATC':'Ile', 'ATA':'Ile', 'ATG':'Met',
    'ACT':'Thr', 'ACC':'Thr', 'ACA':'Thr', 'ACG':'Thr',
    'AAT':'Asn', 'AAC':'Asn', 'AAA':'Lys', 'AAG':'Lys',
    'AGT':'Ser', 'AGC':'Ser', 'AGA':'Arg', 'AGG':'Arg',
    'GTT':'Val', 'GTC':'Val', 'GTA':'Val', 'GTG':'Val',
    'GCT':'Ala', 'GCC':'Ala', 'GCA':'Ala', 'GCG':'Ala',
    'GAT':'Asp', 'GAC':'Asp', 'GAA':'Glu', 'GAG':'Glu',
    'GGT':'Gly', 'GGC':'Gly', 'GGA':'Gly', 'GGG':'Gly'
}

AA_ORDER = [
    'Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'End', 'Gln', 'Glu',
    'Gly', 'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro',
    'Ser', 'Thr', 'Trp', 'Tyr', 'Val'
]

AA_TO_CODONS = defaultdict(list)
for codon, aa in CODON_TO_AA.items():
    AA_TO_CODONS[aa].append(codon)
for aa in AA_TO_CODONS:
    AA_TO_CODONS[aa] = sorted(AA_TO_CODONS[aa])

THIRD_BASE_COLORS = {
    "A": "#ff5a3d",
    "T": "#e3bd7f",
    "C": "#71c7c7",
    "G": "#bfe3c0",
}

# =========================================================
# FUNCTIONS
# =========================================================

def clean_sequence(seq):
    seq = str(seq).upper().replace("U", "T")
    seq = re.sub(r"[^ATGC]", "", seq)
    trim = len(seq) - (len(seq) % 3)
    return seq[:trim]


def count_codons_strict(filepath):
    codon_counts = Counter()
    invalid_records = []

    for record in SeqIO.parse(filepath, "fasta"):
        seq = clean_sequence(record.seq)

        if len(seq) == 0 or len(seq) % 3 != 0:
            invalid_records.append((record.description, "length_not_multiple_of_3_or_empty"))
            continue

        codons = [seq[i:i+3] for i in range(0, len(seq), 3)]

        # Remove terminal stop codon from coding statistics, if present.
        if codons and CODON_TO_AA.get(codons[-1]) == "End":
            codons = codons[:-1]

        internal_stops = sum(1 for codon in codons if CODON_TO_AA.get(codon) == "End")
        if EXCLUDE_INTERNAL_STOPS and internal_stops > 0:
            invalid_records.append((record.description, f"internal_stop_codons:{internal_stops}"))
            continue

        valid = True
        for codon in codons:
            if codon not in CODON_TO_AA:
                invalid_records.append((record.description, f"unknown_or_ambiguous_codon:{codon}"))
                valid = False
                break

        if valid:
            codon_counts.update(codons)

    return codon_counts, invalid_records


def calculate_rscu(codon_counts, sample_name):
    rows = []
    grand_total = sum(codon_counts.values())

    for aa in AA_ORDER:
        codons = AA_TO_CODONS[aa]
        total = sum(codon_counts[c] for c in codons)
        n = len(codons)

        for codon in codons:
            count = codon_counts[codon]
            if total == 0:
                rscu = 0.0
            else:
                expected = total / n
                rscu = count / expected if expected > 0 else 0.0

            freq = (count / grand_total * 1000) if grand_total > 0 else 0.0

            rows.append({
                "Sample": sample_name,
                "AA": aa,
                "Codon_DNA": codon,
                "Codon_RNA": codon.replace("T", "U"),
                "Third_base": codon[2],
                "Count": count,
                "RSCU": round(rscu, 4),
                "Freq_per_1000": round(freq, 4),
            })

    return pd.DataFrame(rows)


def plot_rscu_grouped(df, output_png):
    """
    Plot individual RSCU values per codon.
    This avoids misleading stacked bars, since RSCU sums to the number of synonymous codons.
    """
    samples = list(df["Sample"].unique())
    fig, axes = plt.subplots(len(samples), 1, figsize=(18, 5 * len(samples)), sharex=True)
    if len(samples) == 1:
        axes = [axes]

    for ax, sample in zip(axes, samples):
        sub_sample = df[df["Sample"] == sample].copy()
        x_positions = []
        labels = []
        values = []
        colors = []
        x = 0

        for aa in AA_ORDER:
            sub_aa = sub_sample[sub_sample["AA"] == aa].sort_values("Codon_DNA")
            for _, row in sub_aa.iterrows():
                x_positions.append(x)
                labels.append(row["Codon_RNA"])
                values.append(row["RSCU"])
                colors.append(THIRD_BASE_COLORS[row["Third_base"]])
                x += 1
            x += 1  # space between amino acid groups

        ax.bar(x_positions, values, color=colors, edgecolor="white")
        ax.axhline(1.0, color="black", linestyle="--", linewidth=1)
        ax.set_ylabel("RSCU", fontsize=12)
        ax.set_title(f"RSCU codon usage - {sample}", fontsize=14, fontweight="bold")
        ax.set_xticks(x_positions)
        ax.set_xticklabels(labels, rotation=90, fontsize=8)

    plt.tight_layout()
    plt.savefig(output_png, dpi=300, bbox_inches="tight")
    plt.show()


# =========================================================
# MAIN
# =========================================================

all_rscu = []
all_invalid = []

for sample_name, fasta_file in FASTA_FILES.items():
    print(f"Processing: {sample_name} | {fasta_file}")
    codon_counts, invalid_records = count_codons_strict(fasta_file)
    rscu_df = calculate_rscu(codon_counts, sample_name)
    all_rscu.append(rscu_df)

    rscu_df.to_csv(f"{sample_name}_rscu.tsv", sep="\t", index=False)

    for header, reason in invalid_records:
        all_invalid.append({"Sample": sample_name, "Header": header, "Reason": reason})

final_df = pd.concat(all_rscu, ignore_index=True)
final_df.to_csv(f"{OUTPUT_PREFIX}_combined_rscu.tsv", sep="\t", index=False)

pd.DataFrame(all_invalid).to_csv(f"{OUTPUT_PREFIX}_invalid_sequences.tsv", sep="\t", index=False)

plot_rscu_grouped(final_df, f"{OUTPUT_PREFIX}_rscu_plot.png")

print("Generated files:")
print(f"- {OUTPUT_PREFIX}_combined_rscu.tsv")
print(f"- {OUTPUT_PREFIX}_invalid_sequences.tsv")
print(f"- {OUTPUT_PREFIX}_rscu_plot.png")

files.download(f"{OUTPUT_PREFIX}_combined_rscu.tsv")
files.download(f"{OUTPUT_PREFIX}_invalid_sequences.tsv")
files.download(f"{OUTPUT_PREFIX}_rscu_plot.png")
