# =========================================================
# PR2 bias analysis for organellar CDS FASTA files
# Colab-friendly script
# =========================================================

# If running in Google Colab, keep the next two lines.
# If running locally, comment them and place your FASTA files in the same folder.
!pip install biopython -q

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from google.colab import files

# Upload your CDS FASTA files in Colab
uploaded = files.upload()

# =========================================================
# USER CONFIGURATION
# Replace the example file names below with your own generic/real file names.
# Each FASTA must contain protein-coding DNA sequences in frame.
# =========================================================

FASTA_FILES = {
    ("Species_1", "Plastome"): "species1_plastome_cds.fasta",
    ("Species_1", "Mitogenome"): "species1_mitogenome_cds.fasta",
    ("Species_2", "Plastome"): "species2_plastome_cds.fasta",
    ("Species_2", "Mitogenome"): "species2_mitogenome_cds.fasta",
}

OUTPUT_TABLE = "PR2_bias_results.tsv"
OUTPUT_FIGURE = "PR2_bias_plot.png"

# Plot options
LABEL_EXTREME_GENES = True       # True = label only the most biased genes
N_LABELS_PER_PANEL = 8           # number of gene labels per panel
LABEL_ALL_GENES = False          # True = label all genes, can make the figure crowded
EXCLUDE_INTERNAL_STOPS = False   # True = remove CDS with internal stop codons from the plot/table
MIN_CDS_LENGTH = 300             # minimum CDS length in bp

# =========================================================
# FUNCTIONS
# =========================================================

def clean_sequence(seq):
    seq = str(seq).upper().replace("U", "T")
    seq = re.sub(r"[^ATGC]", "", seq)
    trim = len(seq) - (len(seq) % 3)
    return seq[:trim]


def get_gene_name(header):
    """Extract a readable gene name from common FASTA header formats."""
    header = str(header)
    header = header.replace("gene=", "gene:")

    patterns = [
        r"\[gene=([^\]]+)\]",
        r"gene:([A-Za-z0-9_.-]+)",
        r"gene-([A-Za-z0-9_.-]+)",
        r"cds-([A-Za-z0-9_.-]+)",
    ]

    for pattern in patterns:
        match = re.search(pattern, header)
        if match:
            gene = match.group(1)
            break
    else:
        gene = header.split()[0]

    # Remove common annotation prefixes/suffixes
    for prefix in ["chloe_", "blatx_", "blastx_", "cds-blatx_", "cds-blastx_"]:
        gene = gene.replace(prefix, "")

    gene = re.sub(r"_fragment$", "", gene)
    gene = re.sub(r"_fragment_[0-9]+$", "", gene)
    gene = re.sub(r"_[0-9]+_[0-9]+$", "", gene)
    gene = re.sub(r"_[0-9]+$", "", gene)

    return gene


def count_internal_stops(seq):
    stop_codons = {"TAA", "TAG", "TGA"}
    codons = [seq[i:i+3] for i in range(0, len(seq), 3) if len(seq[i:i+3]) == 3]

    # terminal stop is allowed and ignored
    if codons and codons[-1] in stop_codons:
        codons = codons[:-1]

    return sum(1 for codon in codons if codon in stop_codons)


def calculate_pr2(seq):
    third_bases = []

    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        if len(codon) != 3:
            continue
        if codon in {"TAA", "TAG", "TGA"}:
            continue
        third_bases.append(codon[2])

    A3 = third_bases.count("A")
    T3 = third_bases.count("T")
    G3 = third_bases.count("G")
    C3 = third_bases.count("C")

    if (G3 + C3) == 0 or (A3 + T3) == 0:
        return np.nan, np.nan, A3, T3, G3, C3

    x = G3 / (G3 + C3)
    y = A3 / (A3 + T3)

    return x, y, A3, T3, G3, C3


def process_fasta(filepath, species, genome_type):
    rows = []

    for record in SeqIO.parse(filepath, "fasta"):
        seq = clean_sequence(record.seq)

        if len(seq) < MIN_CDS_LENGTH:
            continue

        gene = get_gene_name(record.description)
        internal_stops = count_internal_stops(seq)

        if EXCLUDE_INTERNAL_STOPS and internal_stops > 0:
            continue

        x, y, A3, T3, G3, C3 = calculate_pr2(seq)

        if np.isnan(x) or np.isnan(y):
            continue

        rows.append({
            "Species": species,
            "Genome": genome_type,
            "Gene": gene,
            "Header": record.description,
            "Length_bp": len(seq),
            "Internal_stops": internal_stops,
            "A3": A3,
            "T3": T3,
            "G3": G3,
            "C3": C3,
            "G3_over_G3plusC3": x,
            "A3_over_A3plusT3": y,
        })

    return rows


def plot_pr2(df, output):
    fig, axes = plt.subplots(2, 2, figsize=(13, 12), sharex=True, sharey=True)

    panel_order = list(FASTA_FILES.keys())
    panel_labels = ["(a)", "(b)", "(c)", "(d)"]

    for ax, (species, genome_type), label in zip(axes.flatten(), panel_order, panel_labels):
        sub = df[(df["Species"] == species) & (df["Genome"] == genome_type)].copy()

        ax.scatter(
            sub["G3_over_G3plusC3"],
            sub["A3_over_A3plusT3"],
            s=38,
            alpha=0.85,
        )

        ax.axhline(0.5, color="black", linestyle="--", linewidth=1)
        ax.axvline(0.5, color="black", linestyle="--", linewidth=1)

        if len(sub) > 0:
            sub["Distance_from_center"] = np.sqrt(
                (sub["G3_over_G3plusC3"] - 0.5) ** 2
                + (sub["A3_over_A3plusT3"] - 0.5) ** 2
            )

            if LABEL_ALL_GENES:
                label_df = sub.copy()
            elif LABEL_EXTREME_GENES:
                label_df = sub.nlargest(N_LABELS_PER_PANEL, "Distance_from_center").copy()
            else:
                label_df = pd.DataFrame()

            offsets = [(25, 15), (30, -15), (-40, 20), (-45, -20),
                       (35, 25), (35, -25), (-50, 30), (-50, -30)]

            for i, (_, row) in enumerate(label_df.iterrows()):
                dx, dy = offsets[i % len(offsets)]
                ax.annotate(
                    row["Gene"],
                    xy=(row["G3_over_G3plusC3"], row["A3_over_A3plusT3"]),
                    xytext=(dx, dy),
                    textcoords="offset points",
                    fontsize=8,
                    fontweight="bold",
                    arrowprops=dict(
                        arrowstyle="-",
                        color="black",
                        linewidth=0.8,
                        shrinkA=0,
                        shrinkB=2,
                    ),
                    ha="center",
                    va="center",
                )

        ax.set_title(f"{label} {species} - {genome_type}", fontsize=13, fontweight="bold")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xlabel("G3 / (G3 + C3)", fontsize=11)
        ax.set_ylabel("A3 / (A3 + T3)", fontsize=11)
        ax.tick_params(labelsize=10)

        for spine in ax.spines.values():
            spine.set_linewidth(1.2)

    plt.tight_layout()
    plt.savefig(output, dpi=300, bbox_inches="tight")
    plt.show()


# =========================================================
# MAIN
# =========================================================

all_rows = []
for (species, genome_type), fasta in FASTA_FILES.items():
    print(f"Processing: {species} | {genome_type} | {fasta}")
    all_rows.extend(process_fasta(fasta, species, genome_type))

df = pd.DataFrame(all_rows)
df.to_csv(OUTPUT_TABLE, sep="\t", index=False)

print("Output table:", OUTPUT_TABLE)
print(df.head())

plot_pr2(df, OUTPUT_FIGURE)

files.download(OUTPUT_TABLE)
files.download(OUTPUT_FIGURE)
