# =========================================================
# ENC-GC3s analysis for organellar CDS FASTA files
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

OUTPUT_TABLE = "ENC_GC3s_results.tsv"
OUTPUT_FIGURE = "ENC_GC3s_plot.png"

# Plot/filter options
EXCLUDE_INTERNAL_STOPS = False   # True = remove CDS with internal stop codons
MIN_CDS_LENGTH = 300             # minimum CDS length in bp
N_LABELS_PER_PANEL = 15          # number of gene labels per panel

# =========================================================
# CODON TABLE: standard genetic code
# =========================================================

CODON_TO_AA = {
    "TTT":"Phe", "TTC":"Phe", "TTA":"Leu", "TTG":"Leu",
    "TCT":"Ser", "TCC":"Ser", "TCA":"Ser", "TCG":"Ser",
    "TAT":"Tyr", "TAC":"Tyr", "TAA":"Stop", "TAG":"Stop",
    "TGT":"Cys", "TGC":"Cys", "TGA":"Stop", "TGG":"Trp",
    "CTT":"Leu", "CTC":"Leu", "CTA":"Leu", "CTG":"Leu",
    "CCT":"Pro", "CCC":"Pro", "CCA":"Pro", "CCG":"Pro",
    "CAT":"His", "CAC":"His", "CAA":"Gln", "CAG":"Gln",
    "CGT":"Arg", "CGC":"Arg", "CGA":"Arg", "CGG":"Arg",
    "ATT":"Ile", "ATC":"Ile", "ATA":"Ile", "ATG":"Met",
    "ACT":"Thr", "ACC":"Thr", "ACA":"Thr", "ACG":"Thr",
    "AAT":"Asn", "AAC":"Asn", "AAA":"Lys", "AAG":"Lys",
    "AGT":"Ser", "AGC":"Ser", "AGA":"Arg", "AGG":"Arg",
    "GTT":"Val", "GTC":"Val", "GTA":"Val", "GTG":"Val",
    "GCT":"Ala", "GCC":"Ala", "GCA":"Ala", "GCG":"Ala",
    "GAT":"Asp", "GAC":"Asp", "GAA":"Glu", "GAG":"Glu",
    "GGT":"Gly", "GGC":"Gly", "GGA":"Gly", "GGG":"Gly",
}

AA_TO_CODONS = {}
for codon, aa in CODON_TO_AA.items():
    if aa != "Stop":
        AA_TO_CODONS.setdefault(aa, []).append(codon)

# Met and Trp are excluded from ENC calculation because they have only one codon.
AA_TO_CODONS_ENC = {aa: codons for aa, codons in AA_TO_CODONS.items() if len(codons) > 1}

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
    header = str(header).replace("gene=", "gene:")
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

    for prefix in ["chloe_", "blatx_", "blastx_", "cds-blatx_", "cds-blastx_"]:
        gene = gene.replace(prefix, "")
    gene = re.sub(r"_fragment$", "", gene)
    gene = re.sub(r"_fragment_[0-9]+$", "", gene)
    gene = re.sub(r"_[0-9]+_[0-9]+$", "", gene)
    gene = re.sub(r"_[0-9]+$", "", gene)
    return gene


def codon_counts(seq):
    counts = {codon: 0 for codon, aa in CODON_TO_AA.items() if aa != "Stop"}
    total_codons = 0
    internal_stops = 0

    codons = [seq[i:i+3] for i in range(0, len(seq), 3) if len(seq[i:i+3]) == 3]
    if codons and CODON_TO_AA.get(codons[-1]) == "Stop":
        codons = codons[:-1]

    for codon in codons:
        aa = CODON_TO_AA.get(codon)
        if aa is None:
            continue
        if aa == "Stop":
            internal_stops += 1
            continue
        counts[codon] += 1
        total_codons += 1

    return counts, total_codons, internal_stops


def calculate_gc3s(seq):
    codons = [seq[i:i+3] for i in range(0, len(seq), 3) if len(seq[i:i+3]) == 3]
    if codons and CODON_TO_AA.get(codons[-1]) == "Stop":
        codons = codons[:-1]

    third_bases = [c[2] for c in codons if c in CODON_TO_AA and CODON_TO_AA[c] != "Stop"]
    if not third_bases:
        return np.nan
    return sum(1 for b in third_bases if b in ["G", "C"]) / len(third_bases)


def calculate_enc(counts):
    """
    Approximate ENC using Wright's classic formula:
    ENC = 2 + 9/F2 + 1/F3 + 5/F4 + 3/F6
    """
    f_by_k = {2: [], 3: [], 4: [], 6: []}

    for aa, codons in AA_TO_CODONS_ENC.items():
        k = len(codons)
        n = sum(counts[c] for c in codons)
        if n <= 1:
            continue

        sum_pi2 = sum((counts[c] / n) ** 2 for c in codons)
        f_value = (n * sum_pi2 - 1) / (n - 1)
        if k in f_by_k:
            f_by_k[k].append(f_value)

    mean_f = {k: np.mean(v) if v else np.nan for k, v in f_by_k.items()}
    if any(np.isnan(mean_f[k]) or mean_f[k] <= 0 for k in [2, 3, 4, 6]):
        return np.nan

    enc = 2 + 9 / mean_f[2] + 1 / mean_f[3] + 5 / mean_f[4] + 3 / mean_f[6]
    return max(20, min(61, enc))


def expected_enc(gc3s):
    if pd.isna(gc3s):
        return np.nan
    denominator = gc3s**2 + (1 - gc3s)**2
    if denominator == 0:
        return np.nan
    return 2 + gc3s + 29 / denominator


def process_fasta(filepath, species, genome_type):
    rows = []

    for record in SeqIO.parse(filepath, "fasta"):
        seq = clean_sequence(record.seq)
        if len(seq) < MIN_CDS_LENGTH:
            continue

        gene = get_gene_name(record.description)
        counts, total_codons, internal_stops = codon_counts(seq)
        if EXCLUDE_INTERNAL_STOPS and internal_stops > 0:
            continue

        gc3s = calculate_gc3s(seq)
        enc_obs = calculate_enc(counts)
        enc_exp = expected_enc(gc3s)
        if pd.isna(gc3s) or pd.isna(enc_obs):
            continue

        enc_ratio = (enc_exp - enc_obs) / enc_exp if enc_exp else np.nan

        rows.append({
            "Species": species,
            "Genome": genome_type,
            "Gene": gene,
            "Length_bp": len(seq),
            "Codons": total_codons,
            "GC3s": gc3s,
            "ENC_obs": enc_obs,
            "ENC_exp": enc_exp,
            "ENC_ratio": enc_ratio,
            "Internal_stops": internal_stops,
        })

    return rows


def plot_enc_gc3s(df, output):
    fig, axes = plt.subplots(2, 2, figsize=(13, 10), sharex=True, sharey=True)
    panel_order = list(FASTA_FILES.keys())
    panel_labels = ["(a)", "(b)", "(c)", "(d)"]
    x = np.linspace(0.001, 0.999, 500)
    y = [expected_enc(i) for i in x]

    for ax, (species, genome_type), label in zip(axes.flatten(), panel_order, panel_labels):
        sub = df[(df["Species"] == species) & (df["Genome"] == genome_type)].copy()

        ax.plot(x, y, color="black", linewidth=1.5)
        ax.scatter(sub["GC3s"], sub["ENC_obs"], s=35, alpha=0.85)

        if len(sub) > 0:
            label_df = pd.concat([
                sub.nsmallest(5, "ENC_obs"),
                sub.nlargest(5, "ENC_obs"),
                sub.nlargest(5, "ENC_ratio"),
            ]).drop_duplicates()

            for _, row in label_df.iterrows():
                ax.text(row["GC3s"] + 0.006, row["ENC_obs"] + 0.4,
                        row["Gene"], fontsize=8, fontweight="bold")

        ax.set_title(f"{label} {species} - {genome_type}", fontsize=13, fontweight="bold")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 65)
        ax.set_xlabel("GC3s")
        ax.set_ylabel("ENC")
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

plot_enc_gc3s(df, OUTPUT_FIGURE)

files.download(OUTPUT_TABLE)
files.download(OUTPUT_FIGURE)
