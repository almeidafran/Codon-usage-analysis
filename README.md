# Codon-usage-analysis
This repository contains Python scripts for codon usage analyses of plastid and mitochondrial protein-coding genes (PCGs). The workflows were designed for organelle genome datasets and can be easily adapted to other species by providing FASTA files containing coding sequences (CDS).

The repository includes scripts for:

- Relative synonymous codon usage (RSCU) analysis
- Effective number of codons (ENC) analysis
- ENC-GC3s plots
- PR2-bias plots
- Detection of internal stop codons
- Generation of publication-quality figures and summary tables

---

# Repository contents

## `codon_usage_general.py`

Performs codon usage and RSCU analyses from CDS FASTA files.

### Outputs
- Codon frequency table
- RSCU table
- Invalid sequence report
- Stacked codon usage plot

### Main analyses
- Codon counting
- RSCU calculation
- Codon frequency estimation
- Detection of problematic sequences

---

## `enc_analysis_general.py`

Performs ENC and GC3s analyses and generates ENC-GC3s plots.

### Outputs
- ENC/GC3s summary table
- ENC-GC3s figure
- Detection of internal stop codons

### Main analyses
- GC3s calculation
- Observed ENC estimation
- Expected ENC estimation
- ENC ratio calculation
- Wright’s theoretical ENC curve plotting

The script allows filtering genes containing internal stop codons before plotting.

---

## `pr2_analysis_general.py`

Performs PR2-bias analysis and generates PR2 plots.

### Outputs
- PR2 summary table
- PR2-bias figure

### Main analyses
- A3/(A3+T3) calculation
- G3/(G3+C3) calculation
- Detection of parity-rule deviations
- Identification of genes showing stronger codon bias

The script automatically extracts gene names from common annotation header formats.

---

# Input requirements

All scripts require CDS sequences in FASTA format.

Example:

```fasta
>gene-psbA
ATGGCT...
```

or

```fasta
>cds-rps14
ATGAAA...
```

Sequences should:

- contain only CDS regions
- be oriented in the correct strand
- preferably start with a start codon and end with a stop codon
- have lengths divisible by three

---

# Dependencies

The scripts were developed for Python 3 and require:

```bash
pip install biopython pandas matplotlib numpy
```

---

# Running the scripts

Example:

```bash
python codon_usage_general.py
python enc_analysis_general.py
python pr2_analysis_general.py
```

Input FASTA filenames can be edited directly in the scripts under:

```python
FASTA_FILES = {
    ...
}
```

---

# Applications

These scripts can be applied to:

- plastome analyses
- mitogenome analyses
- comparative codon usage studies
- codon bias analyses
- organelle genome evolution studies
- synonymous codon preference analyses

---

# Notes

- The scripts were optimized for organelle genomes but can also be adapted for nuclear CDS datasets.
- Genes containing internal stop codons can be automatically detected and excluded from downstream analyses.
- Figures are generated in high resolution (`300 dpi`) and are suitable for publication.

---

# Citation

If you use these scripts in your research, please cite the corresponding publication or acknowledge this repository.
