# syn_sgFinder

This repository is a code hub for identifying sgRNA synonymous mutational consequences and analyzing sequencing data.

---

## 🧬 Overview

The project consists of two main components:

1. **sg_Finder** – Identify sgRNAs that produce only **synonymous mutations** given a list of genes. This pipeline first uses [CHOPCHOP](https://github.com/ryy1221/synSg/blob/dev/chopchop) for sgRNA selection. Then use ''sg_Finder.py'' to identify sgRNAs that only produce synonymous sgRNAs
2. **Analysis** – Downstream analysis of sequencing data, including processing, visualization, and figure generation for publication.

---

## 🔧 Prerequisites

- Python ≥ 3.8
- [CHOPCHOP command-line tool](https://bitbucket.org/valenlab/chopchop/src/master/) [Setup Instructures](https://github.com/ryy1221/synSg/blob/dev/chopchop)
- Conda (recommended) or pip for installing dependencies

Install dependencies:
```bash
conda env create -f environment.yml
conda activate synSg
```
---

## Folder Structure
```
synSg/
│
├── notebooks/               # Jupyter notebooks
│   ├── figures/             # Notebooks for generating publication figures
│   ├── sgFinder/            # sgRNA identification pipeline (15 bp window)
│   ├── analysis/            # Sequencing analysis and plots
│   ├── slurm_scripts/       # SLURM job scripts
│   └── legacy/              # Deprecated or old notebooks
│
├── src/                     # Python source code (core logic and utilities)
│
├── data/                    # Raw and processed data
│   ├── output/              # DESeq2, read depth, other processed results
│   ├── sg_Finder/           # Outputs from sgRNA identification
│   ├── sg_out/              # CHOPCHOP output files (by gene)
│   └── sequencing/          # Raw and aligned sequencing data
│
├── ref_fasta/               # Reference genome and annotation files
├── .gitignore               # Files to be excluded from version control
├── README.md                # Project overview and instructions
├── requirements.txt         # Python dependencies (if using pip)
└── environment.yml          # Conda environment configuration
```
