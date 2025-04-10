# syn_sgFinder
This repository is a code hub for identifying sgRNA mutational consequences and notebook analysis

Prerequisites for sgRNA identification:
sgRNAs are identified via CHOPCHOP. Instructions 
1) sg_Finder: Identify sgRNAs that only have synonymous consequences given a list of genes
2) Analysis: Sequencing data analysis for manuscript

Folder Structure:
synSg/
│
├── notebooks/               # Jupyter notebooks
│   ├── figures/             # Notebooks for generating figures
│   ├── sgFinder             # Pipeline for identifying qualifying sgRNAs in 15bp window
│   ├── analysis/            # Sequencing analysis and plots
│   ├── slurm_scripts/       # Slurm job scripts
│   └── legacy/              # Deprecated or old notebooks
│
├── src/                     # All Python scripts
│
├── data/                    # Raw or intermediate data (can be .gitignored)
│   ├── output/              # Library sequencing analysis outputs (DESeq2, RD) 
│   ├── sg_Finder/           # Output for identifying qualifying sgRNAs in 15bp window
│   └── sg_out/              # Chopchop output files by gene
│   └── sequencing/          # Sequencing data 
│
├── ref_fasta/               # Reference sequences
├── README.md
├── .gitignore
├── requirements.txt     
└── environment.yml      
