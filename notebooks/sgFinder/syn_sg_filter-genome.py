### This notebook contains code that filter sgRNAs that only edit synonymous mutations
### The main difference from previous version is it fetch the BE window sequences from annotation file directly
import os
import pandas as pd
from tqdm import tqdm
from itertools import product
from os.path import exists, join
from os import makedirs
from pandas.errors import EmptyDataError
from multiprocessing import Pool, Manager, Lock
import pickle
from Bio import SeqIO
import re

import sys
sys.path.append('../')
from src import syn_sgFinder
from itertools import chain,product

out_path = '../../data/output_genome/sg_Finder'
mrna_path = '../../data/MANE'

# Load all expressed transcripts
# genelist = pd.read_csv('../chopchop/mart_export.txt', header =None)
list_exp_trsp = pickle.load(open('../../data/22Q2_gene_effect/expressed_transcripts.pkl','rb'))
# Parse mRNA sequence dictionary, for use when sequence aligned to intron-exon junctions
record_dict = SeqIO.to_dict(SeqIO.parse(join(mrna_path,'MANE.GRCh38.v1.0.refseq_rna.fna'), "fasta"))
lock = Lock()

def process_gene(gene, length, tag, sg_p, out_path, record_dict, list_exp_trsp, shared_empty_genes):
    df_synsg_abe = pd.DataFrame()
    df_synsg_cbe = pd.DataFrame()
    """Process a single gene safely, skipping errors and returning DataFrames."""
    try:
        if not exists(join(out_path, gene, tag)):
            makedirs(join(out_path, gene, tag))

        # Get gene transcript, sequence information
        sg_gene = syn_sgFinder.sgFinder_gene(gene, sg_path=sg_p)
        sg_gene.gene_cds_proc(record_dict, list_exp_trsp)  # Filter non-expressed transcripts

        # Initialize DataFrames
        df_synsg_abe_dt = pd.DataFrame()
        df_synsg_cbe_dt = pd.DataFrame()

        for idx, transcript in product(sg_gene.sg_file.index, sg_gene.dict_exons.keys()):
            sg_trsp = syn_sgFinder.sgFinder_transcript(sg_gene.gene, idx, transcript)
            check_point = sg_trsp.process_BE_window(sg_gene, {}, window_length=length)

            if check_point:
                sg_trsp.all_pos_consq(sg_gene)
                df_synsg_abe = df_synsg_abe.append(sg_trsp.df_abe,ignore_index = True)
                df_synsg_cbe = df_synsg_cbe.append(sg_trsp.df_cbe,ignore_index = True)
                
                df_synsg_abe_dt = df_synsg_abe_dt.append(sg_trsp.df_abe_ind,ignore_index = True)
                df_synsg_cbe_dt = df_synsg_cbe_dt.append(sg_trsp.df_cbe_ind,ignore_index = True)

        df_synsg_abe_dt.to_csv(join(out_path,gene,tag,f'df_abe_detail.csv'))
        df_synsg_cbe_dt.to_csv(join(out_path,gene,tag,f'df_cbe_detail.csv'))

        return df_synsg_abe, df_synsg_cbe
    
    except Exception as e:
        with lock:
            print(f"Skipping gene {gene} due to error: {e}")
        shared_empty_genes.append(gene)
        return pd.DataFrame(), pd.DataFrame()  # Return empty DataFrames


def run_sgFinder(input_genes, length, tag='ess', sg_p='../data/output_genome/sg_out'):
    """Parallelized sgFinder with error handling and tqdm progress bar."""
    out_path = '../../data/output_genome/sg_Finder'
    
    manager = Manager()
    shared_empty_genes = manager.list()
    shared_dict_filter = manager.dict()  # Shared dictionary for filtering sgRNAs

    # Parallel execution
    with Pool(processes=os.cpu_count()) as pool:
        results = pool.starmap(
            process_gene,
            [(gene, length, tag, sg_p, out_path, record_dict,list_exp_trsp, shared_empty_genes) for gene in input_genes]
        )
    return (results)

# for name in ['aa','ab','ac','ad','ae','af','ag','ah','ai','aj','ak','al']
for name in ['ae','af','ag','ah','ai','aj','ak','al']:
    with open(f'../../chopchop/whole_genome_genes/mart_export_{name}', "r") as f:
        genes = [line.strip() for line in f if line.strip()]
    res = run_sgFinder(genes, length = 15, tag = 'ess_15' )
    df_synsg_abe = pd.concat([res[i][0] for i in range(0, len(res), 2) ], ignore_index=True)
    df_synsg_cbe = pd.concat([res[i+1][0] for i in range(0, len(res)-1, 2)], ignore_index=True)
    ABE_synsg = df_synsg_abe[df_synsg_abe['syn_or_not'].apply(lambda x: isinstance(x, list) and all(x) & any(x))]
    CBE_synsg = df_synsg_cbe[df_synsg_cbe['syn_or_not'].apply(lambda x: isinstance(x, list) and all(x) & any(x))]
    ABE_synsg.to_csv(join(out_path, f'ABE_synsg_{name}.csv'), index=False)
    CBE_synsg.to_csv(join(out_path, f'CBE_synsg_{name}.csv'), index=False)