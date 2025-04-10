import pandas as pd
import re
from Bio.Seq import Seq
import os
from Bio import SeqIO
from be_predict_bystander import predict as bystander_model
from be_predict_efficiency import predict as be_efficiency_model
from tqdm import tqdm

'''
This python script is the wrapper for running be-hive on sgRNA candidates
Input: sgRNA
Returns: efficiency(logit score based on default 30%), predicted mutation (gDNA position)
'''

chop_path = '/storage/group/epo2/default/yur97/github/synSg/data/sg_out'
finder_path = '/storage/group/epo2/default/yur97/github/synSg/data/sg_Finder'

# Read ABE and CBE complete set(in sequencing)
df_ABE = pd.read_csv('/storage/group/epo2/default/yur97/github/synSg/data/complete_ABE_df.csv')
df_CBE = pd.read_csv('/storage/group/epo2/default/yur97/github/synSg/data/complete_CBE_df.csv')

# Read and process MANE file
mane = pd.read_csv('/storage/group/epo2/default/yur97/synSg/MANE/MANE.GRCh38.v1.0.refseq_genomic.gtf', sep = '\t', header = None).dropna()
mane['gene_name'] = mane[8].str.split(';', expand = True)[0].apply(
    lambda x: re.search(r'gene_id "(.*?)"', x).group(1))


def find_sg(x,sg):
    '''
    Map the sgRNA sequence with PAM to the sgRNA sequence
    '''
    if x[0:-3] == sg:
        return(True)
    else:
        return(False)
    
def get_mane_gene_df(gene_name):
    '''
    Filter the MANE file by gene, then return the dataframe and gene strand
    '''
    global mane
    
    mane_gene = mane[mane['gene_name'] == gene_name]
    gene_strand = mane_gene[6].unique()[0]
    
    return(mane_gene, gene_strand)

def init_models(BE, celltype='mES'):
    """
    Initialize the models for base editing efficiency and bystander effects.
    Args:
        BE (str): Base editor type (e.g., 'ABE', 'CBE').
        celltype (str): Cell type for the models.
    """
    bystander_model.init_model(base_editor=BE, celltype=celltype)
    be_efficiency_model.init_model(base_editor=BE, celltype=celltype)
    
def predict_efficiency(x):
    """
    Predicts base editing efficiency using the loaded model.
    
    Args:
        x (pd.Series): A row of the DataFrame containing sequences and related information.
    
    Returns:
        tuple: Predicted logit score, predicted read fraction, and input sequence.
    """
    before = x['before']
    after = x['after']
    PAM_sg = x['sg_PAM']
    sg_strand = x['sg_strand']
    gene_strand = x['gene_strand']
    
    output_seq = before[40:] + PAM_sg + after[:7]
    pred_d = be_efficiency_model.predict(output_seq, mean=0.3)
    
    return pred_d['Predicted logit score'], pred_d['Predicted fraction of sequenced reads with base editing activity'], output_seq

def process_sequence_files(df, mutation_category):
    global chop_path
    """
    Processes sgRNA sequences and metadata, filtering for the specified mutation category.
    
    Args:
        df (pd.DataFrame): DataFrame containing mutation data.
        mutation_category (str): Mutation category to filter by ('MIS' or 'SYN').
        chop_path (str): Path to the sequence data.
    
    Returns:
        pd.DataFrame: Processed DataFrame containing sgRNA sequences and metadata.
    """
    l_sg_df = []
    
    for i in tqdm(df[df['categ'].isin(mutation_category)].index):
        sg = df.loc[i, 'sgRNA']
        gene = df.loc[i, 'gene']
        
        if sg.startswith('g') and len(sg) > 20:
            sg = sg[1:]
        
        gene_file = [i for i in SeqIO.parse(os.path.join(chop_path, gene, 'sequence.fa'), 'fasta')]
        sg_file = pd.read_csv(os.path.join(chop_path, gene, f'{gene}.txt'), sep='\t')
        
        sg_df = sg_file[sg_file['Target sequence'].apply(lambda x: find_sg(x, sg))][['Target sequence', 'Genomic location', 'Strand']]
        sg_df.columns = ['sg_PAM', 'Genomic location', 'sg_strand']
        sg_df['gene'] = gene
        
        mane_gene, gene_strand = get_mane_gene_df(gene)
        sg_df['gene_strand'] = gene_strand
        sg_df['converted_sgPMA'] = sg_df.apply(lambda x: str(Seq(x['sg_PAM']).reverse_complement()), axis=1)
        
        dict_sg = []
        for i in gene_file:
            before, after, sequence_strand, PAM_sg = i.name.split(':')[3:7]
            if PAM_sg in sg_df['converted_sgPMA'].tolist():
                sg_idx = sg_df.index[sg_df['converted_sgPMA'] == PAM_sg][0]
                dict_ = {'index': sg_idx, 'before': before, 'PAM_sg_sequence': PAM_sg, 'after': after, 'sequence_strand': sequence_strand}
                dict_sg.append(dict_)
        
        sg_search = pd.DataFrame(dict_sg)
        sg_search = sg_search[~sg_search.duplicated()].set_index('index')
        
        sg_df = pd.concat([sg_df, sg_search], axis=1)
        l_sg_df.append(sg_df)
    
    return pd.concat(l_sg_df)

def merge_predictions(df, mutation_category, BE):
    global chop_path
    """
    Merges sgRNA context with predicted editing efficiency.
    
    Args:
        df (pd.DataFrame): DataFrame containing mutation data.
        mutation_category (str): Category of mutation to filter.
        BE (str): Base editor type (e.g., 'ABE', 'CBE').
        chop_path (str): Path to the sequence data.
    
    Returns:
        pd.DataFrame: DataFrame with predictions merged.
    """
    l_sg_df = process_sequence_files(df, mutation_category)
    results_series = l_sg_df.apply(lambda x: predict_efficiency(x), axis=1)
    
    results_df = pd.DataFrame(results_series.tolist(), index=l_sg_df.index)
    results_df.columns = ['predicted logit score', 'predicted read fraction', '50bp_input']
    
    l_sg_df = pd.concat([l_sg_df, results_df], axis=1)
    return l_sg_df

def strip_g(x):
    """
    Removes the extra 'g' at the start of the sgRNA sequence if present.
    
    Args:
        x (pd.Series): A row containing the 'sgRNA' sequence.
    
    Returns:
        str: The sgRNA sequence without the extra 'g'.
    """
    if x['sgRNA'].startswith('g'):
        return x['sgRNA'].strip('g')
    return x['sgRNA']

def get_df_fraction(l_sg_df, df, mutation_category, df_reads_path, column_names):
    """
    Combines the predicted read fraction with the actual read counts from sequencing.
    
    Args:
        l_sg_df (pd.DataFrame): DataFrame containing sgRNA data with predicted read fraction.
        df (pd.DataFrame): DataFrame containing sgRNA sequences.
        mutation_category (str): The mutation category ('MIS' or 'SYN').
        df_reads_path (str): Path to the sequencing data CSV.
    
    Returns:
        pd.DataFrame: Combined DataFrame with adjusted read counts.
    """
    l_sg_df['sg'] = l_sg_df.apply(lambda x: x['sg_PAM'][:20], axis=1)
    df_fraction = l_sg_df[['sg', 'predicted read fraction']].set_index('sg')
    
    df['sg_sequence'] = df[['name', 'sgRNA']].apply(lambda x: strip_g(x), axis=1)
    df = df[df['categ'] == mutation_category][['name', 'sg_sequence','categ']].set_index('sg_sequence')
    
    df_reads = pd.read_csv(df_reads_path, index_col=0)
    df_fraction_ = pd.concat([df_fraction, df], axis=1).set_index('name')
    df_fraction_.columns = ['frac']
    df_reads_frac = pd.concat([df_fraction_, df_reads], axis=1).dropna()
    
    # Adjusted read counts
    df_reads_frac['frac_1'] = df_reads_frac[column_names[0]] * df_reads_frac['frac']
    df_reads_frac['frac_2'] = df_reads_frac[column_names[1]] * df_reads_frac['frac']
    df_reads_frac['sgRNA'] = df_reads_frac.apply(lambda x: strip_g(x), axis=1)
    
    return df_reads_frac

def get_mutation_count(df, mutation_category, finder_path, BE):
    """
    Retrieves mutation count from gene sequence files for the specified mutation category.
    
    Args:
        df (pd.DataFrame): DataFrame with sgRNA sequences.
        mutation_category (str): The mutation category ('MIS' or 'SYN').
        finder_path (str): Path to the finder data files.
        BE (str): Base editor type (e.g., 'ABE', 'CBE').
    
    Returns:
        pd.DataFrame: DataFrame with sgRNA sequences and mutation counts.
    """
    l_sg_out = []
    
    for i in tqdm(df[df['categ'] == mutation_category].index):
        sg = df.loc[i, 'sgRNA']
        gene = df.loc[i, 'gene']
        if sg.startswith('g') and len(sg) > 20:
            sg = sg[1:]
        finder = pd.read_csv(os.path.join(finder_path, gene, 'ess_15', f'df_{BE.lower()}_detail.csv'), index_col=0)
        l_sg_out.append(finder[finder['sgRNA'] == sg])
    
    df_sg_out = pd.concat(l_sg_out)
    df_mutation_count = df_sg_out.groupby('sgRNA').size().reset_index()
    df_mutation_count.columns = ['sgRNA', 'mut_n']
    
    return df_mutation_count

def merge_dataframes(df_reads_frac, df_mutation_count):
    """
    Merges the read fraction DataFrame with mutation counts.
    
    Args:
        df_reads_frac (pd.DataFrame): DataFrame with predicted and actual read fractions.
        df_mutation_count (pd.DataFrame): DataFrame with sgRNA mutation counts.
    
    Returns:
        pd.DataFrame: Merged DataFrame with sgRNA, read fractions, and mutation counts.
    """
    df_final = pd.merge(df_reads_frac, df_mutation_count, on='sgRNA')
    return df_final

def filter_within_std(df_final, col_name):
    """
    Filters rows in the DataFrame where the values fall within mean ± 2 std.
    
    Args:
        df_final (pd.DataFrame): DataFrame to filter.
        col_name (str): Column name to apply the filtering.
    
    Returns:
        pd.Index: Index of rows that satisfy the filtering condition.
    """
    m_3std = (df_final[col_name].mean() - 2 * df_final[col_name].std(), df_final[col_name].mean() + 2 * df_final[col_name].std())
    filter_idx = df_final[(df_final[col_name] > m_3std[0])].index
    return filter_idx

def estimate_mutation_number(df_final):
    """
    Estimates the total number of mutations by limiting mutations to 2 when higher.
    
    Args:
        df_final (pd.DataFrame): DataFrame with mutation counts.
    
    Returns:
        int: Total mutation number.
    """
    mutation_n = 0
    for i in df_final['mut_n']:
        if i <= 2:
            mutation_n += i
        else:
            mutation_n += 2
    return mutation_n


    