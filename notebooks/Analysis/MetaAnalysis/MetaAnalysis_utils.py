import pandas as pd
from os.path import join
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import matplotlib.pyplot as plt
import numpy as np

DATA_PATH = '../../data/MetaAnalysis'

def make_subdf(df, hannah_annot, cell_line):
    """
    Create a subset of the dataframe by merging relevant columns for the specified cell line.
    
    Parameters:
        df (pd.DataFrame): Main dataframe containing sgRNA and replicate data.
        hannah_annot (pd.DataFrame): Annotation dataframe containing sgRNA details.
        cell_line (str): The cell line for which to subset the dataframe.
        
    Returns:
        pd.DataFrame: Merged dataframe containing sgRNA counts and annotation.
    """
    sub_df = pd.concat([
        df[['sgRNA', 'pDNA', f'{cell_line}_RepA', f'{cell_line}_RepB']].set_index('sgRNA'),
        hannah_annot[['sgRNA sequence', 'Gene symbol', '# edits', '# silent edits', 'Mutation category']].set_index('sgRNA sequence')
    ], axis=1)
    
    return sub_df


def sig_perc(df, col_name, category, cut_off=0.05, category_col = None):
    """
    Calculate the number and percentage of significant hits within a given category.
    
    Parameters:
        df (pd.DataFrame): Dataframe containing DESeq results and categories.
        category (str): The mutation category to filter by.
        cut_off (float): The p-adjusted value cutoff for significance.
        
    Returns:
        tuple: A tuple containing the number of significant hits and the total number in the category.
    """
    if category == 'all':
        total_count = len(df)
        sig_count = len(df[df[col_name]<cut_off])
        return sig_count, total_count
    else:
        df_category = df[df[category_col] == category]
        total_count = len(df_category)
        significant_count = len(df_category[df_category[col_name] < cut_off])
        return significant_count, total_count

def sgRNA_categ(mutation_category):
    """
    Categorize the mutation type based on its severity.
    
    Parameters:
        mutation_category (str): The mutation category string for the sgRNA.
        
    Returns:
        str: The mutation category, ordered by severity, or 'No edits' if missing.
    """
    if pd.isna(mutation_category):  # Handle missing data
        return 'No edits'
    
    # Priority-based mutation category assignment
    if 'Nonsense' in mutation_category:
        return 'Nonsense'
    elif 'Splice-donor' in mutation_category or 'Splice-acceptor' in mutation_category or 'Splice site' in mutation_category:
        return 'Splice site'
    elif 'Missense' in mutation_category:
        return 'Missense'
    elif 'Intron' in mutation_category:
        return 'Intron'
    elif 'Silent' in mutation_category:
        return 'Synonymous'
    elif 'UTR' in mutation_category:
        return 'UTR'
    else:
        return 'No edits'

def analyze_DESeq(cell_line, df):
    """
    Perform DESeq2 analysis for a given cell line.
    
    Parameters:
        cell_line (str): The cell line for which to perform analysis.
        df (pd.DataFrame): The dataframe containing sgRNA count data.
        
    Returns:
        pd.DataFrame: DESeq2 results dataframe.
    """
    cols = ['pDNA', f'{cell_line}_RepA', f'{cell_line}_RepB']
    
    # Create metadata for DESeq analysis
    clinical_df = pd.DataFrame(index=cols, columns=['Day'])
    clinical_df['Condition'] = ['Plasmid', 'Assay', 'Assay']
    
    try:
        dds = DeseqDataSet(
            counts=df.T.loc[cols, :],
            metadata=clinical_df.loc[cols, :],
            design_factors='Condition',
            refit_cooks=True
        )
        dds.deseq2()  # Fit dispersion and normalization
        
        # Perform DESeq analysis and summarize results
        stat_res = DeseqStats(dds, contrast=('Condition', 'Assay', 'Plasmid'))
        stat_res.summary()
        
        return stat_res.results_df
    except Exception as e:
        print(f"Error in DESeq analysis: {e}")
        return None


def plot_all_cell_lines(data):
    cell_lines = list(data.keys())
    categories = list(data[cell_lines[0]].keys())
    bar_width = 0.15  # Width of each bar
    index = np.arange(len(categories))  # Base x locations for each category

    plt.figure(figsize=(12, 8))

    for i, cell_line in enumerate(cell_lines):
        percentages = [(data[cell_line][category]['sig_n'] / data[cell_line][category]['tot_n']) * 100 for category in categories]
        plt.bar(index + i * bar_width, percentages, bar_width, label=cell_line)
    ax = plt.gca()
    plt.set_ylim(0,100)
    plt.title('Percentage of Significant Hits by sgRNA Category for All Cell Lines')
    plt.ylabel('Percentage of Significant Hits (%)')
    plt.xticks(index + bar_width * (len(cell_lines) - 1) / 2, categories, rotation=45)
    plt.legend(title="Cell Lines")
    plt.tight_layout()
    
def plot_ratio(dict_perc, labels, label_index, denominator_index ):
    bar_width = 0.15  # Width of each bar
    index = np.arange(len(labels))  # Base x locations for each category

    plt.figure(figsize=(8, 6))

    for i, cell_line in enumerate(cell_lines):
        denominator = dict_perc[cell_line][denominator_index]
        edit_ratios = [dict_perc[cell_line][idx] / denominator for idx in label_index]  # Misssense, Synonymous, Nonsense categories
        bars = plt.bar(index + i * bar_width, edit_ratios, bar_width, label=cell_line)

        for bar, ratio in zip(bars, edit_ratios):
            height = bar.get_height()
            plt.text(bar.get_x() + bar.get_width() / 2, height, f'{ratio:.2f}', ha='center', va='bottom')

    ax = plt.gca()
    ax.set_ylim(0,10)
    plt.title('Significant Hits Ratio to No Edits sgRNAs')
    plt.ylabel('Significant Hits Ratio')
    plt.xticks(index + bar_width * (len(cell_lines) - 1) / 2, labels, rotation=45)
    plt.legend(title="Cell Lines")
    plt.tight_layout()
    
def calc_L2FC(df, cell_line):
    """
    Calculates the log2 fold change (L2FC) between the cell line replicates and pDNA.
    Filters out sgRNAs where log-normalized reads per million (RPM) of pDNA is 
    greater than 3 standard deviations from the mean.

    Parameters:
        df (pd.DataFrame): The dataframe containing the data.
        cell_line (str): The cell line name for which to calculate L2FC.

    Returns:
        np.array: An array containing the average L2FC for replicates A and B.
    """
    
    # Log-normalize the reads per million (RPM) for pDNA
    l2rd0 = np.log2((df['pDNA'] / df['pDNA'].sum() * 1e6) + 1)
    
    # Calculate mean and standard deviation of the log-normalized pDNA
    mean_l2rd0 = l2rd0.mean()
    std_l2rd0 = l2rd0.std()

    # Filter out sgRNAs where pDNA log-normalized RPM is > 3 standard deviations from the mean
    filter_condition = np.abs(l2rd0 - mean_l2rd0) <= 3 * std_l2rd0
    filtered_df = df[filter_condition]

    # Log-normalize the reads per million (RPM) for the replicates
    l2rdA = np.log2((filtered_df[f'{cell_line}_RepA'] / filtered_df[f'{cell_line}_RepA'].sum() * 1e6) + 1)
    l2rdB = np.log2((filtered_df[f'{cell_line}_RepB'] / filtered_df[f'{cell_line}_RepB'].sum() * 1e6) + 1)

    # Recalculate l2rd0 for the filtered dataframe
    l2rd0_filtered = np.log2((filtered_df['pDNA'] / filtered_df['pDNA'].sum() * 1e6) + 1)

    # Calculate L2FC for replicates A and B
    l2fcA = l2rdA - l2rd0_filtered
    l2fcB = l2rdB - l2rd0_filtered
    
    filtered_df[f'{cell_line}_L2FC'] = np.mean(np.array([l2fcA, l2fcB]), axis=0)

    # Return the average of L2FC for replicates A and B
    return filtered_df

def calc_z_score(df, cell_line, _controls):
    # Filter the dataframe to find the entries where the 'Gene symbol' is 'negative control'
    targeting_control_df = df[df['Gene symbol'].isin(_controls)]

    # Calculate mean and standard deviation of the L2FC values for the negative controls
    mean_control = targeting_control_df[f'{cell_line}_L2FC'].mean()
    std_control = targeting_control_df[f'{cell_line}_L2FC'].std()

    # Now calculate the Z-score for the entire dataframe based on the control mean and std
    df['Z_score'] = (df[f'{cell_line}_L2FC'] - mean_control) / std_control
    return(df)

# Process dataframe so that there is a category column
def split_index(x):
    if x.split('_')[1] == 'NT':
        return('Non-targeting')
    else:
        return(x.split('_')[0])
    
def process_df(df):
    # Ensure the index is reset to make the operation clear
    df = df.reset_index()
    # Create the 'categ' column based on the split index for both DataFrames
    df['categ'] = df['index'].apply(lambda x: split_index(x))
    # If needed, set the index back to its original state
    df = df.set_index('index')
    
    return(df)
    
