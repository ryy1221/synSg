### In this notebook we do whole genome sequencing analysis 
import pandas as pd
from os.path import join
data_path = '../../data/sequencing/WGS/03.Result_X202SC24094080-Z01-F001_cancer/result/Mutation/SNP/Annotation'
output_path = '../Analysis_output'
f_prefix = 'WGS_analysis'

# Function for calculate the allele frequency of heterozygous mutation 
def calc_allele_freq(df_, colname):
    # Process the allele frequency info of the dataframe
    split_df = df_[colname].str.split(':', expand = True)
    split_df.columns = ['Genotype', 'Allele Count', 'Read Depth', 'GQ', 'PL']

    # Operate only on heterozygous mutations
    het_idx = split_df[split_df['Genotype'] == '0/1'].index
    allele_count =pd.to_numeric(split_df.loc[het_idx,'Allele Count'].str.split(',').str[1])
    allele_depth = pd.to_numeric(split_df.loc[het_idx,'Read Depth'])
    return(allele_count.div(allele_depth))

merge_df = pd.DataFrame(columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'GeneName', 
            'Func', 'Gene', 'ExonicFunc', 'AAChange']) 

for file_name in ['C2_S1','C2_S2','C2_E1','C2_E2']: # Loop through the file names
    # Define the fields that we need so that we can locate the mutations
    print(f'Start processing {file_name}...')

    fields_e1 = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'GeneName', 
          'Func', 'Gene', 'ExonicFunc', 'AAChange', 'INFO', 'FORMAT', file_name]
    df = pd.read_csv(join(data_path, f'{file_name}.GATK.snp.annovar.hg38_multianno.xls'),\
            usecols=fields_e1, sep = '\t')
    mut_freq = calc_allele_freq(df, file_name)
    assert len(df) > len(mut_freq) # Make sure the function return same length
    df.loc[mut_freq.index, f'Allele Freq {file_name}'] = mut_freq
    assert len(df[df['GeneName'] == 'PSMB5'])>0 # make sure that the PSMB5 have mutation

    # Merge the dataframe so that we can compare the allele frequency between 2 samples
    merge_df.merge(df,
        how='left',
        on=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'GeneName', 
            'Func', 'Gene', 'ExonicFunc', 'AAChange']) 
    merge_df = merge_df.dropna()
    print(f'Finish processing {file_name}...')

# Save merged dataframe
print('Saving files...')
merge_df.to_csv(join(output_path,f'{f_prefix}_merged_AF.csv'))







