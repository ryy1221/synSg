import os
import pandas as pd
import itertools
from os.path import join
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc

sns.set_theme(style="white", font_scale=2,palette='viridis')
font = {'weight' : 'bold', 'family':'Nimbus Sans'}
rc('font', **font) 


# Function to filter sequences with significant read counts
def filter_sg_df(out_path, f):
    df = pd.read_csv(join(out_path, f), sep='\t')
    return df[(df['%Reads'] >= 1)]  # Only keep sequences >1% mapped


def find_modified_positions(aligned, reference):
    return [i for i, (a, r) in enumerate(zip(aligned, reference)) if a != r]  # Faster list comprehension


def get_sg_edit_info(sgRNA, df_summ):
    row = df_summ[df_summ['sgRNA'] == sgRNA]
    if row.empty:
        return [], None, None  # Handle missing sgRNA cases
    return row['edit_genome_pos'].values[0], row['gene_strand'].values[0], row['sg_strand'].values[0]


def confirm_edits(x, list_pred):
    return any(element in list_pred for element in x) if len(x) > 0 else True


def process_sample(name, day, replicate, info_dict, full_len_af_outpath, sg_af_outpath, df_summ):
    sgRNA, amplicon_seq, amplicon_start = info_dict[name].values()
    sample_idx = f'{name}_D{day}_{replicate}'

    # Locate files
    full_len_af_files = [i for i in os.listdir(full_len_af_outpath) if sample_idx in i]
    sg_af_files = [i for i in os.listdir(sg_af_outpath) if i.startswith(sample_idx)]

    # If there is not file then skip the sample
    if len(full_len_af_files) != 1 or len(sg_af_files) != 1:
        print(f"Skipping {sample_idx}: Missing files")
        return None

    df_full_len_af = pd.read_csv(join(full_len_af_outpath, full_len_af_files[0]), sep='\t')
    if day == 5:
        df_sg_af = filter_sg_df(sg_af_outpath, sg_af_files[0])
    else:
        df_sg_af = pd.read_csv(join(sg_af_outpath, sg_af_files[0]), sep = '\t')
        
    # If total read number is small then the file has problem
    if df_sg_af['#Reads'].sum()<100:
        print(f'Skipping {sample_idx}: Alignment low quality')
        return None

    # Alignment position extraction (safer)
    alignments = pairwise2.align.localms(df_full_len_af['Reference_Sequence'][0],
                                         df_sg_af['Reference_Sequence'][0], 2, -1, -0.5, -0.1)
    start_pos = alignments[0].start if alignments else 0

    df_sg_af = df_sg_af.sort_values('%Reads', ascending=False)
    sequences = df_sg_af['Aligned_Sequence']

    # Mutation analysis
    mutpos_mod = df_sg_af.apply(lambda row: find_modified_positions(row['Aligned_Sequence'], row['Reference_Sequence']), axis=1)
    gpos_mod = mutpos_mod.apply(lambda x: [i + amplicon_start + start_pos for i in x])

    pred_edits, gene_strand, sg_strand = get_sg_edit_info(sgRNA, df_summ)

    # Confirm edits
    edit_res = gpos_mod.apply(lambda x: confirm_edits(x, pred_edits))
    if all(edit_res):
        print(f'All edits in {sample_idx} are confirmed')

    return {
        'sequences': sequences.tolist(),
        '%reads': df_sg_af['%Reads'].tolist(),
        'mutation_pos': mutpos_mod.tolist(),
        'genome_pos':gpos_mod.tolist(),
        'predicted_pos':pred_edits
    }


def validate_edits(info_file, full_len_af_outpath, sg_af_outpath, df_summ):
    sample_list = info_file['Name'].unique()
        # sample_list = [i for i in sample_list if not i.startswith('CT')]
    print(f'Processing {len(sample_list)} samples...')

    days, reps = [5, 17], [1, 2, 3]
    info_dict = info_file[['Name', 'sgRNA', 'Amplicon', 'Amplicon_start_gPos']].drop_duplicates().set_index('Name').to_dict('index')

    dict_d5, dict_d17 = {}, {}

    for name in info_dict.keys():
        if name in sample_list:
            print(f'Processing {name}')
            dict_d5[name], dict_d17[name] = {}, {}

            for day, replicate in itertools.product(days, reps):
                try:
                    result = process_sample(name, day, replicate, info_dict, full_len_af_outpath, sg_af_outpath, df_summ)
                    if result is not None:  # Ensure result is valid before storing
                        (dict_d5 if day == 5 else dict_d17)[name][replicate] = result
                except KeyError:
                    print(f'Skipped {name} {replicate} due to missing data.')

    return dict_d5, dict_d17, sample_list

def filter_low_replica_sequences(data, threshold_ratio=0.5):
    """
    Removes replicates that have significantly fewer sequences than others.
    
    :param data: Dictionary containing sequence data per sample and replicate.
    :param threshold_ratio: The minimum fraction of sequences compared to the max replicate count to keep.
    :return: Filtered dictionary with only valid replicates.
    """
    filtered_data = {}

    for sample, replicates in data.items():
        # Count number of sequences in each replicate
        seq_counts = {rep: len(info['sequences']) for rep, info in replicates.items()}

        # Determine the max sequence count among replicates
        max_seq_count = max(seq_counts.values())

        # Identify replicates that meet the threshold criteria
        valid_replicates = {
            rep: info for rep, info in replicates.items()
            if seq_counts[rep] >= max_seq_count * threshold_ratio
        }

        # Store only if at least one valid replicate remains
        if valid_replicates:
            filtered_data[sample] = valid_replicates

    return filtered_data

def get_intersection_of_replicates(sequences_dict):
    # Get intersection of sequences across all replicates
    replicates_sequences = [set(sequences_dict[rep]['sequences']) for rep in sequences_dict.keys()]
    return set.intersection(*replicates_sequences)

def sequences_same_across_replicates(sequences_dict):
    # Extract sequences for each replicate and check if they are the same across all replicates
    replicates_sequences = [set(sequences_dict[rep]['sequences']) for rep in sequences_dict.keys()]
    return all(seq == replicates_sequences[0] for seq in replicates_sequences)

def filter_dict(dict_d5, dict_d17):
    dict_d5 = filter_low_replica_sequences(dict_d5)
    dict_d17 = filter_low_replica_sequences(dict_d17)

    filtered_dict_d5 = {}
    filtered_dict_d17 = {}

    for name in dict_d5.keys():
        if sequences_same_across_replicates(dict_d5[name]):
            filtered_dict_d5[name] = dict_d5[name]
        else:
            # Find intersection of sequences across all replicates in d5
            intersected_sequences = get_intersection_of_replicates(dict_d5[name])
            if intersected_sequences:
                filtered_dict_d5[name] = {rep: {
                    'sequences': [seq for seq in dict_d5[name][rep]['sequences'] if seq in intersected_sequences],
                    '%reads': [dict_d5[name][rep]['%reads'][i] for i, seq in enumerate(dict_d5[name][rep]['sequences']) if seq in intersected_sequences],
                    'mutation_pos': [dict_d5[name][rep]['mutation_pos'][i] for i, seq in enumerate(dict_d5[name][rep]['sequences']) if seq in intersected_sequences],
                    'genome_pos': [dict_d5[name][rep]['genome_pos'][i] for i, seq in enumerate(dict_d5[name][rep]['sequences']) if seq in intersected_sequences],
                    'predicted_pos':dict_d5[name][rep]['predicted_pos']
                    
                } for rep in dict_d5[name].keys()}

    # Filter dict_d17 based on the filtered_dict_d5
    
    for name in dict_d17.keys():
        if name in filtered_dict_d5:
            sequences_d5_all_replicates = list(dict.fromkeys(
        seq for replicate in filtered_dict_d5[name] for seq in filtered_dict_d5[name][replicate]['sequences']
    ))

            filtered_dict_d17[name] = {}

            for replicate in dict_d17[name].keys():
                sequences_d17 = dict_d17[name][replicate]['sequences']

                # Ensure that sequences in d17 follow the order of sequences in d5
                ordered_sequences_d5 = [seq for seq in sequences_d5_all_replicates if seq in sequences_d17]

                filtered_dict_d17[name][replicate] = {
                    'sequences': ordered_sequences_d5,
                    '%reads': [dict_d17[name][replicate]['%reads'][sequences_d17.index(seq)] for seq in ordered_sequences_d5],
                    'mutation_pos': [dict_d17[name][replicate]['mutation_pos'][sequences_d17.index(seq)] for seq in ordered_sequences_d5],
                    'genome_pos': [dict_d17[name][replicate]['genome_pos'][sequences_d17.index(seq)] for seq in ordered_sequences_d5],
                    'predicted_pos':dict_d17[name][replicate]['predicted_pos']
                }
    return(filtered_dict_d5, filtered_dict_d17)

def plot_sequences_with_highlights_and_reads(sequences, highlights, day5_reads, day17_reads):
    """
    Plot sequences with highlighted positions and % reads values for Day 5 and Day 17.

    Parameters:
    sequences (list of str): List of sequences to plot.
    highlights (list of list of int): List of lists containing positions to highlight for each sequence.
    day5_reads (list of list of float): List of lists containing % reads values for Day 5 for each sequence.
    day17_reads (list of list of float): List of lists containing % reads values for Day 17 for each sequence.
    """
    # Ensure all sequences are of equal length
    sequence_length = len(sequences[0])
    for seq in sequences:
        assert len(seq) == sequence_length, "All sequences must be of equal length"
    
    # Number of sequences
    num_sequences = len(sequences)
    
    # Create subplots
    fig, axes = plt.subplots(num_sequences, 1, figsize=(10, 0.5 * num_sequences), constrained_layout=True)
    
    if num_sequences == 1:
        axes = [axes]

    for i, (seq, positions, d5_reads, d17_reads) in enumerate(zip(sequences, highlights, day5_reads, day17_reads)):
        ax = axes[i]
        ax.set_xlim(-0.5, sequence_length + 8)  # Adjusted x limit to make space for reads
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.set_xticks([])  # Remove x-ticks

        # Remove spines
        for spine in ax.spines.values():
            spine.set_visible(False)

        # Highlight the positions with blue color
        for pos in positions:
            ax.axvspan(pos - 0.5, pos + 0.5, color='orange', alpha=0.5)

        # Plot the sequence letters
        for j, char in enumerate(seq):
            ax.text(j, 0.5, char, ha='center', va='center', fontsize=12, fontweight='bold')

        # Add grid lines
        ax.grid(True, axis='x', linestyle='--', color='gray', alpha=0.5)

        # Plot % reads values horizontally
        reads_text = f'D5: {", ".join(f"{d:.1f}%" for d in d5_reads)} | D17: {", ".join(f"{d:.1f}%" for d in d17_reads)}'
        ax.text(sequence_length + 1, 0.5, reads_text, ha='left', va='center', fontsize=12, color='black')

    return(fig)

def get_edit_df(data):
    rows = []

    for feature, replicates in data.items():
        for replicate, values in replicates.items():
            for i in range(len(values['sequences'])):
                row = {
                    'Feature': feature,
                    'Replicate': replicate,
                    'Sequence': values['sequences'][i],
                    '% Reads': values['%reads'][i],
                    'Mutation Positions': values['mutation_pos'][i],
                    'Genome Positions': values['genome_pos'][i] if 'genome_pos' in values else None,
                    'Predicted Edit': values.get('predicted_pos')
                }
                rows.append(row)
    df = pd.DataFrame(rows)
    return(df)


