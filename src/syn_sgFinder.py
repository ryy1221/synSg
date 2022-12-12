# Some simple functions that I use when processing gene lists and sgRNAs
# Author: Yiyun
# Date: Sep, 2022

import pandas as pd
from os import listdir,makedirs
from os.path import join,exists
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
import pickle
import sys
sys.path.append('../src')
from utils import *

 # the output directory for sgRNA result from chopchop
chop_path = '../../chopchop' # There should always be a parallel chopchop repository
out_path = '../data/output/sg_Finder'

if not exists(out_path):
    makedirs(out_path)

# Read gene_table file, only gene_table file is common to every gene
df_gene_tab = pd.read_csv(join(chop_path,'gene_table','hg38.gene_table'), sep = '\t')
# clean genes on alternative chromosomes
df_gene_tab = df_gene_tab[~df_gene_tab.chrom.str.contains('_')]

class sgFinder_gene(): # the sgFinder class for individual genes

    def __init__(self, gene,sg_path = '../data/output/sg_out_enrich'):
        self.gene = gene
        self.sg_path = sg_path
        # The gene class contains all the sgRNAs
        sg_file = pd.read_csv(join(sg_path,gene,f'{gene}.txt'), sep = '\t', index_col = 0)
        self.sg_file = sg_file[(sg_file['MM0']==0)&(sg_file['MM1']==0)]
        print(f'Processing {gene}...Keep {len(self.sg_file)}/{len(sg_file)} low off-target sgRNAs...')

        if not exists(join(out_path,gene)):
            makedirs(join(out_path,gene))

    # Returns all transcripts and exon regions
    # All positions added 1 because the file has different counting method than CHOPCHOP
    def gene_cds_proc(self):
        global df_gene_tab
        df_sel_genetab = df_gene_tab[df_gene_tab.name2 == self.gene]
        self.df_sel_genetab = df_sel_genetab

        # The transcription strand is also very important
        strand = df_sel_genetab.strand.unique()
        raise_error(len(strand)==1, f"Not unique strand for {self.gene} , {strand}")

        self.gene_strand = strand[0]
        self.chrom = df_sel_genetab.chrom.unique()[0]

        # The sgRNAs are intersection of all isoforms
        # This means that they don't necessarily cause the same consequence in all transcripts
        # Therefore we need to make a dictionary to store all transcripts and their exons
        # Negative strans genes information are processed
        dict_gene_exons = {}
        dict_gene_cdst = {};dict_gene_cded = {}
        dict_gene_exonFrame = {}
        for rows in df_sel_genetab.index.values:
            df_sel_trpt = df_sel_genetab.loc[rows,:]
            trpt = df_sel_trpt['name']
            if self.gene_strand == '+':
                cds_st = df_sel_trpt.cdsStart+1
                cds_ed = df_sel_trpt.cdsEnd
                exon_st = [int(i)+1 for i in df_sel_trpt.exonStarts.strip(' ').split(',') if len(i) > 0]
                exon_ed = [int(i) for i in df_sel_trpt.exonEnds.strip(' ').split(',') if len(i) > 0]
                exon_frame = [int(i) for i in df_sel_trpt.exonFrames.strip(' ').split(',') if len(i) > 0]

                if not len(exon_st) == len(exon_ed) and len(exon_st) == len(exon_frame):
                    exit() # raise an error, define later
                exon_regions = list(zip(exon_st,exon_ed))
            elif self.gene_strand == '-': # the 1 position is the start position
                cds_st = df_sel_trpt.cdsEnd
                cds_ed = df_sel_trpt.cdsStart+1
                exon_st = [int(i) for i in df_sel_trpt.exonEnds.strip(' ').split(',') if len(i) > 0]
                exon_ed = [int(i)+1 for i in df_sel_trpt.exonStarts.strip(' ').split(',') if len(i) > 0]
                exon_frame = [int(i) for i in df_sel_trpt.exonFrames.strip(' ').split(',') if len(i) > 0]
                if not len(exon_st) == len(exon_ed) and len(exon_st) == len(exon_frame):
                    exit() # raise an error, define later
                exon_regions = list(zip(exon_st,exon_ed))[::-1]
                exon_frame = exon_frame[::-1]

            dict_gene_exons.update({trpt:exon_regions})
            dict_gene_cdst.update({trpt:cds_st})
            dict_gene_cded.update({trpt:cds_ed})
            dict_gene_exonFrame.update({trpt:exon_frame})

        self.dict_exons = dict_gene_exons
        self.dict_cdst = dict_gene_cdst;self.dict_cded = dict_gene_cded
        self.dict_exonframe =dict_gene_exonFrame

        pickle.dump(self, open(join(out_path,self.gene,f'{self.gene}.pkl'),'wb'))

class sgFinder_sg():
    # This class depends on sgRNA index in the CHOPCHOP output file
    def __init__(self, gene, sg_idx):
        self.sg_idx = sg_idx
        sgFinder_gene = pickle.load(open(join(out_path,gene,f'{gene}.pkl'),'rb'))

    # This function calcualte the BE window genomic position, only need the index of the sg dataframe
    def process_BE_window(self,sgFinder_gene):
        sg_strand = sgFinder_gene.sg_file.loc[self.sg_idx,'Strand']
        gloc = int(sgFinder_gene.sg_file.loc[self.sg_idx,'Genomic location'].split(':')[-1])
        target_seq = Seq(sgFinder_gene.sg_file.loc[self.sg_idx,'Target sequence'])
        sg_seq = target_seq[:-3] # exclude PAM sequence

        # calculate the genome position of the editing window, return the sequence as well
        be_window_sgseq = sg_seq[-20:-9] # The window region is always 5 nt, but position is not
        self.sg_seq = sg_seq
        self.be_window_sgseq = be_window_sgseq
        self.sg_strand = sg_strand
        self.gloc = gloc # All sgRNA start position are from the small end

        self.list_possib_pos_cbe, self.list_possib_pos_abe, if_extend = besite_gloc(self.be_window_sgseq,\
        self.sg_strand, self.gloc, sgFinder_gene.gene_strand,list(sgFinder_gene.dict_cdst.values()),list(sgFinder_gene.dict_cded.values()))

        # Only when edit base is at the start of sgRNA and extra sequence needed
        if if_extend:
            try:
                self.extend_sg, self.intron_idx = extend_seq(sgFinder_gene.gene,self.sg_strand,self.sg_seq, sgFinder_gene.sg_path)
            except ValueError as e:
                print(e)
            #     print(f'{sgFinder_gene.gene} sgRNA {self.sg_seq} has trouble aligning to fasta file')
                return(False)
        pickle.dump(self,open(join(out_path,sgFinder_gene.gene,f'{self.sg_idx}.pkl'),'wb'))
        return(True)

class sgFinder_transcript():
    # This class depends on sgRNA, and transcript ID
    def __init__(self, gene, sg_idx, trpt_ID):
        self.transcript_ID = trpt_ID
        sgFinder_gene = pickle.load(open(join(out_path,gene,f'{gene}.pkl'),'rb'))
        sgFinder_sg = pickle.load(open(join(out_path,gene,f'{sg_idx}.pkl'),'rb'))

    # This function calculate the distance between cds start site and the potential be site
    # Returns 1, 2, 3 meand 1st, 2nd, 3rd position in an amino acid
    # Then find the amino acids and it's mutation according to the position
    # Return a True/False indication synonymous or not
    def get_aa_consq(self, be_gpos, sgFinder_gene, sgFinder_sg, record_dict):
        cds_st = sgFinder_gene.dict_cdst[self.transcript_ID]
        exon_regions = sgFinder_gene.dict_exons[self.transcript_ID]
        exon_frame = sgFinder_gene.dict_exonframe[self.transcript_ID]
        # We just infer the amino acid position from the exon frame
        cut_site_exon_idx = find_be_exon(exon_regions,sgFinder_gene.gene_strand,be_gpos)
        cdst_exon_idx = find_be_exon(exon_regions,sgFinder_gene.gene_strand,cds_st)
        exonFrame = exon_frame[cut_site_exon_idx]

        if cut_site_exon_idx == cdst_exon_idx:
            start = cds_st
        else: start = exon_regions[cut_site_exon_idx][0]

        if sgFinder_gene.gene_strand == '+':
            pos = (be_gpos - start + 1 + exonFrame )%3
        else:
            pos = (start - be_gpos + exonFrame + 1 )%3

        if pos == 0: pos =3
        pos = int(pos)

        # Get the original and mutated amino acid
        if hasattr(sgFinder_sg,'extend_sg'):
            ext_sg = True
        else:ext_sg=False

        aa_gpos,aa_sgpos = find_be_sgpos(sgFinder_gene.gene_strand,\
        sgFinder_sg.sg_strand,be_gpos,pos,sgFinder_sg.gloc, extend_sg = ext_sg)
        aa_st_gpos = aa_gpos[0];aa_ed_gpos = aa_gpos[1];aa_st_sgpos = aa_sgpos[0];aa_ed_sgpos = aa_sgpos[1]


        if ext_sg:
            if (len(sgFinder_sg.intron_idx)>0) and (sgFinder_sg.extend_sg!=False):
                self.junction_flag = True
                if not all(i.isupper() for i in sgFinder_sg.extend_sg):
                    print('Finding mrna sequence in another exon...')
                    for keys in record_dict:
                        if keys.startswith(self.transcript_ID+'.'):
                            sgFinder_sg.extend_sg = find_junction_nuc(record_dict[keys].seq, \
                            sgFinder_sg.intron_idx, sgFinder_sg.extend_sg, sgFinder_sg.sg_strand)
            if sgFinder_sg.extend_sg != False:
                target_aa = sgFinder_sg.extend_sg[aa_st_sgpos:aa_ed_sgpos+1]
            else:return(False)
        else:
            target_aa = sgFinder_sg.sg_seq[aa_st_sgpos:aa_ed_sgpos+1]

        if sgFinder_gene.gene_strand == '+':
            if sgFinder_sg.sg_strand == '+':
                mut_aa = MutableSeq(target_aa)
                mut_aa[pos-1] = BE_mutate(mut_aa[pos-1])
            else:
                target_aa = target_aa.reverse_complement()
                mut_aa = MutableSeq(target_aa)
                mut_aa[pos-1] = BE_mutate(mut_aa[pos-1])

        elif sgFinder_gene.gene_strand == '-':
            if sgFinder_sg.sg_strand == '+':
                mut_aa = MutableSeq(target_aa)
                mut_aa[3-pos] = BE_mutate(mut_aa[3-pos])
                target_aa = target_aa.reverse_complement()
                mut_aa = Seq(mut_aa).reverse_complement()
            else:
                mut_aa = MutableSeq(target_aa)
                mut_aa[pos-1] = BE_mutate(mut_aa[pos-1])

        # For checking
        # print([be_gpos, pos,sgFinder_sg.sg_seq, (aa_st_gpos,aa_ed_gpos), (aa_st_sgpos,aa_ed_sgpos),target_aa,mut_aa])

        if str(mut_aa.translate()) == str(target_aa.translate()):
            return(True,f'{target_aa}>{mut_aa}',f'{str(target_aa.translate())}>{str(mut_aa.translate())}',f'{aa_st_gpos}-{aa_ed_gpos}')
        else:
            return(False,f'{target_aa}>{mut_aa}',f'{str(target_aa.translate())}>{str(mut_aa.translate())}',f'{aa_st_gpos}-{aa_ed_gpos}')

    def write_res_trsp(self, gene_class,sg_class, res, edit_pos):
        dict_res_sg_trsp = {'gene':gene_class.gene,'transcript':self.transcript_ID,'chrom':gene_class.chrom,\
        'gene_strand':gene_class.gene_strand,\
        'sg_strand':sg_class.sg_strand,'sgRNA':str(sg_class.sg_seq),'edit_window':str(sg_class.be_window_sgseq),\
        'edit_genome_pos':edit_pos,'Synonymous':str(res[0]),'Codon_Change':str(res[1]),\
        'AA_change':str(res[2]),'AA_pos':str(res[3]), 'junction':str(self.junction_flag)}
        df_res_sg_trsp = pd.DataFrame(dict_res_sg_trsp,index=[0])
        return(df_res_sg_trsp)

    def write_res_summ(self, gene_class,sg_class, lsyn, BE_type):
        if BE_type == 'ABE':
            dictionary = {'gene':gene_class.gene,'transcript':self.transcript_ID,\
            'sgRNA':str(sg_class.sg_seq),'ABE_locs':str(sg_class.list_possib_pos_abe),\
            'syn_or_not':lsyn}
            df_summ=pd.DataFrame.from_dict(dictionary,orient='index').transpose()
        elif BE_type =='CBE':
            dictionary = {'gene':gene_class.gene,'transcript':self.transcript_ID,\
            'sgRNA':str(sg_class.sg_seq),'CBE_locs':str(sg_class.list_possib_pos_cbe),\
            'syn_or_not':lsyn}
            df_summ=pd.DataFrame.from_dict(dictionary,orient='index').transpose()
        return(df_summ)

    def all_pos_consq(self, sgFinder_gene, sgFinder_sg, record):
        # ```
        # Create dataframe for CBE and ABE
        # One summary frame contains gene, sg, transcript, possib locations, syn_or_not(as list)
        # Another individual frame contains gene, sg, transcript, location, st en pos, amino acid
        # Loop through all possible mutate positions
        # For each position, run get_aa_consq
        # return True/False and save sg into ABE/CBE edit
        #
        self.junction_flag = False

        # First identify and save controls
        df_abe_control,df_cbe_control = write_CT_df(sgFinder_gene,sgFinder_sg)
        self.df_abe_control = df_abe_control;self.df_cbe_control = df_cbe_control

        list_syn_or_not_cbe = []
        dict_res_sg_trsp_cbe=pd.DataFrame()
        for idx, be_pos in enumerate(sgFinder_sg.list_possib_pos_cbe):
            syn_or_not = self.get_aa_consq(be_pos, sgFinder_gene, sgFinder_sg,record)
            if syn_or_not!=False:
                list_syn_or_not_cbe.append(syn_or_not[0])
                # If this is a synonymous mutation3
                df_res_sg_trsp = self.write_res_trsp(sgFinder_gene,sgFinder_sg, syn_or_not, be_pos)
                dict_res_sg_trsp_cbe = dict_res_sg_trsp_cbe.append(df_res_sg_trsp,ignore_index = True)
            else:
                break
        df_cbe_summ = self.write_res_summ(sgFinder_gene,sgFinder_sg, list_syn_or_not_cbe,'CBE')

        list_syn_or_not_abe = []
        dict_res_sg_trsp_abe=pd.DataFrame()
        for idx, be_pos in enumerate(sgFinder_sg.list_possib_pos_abe):
            syn_or_not = self.get_aa_consq(be_pos,sgFinder_gene, sgFinder_sg, record)
            if syn_or_not!=False:
                list_syn_or_not_abe.append(syn_or_not[0])
                # If this is a synonymous mutation
                df_res_sg_trsp = self.write_res_trsp(sgFinder_gene,sgFinder_sg, syn_or_not, be_pos)
                dict_res_sg_trsp_abe=dict_res_sg_trsp_abe.append(df_res_sg_trsp,ignore_index = True)
            else:
                break
        df_abe_summ = self.write_res_summ(sgFinder_gene,sgFinder_sg, list_syn_or_not_abe, 'ABE')

        self.df_abe = df_abe_summ
        self.df_cbe = df_cbe_summ
        self.df_cbe_ind = dict_res_sg_trsp_cbe
        self.df_abe_ind = dict_res_sg_trsp_abe
