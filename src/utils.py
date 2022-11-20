# This is the util functions for syn_sgFinder.py
from Bio import SeqIO
from Bio.Seq import Seq
from os.path import join
from os.path import join
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import random
import primer3 as p3
import numpy as np
from tqdm import tqdm


ivan_path = '../data/fromIvan'
primer_path = '../data/output/primers/'



def BE_mutate(nuc): # Input is the nucleotide
    if nuc == 'C':
        mut = "T"
    elif nuc == "A":
        mut = 'G'
    elif nuc == 'T':
        mut = 'C'
    elif nuc == 'G':
        mut = 'A'
    return(mut)

def raise_error(condition, msg, error_type=ValueError):
    if not condition:
        raise error_type(msg)

def gsite_sted(dict_cdst, dict_cded, gene_strand):
    if gene_strand == '+':
        start_codon_sites = list(set([i for i in dict_cdst] + [i+1 for i in dict_cdst] + [i+2 for i in dict_cdst]))
        end_codon_sites = list(set([i for i in dict_cded] + [i-1 for i in dict_cded] + [i-2 for i in dict_cded]))
    else:
        start_codon_sites = list(set([i for i in dict_cdst] + [i-1 for i in dict_cdst] + [i-2 for i in dict_cdst]))
        end_codon_sites = list(set([i for i in dict_cded] + [i+1 for i in dict_cded] + [i+2 for i in dict_cded]))
    return(start_codon_sites+end_codon_sites)

# we only extend those sgRNAs that has edit site at the start position(0 or 1 distance)
def besite_gloc(sg_seq, sg_strand, sg_start_gloc, gene_strand, dict_cdst, dict_cded):
    pos_cbe = []
    pos_abe = []
    extend_seq = False
    sted_gsites = gsite_sted(dict_cdst,dict_cded,gene_strand)

    for idx, letter in enumerate(sg_seq):
        if letter =='C':
            if sg_strand == '-':
                possib_edit_gsite = sg_start_gloc + 22 - idx
            elif sg_strand == '+':
                possib_edit_gsite = sg_start_gloc + idx
            pos_cbe.append(possib_edit_gsite)
            if idx in [0,1]:
                extend_seq = True

        if letter =='A':
            if sg_strand == '-':
                possib_edit_gsite = sg_start_gloc + 22 - idx
            elif sg_strand == '+':
                possib_edit_gsite = sg_start_gloc + idx
            pos_abe.append(possib_edit_gsite)
            if idx in [0,1]:
                extend_seq = True

    if any([i for i in pos_cbe if i in sted_gsites]):
        pos_cbe = []
        print(f'{possib_edit_gsite} is in start/stop codon ...skipping{sg_seq}')
    if any([i for i in pos_abe if i in sted_gsites]):
        pos_abe = []
        print(f'{possib_edit_gsite} is in start/stop codon ...skipping{sg_seq}')

    return(pos_cbe, pos_abe, extend_seq)

def find_junction_nuc(record, intron_index, extend_sgseq, sg_strand):
    max_idx = max(intron_index)
    sg_strand_upper = extend_sgseq[max_idx+1::]
    if sg_strand == '-':
        start_idx = record.find(sg_strand_upper.reverse_complement())
        extend_sg = record[start_idx:start_idx+22]
        extend_sg = extend_sg.reverse_complement()
    else:
        start_idx = record.find(sg_strand_upper)
        extend_sg = record[start_idx-max_idx-1:start_idx+20+max_idx]

    if not start_idx>=0:
        print('Finding extra sequence in mRNA failed, skip...')
        return(False)
    else:return(extend_sg)

def extend_seq(gene,sg_strand,sg_seq, sg_path):
    for record in SeqIO.parse(join(sg_path,gene,'gene_file.fa'), "fasta"):
        if sg_strand == '-':
            start_idx = record.seq.find(sg_seq.reverse_complement())
            extend_sg = record.seq[start_idx:start_idx+22]
            extend_sg = extend_sg.reverse_complement()
        else:
            start_idx = record.seq.find(sg_seq)
            extend_sg = record.seq[start_idx-2:start_idx+20]

        # For those aligned to intron/exon junctions:
        intron_idx = []
        for idx, i in enumerate(extend_sg):
            if not i.isupper():intron_idx.append(idx)

        # Below are a few check points that things can go wrong
        raise_error(start_idx>=0, \
        f'{sg_seq} in gene {gene} not aligned, skip...')
        raise_error(len(extend_sg)==22, \
        f'{sg_seq} in gene {gene} need sequence not included in fasta file...' )

    return(extend_sg, intron_idx)

def find_be_exon(dict_exon, gene_strand,be_pos):
    for idx, exon_reg in enumerate(dict_exon):
        if gene_strand == '+':
            if be_pos >= exon_reg[0] and be_pos <= exon_reg[1]:
                be_exon_idx = idx
        else:
            if be_pos >= exon_reg[1] and be_pos <= exon_reg[0]:
                be_exon_idx = idx
    return(be_exon_idx)

def find_be_sgpos(gene_strand,sg_strand,be_pos,aa_pos,sg_gloc,extend_sg = False):
    if extend_sg:
        sg_len = 24
    else:sg_len = 22

    if gene_strand == '+':
        aa_st_gpos = be_pos-aa_pos+1
        aa_ed_gpos = be_pos+3-aa_pos
        if sg_strand == '+':
            aa_st_sgpos = aa_st_gpos-sg_gloc+(sg_len-22)
            aa_ed_sgpos = aa_ed_gpos-sg_gloc+(sg_len-22)
        elif sg_strand == '-':
            aa_st_sgpos = sg_len-(aa_ed_gpos-sg_gloc)
            aa_ed_sgpos = sg_len-(aa_st_gpos-sg_gloc)

    elif gene_strand == '-':
        aa_st_gpos = be_pos+aa_pos-1
        aa_ed_gpos = be_pos-3+aa_pos
        if sg_strand == '+':
            aa_st_sgpos = aa_ed_gpos-sg_gloc+2
            aa_ed_sgpos = aa_st_gpos-sg_gloc+2
        elif sg_strand == '-':
            aa_st_sgpos = sg_len-(aa_st_gpos-sg_gloc)
            aa_ed_sgpos = sg_len-(aa_ed_gpos-sg_gloc)

    return([aa_st_gpos,aa_ed_gpos],[aa_st_sgpos,aa_ed_sgpos])

def make_CT_df(BE_type,gene_class,sg_class):
    df_control_gene = pd.DataFrame({'gene':gene_class.gene,\
                        'chrom':gene_class.chrom,\
                         'gene_strand':str(gene_class.gene_strand),\
                             'sg_strand':sg_class.sg_strand,\
                        'sgRNA':str(sg_class.sg_seq),\
                       'sgRNA_start':str(sg_class.gloc),\
                        'BE_type':BE_type}, index = [0])
    return(df_control_gene)

def write_CT_df(gene_class,sg_class):
    df_control_ABE = pd.DataFrame()
    df_control_CBE = pd.DataFrame()
    if all([i!='A' for i in sg_class.be_window_sgseq]):
        df_control_gene = make_CT_df('ABE', gene_class,sg_class)
        df_control_ABE = df_control_ABE.append(df_control_gene,ignore_index = True)
    elif all([i!='C' for i in sg_class.be_window_sgseq]):
        df_control_gene = make_CT_df('CBE',gene_class,sg_class)
        df_control_CBE = df_control_CBE.append(df_control_gene,ignore_index = True)
    return(df_control_ABE,df_control_CBE)


RE = Seq('CGTCTC')
RE_l = Seq('CGTCTCACACC')
RE_r = Seq('GTTTCGAGACG')

BadSites = ['CGTCTC', 'GAGACG', 'AAAA', 'CCCC', 'TTTT', 'GGGG']

def attach_RE(x):
    new_x = RE_l+x+RE_r if x[0] in ['A','G'] else RE_l+'G'+x+RE_r
    return(new_x)

def count_RE(input_seq):
    global RE
    if (input_seq.count(RE) ==1) and (input_seq.count(RE.reverse_complement()) == 1):
        return(True)
    else:return(False)

def readnCT(fname, sel = True):
    lsg = []
    with open(join(ivan_path,fname),'r') as f:
        for lines in f:
            sg = lines.split('\n')[0]
            lsg.append(sg)
    df = pd.DataFrame()
    df['sgRNA'] = lsg
    df['gene'] = fname.split('_')[1].split('.')[0][2:]

    df = df.reset_index()
    df['index'] = df['index'].astype(str)
    df['ID'] = df['gene']+'_'+df['index']
    df = df[['sgRNA','ID','gene']]
    if sel:
        df = df.loc[random.sample(list(df.index),40)]
    return(df)

def make_sg_df(lsg,name):
    df = pd.DataFrame()
    df['sgRNA'] = lsg
    df['gene'] = name
    df = df.reset_index()
    df['index'] = df['index'].astype(str)
    df['ID'] = df['gene']+'_'+df['index']
    df = df[['sgRNA','ID','gene']]
    return(df)

def getpCT(fname,sel_range,label_tab):
    K562_Fit = pd.read_csv(join(ivan_path,fname), sep='\t',index_col = 0)
    df_K562_topess = K562_Fit.sort_values('neg|lfc')[sel_range[0]:sel_range[1]]
    df_K562_topenr = K562_Fit.sort_values('neg|lfc')[-sel_range[1]:-sel_range[0]]
    K562_topess = df_K562_topess.index
    K562_topenr = df_K562_topenr.index

    K562_topess = [i for i in K562_topess if 'Control' not in i]
    K562_topenr = [i for i in K562_topenr if 'Control' not in i]

    list_posCT = label_tab.loc[K562_topenr]['sgRNA Seq'].to_list() + label_tab.loc[K562_topess]['sgRNA Seq'].to_list()
    return(list_posCT)

def attach_group_primer(l,F,R):
    l= list(zip(l[::2], l[1::2]))
    new_l = []
    for i in l:
        oligo1 = attach_RE(i[0])
        oligo2 = attach_RE(i[1])
        new_l.append(F+oligo1+oligo2+R)
    return(new_l)

def get_noprimer_list(BE_list):
    # Does strand matter here? Add RE site
    RE_BE_list = [attach_RE(x) for x in BE_list]
    BE_pass_list = []
    for idx,oligo in enumerate(RE_BE_list):
        if count_RE(oligo):
            BE_pass_list.append(BE_list[idx])
    rBE_pass_list=random.sample(list(BE_pass_list), len(BE_pass_list))
    rBE_noprimer_list = attach_group_primer(rBE_pass_list,'','')
    return(rBE_noprimer_list,rBE_pass_list)

def write_noprime_list(fname,l):
    with open(join(primer_path,fname), 'w') as f:
        for i in range(0,len(l)):
            f.write('>'+ 'Seq' + str(i))
            f.write('\n')
            f.write(str(l[i]))
            f.write('\n')

def read_blast_res(fname):
    blastf= open(join(primer_path,fname), "r").read().splitlines()
    primers = []
    for i in range(len(blastf)):
        if 'Query=' in blastf[i]:
            Qseq = blastf[i][-20:]
            if 'No hits found' in blastf[i+5]: #if no matches
                primers.append(Qseq)
                continue
    return(primers)

def calc_TM(primers,noprime_list):
    PrimerHeteroTM = []
    for i in tqdm(primers):
        Hetero_List_S = []
        Hetero_List_A = []
        for oligo in noprime_list :
            oligo = str(oligo)
            senseTM = p3.calcHeterodimer(i, oligo).tm
            antiTM = p3.calcHeterodimer(i, str(Seq(oligo).reverse_complement())).tm
            Hetero_List_S.append(senseTM)
            Hetero_List_A.append(antiTM)
        PrimerHeteroTM.append([i, max(Hetero_List_S), max(Hetero_List_A)])
    return(PrimerHeteroTM)

def sel_primers(thrd, heteroTM_list):
    GoodPrimer = []
    for i,m,n in heteroTM_list:
        if m < thrd:
            if n < thrd:
                GoodPrimer.append(i)

    GoodPrimerCheck = []
    for i in GoodPrimer:
        PrimeA = p3.calcTm(i)
        for m in GoodPrimer:
            PrimeB = p3.calcTm(m)
            HD_S = p3.calcHeterodimer(i, m).tm
            HD_A = p3.calcHeterodimer(i, str(Seq(m).reverse_complement())).tm
            HD_sum = (HD_S + HD_A)
            Tm_diff = abs((PrimeA-PrimeB))
            GoodPrimerCheck.append([i,PrimeA,m,PrimeB,HD_S, HD_A, HD_sum, Tm_diff])

    coln = ['Primer A', 'Tm A', 'Primer B', 'Tm B', 'Sense Tm', 'AntiSense B Tm', 'SenseAnti Sum', 'Tm Dif']
    GoodPrimer_DF = pd.DataFrame(GoodPrimerCheck, columns=coln)
    GoodPrimer_DF = GoodPrimer_DF.sort_values(by='SenseAnti Sum')

    return(GoodPrimer_DF)

def org_sgdf(p_df):
    input_df = pd.read_csv(p_df, index_col = 0)
    input_df = input_df.reset_index()
    input_df['index'] = input_df['index'].astype(str)
    input_df['ID'] = input_df['gene']+'_'+input_df['index']
    input_df = input_df[['sgRNA','ID','gene']]
    return input_df

def sel_empty_window(df):
    sel_df = pd.DataFrame()
    for genes in list(df.gene.unique()):
        CTsg = df[df.gene == genes].head(1)
        sel_df = sel_df.append(CTsg, ignore_index = True)
    if len(sel_df) < 120:
        remain_df = df.loc[random.sample(list(df.index), 250-len(sel_df))]
    sel_df = sel_df.append(remain_df,ignore_index = True)
    sel_df = sel_df[~sel_df.duplicated()]
    return(sel_df)

def calc_best_primerTM(GoodPrimer_DF, noprimer_list):
    F1 = GoodPrimer_DF.head(1)['Primer A'].values[0]
    R1 = GoodPrimer_DF.head(1)['Primer B'].values[0]
    LikestPrimer = [F1, R1]
    PrimerHeteroTM2 = []
    for i in LikestPrimer:
        Hetero_List = []
        for oligo in noprimer_list:
            oligo = str(oligo)
            Hetero_List.append(p3.calcHeterodimer(i, oligo).tm)
        PrimerHeteroTM2.append([i, max(Hetero_List)])

    return(F1,R1,PrimerHeteroTM2,p3.calcHeterodimerTm(F1,str(Seq(R1).reverse_complement())))
