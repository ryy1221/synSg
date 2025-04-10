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

def RevComplment(s): #Case-senstive reverse complmnentarty conversion
    comp = []
    for i in s[::-1]:
        if i == 'a':
            comp.append('t')
        if i == 'A':
            comp.append('T')
        if i == 't':
            comp.append('a')
        if i == 'T':
            comp.append('A')
        if i == 'c':
            comp.append('g')
        if i == 'C':
            comp.append('G')
        if i == 'g':
            comp.append('c')
        if i == 'G':
            comp.append('C')
    return ''.join(comp)



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

def exon_boundary(dict_exons,cdst,cded,gene_strand):
    exon_neighbor_pos = []
    if gene_strand == '+':
        for i in dict_exons:
            exon_neighbor_pos.extend([i[0],i[0]+1,i[0]+2,i[1],i[1]-1,i[1]-2])
            exon_neighbor_pos.extend([cdst,cdst+1,cdst+2])
            exon_neighbor_pos.extend([cded,cded-1,cded-2])
    else:
        for i in dict_exons:
            exon_neighbor_pos.extend([i[0],i[0]-1,i[0]-2,i[1],i[1]+1,i[1]+2])
            exon_neighbor_pos.extend([cdst,cdst-1,cdst-2])
            exon_neighbor_pos.extend([cded,cded+1,cded+2])
    return(exon_neighbor_pos)

def filter_sg_boundary(boundary_pos, sg_start_gloc, sg_strand):
    if sg_strand == '-':
        map_boundary_site = sg_start_gloc + 22
    elif sg_strand == '+':
        map_boundary_site = sg_start_gloc
#     print(sg_start_gloc, map_boundary_site, boundary_pos)

    if map_boundary_site in boundary_pos:
#         print('sgRNA is at boundary')
#         print(sg_start_gloc, map_boundary_site, boundary_pos)
        return(True)
    else:
        return(False)
    
def read_summ_files(name):
    name_abe = 'ABE_synsg.csv';name_cbe = 'CBE_synsg.csv'
    name_CTabe = 'ABE_CT.csv';name_CTcbe = 'CBE_CT.csv'
    name_abe_detail = 'df_abe_detail.csv';name_cbe_detail = 'df_cbe_detail.csv'
    out_path = '../data/output/sg_Finder'

    f_abe = pd.read_csv(join(out_path,name,name_abe), index_col = 0)
    f_cbe = pd.read_csv(join(out_path,name,name_cbe), index_col = 0)
    f_CTabe = pd.read_csv(join(out_path,name,name_CTabe), index_col = 0)
    f_CTcbe = pd.read_csv(join(out_path,name,name_CTcbe), index_col =0)
    return(f_abe,f_cbe,f_CTabe,f_CTcbe)

def is_stop(x):
    if all(i.split('>')[-1]=='*' for i in x) and any(i.split('>')[-1]=='*' for i in x):
        return(True)
    else:return(False)

def is_missense(x):
    if all(i.split('>')[0]!=i.split('>')[1] for i in x) and any(i.split('>')[0]!=i.split('>')[1] for i in x):
        if not any('*' in i for i in x):
            return(True)
        else:return(False)
    else:return(False)

def gene_name(x):
    if len(set(x)) == 1:
        return(list(x)[0])
    else:raise ValueError
    

def get_cons(df):
    df_group = df.groupby('sgRNA').agg({'gene':lambda x:gene_name(x),'Codon_Change':list,'AA_change':list,'AA_pos':list,'Synonymous':list})
    syn_sg = df_group[df_group['Synonymous'].apply(lambda x: all(x)&any(x))]
    stop_sg = df_group[df_group['AA_change'].apply(lambda x: is_stop(x))]
    mis_sg = df_group[df_group['AA_change'].apply(lambda x: is_missense(x))]
    return(syn_sg,stop_sg,mis_sg)


def sg_consequence(name,list_gene,tag,name_detail):
    list_df_stop = [];list_df_mis = [];list_df_syn = []
    out_path = '../data/output/sg_Finder'
    for gene in list_gene:
        df = pd.read_csv(join(out_path,gene, name, name_detail), index_col = 0)
        if not df.empty:
            syn_df,stop_df,mis_df = get_cons(df)
            if not stop_df.empty:list_df_stop.append(stop_df)
            if not mis_df.empty: list_df_mis.append(mis_df)
            if not syn_df.empty:list_df_syn.append(syn_df)
    if len(list_df_stop) > 0:
        df_stop = pd.concat(list_df_stop,axis = 0)
    else: df_stop = pd.DataFrame()
    if len(list_df_mis)>0:
        df_mis = pd.concat(list_df_mis,axis = 0)
    else:df_mis = pd.DataFrame()
    if len(list_df_syn)>0:
        df_syn = pd.concat(list_df_syn,axis = 0)
    else:df_syn = pd.DataFrame()
    
    return(df_stop, df_mis,df_syn)
            


def make_analyze_df(list_df):
    df_alys = pd.concat(list_df, axis = 1).fillna(0)
    df_alys['sum_11'] = df_alys['non_CT11']+df_alys['CT11']
    df_alys['sum_13'] = df_alys['non_CT13']+df_alys['CT13']
    df_alys['sum_15'] = df_alys['non_CT15']+df_alys['CT15']
    return(df_alys)


def get_control_sg(dfabeCT,dfcbeCT,tag):
    tot_abe_CT = dfabeCT.groupby('gene')['sgRNA'].apply(lambda x: len(set(x)))
    tot_abe_CT.name='CT'+tag
    tot_cbe_CT = dfcbeCT.groupby('gene')['sgRNA'].apply(lambda x: len(set(x)))
    tot_cbe_CT.name='CT'+tag
    return(tot_abe_CT,tot_cbe_CT)
    

def get_genome_seq(sequence,list_index, sg_gloc,window_len, gene_strand, sgRna_strand):
    sequence_idx = list_index.index(sg_gloc)
    shift = int((window_len - 11)/2)
    center_pos = int((window_len-1)/2)

    if gene_strand == '+':
        if sgRna_strand == '+':
            seq = sequence[sequence_idx-shift:sequence_idx+window_len-shift]
        else:
            seq = sequence[sequence_idx+12-shift:sequence_idx+12+window_len-shift].reverse_complement()
    else:
        if sgRna_strand == '+':
            seq = sequence[sequence_idx-window_len+1+shift:sequence_idx+1+shift].reverse_complement()
        else:seq = sequence[sequence_idx-17-center_pos:sequence_idx-17+center_pos+1]
    return(seq)

# we only extend those sgRNAs that has edit site at the start position(0 or 1 distance)
def besite_gloc(sg_seq, sg_strand, sg_start_gloc, gene_strand, shift):
    pos_cbe = []
    pos_abe = []
    extend_seq = False
#     sted_gsites = gsite_sted(dict_cdst,dict_cded,gene_strand)

    for idx, letter in enumerate(sg_seq):
        if letter =='C':
            if sg_strand == '-':
                possib_edit_gsite = sg_start_gloc + 22 - idx + shift
            elif sg_strand == '+':
                possib_edit_gsite = sg_start_gloc + idx - shift
            pos_cbe.append(possib_edit_gsite)
            if idx in [0,1,len(sg_seq)-1,len(sg_seq)]:
                extend_seq = True

        if letter =='A':
            if sg_strand == '-':
                possib_edit_gsite = sg_start_gloc + 22 - idx + shift
            elif sg_strand == '+':
                possib_edit_gsite = sg_start_gloc + idx - shift
            pos_abe.append(possib_edit_gsite)
            if idx in [0,1,len(sg_seq)-1,len(sg_seq)]:
                extend_seq = True
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
        new_l.append(F+oligo1+oligo2+RevComplment(R))
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

def sel_CT(df, number):
    sel_df = pd.DataFrame()
    for genes in list(df.gene.unique()):
        CTsg = df[df.gene == genes].head(1)
        sel_df = sel_df.append(CTsg, ignore_index = True)
    if len(sel_df) < number:
        remain_df = df.loc[random.sample(list(df.index),  number-len(sel_df))]
    sel_df = pd.concat([sel_df,remain_df], axis = 0)
    sel_df = sel_df[~sel_df['sgRNA'].duplicated()]
    return(sel_df)

def calc_best_primerTM(F1,R1, noprimer_list):
    LikestPrimer = [F1, R1]
    PrimerHeteroTM2 = []
    for i in LikestPrimer:
        Hetero_List = []
        for oligo in noprimer_list:
            oligo = str(oligo)
            Hetero_List.append(p3.calcHeterodimer(i, oligo).tm)
        PrimerHeteroTM2.append([i, max(Hetero_List)])

    return(F1,R1,PrimerHeteroTM2,p3.calcHeterodimerTm(F1,str(Seq(R1).reverse_complement())))
