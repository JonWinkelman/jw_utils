#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 14:39:32 2022

@author: jonwinkelman
"""

from jw_utils import parse_gff as pgf
import pandas as pd
from jw_utils import parse_fasta as pfa
import numpy as np
from jw_utils import dna_utils as du
from plotly import graph_objects as go

def make_arrow_trace(protein, par, strand, prod, start, end, br,**kwargs): #color=None,
    if kwargs['color']:
        color = kwargs['color']
    else:
        color = 'rgb(100,100,100)'

    'return a trace with an arrow drawn for the input feature'
    
    trace = go.Scatter(
            x = [start, start, br, br, end,br,br,start,start],
            y = [0, 0.25,0.25,0.5,0,-0.5,-0.25,-0.25,0],
            mode = 'lines',
            fill='toself',
            line = {'width': 1.5,
                    'color':color},
            name = protein,
            text = f'{protein}<br><br>{par}<br>strand: {strand}<br>{prod}<br>start: {start}<br>end: {end}',
            hoverinfo = 'text')
    return trace


def flip_block(df):
    """
    Flip a genomic block (reverse coordinates and strand)
    for a dataframe containing columns:
    ['ID', 'start', 'end', 'strand', 'product', 'gene_ID', 'HOG']
    """
    
    # Determine boundaries of the block
    block_min = df[['start', 'end']].min().min()
    block_max = df[['start', 'end']].max().max()

    df = df.copy()  # to avoid modifying original dataframe

    # Flip coordinates
    df['new_starts'] = block_max - df['end']
    df['new_ends']   = block_max - df['start']

    # Flip strand
    df['new_strand'] = df['strand'].map({'+': '-', '-': '+'}).fillna(df['strand'])

    # Replace old columns with new ones
    df['start'] = df['new_starts']
    df['end']   = df['new_ends']
    df['strand'] = df['new_strand']

    # Drop temp columns
    df = df.drop(columns=['new_starts', 'new_ends', 'new_strand'])

    # Sort by the new starts coordinate
    df = df.sort_values('start')

    return df

def get_df_for_feature(gff_fp, feature, fts=15, feature_type='cds'):
    """Return a local df containing features near the passed feature
    
    feature (str): protein ID or gene locus tag (no "cds-" or "gene-" prefix btw...)
    fts (int): num features to show in map on each side of of feature of interest 
    """

 
    df = pgf.make_simple_annot_df(gff_fp, start_end=True, contig=True)
    df.index=df.index.str.replace('gene-', '')
    df['protein_ID']=df['protein_ID'].str.replace('cds-', '')
    if feature_type == 'cds':
        df = df.reset_index().set_index('protein_ID')
    if feature in df.index:  
        f_index = df.index.get_loc(feature)
        if f_index>=fts and f_index<= (len(df.index)-fts):
            trimmed_df = df.iloc[(f_index-fts):(f_index+fts),:]
        
        elif f_index-fts<0:
            trimmed_df = df.iloc[:(f_index+fts),:]
        elif f_index + fts > df.shape[0]:
            trimmed_df = df.iloc[(f_index-fts):,:]
    else:
        print(f'{feature} not in in the dataframe index, e.g. index = {df.index[0]}')
        trimmed_df = df.iloc[:10,:]
    if trimmed_df.loc[feature, 'strand'] == '-': 
        trimmed_df=flip_block(trimmed_df)
    return trimmed_df


def make_map(trimmed_df, feature_name = None, height = 150, yrange = [-1,1], x_spread = 10000,label_feature=False ):
    
    """
    x_spread (int): number of nts on each side of feature ID. 
    """
    traces = []
    xrange = [0,0]
    #make traces for each feature

    comon_names = []
    for i,protein in enumerate(trimmed_df.index):
        #hoverinfo variables
        gene =  trimmed_df.loc[protein,'gene_ID']
        prod = trimmed_df.loc[protein,'product']
        start = trimmed_df.loc[protein, 'start']
        end = trimmed_df.loc[protein, 'end']
        strand = trimmed_df.loc[protein,'strand']
        comon_name = trimmed_df.loc[protein,'common_name']
        comon_names.append(comon_name)
        #make backbone trace between features
        if i< (len(trimmed_df.index)-1):
            next_orf_start = trimmed_df.iloc[(i+1),0]
            traces.append(go.Scatter(x = [end, next_orf_start],
                          y = [0,0],
                          mode = 'lines',
                          line = {'width': 3,
                                  'color':'black'},
                          showlegend = False,
                          hoverinfo = None))
        if strand == '-':
            start =  end
            end = trimmed_df.loc[protein, 'start']   
        l = (end-start)     #lenght of the arrow 
        br = start+(0.65*l) #defines where head of arrow is placed
         #make feature traces, highlighting the feature of interest
    
        if protein == feature_name:
            arrow_trace = make_arrow_trace(protein, gene, strand, prod, start, end, br, color='red')
            traces.append(arrow_trace)
            if label_feature:
                traces.append(make_annot_trace(x=start,y=0.1,text=comon_name))
            xrange  =[start-x_spread, end+x_spread]
        else:
            arrow_trace = make_arrow_trace(protein, gene, strand, prod, start, end, br,color='grey')
            traces.append(arrow_trace)
            if label_feature:
                traces.append(make_annot_trace(x=start,y=0.1,text=comon_name))
    dl = {'data':traces,
            'layout':go.Layout(#title = f'local genome map around {feature_name}',
                               paper_bgcolor='rgb(255,255,255)',
                               plot_bgcolor='rgb(255,255,255)',
                               width = 700,
                               height = height,
                               margin={'t':0, 'b':0, 'l':0, 'r':0},
                               showlegend = False,
            )} 
 
    fig = go.Figure(dl)  
    fig.update_yaxes(showgrid=False, showticklabels=False)
    fig.update_yaxes(range=yrange)
    fig.update_xaxes(range=xrange)
    return fig, comon_names


def make_annot_trace(x,y,text, fontsize=4):  
    return go.Scatter(
        x=[x],
        y=[y],
        text=[text],
        mode="text",
        textposition="top right",
        showlegend=False,
        textfont=dict(
            size=fontsize,
            color="black",
            family="Arial"
        ),
    )

def get_genes_wo_neighbors(path_to_gff, clearance=-5000, return_df=False):
    ''''
    get genes that have > clearance separation from neighbors on same strand
    
    parameters:
        path_to_gff (str): path to the gff3 genome annotation file
        clearance (int): number of nt that must exist beteween gene and 
            neighbor before start and after stop codon (on the same dna strand)
        return_df (bool): if False, returns list of genes that meet clearance
            criteria. if True, data-rich df is returned
        
        returns (list or df): see return_df parameter   
    '''
    
    gff_dict = pgf.make_seq_object_dict(path_to_gff,feature_type='CDS')
    strands = ['+', '-']
    df_lst = []
    for strand in strands:
        starts = []
        ends = []
        genes = []
        phase = []
        gene_df = pd.DataFrame()
        for gene in gff_dict.keys():
            gene = gff_dict[gene]
            if gene.strand == strand:
                starts.append(gene.start)
                ends.append(gene.end)
                genes.append(gene.ID)
                phase.append(gene.strand)
        gene_df['feature_ID'] = genes    
        gene_df['start'] = starts 
        gene_df['end'] = ends 
        gene_df['phase'] = phase 
        gene_df = gene_df.sort_values('start', ascending=True).reset_index(drop=True)
        gene_df['start_clearance'] = gene_df['start']         - gene_df['end'].shift()
        gene_df['end_clearance'] = gene_df['start'].shift(-1) - gene_df['end']
        gene_df = gene_df.fillna(np.inf)
        filt = (gene_df.loc[:,'start_clearance']>clearance) & (gene_df.loc[:,'end_clearance'] >clearance)
        gene_df = gene_df.loc[filt,:]
        df_lst.append(gene_df)
    if return_df:
        return  pd.concat(df_lst).set_index('feature_ID')
    else:
        return list(df_lst[0]['feature_ID']) + list(df_lst[1]['feature_ID'])



def get_genes_wo_neighbors_df(path_to_gff, clearance=-5000):
    ''''
    get genes that have > clearance separation from neighbors on same strand
    
    parameters:
        path_to_gff (str): path to the gff3 genome annotation file
        clearance (int): number of nt that must exist beteween gene and 
            neighbor before start and after stop codon (on the same dna strand)
            
        returns (df): Data-rich df is returned
    '''
    return get_genes_wo_neighbors(path_to_gff, clearance=clearance, return_df=True)
    


def rev_comp(seq):
    ''' 
    return the reverse complement of a dna sequences

    This is not optimized for large sequences, but code is simple
    
    Parameters:
    
    seq (str): can be a string of upper or lower case DNA bases: (a,t,c,g)
    
    
    '''
    base_pair_dict = {
        'A':'T',
        'T':'A',
        'C':'G',
        'G':'C',
        'a':'t',
        't':'a',
        'c':'g',
        'g':'c',
        'n':'n',
        'N':'N'
        }
    rev_comp_seq = ''
    for base in seq:
        rev_comp_seq = base_pair_dict[base] + rev_comp_seq
    return rev_comp_seq





def _get_UTR(seq_d, UTR_coords, contig, strand, b4_start,after_start):
    """return dna sequence in fasta genome around start codon. 
    
    seq_d (dict): dict with fasta genome where the ID is the contig name
    UTR_coords (list): [start, end], where start, end are top strand gene coords, not start codon
    contig (str): contig name, corresponding to a fasta id line
    strand (str): + or -, strand of DNA which encodes gene
    b4_start (int): how much sequence upstream of start codon to return
    after_start (int): how much sequence downstream of start codon to return
    """
    if UTR_coords[0] <0:
        UTR_coords[0] = 0
    if UTR_coords[1] > len(seq_d[contig]) - b4_start:
        UTR_coords[1] = len(seq_d[contig])
    seq = seq_d[contig][UTR_coords[0]:UTR_coords[1]]
    if strand == '-':
        seq = du.rev_comp(seq)
    return seq



def get_sequence_around_start(start, end, strand, contig, path_to_fasta_genome, b4_start, after_start):
    """return dna sequence in fasta genome around start codon. 
    
    seq_d (dict): dict with fasta genome where the ID is the contig name
    UTR_coords (list): [start, end], where start, end are top strand gene coords, not start codon
    contig (str): contig name, corresponding to a fasta id line
    strand (str): + or -, strand of DNA which encodes gene
    b4_start (int): how much sequence upstream of start codon to return
    after_start (int): how much sequence downstream of start codon to return
    """
    seq_d = pfa.get_seq_dict(path_to_fasta_genome)
    if strand == '+':
        UTR_coords = [start-(b4_start+1),start+(after_start-1)]
        seq = _get_UTR(seq_d, UTR_coords, contig, strand,b4_start,after_start)
    if strand == '-':
        UTR_coords = [end-after_start, end+b4_start]
        seq = _get_UTR(seq_d, UTR_coords, contig, strand, b4_start,after_start)
    return seq

    
    


        
    
            
            
            
    
