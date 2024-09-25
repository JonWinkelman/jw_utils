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

    
    


        
    
            
            
            
    
