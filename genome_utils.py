#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 14:39:32 2022

@author: jonwinkelman
"""

from jw_utils import parse_gff as pgf
import pandas as pd
from jw_utils import parse_fasta as pf
import numpy as np


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



def get_feature_seq(path_to_fasta_genome, path_to_annot_file,feature_ID):   
       
    
    fasta_dict = pf.get_seq_dict(path_to_fasta_genome) 


# def get_feature_dna_seq(path_to_fasta_genome,  path_to_annot_file, feature_ID, feature_type='gene'):
#     """returns the coding strand of the feature"""
#     gff_geneObject_dict = pgf.make
#     polarity = self.gff_geneObject_dict[gene].strand
#     gene_start = self.gff_geneObject_dict[gene].start
#     gene_end = self.gff_geneObject_dict[gene].end
#     gene_seq = self.genome_seq_dict[self.chromosome_names[0]][(gene_start-1):gene_end]
#     if polarity == '-':
#         start_codon = gene_end
#         stop_codon = gene_start
#         gene_seq = rev_comp(gene_seq)
#     return gene_seq

    
    
class GenomeUtils:
    def __init__(self, path_to_fasta_genome, path_to_gff):
         
        self.path_to_fasta_genome = path_to_fasta_genome
        self.genome_seq_dict = pf.get_seq_dict(self.path_to_fasta_genome) 
        self.path_to_gff = path_to_gff
        self.fasta_dict = pf.get_seq_dict(self.path_to_fasta_genome)
        self.gff_cdsObject_dict = pgf.make_seq_object_dict(
                                            self.path_to_gff, 
                                            feature_type='CDS'
                                            ) 
        self.gff_geneObject_dict = pgf.make_seq_object_dict(
                                            self.path_to_gff, 
                                            feature_type='gene'
                                            ) 
        self.chromosome_names = list(self.genome_seq_dict.keys())

        
        
    def get_gene_dna_seq(self, gene):
        """returns the coding strand of the feature"""
        polarity = self.gff_geneObject_dict[gene].strand
        gene_start = self.gff_geneObject_dict[gene].start
        gene_end = self.gff_geneObject_dict[gene].end
        gene_seq = self.genome_seq_dict[self.chromosome_names[0]][(gene_start-1):gene_end]
        if polarity == '-':
            start_codon = gene_end
            stop_codon = gene_start
            gene_seq = rev_comp(gene_seq)
        return gene_seq

    def get_feature_dna_seq(self, feature_id, feature_type, contig):
        """Return the sequence of any given feature type within the annotation file"""
        
        feat_obj_dict = pgf.make_seq_object_dict(
                                            self.path_to_gff, 
                                            feature_type=feature_type
                                            )   
        polarity = feat_obj_dict[feature_id].strand
        feature_start = feat_obj_dict[feature_id].start
        feature_end = feat_obj_dict[feature_id].end
        feature_seq = self.genome_seq_dict[contig][(feature_start-1):feature_end]
        if polarity == '-':
            start_codon = feature_end
            stop_codon = feature_start
            feature_seq = rev_comp(feature_seq)
        return feature_seq


    

    def get_upstream_sequence(self, length_before_start,length_after_start=0,
                            gene_list=None, feature_type='gene', prefix_format=True):
        """Return the desired length of DNA sequence before the start codon.

        **Sequence returned is on the same strand and orientaion as the mRNA

        Parameters:
        length_before_start (int): length of sequence to return before start codon.
        length_after_start (int): length of sequence to return after start codon.
        gene_list=None (list): list of genes. Preferably with the prefix 'gene-' or 'cds-' on front, else
                                funtion will add this prefix. 
        feature_type (str): can be 'gene', 'CDS',.
        prefix_format (bool): if True, return dictionary keys with 'gene-' or 'cds-' as prefix.

        Return (dict):
        dictionary with gene protein ID ('gene-...' or 'cds-...) as the key and the upstream
        sequence as the value
        
        """
        
        avail_feature_types = ['CDS', 'gene']
        if feature_type not in avail_feature_types:
            raise TypeError(f'feature_type {feature_type} not in available feature types.\n Try one of the following: {avail_feature_types}') 
 
        contigs = list(self.genome_seq_dict.keys())
        genome_seq = self.genome_seq_dict[contigs[0]]
        if len(contigs) > 1: 
            print(f'There are multiple contigs in the fasta: {self.path_to_fasta_genome}, \n attempting to use {contigs[0]}' )
        if feature_type=='CDS': 
            prefix = 'cds-'
            gff_obj_dict = self.gff_cdsObject_dict
        else:
            prefix = 'gene-'
            gff_obj_dict = self.gff_geneObject_dict
        if type(gene_list) != list:
            gene_list=list(gff_obj_dict.keys())
        else:
            #check for "gene-" or "CDS-" on front of feature ID
            for i, gene in enumerate(gene_list):
                gene = str(gene)
                if not gene.startswith(prefix):
                    gene = prefix + gene
                    gene_list[i] = gene

        strand= [gff_obj_dict[gene].strand if gff_obj_dict.get(gene) else None for gene in gene_list]
        upstream_regions = {}
        i=0
        for gene in gene_list:
            if gff_obj_dict.get(gene):
                if gff_obj_dict[gene].strand == '-':
                    region_end  = gff_obj_dict[gene].end + length_before_start
                    region_beginning = gff_obj_dict[gene].end - length_after_start
                    region = genome_seq[region_beginning:region_end]
                    if prefix_format:
                        upstream_regions[gene] = rev_comp(region)
                    else:
                        upstream_regions[gene.replace(prefix, '')] = rev_comp(region)
                elif gff_obj_dict[gene].strand == '+':
                    region_beginning = (gff_obj_dict[gene].start-length_before_start) - 1
                    region_end = (gff_obj_dict[gene].start + length_after_start) - 1
                    region = genome_seq[region_beginning:region_end]
                    if prefix_format:
                        upstream_regions[gene] = region
                    else:
                        upstream_regions[gene.replace(prefix, '')] = region
            else:
                upstream_regions[gene] = None

        return upstream_regions

        
        

    def get_downstream_sequence(self, length_before_stop,length_after_stop=0,
                            gene_list=None, feature_type='gene', prefix_format=True):
        """Return the desired length of DNA sequence after the stop codon.

        **Sequence returned is on the same strand and orientaion as the mRNA

        Parameters:
        length_before_stop (int): length of sequence to return before the end of the stop codon,
                                  e.g. if 0, then sequence starts on first nt after stop codon.
        length_after_stop (int): length of sequence to return downstream of the coding sequence.
        gene_list=None (list): list of genes. Preferably with the prefix 'gene-' or 'cds-' on front, else
                                funtion will add this prefix. 
        feature_type (str): can be 'gene', 'CDS',.
        prefix_format (bool): if True, return dictionary keys with 'gene-' or 'cds-' as prefix.

        Return (dict):
        dictionary with gene protein ID ('gene-...' or 'cds-...) as the key and the upstream
        sequence as the value
        
        """
        avail_feature_types = ['CDS', 'gene']
        if feature_type not in avail_feature_types:
            raise TypeError(f'feature_type {feature_type} not in available feature types.\n Try one of the following: {avail_feature_types}') 
            
        contigs = list(self.genome_seq_dict.keys())
        genome_seq = self.genome_seq_dict[contigs[0]]
        if len(contigs) > 1: 
            print(f'There are multiple contigs in the fasta: {self.path_to_fasta_genome}, \n attempting to use {contigs[0]}' )
        if feature_type=='CDS': 
            prefix = 'cds-'
            gff_obj_dict = self.gff_cdsObject_dict
        else:
            prefix = 'gene-'
            gff_obj_dict = self.gff_geneObject_dict
        if type(gene_list) != list:
            gene_list=list(gff_obj_dict.keys())
        else:
            #check for "gene-" or "CDS-" on front of feature ID
            for i, gene in enumerate(gene_list):
                gene = str(gene)
                if not gene.startswith(prefix):
                    gene = prefix + gene
                    gene_list[i] = gene

        # strand= [gff_obj_dict[gene].strand if gff_obj_dict.get(gene) else None for gene in gene_list]
        upstream_regions = {}
        i=0
        for gene in gene_list:
            if gff_obj_dict.get(gene):
                if gff_obj_dict[gene].strand == '-':
                    region_end  = (gff_obj_dict[gene].start + length_before_stop)-1
                    region_beginning = (gff_obj_dict[gene].start - length_after_stop)-1
                    region = genome_seq[region_beginning:region_end]
                    if prefix_format:
                        upstream_regions[gene] = rev_comp(region)
                    else:
                        upstream_regions[gene.replace(prefix, '')] = rev_comp(region)
                elif gff_obj_dict[gene].strand == '+':
                    region_beginning = gff_obj_dict[gene].end-length_before_stop
                    region_end = gff_obj_dict[gene].end + length_after_stop
                    region = genome_seq[region_beginning:region_end]
                    if prefix_format:
                        upstream_regions[gene] = region
                    else:
                        upstream_regions[gene.replace(prefix, '')] = region
            else:
                upstream_regions[gene] = None

        return upstream_regions    

    
    
    
    
class ProteomeUtils:
    def __init__(self, path_to_proteome, path_to_gff):  
        self.path_to_proteome = path_to_proteome
        self.path_to_gff = path_to_gff
        self.protein_seq_dict =  pf.get_seq_dict(self.path_to_proteome)
     
    
    def get_AA_seq_for_protein(self, protein):
        return self.protein_seq_dict[protein]
        
    def get_get_AA_seq_for_proteins(self, protein_list):
        AA_seqs = [self.protein_seq_dict.get(protein.replace('cds-','')) for protein in protein_list]
        return pd.DataFrame({'protein ID':protein_list, 'AA_seq':AA_seqs}).set_index('protein ID')
            
            
    def get_DNA_seq_for_proteins(self, path_fasta_genome, protein_list):
        '''
        returns the DNA sequences for the input proteins
        parameters]
    
        '''
        genome_obj = GenomeUtils(path_fasta_genome, self.path_to_gff)
        prot2gene_dict = pgf.make_prot2gene_dict(self.path_to_gff)
        seqs = []
        genes = []
        proteins = []
        for protein in protein_list:
            gene = prot2gene_dict.get(protein)
            if genome_obj.gff_geneObject_dict.get(gene):
                seq = genome_obj.get_gene_dna_seq(gene)
            else:
                seq = ''
            seqs.append(seq)
            genes.append(gene)
            proteins.append(protein)
        df = pd.DataFrame({'protein ID':proteins,'gene ID':genes, 'sequence':seqs})
        return df.set_index('protein ID')

        
    
            
            
            
    
