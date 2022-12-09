#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 22 09:02:17 2022

@author: jonwinkelman
"""
from Bio import SeqIO
import pandas as pd


path_to_genbank = './data/references/Reference_FERM_BP3421.gbk'
bam_filepath = './RNAseq_bam_files/BAN-1.bam'
print('hello')


def build_genbank_dict(path_to_genbank):
    """
    Builds dict {chromosomeID:{featureID:gbk_feature_object}} from a genbank annotation file.
    
    """

    chromosome_dict = {}
    gb_obj = SeqIO.parse(path_to_genbank, "genbank")
    for chromosome in gb_obj:  
        for feature in chromosome.features:
            if feature.type == 'CDS':
                if chromosome_dict.get(chromosome.name): 
                    chromosome_dict[chromosome.name].update({feature.qualifiers['locus_tag'][0]:feature})
                else:
                    chromosome_dict[chromosome.name] = {feature.qualifiers['locus_tag'][0]:feature}
    return chromosome_dict


def make_simple_annot_df(path_to_genbank, contig_name, start_end=False):
    """
    Return df derived from genbank file.

    Arguments:
    path_to_genbank (str): path to genbank File
    contig_name (str): name of the contig/plasmid/chromosome to be parsed
    start_end (bool): If True, adds gene start, end, and strand columns to df
    """
    gene_obj_dict = build_genbank_dict(path_to_genbank)[contig_name]
    columns = ['gene_ID', 'protein_ID', 'common_name', 'product']
    df = pd.DataFrame()
    df[columns[0]] = [gene_obj_dict[gene_ID].qualifiers['locus_tag'][0]
                      for gene_ID in gene_obj_dict.keys()]
    df[columns[1]] = [gene_obj_dict[gene_ID].qualifiers['protein_id'][0] if gene_obj_dict[gene_ID].qualifiers.get(
        'protein_id') else gene_ID for gene_ID in gene_obj_dict.keys()]
    df[columns[2]] = [gene_obj_dict[gene_ID].qualifiers.get(
        'name') for gene_ID in gene_obj_dict.keys()]
    df[columns[3]] = [gene_obj_dict[gene_ID].qualifiers['product'][0]
                      for gene_ID in gene_obj_dict.keys()]
    if start_end:
        columns = ['gene_ID', 'protein_ID', 'common_name',
                   'product', 'start', 'end', 'strand']
        df[columns[4]] = [
            gene_obj_dict[gene_ID].location.start for gene_ID in gene_obj_dict.keys()]
        df[columns[5]] = [
            gene_obj_dict[gene_ID].location.end for gene_ID in gene_obj_dict.keys()]
        df[columns[6]] = ['+' if gene_obj_dict[gene_ID].location.strand ==
                          1 else '-' for gene_ID in gene_obj_dict.keys()]
    df = df.set_index('gene_ID')
    return df
