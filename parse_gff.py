#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 11 13:34:22 2022

@author: jonwinkelman
"""
import pandas as pd

class _seq_attributes:
    'make a sequence object from a partially processed line in a GFF file'
    def __init__(self, line, anot_lst, attributes_dict, *args, **kwargs):
        """ make a sequence object for a each feature of a given type
        
        parameters
        line (str): the entire feature line in gff file
        annot_lst (list): beginning of feature line, tab-delimited
        attribute_dict (dict): second part of feature line. Field are variable depending
                            how it was annotated

        """
        self.line = line
        self.ID = attributes_dict['ID']
        self.Parent = attributes_dict.get('Parent')
        self.Dbxref = attributes_dict.get('Dbxref')
        self.Name = attributes_dict.get('Name')
        self.gbkey = attributes_dict.get('gbkey')
        self.inference = attributes_dict.get('inference')
        self.locus_tag = attributes_dict.get('locus_tag')
        self.product = attributes_dict.get('product')
        self.protein_id = attributes_dict.get('protein_id')
        self.transl_table = attributes_dict.get('transl_table')
        self.gene_biotype = attributes_dict.get('gene_biotype')
        self.Note = attributes_dict.get('Note')
        self.gene = attributes_dict.get('gene')


        self.chromosome = anot_lst[0]
        self.source = anot_lst[1]
        self.feature_type = anot_lst[2]
        self.start = int(anot_lst[3])
        self.end = int(anot_lst[4])
        self.score = anot_lst[5]
        self.strand = anot_lst[6]   
        self.phase = anot_lst[7]
        
        

def make_seq_object_dict(path_to_gff3, feature_type = 'gene'):               
    """
    Parse GFF and make dict containing feature objects.
    
    feature_types = 'gene', 'CDS', 'tRNA', 'exon', 'pseudogene', 'SRP_RNA', 'sequence_feature',...
    
    Parameters:
    path_to_gff3 (str): path to gff3 file
    feature_type (str): 
    
    Returns:
    dict: dictionary key is the feature name, values are feature objects that contain
    various attributes such as start and end position in the genome, product description,
    what strand of DNA feature is encoded on, parent gene name, 
    """
    with open(path_to_gff3, 'r') as f:
        seq_dict = {}     
        for i,line in enumerate(f):
            #if i>6: break
            if line[0] != '#' and line[0] != ' ':
                anot_lst = line.split('\t')
                if len(anot_lst) <6:
                    raise Exception(f'line {i+1} in {path_to_gff3} does not contain all annotation fields')
                #get only protein homology lines
                if anot_lst[2] == feature_type: #make this a passed variable
                    if anot_lst[-1].find('ID=') !=-1:
                        temp = [key_val.split('=') for key_val in anot_lst[-1].split(';')]
                        attributes_dict = {key_val[0]:key_val[1] for key_val in temp}
                        ID = attributes_dict['ID']
                        seq_dict[ID] = _seq_attributes(line, anot_lst, attributes_dict)
    return seq_dict   
    
    
def make_gene2prot_dict(path_to_gff3):               
    """
    Return a dictionary {geneID:proteinID}
    
    
    Parameters:
    path_to_gff3 (str): path to gff3 file
    
    Returns:
    dict
    """
    seq_dict = make_seq_object_dict(path_to_gff3, feature_type = 'CDS')
    return {seq_dict[protein].Parent:protein for protein in seq_dict.keys()}
   

def make_prot2gene_dict(path_to_gff3):
    """
    Return a dictionary {proteinID:geneID}
    
    
    Parameters:
    path_to_gff3 (str): path to gff3 file
    
    Returns:
    dict
    """
    seq_dict = make_seq_object_dict(path_to_gff3, feature_type = 'CDS')
    return {protein:seq_dict[protein].Parent for protein in seq_dict.keys()}


def write_simple_annot_file(path_to_gff, filepath):
    '''
    write tsv file with protein accession, parent gene, common name and function in each row
    
    parameters:
    
        path_to_gff (str): path to gff3 annotation file.
        
        filename/path (str): name you want to give file with ".txt" extension added. If
        path is given, then file will be written to that path, else it will be written in 
        current working directory.
    '''
    jb_bug = make_seq_object_dict(path_to_gff, feature_type='CDS')
    jb_bug_gene = make_seq_object_dict(path_to_gff, feature_type='gene')
    with open(f'{filepath}.txt', 'w') as f:
        for prot in jb_bug.keys():
            gene = jb_bug[prot].Parent.strip()
            if jb_bug_gene.get(gene):
                common = jb_bug_gene[gene].Name
            else: common = None
            product = jb_bug[prot].product
            f.write(f'{prot}\t{gene}\t{common}\t{product}\n')
            
            
# def get_gff_summary(path_to_gff):
#     """return summary of genome assembly derived from gff
    
#     parameters
#     path_to_gff (str): 
#     """ 
#     summmary = {}
#     with open(path_to_gff, 'w') as f:
#         for line in f:
#             line.split()
#             if line
                     

def make_simple_annot_df(path_to_gff, start_end = False):
    """generate simple dataframe from gff

    index: gene_ID
    columns: 'protein_ID', 'common_name', 'product'

    parameters:
    path_to_gff (str): path to the gff file
    start_end (bool): include start and end of gene as individual columns
    
    """
    jb_bug = make_seq_object_dict(path_to_gff, feature_type='CDS')
    jb_bug_gene = make_seq_object_dict(path_to_gff, feature_type='gene')
    genes = []
    proteins = []
    common_name = []
    product = []
    starts = []
    end = []
    strand = []
    for prot in jb_bug.keys():
        proteins.append(prot)
        gene = jb_bug[prot].Parent.strip()
        genes.append(gene)
        if start_end:
            starts.append(jb_bug[prot].start)
            end.append(jb_bug[prot].end)
            strand.append(jb_bug[prot].strand)
        if jb_bug_gene.get(gene):
            common_name.append(jb_bug_gene[gene].Name)
        else: 
            common_name.append('')
        product.append(jb_bug[prot].product)
    df = pd.DataFrame()
    df['gene_ID'] = genes
    df['protein_ID'] = proteins
    df['common_name'] = common_name
    df['product'] = product
    if start_end:
        df['start'] = starts
        df['end'] = end
        df['strand'] = strand
        
    df = df.set_index('gene_ID')
    return df


if __name__ == '__main__':
    gff_obj = make_seq_object_dict('/Users/jonwinkelman/Dropbox/Trestle_projects/Mukherjee_lab/ncbi_dataset/ncbi_dataset/data/GCA_000014625.1/genomic.gff')
    print(gff_obj.keys())
    a = gff_obj['gene-PA14_00060']
    print(a.Note)