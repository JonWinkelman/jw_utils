#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 14:08:25 2022

@author: jonwinkelman

"""

def get_seq_dict(path_to_fasta):
    '''
    parse fast file and return as a dictionary
    
    parameters:
        path_to_fasta (str): path to a fasta proteome
        
        return (dict): eqID as key and sequence as value
    '''
    seq = ''
    seq_dict = {}
    prot_id = None
    with open(path_to_fasta, 'r') as f:
        for line in f:
            if line[0] == '>':
                if prot_id:
                    seq_dict[prot_id] = seq
                prot_id = line.split(' ')[0][1:].strip()
                seq = ''
            else:
                seq = seq + line.strip()
        seq_dict[prot_id] = seq
            
    return seq_dict



def get_seq_dict_fulldescript(path_to_fasta):
    '''
    parse fast file and return as a dictionary
    
    parameters:
        path_to_fasta (str): path to a fasta proteome
        
        return (dict): eqID as key and sequence as value
    '''
    seq = ''
    seq_dict = {}
    prot_id = None
    with open(path_to_fasta, 'r') as f:
        for line in f:
            if line[0] == '>':
                if prot_id:
                    seq_dict[prot_id] = seq
                line_list = line.split(' ')  
                d1 = line_list[0][1:].strip()
                d2 = line_list[1].strip()
                prot_id = f'{d1}_{d2}'
                seq = ''
            else:
                seq = seq + line.strip()
        seq_dict[prot_id] = seq
            
    return seq_dict



def get_protein_subset(path_to_fasta, seq_ids):
    '''
    return a dict with seqID as key and sequence as value
    
    parameters:
        path_to_fasta (str): path to a fasta proteome
    '''
    
    seq_dict_subset = {}
    seq_dict = get_seq_dict(path_to_fasta)
    for i in seq_ids:
        seq = seq_dict.get(i)
        if seq:
            seq_dict_subset[i] = seq
    return seq_dict_subset


def write_to_fasta(seq_dict, path):
    with open(path, 'w') as f:
        for name, seq in seq_dict.items():
            f.write(f'>{name}\n{seq}\n')
            
            
            
            
def fasta_to_phyl(fasta_aln_path, output_path):
    """
    turn fasta format into simple phylip format
    Note: other info in fasta will be lost in this phylip format
    """

    fasta_dict = get_seq_dict(fasta_aln_path)
    num_algnments = len(list(fasta_dict.values()))
    aln_length = len(list(fasta_dict.values())[0])
    with open(output_path, 'w') as f:
        f.write(str(num_algnments)+' ' + str(aln_length) + '\n')
        for protID, seq in fasta_dict.items():
            f.write(f'{protID} {seq}\n')
            
            
        
        
        
            