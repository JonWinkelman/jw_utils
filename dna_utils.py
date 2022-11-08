#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 26 09:12:04 2022

@author: jonwinkelman
"""
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
        }
    rev_comp_seq = ''
    for base in seq:
        rev_comp_seq = base_pair_dict[base] + rev_comp_seq
    return rev_comp_seq
    