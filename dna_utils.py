#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 26 09:12:04 2022

@author: jonwinkelman
"""
from collections import defaultdict, Counter
import itertools


def rna_to_dna(rna_sequence):
    # Replace every 'U' in the RNA sequence with 'T'
    rna_sequence = rna_sequence.upper().strip()
    dna_sequence = rna_sequence.replace('U', 'T')
    return dna_sequence


def dna_to_rna(dna_sequence):
    # Replace every 'U' in the RNA sequence with 'T'
    dna_sequence = dna_sequence.upper().strip()
    rna_sequence = dna_sequence.replace('T', 'U')
    return rna_sequence


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
        'Y':'R',
        'R':'Y',
        'S':'S',
        'W':'W',
        'M':'K',
        'K':'M',
        'V':'B',
        'B':'V',
        'D':'H',
        'H':'D',
        'N':'N'
        }
    rev_comp_seq = ''
    for base in seq:
        base = base.upper()
        rev_comp_seq = base_pair_dict[base] + rev_comp_seq
    return rev_comp_seq


def generate_kmers(k, alphabet = 'ACGU'):
    """Generates all possible k-mers of length k."""
    return [''.join(p) for p in itertools.product(alphabet, repeat=k)]

# def count_kmers(sequences, k):
#     """Counts occurrences of each k-mer in the sequences."""
#     kmer_counts = Counter()
#     for seq in sequences:
#         for i in range(len(seq) - k + 1):
#             kmer = seq[i:i+k]
#             kmer_counts[kmer] += 1
#     return kmer_counts



def count_kmers(sequences, k):
    """Counts occurrences of each k-mer in the sequences, including those with zero counts."""
    # Generate all possible k-mers
    all_kmers = generate_kmers(k)
    
    # Initialize kmer_counts with all possible k-mers set to 0
    kmer_counts = Counter({kmer: 0 for kmer in all_kmers})
    
    # Update counts based on the sequences
    for seq in sequences:
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            kmer_counts[kmer] += 1
    
    return kmer_counts

def calculate_kmer_enrichment(kmer_counts, total_kmers, k):
    """Calculates enrichment of each k-mer.
    *** Assumes equal distribution of A, C, G, T"""
    
    background_freq = 1 / (4 ** k)  # Assuming equal distribution of A, C, G, T
    enrichments = {}
    for kmer, count in kmer_counts.items():
        expected_count = background_freq * total_kmers
        enrichments[kmer] = count / expected_count if expected_count > 0 else float('inf')
    return enrichments

import math

import math

def calculate_kmer_enrichment_log(kmer_counts, total_kmers, k):
    """Calculates log2 enrichment of each k-mer.
    Assumes equal distribution of A, C, G, T."""
    
    background_freq = 1 / (4 ** k)  # Assuming equal distribution of A, C, G, T
    enrichments = {}
    for kmer, count in kmer_counts.items():
        expected_count = background_freq * total_kmers
        if expected_count > 0:
            ratio = count / expected_count
            if ratio > 0:
                enrichments[kmer] = math.log2(ratio)
            else:
                enrichments[kmer] = -float('inf')  # or another small value, depending on your needs
        else:
            enrichments[kmer] = float('inf') if count > 0 else float('-inf')
    return enrichments



def find_enriched_kmers(sequences, k_values=[3, 4, 5], alphabet = 'ACGU', log=False):
    """Finds and returns k-mers and their enrichment/depetion for specified k values.
    *** Assumes equal distribution of A, C, G, T
    """
    all_enriched_kmers = {}
    
    for k in k_values:
        kmers = generate_kmers(k, alphabet)
        kmer_counts = count_kmers(sequences, k)
        total_kmers = sum(kmer_counts.values())
        if log:
            enrichments = calculate_kmer_enrichment_log(kmer_counts, total_kmers, k)
        else:
            enrichments = calculate_kmer_enrichment(kmer_counts, total_kmers, k)
        
        all_enriched_kmers[k] = {kmer: enrichment for kmer, enrichment in enrichments.items()}
    
    return all_enriched_kmers


def find_alphabet(sequences):
    "Returns the set of characters from a list of sequences"
    # Use a set to store unique characters
    alphabet = set()
    
    # Iterate over each sequence in the list
    for seq in sequences:
        # Add each character in the sequence to the set
        alphabet.update(seq)
    alphabet = ''.join(alphabet)
    return alphabet
    