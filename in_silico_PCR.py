from jw_utils import parse_fasta as pfa
import re
import os
import pandas as pd
from Bio.Seq import Seq
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import pairwise2


# IUPAC degenerate nucleotide codes and their regex equivalents
IUPAC_CODES = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]"
}

IUPAC_CODES = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT",
    "K": "GT", "M": "AC", "B": "CGT", "D": "AGT",
    "H": "ACT", "V": "ACG", "N": "ACGT"
}

def expand_degenerate_sequence(seq):
    """Generate all possible sequences from a degenerate primer."""
    bases = [IUPAC_CODES[base] for base in seq.upper()]
    return ["".join(p) for p in itertools.product(*bases)]
    
def degenerate_to_regex(primer):
    """Convert a degenerate primer sequence into a regex pattern."""
    return "".join(IUPAC_CODES[base] for base in primer.upper())

def in_silico_pcr(forward_primer, reverse_primer, sequence, min_amplicon=50, max_amplicon=2000):
    """
    Perform in silico PCR using degenerate primers on a given DNA sequence.

    Parameters:
    - forward_primer (str): Forward degenerate primer sequence.
    - reverse_primer (str): Reverse degenerate primer sequence.
    - sequence (str): DNA sequence to search.
    - min_amplicon (int): Minimum amplicon length.
    - max_amplicon (int): Maximum amplicon length.

    Returns:
    - List of amplified sequences (if found), else an empty list.
    """
    # Convert primers to regex patterns
    fwd_regex = degenerate_to_regex(forward_primer)
    rev_regex = degenerate_to_regex(str(Seq(reverse_primer).reverse_complement()))
    
    # Search for forward primer
    fwd_matches = [m.start() for m in re.finditer(fwd_regex, sequence, re.IGNORECASE)]
    
    # Search for reverse primer (on reverse complement strand)
    rev_matches = [m.start() for m in re.finditer(rev_regex, sequence, re.IGNORECASE)]
    
    amplified_products = []
    
    # Find valid amplicons
    for fwd in fwd_matches:
        for rev in rev_matches:
            if min_amplicon <= ((rev+len(reverse_primer)) - fwd) <= max_amplicon:
                amplified_products.append(sequence[fwd:rev + len(reverse_primer)])

    return amplified_products





def find_local_alignment(primer, sequence):
    """
    Find local alignment of primer within a sequence using Smith-Waterman.

    Parameters:
    - primer (str): Primer sequence.
    - sequence (str): DNA sequence.

    Returns:
    - Tuple: (alignment_start, aligned_sequence, alignment_score) or None if no match found.
    """
    alignments = pairwise2.align.localms(sequence, primer, match=2, mismatch=-1, open=-2, extend=-1)

    if alignments:
        best = alignments[0]  # Best alignment
        alignment_start = best[3]  # Start index of alignment in the sequence
        aligned_sequence = best[0][alignment_start:alignment_start + len(primer)]
        alignment_score = best[2]  # Alignment score

        return alignment_start, aligned_sequence, alignment_score
    return None



