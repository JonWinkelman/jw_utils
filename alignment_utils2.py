import pandas as pd
from collections import Counter
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceMatrix
from Bio import AlignIO
import numpy as np
import math
import subprocess
import tempfile
import os
from jw_utils import parse_fasta as pfa


def map_unaligned_to_aligned(unaligned_seq, aligned_seq):
    """
    Maps each unaligned (raw) position to its aligned (MSA) position.
    
    Parameters:
    - unaligned_seq: str, raw protein sequence
    - aligned_seq: str, aligned version (same sequence with gaps)
    
    **Assumptions:
    1. aligned_seq is derived from unaligned_seq and preserves residue order
    2. gaps (-) are in aligned_seq only
    3. all residues in unaligned_seq appear in aligned_seq (no mismatches)
    
    Returns:
    - Dictionary: {unaligned_index: aligned_index}
    """
    aligned_seq = aligned_seq.replace('*','')
    unaligned_seq = unaligned_seq.replace('*','')
    if aligned_seq.replace('-', '').upper() != unaligned_seq.upper():
        raise ValueError('the aligned sequence and unaligned seq are not the same when "-" are removed')
    mapping = {}
    unaligned_pos = 0

    for aligned_pos, aa in enumerate(aligned_seq):
        if aa != '-':
            if unaligned_pos >= len(unaligned_seq):
                raise ValueError("More residues in aligned_seq than in unaligned_seq.")
            mapping[unaligned_pos] = aligned_pos
            unaligned_pos += 1

    if unaligned_pos != len(unaligned_seq):
        raise ValueError("unaligned_seq not fully accounted for in aligned_seq.")

    return mapping



def run_muscle(inp, output_fp, mode='align', muscle_path='muscle', **kwargs):
    """
    Run MUSCLE v5 alignment using either -align or -super5.

    Parameters:
    - inp (str): Path to input FASTA file or a dict {seqname1:seq1, seqname2:seq2,...}
    - output_fp (str): Path to output aligned FASTA file
    - mode (str): 'align' or 'super5'
    - muscle_path (str): Path to the MUSCLE binary (default: 'muscle')
    - **kwargs: Additional MUSCLE flags (e.g., -perm, -perturb)

    usage example 1: run_muscle("input.fa", "output.aln", mode='super5')
    usage example 2: run_muscle("input.fa", "output.aln", mode='align', **{'-perm': 'abc', '-perturb': '1'})
    usage example 3: run_muscle(None, None, **{'-h': None})
        If '-h' is passed as a keyword, help is printed and the function exits.
    

    """
    if type(inp) != dict:
        inp = str(inp)
        output_fp = str(output_fp)
        tmp_path = None
    else:
        with tempfile.NamedTemporaryFile(delete=False) as tmp:
            tmp_path = tmp.name   
        pfa.write_to_fasta(inp, tmp_path)
        inp = tmp_path
    try:  
        if '-h' in kwargs:
            subprocess.run([muscle_path, '-h'])
            return
    
        if mode not in ['align', 'super5']:
            raise ValueError("Mode must be either 'align' or 'super5'")
    
        cmd = [muscle_path, f'-{mode}', inp, '-output', output_fp]
    
        for key, value in kwargs.items():
            if key != '-h':  # already handled
                cmd.append(str(key))
                if value is not None:
                    cmd.append(str(value))

        print("Running command:", " ".join(cmd))
        subprocess.run(cmd, check=True)
        aln_seq_d = pfa.get_seq_dict(output_fp)

    finally:
        # Clean up temp file
        if tmp_path and os.path.exists(tmp_path):
            os.remove(tmp_path)
            
    return aln_seq_d


def sort_by_column_character_frequency(seq_dict, column):
    """
    Sorts a dictionary of aligned sequences based on the most common character 
    in a specified column. Sequences containing the most common character in 
    that column appear at the top, and those with the least common character appear at the bottom.

    Parameters:
    seq_dict (dict): Dictionary with {seq_id: aligned_sequence} pairs.
    column (int): Column index to analyze.

    Returns:
    dict: Sorted dictionary with sequences ordered by column character frequency.
    """

    # Extract characters at the given column for all sequences
    column_chars = [seq[column] for seq in seq_dict.values() if len(seq) > column]

    # Count occurrences of each character
    char_counts = Counter(column_chars)

    # Sort sequences by frequency of their character in the given column
    sorted_items = sorted(seq_dict.items(), key=lambda x: char_counts.get(x[1][column], 0), reverse=True)

    # Return sorted dictionary
    return dict(sorted_items)


def calculate_entropy(alignment_input, file_format="fasta", ignore_gaps=True):
    """
    Calculate Shannon entropy at each position in a multiple sequence alignment.
    
    Parameters:
    ----------
    alignment_input : str or dict
        - If str: Path to the alignment file (FASTA, Clustal, etc.).
        - If dict: A dictionary where keys are sequence IDs and values are aligned sequences.
    file_format : str, optional
        Format of the alignment file (if a file is provided). Default is "fasta".
    ignore_gaps : bool, optional
        Whether to exclude gaps ('-') from entropy calculation. Default is True.

    Returns:
    -------
    position : Shannon entropy : dict with of float (value) for each alignment position (key) {postion:entropy}
    
    Formula:
    --------
    Shannon entropy is calculated as:

        H = -Σ (p_i * log2(p_i))

    where:
        - p_i is the frequency of each unique character at the given alignment position.
        - The sum is taken over all observed characters (nucleotides/amino acids).

    Example Usage:
    --------------
    entropy_values = calculate_entropy("alignment.fasta")
    entropy_values = calculate_entropy({"seq1": "ATG-", "seq2": "ATGA", "seq3": "ATG-"})
    """

    if isinstance(alignment_input, str):
        # Load the alignment from a file
        alignment = AlignIO.read(alignment_input, file_format)
        sequences = [str(record.seq) for record in alignment]
    elif isinstance(alignment_input, dict):
        # Use the dictionary directly
        sequences = list(alignment_input.values())
    else:
        raise ValueError("alignment_input must be a file path (str) or a dictionary of sequences.")

    num_seqs = len(sequences)
    alignment_length = len(sequences[0])  # Assume all sequences are aligned and of the same length

    entropies = {}

    # Iterate through each column (position) in the alignment
    for i in range(alignment_length):
        column = [seq[i] for seq in sequences]  # Get characters at position i
        unique_chars = set(column)

        # Optionally remove gaps ('-') from calculations
        if ignore_gaps:
            unique_chars.discard('-')

        # Compute frequency of each character
        freqs = {char: column.count(char) / num_seqs for char in unique_chars}

        # Compute Shannon entropy
        entropy = -sum(p * math.log2(p) for p in freqs.values() if p > 0)
        entropies[i] = entropy
    return entropies



def compute_distance_matrix(fasta_file, model="identity"):
    """
    Compute a pairwise distance matrix from a multiple sequence alignment in FASTA format.

    Parameters:
    ----------
    fasta_file : str
        Path to the alignment file (FASTA format).
    model : str, optional
        Distance model to use. Options for nucleotides: "identity", "jc69", "k2p".
        Options for proteins: "identity", "blosum62", "pam250".
        Default is "identity" (simple % identity).

    Returns:
    -------
    distance_matrix_df : pandas.DataFrame
        Distance matrix as a DataFrame.

    Example Usage:
    --------------
    df = compute_distance_matrix("alignment.fasta", model="identity")
    """

    # Load the alignment
    alignment = AlignIO.read(fasta_file, "fasta")

    # Create a DistanceCalculator with the chosen model
    calculator = DistanceCalculator(model)

    # Compute the distance matrix
    distance_matrix = calculator.get_distance(alignment)

    # Convert DistanceMatrix to a DataFrame for easier visualization
    taxa = distance_matrix.names
    matrix_data = [[distance_matrix[t1, t2] for t2 in taxa] for t1 in taxa]
    
    # Create a pandas DataFrame
    distance_matrix_df = pd.DataFrame(matrix_data, index=taxa, columns=taxa)

    return distance_matrix_df



def aln_d_to_dataframe(seq_d):
    """"""
    return pd.DataFrame.from_dict({name: list(val.upper()) for name, val in seq_d.items()}).transpose()

def get_character_frequency(s, chars):
    """"""
    counts = Counter(s)
    char_freq = []
    for char in chars:
        char_freq.append(counts.get(char, 0)/len(s))
    return char_freq



def get_frequency_matrix(seq_d,  molecule_type='prot'):
    """
    seq_d (dict):
    molecule_type (str): 'prot', 'dna', 'rna'

    return (pandas.DataFrame): df columns represent columns in alignment, 0-based, and
                               each row is a character, e.g. if molecule_type='prot' then
                               one of the 20 AAs
    """
    if molecule_type=='dna':
        chars = ['A', 'C', 'T', 'G', '-']
    if molecule_type=='rna':
        chars = ['A', 'C', 'U', 'G', '-']
    if molecule_type=='prot':
        chars = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '-']
        
    aln_df = aln_d_to_dataframe(seq_d)
    cols = []
    for i,col in enumerate(aln_df):
        cols.append(get_character_frequency(aln_df[col], chars))
    freq_df = pd.DataFrame(cols).transpose()
    freq_df.index=chars
    return freq_df




def remove_homogeneous_cols(seq_d, threshold):
    """return aln dict containing columns frequency of any character is not above given threshold
    
    *removes columns where a character frequency in that column is above threshold
    """
    freq_matrix = get_frequency_matrix(seq_d)
    aln_df  = aln_d_to_dataframe(seq_d)
    cols_to_keep=[]
    for col in freq_matrix:
        if sum(freq_matrix[col]>threshold) < 1:
            cols_to_keep.append(col)
    return aln_df[cols_to_keep].apply(lambda x: ''.join(x), axis=1).to_dict()


def remove_gappy_columns(seq_d, threshold=0.8):
    """return aln dict where columns containing gaps ('-') above at frequency abouve threshold have been removed.
    
    """
    freq_matrix = get_frequency_matrix(seq_d)
    aln_df  = aln_d_to_dataframe(seq_d)
    (aln_df == '-').sum(axis=0)/aln_df.shape[0]
    col_gap_freqs = (aln_df == '-').sum(axis=0)/aln_df.shape[0]
    cols_to_keep = list(col_gap_freqs[col_gap_freqs < threshold].index)
    return aln_df[cols_to_keep].apply(lambda x: ''.join(x), axis=1).to_dict()


def align_by_motif(loop_seqs_series, motif='GGA'):
    """aligns first instace of motif and pads left and right with '-'.
    
    *simply pads the sequence with '-' so that the beginning of the motif is aligned
    
    loop_seqs (pd.Series): series containing loop sequences

    return (list):
    
    ['---AUGGAC-----',
     '----UGGA------',
     '----CGGA------',
     '---AAGGAC-----']
    """
    loops_seqs = list(loop_seqs_series)
    lens = list(loop_seqs_series.apply(lambda x: len(x)))
    gga_indeces = list(loop_seqs_series.apply(lambda x: x.find(motif)))
    l_pad = max(gga_indeces)
    l_pads = [l_pad - gga_index for gga_index in gga_indeces]
    l_aligned_loops = [('-'*l_pad)+loop_seq for l_pad, loop_seq in zip(l_pads, loops_seqs) ]
    max_len = max([len(l) for l in l_aligned_loops])
    aligned_loops = [loop+'-'* (max_len-len(loop)) for loop in l_aligned_loops]
    loop_seqs_df = loop_seqs_series.reset_index()
    new_col_name = loop_seqs_df.columns[1] + '_aligned'
    loop_seqs_df[new_col_name] = aligned_loops
    return loop_seqs_df