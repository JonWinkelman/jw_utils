import pandas as pd
from collections import Counter

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