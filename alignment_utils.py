from jw_utils import parse_fasta as pfa
import pandas as pd
import tempfile
from Bio import AlignIO
from Bio.Align import AlignInfo


def get_frequency_dict(aln_dict):
    """retrun a dict with the frequency of each character in each column of alignment
    
    parameters
    aln_dict : dict
        dictionary with sequence id as key and aligned sequence as value

    return
        nested dictionary with main key as aln column,and nested d as
        with character as key and frequency as value
        
        e.g. for the zero-ith column in alignment:
        d[0] = {'A': 0.25,  'C': 0,  'G': 0.5.,  'T': 0.25,  'N': 0, '-': 0}
    """
    
    df = pd.DataFrame.from_dict(aln_dict, orient='index', columns=['sequence'])
    df = df['sequence'].apply(list).apply(pd.Series)
    num_sequences = df.shape[0]
    
    def make_blank_char_count_dict():
        return {
            'A': 0,  # Adenine
            'C': 0,  # Cytosine
            'G': 0,  # Guanine
            'T': 0,  # Thymine
            'N': 0,  # Any base (A or C or G or T)
            '-': 0,  # gap in alignment
        }
        
    char_counts_ineach_column = {}    
    for i,col in enumerate(df):
        character_counts = make_blank_char_count_dict()
        col_lst = list(df[col])
        unique_chars = set(col_lst)
        if len(unique_chars) == 1:
            character_counts[ list(unique_chars)[0] ] = 1
            char_counts_ineach_column[i] = character_counts
        else:
            d = df[col].value_counts().to_dict()
            for base, count in d.items():
                freq = count/num_sequences
                character_counts[base] = freq
                char_counts_ineach_column[i] = character_counts
    return char_counts_ineach_column
    


def positions_to_delete(aln_dict, threshold=1, character_of_interest='-'):
    """"Return a list of columns (0-based) in alignment that are above threshold.

    parameters
    aln_seq_dict (dict)
    threshold (int): if the frequency of value of interest is above the threshold value, 
                    then that alignment position is added to list
    value_of_interest (str): value of interest to return frequency for, e.g. for gaps enter '-',

    return
    list of ints corresponding to columns in alignment where the frequency of the value_of_interest
        exceeds the threshold
    
    """
    columns_to_delete = []
    columns_freq_dict = get_frequency_dict(aln_dict)
    
    for column, freq_d in columns_freq_dict.items():
        for char, freq in freq_d.items():
            if char == character_of_interest:
                if freq >= threshold:
                    columns_to_delete.append(column) 
    return columns_to_delete


def delete_alignment_positions(aln_dict, threshold=0.8, character_of_interest='-' ):
    """delete 0-based position in alignment if character_of_interest's frequency in column is >= threshold.
     **useful for removing columns that contain all or mostly gaps, or removing positions that 
     do not have phylogenetic information
    aln_dict (dict):
    threshold (int): if the frequency of value of interest is above the threshold value, 
            then that alignment position is added to list
    value_of_interest (str): value of interest to obtain frequency for in columns. e.g. to delete columns
            with frequecy of gaps greater than the threshold, enter '-'.


    return
    edited sequence alignment dictionary with offending columns removed.
    """

    positions = positions_to_delete(aln_dict, threshold=threshold, character_of_interest=character_of_interest)
    df = pd.DataFrame.from_dict(aln_dict, orient='index', columns=['sequence'])
    df = df['sequence'].apply(list).apply(pd.Series)
    df = df.drop(columns=positions)
    seq_dict_edit = df.apply(lambda row: ''.join(row), axis=1).to_dict()
    return seq_dict_edit


def collapse_alignment(aln_d):
    """"""
    collapsed_aln_d = {}
    for name, seq in aln_d.items():
        n_seq = seq.replace('-','')
        collapsed_aln_d[name] = n_seq
    return collapsed_aln_d



def alignment_to_matrix(aln_fp, molecule_type='prot', to_type='probability'):
    """Make a matrix from an alignment file.
    
    aln_fp        (str): path to alignment file. 
    molecule_type (str): 'prot', 'dna', 'rna'
    to_type       (str): 'counts', 'probability'

    return (pandas.DataFrame): df where columns represent columns in alignment, 0-based, and
                               each row is a character, e.g. if molecule_type='prot' then
                               one of the 20 AAs. 
    """
    
    alignment = AlignIO.read(aln_fp, "fasta")
    summary_align = AlignInfo.SummaryInfo(alignment)
    consensus = summary_align.dumb_consensus()
    freq_matrix = summary_align.pos_specific_score_matrix(consensus, chars_to_ignore=['-'])
    num_seqs = len(alignment)
    type='prot'
    if molecule_type=='prot':
        chars = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    if molecule_type=='dna':
        chars = ['A', 'C', 'T', 'G']
    if molecule_type=='rna':
        chars = ['A', 'C', 'U', 'G']
    col_names = []
    col = {char:[] for char in chars}
    for position, scores in enumerate(freq_matrix):
        col_names.append(position)
        for base, score in scores.items():
            if to_type == 'probability':
                col[base].append(score/num_seqs)
            elif to_type == 'counts':
                col[base].append(score)
    
    return pd.DataFrame(col).transpose()




def dict_alignment_to_matrix(aln_dict, molecule_type='prot', to_type='probability'):
    """Make a matrix from an alignment dict.
    
    aln_dict      (dict): aligment dict, {id:seq} 
    molecule_type (str): 'prot', 'dna', 'rna'
    to_type       (str): 'counts', 'probability'

    return (pandas.DataFrame): df where columns represent columns in alignment, 0-based, and
                               each row is a character, e.g. if molecule_type='prot' then
                               one of the 20 AAs. 
    """
    with tempfile.NamedTemporaryFile(mode='w+', delete=True) as tmp_file:
        print(f'temp_file created at {tmp_file.name}')
        pfa.write_to_fasta(aln_dict, tmp_file.name)
        tmp_file.seek(0) #moves the file pointer back to start after writing
        return alignment_to_matrix(tmp_file.name, molecule_type, to_type)