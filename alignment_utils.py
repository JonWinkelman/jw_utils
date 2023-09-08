from jw_utils import parse_fasta as pfa


def delete_alignment_positions(seq_dict, positions_to_delete):
    """delete position in alignment, 0-based and return alignment dict
    
    parameters:
    seq_dict (dict): {seq_id:aligned sequence}
    positions_to_delete (list or other iterable): positions/columns (zero-based) in the alignment to delete
    """

    seq_dict_edit = {}
    for name, seq in seq_dict.items():
        seq_edit = ''
        for i, val in enumerate(seq):
            if i not in positions_to_delete:
                seq_edit += val
        
        seq_dict_edit[name] = seq_edit
    return seq_dict_edit



def get_value_frequency(seq_dict, value_of_interest='-'):
    """Return dict of columns and fraction of a values occurrence in that column
    
    parameters
    seq_dict (dict) {seq_id:aligned sequence}
    value_of_interest (str): value of interest to return frequency for, e.g. for gaps enter '-',

    return:
    frequency dict, {
    """
    num_seqs = len(seq_dict)
    aln_length = len(list(seq_dict.values())[0])
    freq_dict = {}
    cols = {i:[] for i in range(aln_length)}
    for name, aln_seq in seq_dict.items():
        for i, col in enumerate(aln_seq):
            cols[i].append(col)
    i=0
    for column, values in cols.items():
        fraction = len([val for val in values if val == value_of_interest])/num_seqs
        freq_dict[i] = fraction
        i+=1
    return freq_dict


def positions_to_delete(aln_seq_dict, threshold, value_of_interest='-'):
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
    positions_to_delete = []
    freq_dict = get_value_frequency(aln_seq_dict, value_of_interest=value_of_interest)
    for pos, freq in freq_dict.items():
        if freq > threshold:
            positions_to_delete.append(pos) 
    return positions_to_delete


def delete_alignment_positions(aln_seq_dict, threshold=0.8, value_of_interest='-' ):
    """delete position in alignment, 0-based

    aln_seq_dict (dict):
    threshold (int): if the frequency of value of interest is above the threshold value, 
            then that alignment position is added to list
    value_of_interest (str): value of interest to obtain frequency for in columns. e.g. to delete columns
            with frequecy of gaps greater than the threshold, enter '-'.


    return
    edited sequence alignment dictionary with offending columns removed.
    """

    positions = positions_to_delete(aln_seq_dict, threshold=threshold, value_of_interest='-')
    seq_dict_edit = {}
    for name, seq in aln_seq_dict.items():
        seq_edit = ''
        for i, val in enumerate(seq):
            if i not in positions:
                seq_edit += val
        
        seq_dict_edit[name] = seq_edit
    return seq_dict_edit