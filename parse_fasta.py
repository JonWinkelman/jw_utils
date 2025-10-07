#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 24 14:08:25 2022

@author: jonwinkelman

"""

# def get_seq_dict(path_to_fasta):
#     '''
#     parse fast file and return as a dictionary
    
#     parameters:
#         path_to_fasta (str): path to a fasta proteome
        
#         return (dict): eqID as key and sequence as value
#     '''
#     seq = ''
#     seq_dict = {}
#     prot_id = None
#     with open(path_to_fasta, 'r') as f:
#         for line in f:
#             if line[0] == '>':
#                 if prot_id:
#                     seq_dict[prot_id] = seq
#                 prot_id = line.split(' ')[0][1:].strip()
#                 seq = ''
#             else:
#                 seq = seq + line.strip()
#         seq_dict[prot_id] = seq
            
#     return seq_dict

import gzip
from typing import Dict, Callable, Optional, Union, TextIO

def get_seq_dict(
    path_to_fasta: str,
    use_full_header: bool = False,
    key_func: Optional[Callable[[str], str]] = None,
) -> Dict[str, str]:
    """
    Parse a FASTA (optionally gzipped) into a dict: {header_key: sequence}.

    Parameters
    ----------
    path_to_fasta : str
        Path to a FASTA file (.fa, .fasta, .faa, with or without .gz).
    use_full_header : bool, default False
        If False (default), use only the first whitespace-delimited token
        from the header (sans '>') as the key. If True, use the entire
        header line (minus the leading '>') as the key.
    key_func : Callable[[str], str], optional
        Optional function to transform the header string into a key.
        Applied after `use_full_header` logic.

    Returns
    -------
    Dict[str, str]
        Mapping from header-derived key to sequence.
    """
    seq_dict: Dict[str, str] = {}
    key: Optional[str] = None
    seq_parts = []
    path_to_fasta =  str(path_to_fasta)
    # Open transparently, depending on file extension
    open_func = gzip.open if path_to_fasta.endswith(".gz") else open
    mode = "rt"  # text mode (handles universal newlines)

    with open_func(path_to_fasta, mode) as fh:  # type: Union[TextIO, gzip.GzipFile]
        for raw in fh:
            line = raw.strip()
            if not line:
                continue  # skip blanks

            if line.startswith(">"):
                # finalize previous record
                if key is not None:
                    seq_dict[key] = "".join(seq_parts)

                # new header
                header = line[1:].strip()
                if use_full_header:
                    next_key = header
                else:
                    next_key = header.split(None, 1)[0]

                if key_func is not None:
                    next_key = key_func(next_key)

                key = next_key
                seq_parts = []
            else:
                seq_parts.append(line.replace(" ", "").replace("\t", "").upper())

    # finalize last record
    if key is not None:
        seq_dict[key] = "".join(seq_parts)
    else:
        raise ValueError(f"No FASTA headers found in {path_to_fasta!r}")

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


# def write_to_fasta(seq_dict, path, line_size=None):
#     """Write sequence dict to a fasta file.
    
#     seq_dict (dict): {name:sequence}
#     path (str): path for saved file
#     line_size (int): length of line to break the full sequence into in fasta file.
#     """
#     with open(path, 'w') as f:
#         for name, seq in seq_dict.items():
#             f.write(f'>{name}\n')
#             if line_size:
#                 for i in range(int(len(seq)/line_size)):
#                     start = (i*line_size)
#                     end   = start+line_size
#                     f.write(seq[start:end] + '\n')
#                 f.write(seq[end:]+ '\n')
#             else:
#                 f.write(f'{seq}'+ '\n')
def write_to_fasta(seq_dict, path, line_size: int = None):
    """
    Write sequence dict to a FASTA file.

    Parameters
    ----------
    seq_dict : dict
        {header: sequence}, where header is written after ">" exactly as provided
    path : str
        Output FASTA file path
    line_size : int, optional
        Wrap sequence lines at this length. If None, writes the full sequence on one line.
    """
    with open(path, "w") as f:
        for header, seq in seq_dict.items():
            f.write(f">{header}\n")  # write the whole key as the header
            if line_size:
                for i in range(0, len(seq), line_size):
                    f.write(seq[i:i+line_size] + "\n")
            else:
                f.write(seq + "\n")

            
            
            
            
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
            
            
        

def concat_fasta_files(fp_list):
    """Takes multiple fasta files and concatenates them into one dictionary
    
    fp_list (list): A list of filepaths to the individual fasta files to be concatenated
    """

    concat_d = {}
    names = []
    for fp in fp_list:
        
        for name, seq in get_seq_dict(fp).items():
            if name not in names:
                names.append(name)
                concat_d[name] = seq
            else:
                raise Exception(f'{name} is present more than once in sequences')

    return concat_d 


def concatenate_fasta_with_prefix(file_paths, output_file, prefix_list):
    """
    Concatenates multiple proteomes into a single file and appends an assembly accession prefix to each protein ID.

    Parameters:
    file_paths (list): List of file paths to the proteome FASTA files.
    output_file (str): Path to the output file where the concatenated proteome will be saved.
    prefix_list (list of str): List of prefixes to append to each protein ID, matching the order of the file paths.
    """
    if len(file_paths) != len(prefix_list):
        raise ValueError("The number of file paths and prefixes must match.")
    
    with open(output_file, 'w') as outfile:
        for file_path, prefix in zip(file_paths, prefix_list):
            with open(file_path, 'r') as infile:
                for line in infile:
                    if line.startswith('>'):
                        # Append the prefix to the protein ID in the header line
                        outfile.write(f">{prefix}_{line[1:]}")
                    else:
                        outfile.write(line)
        
            