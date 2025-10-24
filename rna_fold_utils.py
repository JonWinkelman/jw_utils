
import re
import pandas as pd
import shlex
import subprocess
import os
from jw_utils import parse_fasta as pfa



def fold_RNA(input_fasta_fp, output_file, folding_temp=37, flags = '--noPS'):
    """Takes rna sequences from fasta file and folds at a given temp.
    
    input_fasta_fp (str): fasta file contianing rna sequences
    output_file (str): path of file you would like to output to
    folding_temp (int): temperature to fold RNA. Default=37
    noPS (bool): default True. If False, will output a Postscript file for each rna
    flags (str): e.g. '-d2 --noLP'
    """
    
    if not flags:
        flags=''
    
    if not os.path.exists(input_fasta_fp):
        raise FileNotFoundError(f"{input_fasta_fp} not found")
    
    cmd = f"rnafold -i {input_fasta_fp} -T {folding_temp} {flags}"
    cmd = shlex.split(cmd)
    print(cmd)

    with open(output_file, 'w') as file:
        # Execute the command and redirect stdout to the file
        result = subprocess.run(cmd, stdout=file, text=True)

    if result.returncode != 0:
        print("Error occurred during command execution")
    

def parse_rnafold_fasta(fasta_fp, prot_name='flag'):
    """parse raw RNAfold output fasta into a dataframe with seq, structure, and free energy columns, and indexed by protein ID"""
    d = {}
    with open(fasta_fp, 'r') as f:
        for i, line in enumerate(f):
            if i%3 == 0: 
                l_id = line[1:].strip()
            elif (i-1)%3 == 0:
                seq=line.strip()
            elif (i-2)%3 == 0:
                ln_lst = line.split(' ')
                struct=ln_lst[0]
                energy =  float(ln_lst[-1].replace('(','').replace(')','').strip())
                d[l_id] = [seq, struct, energy]
    return pd.DataFrame(d, index=[f'{prot_name}_rna_seqs', f'{prot_name}_fold_structure', f'{prot_name}_free_energy']).transpose()


    

def find_simple_hairpin_indices(structure, min_loop_nt=3, min_stem_len=3):
    """Return list of tuples of start and end indices for an rna hairpin from a structure string
    
    Finds indices of simple hairpin structures with a specified minimum loop length in an RNA sequence.
    Note: This function identifies hairpins based on the loop length and provides start and end indices
    for each identified hairpin structure.

    Parameters:
    - structure (str): The RNA structure representation as a string.
    - min_loop_nt (int): The minimum number of unpaired nucleotides in the loop part of the hairpin.

    Returns:
    - list of tuples: A list containing tuples of (start, end) 0-index indices for matched simple hairpin structures.

    Limitations:
    - The function uses regular expressions which cannot inherently balance parentheses to match
      exact stem lengths. It identifies patterns that resemble hairpins based on loop length.
    - It is designed for simple hairpin structures and may not accurately capture all valid
      hairpin structures in complex scenarios with nested structures or pseudoknots.
    """
    
    # Constructing the regex pattern for the loop
    loop_pattern = r'\.{' + str(min_loop_nt) + r',}'
    # Constructing a simplified pattern that looks for hairpin structures
    hairpin_pattern = r'\({3,8}' + loop_pattern + r'\){3,8}'
    
    simple_hairpin_indices = []
    for match in re.finditer(hairpin_pattern, structure):
        start, end = match.start(), match.end()-1
        simple_hairpin_indices.append((start, end))
    simple_hairpin_indices
    
    return simple_hairpin_indices



def find_simple_hairpin_loop_indices(structure, min_loop_nt=3):
    """
    Finds indices of simple hairpin structures with a specified minimum loop length in an RNA sequence.
    Note: This function identifies hairpins based on the loop length and provides start and end indices
    for each identified hairpin structure.

    Parameters:
    - structure (str): The RNA structure representation as a string.
    - min_loop_nt (int): The minimum number of unpaired nucleotides in the loop part of the hairpin.

    Returns:
    - list of tuples: A list containing tuples of (start, end) 0-index indices for matched simple hairpin structures.

    Limitations:
    - The function uses regular expressions which cannot inherently balance parentheses to match
      exact stem lengths. It identifies patterns that resemble hairpins based on loop length.
    - It is designed for simple hairpin structures and may not accurately capture all valid
      hairpin structures in complex scenarios with nested structures or pseudoknots.
    """
    
    # Constructing the regex pattern for the loop
    loop_pattern = r'\.{' + str(min_loop_nt) + r',}'
    # Constructing a simplified pattern that looks for hairpin structures
    hairpin_pattern = r'\(\(\(' + loop_pattern + r'\)\)\)'
    
    simple_hairpin_indices = []
    for match in re.finditer(hairpin_pattern, structure):
        start, end = match.start()+3, match.end()-4
        simple_hairpin_indices.append((start, end))
    simple_hairpin_indices
    
    return simple_hairpin_indices




def get_loop_seq(rna_struc,rna_seq, min_loop_nt=3):
    """"""
    simple_hairpin_indices = find_simple_hairpin_loop_indices(rna_struc, min_loop_nt=min_loop_nt)
    loop_seqs = []
    for hairpin_indeces in simple_hairpin_indices:
        loop_seq = rna_seq[hairpin_indeces[0]:hairpin_indeces[1]+1].upper()
        loop_seqs.append(loop_seq)
    return loop_seqs


def get_loop_seqs(rna_strucures,rna_seqs, min_loop_nt=3):
    """"""
    
    loop_seqs_lst = []  
    for rna_struc, rna_seq in zip(rna_strucures, rna_seqs):
        loop_seqs_lst.append( get_loop_seq(rna_struc,rna_seq) )
    return loop_seqs_lst



def get_hairpin_seqAndstruct(rna_struc,rna_seq, min_loop_nt=3):
    """Returns list of hp sequences and corresponding list of associated structures"""
    simple_hairpin_indices = find_simple_hairpin_indices(rna_struc, min_loop_nt=min_loop_nt)
    hp_seqs = []
    hp_structs = []
    for hairpin_indeces in simple_hairpin_indices:
        hp_seq = rna_seq[hairpin_indeces[0]:hairpin_indeces[1]+1].upper()
        hp_struct = rna_struc[hairpin_indeces[0]:hairpin_indeces[1]+1]
        hp_seqs.append(hp_seq)
        hp_structs.append(hp_struct)
    return hp_seqs, hp_structs






