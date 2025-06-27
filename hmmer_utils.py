import re
import pandas as pd
import os
import subprocess
import shlex
import tempfile
from jw_utils import parse_fasta 

def run_simple_search(hmm_profile_path, output_fp, proteome_fp, raw_output_fp="hmmsearch_raw.txt"):
    """
    Run hmmsearch and output results in table format.

    Parameters:
        hmm_profile_path (str): Path to the HMM profile.
        output_fp (str): Path to the summary output table (--tblout).
        proteome_fp (str): Path to the FASTA protein file.
        raw_output_fp (str): Path to full hmmsearch output (-o).
    """
    hmmsearch_cmd = [
        "hmmsearch",
        "--tblout", output_fp,
        "-o", raw_output_fp,
        hmm_profile_path,
        proteome_fp
    ]

    try:
        subprocess.run(hmmsearch_cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running hmmsearch: {e}")
        raise




def get_db_versions(pfam_db_fp):
    """return corresponding db version for each base pfam id {base_name:base_name.version} 

    parameters:
    pfam_db_fp (str) path to the pfam database. Database contains one super large text file with 
                     each hmm profile concatenated together. Each profile name has a version which
                     I assume is updated and changes periodically

    return
    (dict) dict with base name as the key and the base.version as the value.
    """
    version_d = {}
    with open(pfam_db_fp, 'r') as f:
        for line in f:
            if line.startswith('ACC'):
                pfam_id = line.split(' ')[-1].strip()
                pfam_base_name = pfam_id.split('.')[0]
                version_d[pfam_base_name] = pfam_id
    return version_d


def parse_domtblout_to_df(filepath):
    """
    Parse an HMMER --domtblout file into a pandas DataFrame with all columns,
    correctly handling the variable-length 'description' column at the end.
    """
    column_names = [
        "target_name", "target_accession", "tlen",
        "query_name", "query_accession", "qlen",
        "full_seq_Evalue", "full_seq_score", "full_seq_bias",
        "domain_number", "domains_in_target",
        "domain_cEvalue", "domain_iEvalue", "domain_score", "domain_bias",
        "hmm_from", "hmm_to", "ali_from", "ali_to", "env_from", "env_to",
        "acc", "description"
    ]

    records = []
    with open(filepath) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            # Split into exactly len(column_names) parts,
            parts = line.rstrip('\n').split(maxsplit=len(column_names) - 1)
            # In case of missing description, pad with empty string
            if len(parts) < len(column_names):
                parts += [''] * (len(column_names) - len(parts))
            records.append(parts)
    
    df = pd.DataFrame(records, columns=column_names)
    return df


def parse_hmmsearch_output(hmm_results_fp):
    """Return dataframe of hmmsearch '--tblout' results 
    
    parameters:
    hmm_results_fp (str): path to output of hmmsearch result file. This file must have 
                          been produced using the `hmmsearch --tblout` flag.
    """
    
    lines = []
    with open (hmm_results_fp, 'r') as f:
        for line in f:
            if not line.startswith('#'):
                cleaned_line = re.sub(r'\s+', ' ', line)
                cleaned_line_lst = cleaned_line.split(' ')
                targe_desc = ' '.join(cleaned_line_lst[18:])
                cleaned_line_lst = cleaned_line_lst[:18]
                cleaned_line_lst.append(targe_desc)
                lines.append(cleaned_line_lst)
    columns = ['target_name','accession1','query_name','accession',
              'E-value','score','bias','E-value_2','score_2','bias_2',
              'exp','reg','clu','ov','env','dom','rep','inc','description_of_target']
    if len(lines)<1:
        return pd.DataFrame(columns=columns)
    else:
        df = pd.DataFrame(lines, columns=columns)
    # df.columns = ['target_name','accession1','query_name','accession',
    #   'E-value','score','bias','E-value_2','score_2','bias_2',
    #   'exp','reg','clu','ov','env','dom','rep','inc','description_of_target']
    df[['E-value','score','bias','E-value_2','score_2','bias_2', 'exp']] = df[['E-value','score','bias','E-value_2','score_2','bias_2', 'exp']].astype(float)
    df[['reg','clu','ov','env','dom','rep','inc']] = df[['reg','clu','ov','env','dom','rep','inc']].astype(int)
    return df


def run_hmm_fetch_search(hmm_id, db_fp, proteome_fp, output_fp, hmm_file_path=None):
    """
    Runs hmmfetch and hmmsearch commands using the given HMM ID, database, and proteome file paths.
    Optionally uses an HMM file path directly if provided.
    
    Parameters:
    - hmm_id: The ID of the HMM to fetch and search.
    - db_fp: File path to the HMM database.
    - proteome_fp: File path to the proteome file.
    - output_fp: File path to save the hmmsearch output.
    - hmm_file_path: Optional; direct file path to an HMM file.
    
    Returns:
    None
    """
    if hmm_file_path:
        # If an HMM file path is provided, use it directly with hmmsearch
        hmmsearch_cmd = f'hmmsearch --tblout {output_fp} {hmm_file_path} {proteome_fp}'
        subprocess.run(hmmsearch_cmd, shell=True)
    else:
        # Constructing the hmmfetch command to fetch the HMM from the database
        hmmfetch_cmd = f'hmmfetch {db_fp} {hmm_id}'
        # Pipe the output of hmmfetch to hmmsearch
        hmmsearch_cmd = f'hmmsearch --tblout {output_fp} - {proteome_fp}'
        
        # Execute the commands, piping hmmfetch output to hmmsearch
        fetch_proc = subprocess.Popen(hmmfetch_cmd, stdout=subprocess.PIPE, shell=True)
        subprocess.run(hmmsearch_cmd, stdin=fetch_proc.stdout, shell=True);
        fetch_proc.stdout.close()


# def run_simple_search(hmm_profile_path, output_fp, proteome_fp):

#     hmmsearch_cmd = f'hmmsearch --tblout {output_fp} -o temp.txt {hmm_profile_path} {proteome_fp}'
#     subprocess.run(hmmsearch_cmd, shell=True)


    

def run_hmm_alignment(db_fp, seqs_fp, hmm_id, output_fp, 
                      output_format='Stockholm', trim=False):
    """
    Fetches an HMM profile using hmmfetch and aligns sequences using hmmalign,
    allowing selection of the output format.

    Parameters:
    db_fp (str): Filepath to the HMM profile database.
    seqs_fp (str): Filepath to the FASTA file containing the sequences to be aligned.
    hmm_id (str): ID of the HMM profile to fetch from the database.
    output_fp (str): Filepath where the aligned sequences will be saved.
    output_format (str): Format of the output alignment (e.g., 'Stockholm', 'Pfam', 'A2M', 'PSIBLAST').

    Returns:
    None: Writes output to a file specified by output_fp.
    """
    try:
        # Set up the hmmfetch command
        fetch_cmd = ['hmmfetch', db_fp, hmm_id]
        # Set up the hmmalign command with the specified output format
        align_cmd = ['hmmalign', '--outformat', output_format] 
        if trim:
            align_cmd.append('--trim')
        align_cmd = align_cmd + ['-o', output_fp, '-', seqs_fp]

        # Execute hmmfetch
        fetch_process = subprocess.Popen(fetch_cmd, stdout=subprocess.PIPE)
        # Pipe the output of hmmfetch to hmmalign and execute
        align_process = subprocess.Popen(align_cmd, stdin=fetch_process.stdout)
        fetch_process.stdout.close()  # Allow fetch_process to receive a SIGPIPE if align_process exits
        align_process.communicate()  # Wait for hmmalign to complete

        if align_process.returncode != 0:
            print("hmmalign failed to complete successfully.")

    except Exception as e:
        print(f"An error occurred: {e}")



import subprocess

def run_simple_hmm_alignment(seqs, hmm_fp, output_fp,
                             output_format='Pfam', trim=False):
    """
    Aligns sequences using hmmalign.

    Parameters:
    seqs (str or dict): Filepath to the FASTA file **or** a dict {header: seq}.
    hmm_fp (str): Filepath to the HMM profile.
    output_fp (str): Filepath where the aligned sequences will be saved.
    output_format (str): Format of the output alignment ('Stockholm','Pfam','A2M','PSIBLAST').
    trim (bool): Whether to enable trimming in hmmalign.

    Returns:
    None: Writes output to a file specified by output_fp.
    """
    # If seqs is a dict, dump to a temp FASTA
    cleanup_tmp = False
    if isinstance(seqs, dict):
        tmp = tempfile.NamedTemporaryFile(delete=False, suffix='.fasta')
        tmp_fp = tmp.name
        tmp.close()
        pfa.write_to_fasta(seqs, tmp_fp)
        seqs_file = tmp_fp
        cleanup_tmp = True
    else:
        seqs_file = seqs

    # Build the hmmalign command
    align_cmd = ['hmmalign', '--outformat', output_format]
    if trim:
        align_cmd.append('--trim')
    align_cmd += ['-o', output_fp, hmm_fp, seqs_file]

    print(f'Executing command: {" ".join(align_cmd)}')

    try:
        subprocess.run(align_cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running hmmalign: {e}")
    finally:
        # Clean up temporary FASTA if we created one
        if cleanup_tmp:
            try:
                os.unlink(tmp_fp)
            except OSError:
                pass












def clean_up_hmmalignment(hmm_alignment_fp):
    """removes non-sequence lines from hmm stockholm output.
    
    returns (dict): standard sequence alignment dictionary
    """
    new_aln = []
    new_aln_d = {}
    with open(hmm_alignment_fp, 'r') as f:
        for line in f:
            if line[0] != '#':
                if len(line) > 5:
                    if line[0] != '/':
                        line_lst = line.strip().split(' ')
                        seq = line_lst[-1].replace('.', '-')
                        new_aln_d[line_lst[0]] = seq.upper()
    return new_aln_d





