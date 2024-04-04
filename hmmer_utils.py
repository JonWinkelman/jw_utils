import re
import pandas as pd

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
    df = pd.DataFrame(lines)
    df.columns = ['target_name','accession1','query_name','accession',
      'E-value','score','bias','E-value_2','score_2','bias_2',
      'exp','reg','clu','ov','env','dom','rep','inc','description_of_target']
    df[['E-value','score','bias','E-value_2','score_2','bias_2', 'exp']] = df[['E-value','score','bias','E-value_2','score_2','bias_2', 'exp']].astype(float)
    df[['reg','clu','ov','env','dom','rep','inc']] = df[['reg','clu','ov','env','dom','rep','inc']].astype(int)
    return df