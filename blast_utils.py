import subprocess
import pandas as pd



def make_blast_db(input_file, dbtype, database_name):
    """"""
    
    cmd = ['makeblastdb',
           '-in', input_file,
           '-dbtype', dbtype,  
           '-out', database_name]
    subprocess.run(cmd)



def blast_prot_seqs(database_name, input_file, output_file):
    """"""
    
    cmd = ['blastp',
           '-db', database_name,
           '-query', input_file,
            '-outfmt', "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
           '-out', output_file]
    subprocess.run(cmd)
    hits_df = pd.read_csv(output_file, sep='\t', header=None)
    hits_df.columns = ['qseqid', 'sseqid', 'pident', 'length', 
                       'mismatch', 'gapopen','qstart', 'qend', 
                       'sstart', 'send', 'evalue', 'bitscore']

    return hits_df