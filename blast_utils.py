import subprocess
import pandas as pd



def make_blast_db(input_file, dbtype, database_name):
    """
    input_file (str): path to a fasta file
    dbtype (str): nucl, prot
    database_name (str): desired name of the database output file
    """
    
    cmd = ['makeblastdb',
           '-in', input_file,
           '-dbtype', dbtype,  
           '-out', database_name]
    subprocess.run(cmd)



def blast_prot_seqs(database_name, input_file, output_file):
    """
    
    database_name (str): fp to base name of blast db files
    input_file (str): fasta file with one or more sequences
    output_file (str): desired name of output file
    """
    
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
    hits_df.to_csv(output_file)

    return hits_df