#!/usr/bin/env python
import os
from jw_utils import parse_fasta as pfa
from jw_utils import file_utils as fu
import subprocess
import pandas as pd
from pathlib import Path


def full_rblast(proteome_dir, query_file,query_protein_name, 
                reference_proteome_path, blast_db_name=None, 
                proteomes_to_include=None, suffix='.faa'):
    """
    Perform a full reciprocal BLASTp search for orthologous gene detection.

    Args:
        proteome_dir (str): Directory containing individual proteome files.
        query_file (str): Path to the initial query file.
        query_protein_name (str): Name of the query protein.
        reference_proteome_path (str): Path to the reference proteome file.
        blast_db_name (str, optional): Directory to store BLAST databases.
        proteomes_to_include (list, optional): Subset of proteomes to use.
        suffix (str): File extension of proteomes.

    Returns:
        DataFrame: Summary of best hits across proteomes.
    """
    
    if not blast_db_name:
        blast_db_name = './blast_db'
    rblast_results_dir = './rblast_results_' + Path(query_file).stem
    make_blastp_databases(proteome_dir, blast_db_name, proteomes_to_include)
    reciprocal_blastp(blast_db_name, 
                     query_file,
                     rblast_results_dir,
                     proteome_dir, 
                     reference_proteome_path,
                     proteomes_to_include,
                     suffix
                    )
    return get_full_best_hits_dict(rblast_results_dir, query_protein_name, accession_length=15)


def reciprocal_blastp(database_dir, query_file, out_dirname, proteome_dir, reference_proteome_path, accessions_to_include, suffix='.faa'):
    """
    Perform forward and reciprocal BLAST searches.

    Args:
        database_dir (str): Directory with individual BLAST databases.
        query_file (str): Path to initial query file.
        out_dirname (str): Directory to store BLAST results.
        proteome_dir (str): Directory containing individual proteome files.
        reference_proteome_path (str): Path to reference proteome.
        accessions_to_include (list): Subset of proteomes to query.
        suffix (str): Proteome file suffix.

    Raises:
        Exception: If a proteome in 'proteomes_to_include' or 'databases' is missing.
    """
    
    databases  = [db for db in os.listdir(database_dir) if db != '.DS_Store']
    if not accessions_to_include:
        proteomes = [Path(p).stem for p in os.listdir(proteome_dir) if p.endswith(suffix)]
        if len(proteomes) == 0: 
            print(f'no proteomes with suffix {suffix} were found.')
    else:
        proteomes = accessions_to_include
    for proteome in proteomes:
        if proteome not in databases:
            raise Exception(f'{proteome} not present in databases')
    for database in databases:
        if database not in proteomes:
            raise Exception(f'{database} not present in proteomes')
            
    os.makedirs(out_dirname, exist_ok=True)
    forward_results_dir = os.path.join(out_dirname, 'forward_blast_results')
    os.makedirs(forward_results_dir, exist_ok=True)
    print(forward_results_dir)
    multi_database_query(database_dir, query_file, forward_results_dir)

    rblast_queries_dir = os.path.join(out_dirname, 'forward_best_hits')
    write_best_hit_queries(proteome_dir, forward_results_dir, rblast_queries_dir)
    
    
    ref_db_name = reference_proteome_path.replace('.faa', '_blast_db')
    #os.makedirs('reference_db')
    reference_db_path = os.path.join(out_dirname, 'reference_db/reference_db')
    print(reference_db_path)
    create_blast_database(reference_proteome_path, reference_db_path, dbtype='prot')
    recip_blast(rblast_queries_dir, reference_db_path, os.path.join(out_dirname,'reciprocal_blast_results'))


def create_blast_database(proteome_path, out_fp, dbtype='prot'):
    """
    Create a BLAST database from a reference proteome.

    Args:
        proteome_path (str): Path to reference proteome file.
        out_fp (str): Output path for the BLAST database.
        dbtype (str): Type of database ('prot' or 'nucl').
    """
    try:
        # Run the makeblastdb command using subprocess
        subprocess.run(
            ['makeblastdb', '-in', proteome_path, '-dbtype', dbtype, '-out', out_fp],
            check=True  # Raise an exception if the command fails
        )
        print(f"BLAST database created successfully of {proteome_path} at {out_fp} (type: {dbtype})")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred: {e}")





def make_blastp_databases(proteome_dir, out_dir, proteomes_to_include=None, dbtype='prot', proteome_suffix='.faa'):
    """
    Create BLASTP databases for each proteome in a directory.

    Args:
        proteome_dir (str): Directory with individual proteome files.
        out_dir (str): Directory to store created databases.
        proteomes_to_include (list, optional): Subset of proteomes to include.
        dbtype (str): Database type ('prot' or 'nucl').
        proteome_suffix (str): File extension of proteomes.
    """
    
    if not proteomes_to_include:
        proteomes_to_include = [Path(proteome).stem for proteome in os.listdir(proteome_dir)]

    for acc in proteomes_to_include:
        proteome_path = f'{proteome_dir}/{acc}.faa'
        db_name = f"{acc}_db"
        out = os.path.join(out_dir, acc, db_name)

        os.makedirs(os.path.dirname(out), exist_ok=True)

        # Run makeblastdb using subprocess
        create_blast_database(proteome_path, out, dbtype='prot')
        # subprocess.run(
        #     ['makeblastdb', '-in', proteome_path, '-dbtype', 'prot', '-out', out],
        #     check=True
        # )


def multi_database_query(database_dir, query_filepath, out_dir_path=None):
    """
    Query each BLASTP database and store hits in the output directory.

    Args:
        database_dir (str): Directory containing individual BLAST databases.
        query_filepath (str): Path to query file.
        out_dir_path (str): Directory to store query results.
    """
    
    if not out_dir_path:
        out_dir_path = './forward_blast_results'
    os.makedirs(out_dir_path, exist_ok=True)

    for acc in [acc for acc in os.listdir(database_dir) if acc != '.DS_Store']:
        db_filepath = os.path.join(database_dir, f'{acc}/{acc}_db')
        out_filepath = os.path.join(out_dir_path, f'{acc}_blastpOut.txt')

        # Run blastp using subprocess
        subprocess.run(
            ['blastp', '-query', query_filepath, '-db', db_filepath, 
             '-num_descriptions', '5', '-num_alignments', '5', '-out', out_filepath],
            check=True
        )


def recip_blast(rblast_queries_dir, rblast_db_path, out_dir):
    """
    Perform reciprocal BLAST searches for queries against a reference database.

    Args:
        rblast_queries_dir (str): Directory with queries for reciprocal BLAST.
        rblast_db_path (str): Path to reference BLAST database.
        out_dir (str): Directory to store reciprocal BLAST results.
    """
    
    os.makedirs(out_dir, exist_ok=True)

    query_files = [file for file in os.listdir(rblast_queries_dir) if file.endswith('.faa')]

    for query in query_files:
        out_filename = query.replace('.faa', '_blastpOut.txt')
        out_filepath = os.path.join(out_dir, out_filename)
        query_filepath = os.path.join(rblast_queries_dir, query)

        # Run blastp using subprocess
        subprocess.run(
            ['blastp', '-query', query_filepath, '-db', rblast_db_path, 

             '-num_descriptions', '5', '-num_alignments', '5', '-out', out_filepath],
            check=True
        )



def write_best_hit_queries(proteome_dir, bbhits, query_dir=None):
    """
    Write each best hit to a separate FASTA file.

    Args:
        proteome_dir (str): Path to proteome directory.
        bbhits (str): Directory with BLAST output files.
        query_dir (str, optional): Directory to store output FASTA files.

    Returns:
        dict: Nested dictionary of best hit sequences {accession: {protein_id: sequence}}.
    """
    
    best_hit_dict_tot = get_sequence(proteome_dir, bbhits)
    if not query_dir:
        query_dir = './query_fasta_files'
    os.makedirs(query_dir, exist_ok=True)
    for acc, seq_dict in best_hit_dict_tot.items():
        with open(os.path.join(query_dir,f'{acc}.faa'), 'w') as f:
            for prot_id, seq in seq_dict.items():
                f.write(f'>{prot_id}\n')
                f.write(f'{seq}\n')
    return best_hit_dict_tot


def make_besthit_dict(path_to_blastpOut):
    """return a dict {accession:protein_ID} with best blast hit for each genome in the output dir"""
    files = fu.listdir(path_to_blastpOut,end_to_keep='.txt')
    line_of_interest=None
    best_hit_dict = {}
    for file in files:
        path=os.path.join(path_to_blastpOut, file)
        with open(path, 'r') as f:
            accession = file[:15]
            for i, line in enumerate(f):
                if line.startswith('Sequences producing significant alignments:'):
                    line_of_interest = i+2
                if i==line_of_interest:
                    best_hit_dict[accession]=line.split(' ')[0]
    return best_hit_dict


def get_sequence(proteome_dir, path_to_blastpOut):
    """return {accession:{protein_ID:seq}} for each best hit    """
    best_hit_dict = make_besthit_dict(path_to_blastpOut)
    #get protein seq for each best hit
    best_hit_dict_tot = {}
    for acc, prot_id in best_hit_dict.items():
        proteome_path = os.path.join(proteome_dir, f'{acc}.faa')
        d = pfa.get_seq_dict(proteome_path)
        seq = d.get(prot_id, None)
        best_hit_dict_tot[acc] = {prot_id:seq}
        if not seq:
            print(acc, prot_id)
    return best_hit_dict_tot



def best_hit_to_dict(path_to_rblast_results):
    """Return a dict with the accession of query genome and the best hit in the source genome

    These are the results of the reciprical blast. For each genome, there was a best hit in the initial
    blastp. This best hit was reciprical blasted against the genome containing the original query. If best
    hit in this reciprical blastp is the initial query, then these proteins are likely orthologs.
    
    return (dict): {accession:protein_id}  
    """
    paths = fu.get_filepaths_in_dir(path_to_rblast_results) 
    d = {}
    for path in paths:
        line_of_interest = None
        accession = path.split('/')[-1][:15]
        with open(path, 'r') as f:
            for i,line in enumerate(f):
                if line.startswith('>'):
                    prot_id = line.split(' ')[0][1:]
                    d[accession]=prot_id
                    print(prot_id)
    return d



### processing results
def get_blastbesthit_dict(path_to_blast_results, accession_length=13):
    """Return a dict {accession:protein_id} with the accession of genome and the best hit in the source genome"""
    
    files =os.listdir(path_to_blast_results)
    paths = [os.path.join(path_to_blast_results,file) for file in files]
    best_hits_dict = {}
    for file, path in zip(files, paths):
        accession = file[:accession_length]
        with open(path, 'r') as f:
            for i,line in enumerate(f):
                if line.startswith('>'):
                    prot_id = line.split(' ')[0][1:]
                    best_hits_dict[accession]=prot_id
                    break
    return best_hits_dict

def get_full_best_hits_dict(recip_blast_dir, query_protein_name, accession_length=15):
    for_dir = os.path.join(recip_blast_dir, 'forward_blast_results')
    rev_dir = os.path.join(recip_blast_dir, 'reciprocal_blast_results')
    forward_bbhits = get_blastbesthit_dict(for_dir, accession_length=accession_length)
    rev_bbhits = get_blastbesthit_dict(rev_dir, accession_length=accession_length)
    full_df = pd.DataFrame.from_dict(rev_bbhits, orient='index').merge(pd.DataFrame.from_dict(forward_bbhits, orient='index'), left_index=True, right_index=True)
    full_df.columns = ['r_bbhit', 'f_bbhit']
    full_df['forward_query'] = [query_protein_name for row in range(full_df.shape[0])]
    #full_df['recip_best_hit'] = full_df['r_bbhit'].apply(lambda x: True if x==query_protein_name else False)
    full_df['recip_best_hit'] = full_df['forward_query'] == full_df['r_bbhit']
    full_df = full_df[['forward_query','f_bbhit', 'r_bbhit', 'recip_best_hit']]
    full_df.index.name= 'target_proteome'
    return full_df
    