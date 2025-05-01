import subprocess
import os
import zipfile
import json
import pandas as pd
import tempfile
import shutil
from pathlib import Path

def copy_ncbi_files(ncbi_data_dir, dest_dir, suffix = '.faa'):
    """Copies internal files from a ncbi data dir to a given dir
    
   
   Assumptions:
   1) path to proteome is <ncbi_data_dir>/GC#_#########.d/protein.faa
   2) proteome ends in .faa, renames stem to <accession>.faa
   3) internal dirs start with 'GC'
    """
    proteomes_copied = []
    ncbi_data_dir = Path(ncbi_data_dir)
    dest_dir = Path(dest_dir)
    for dir in [f for f in ncbi_data_dir.glob('GC*')]:
        #copy each proteome from dataset to ./data/Proteomes
        suffix = suffix.strip('.')
        src_fp = [f for f in dir.glob(f'*{suffix}')]
        if len(src_fp) > 1:
            raise Exception(f'{len(src_fp)} fasta files with the suffix {suffix} were found')
        src_fp = src_fp[0]
        if src_fp.exists():
            dest_fp = dest_dir / f'{dir.name}.{suffix}'
            shutil.copy(src_fp,dest_fp)
            proteomes_copied.append(dest_fp.name)
        else:
            print(f'source file "{src_fp}" does not exist!')
    return proteomes_copied



def download_genomes_from_accfile(accessions_fp, files_to_include, dataset_fp):
    """Download ncbi datasets genomes from input file containing accessions.

        goes through dehydration, unzipping and rehydration steps.

    accessions_fp: str, the path to the file containing ncbi assembly accessions
    files_to_include: str, the types of files to include (e.g., 'genome,protein,gff3').
            should not contain any spaces.
            * genome:     genomic sequence
            * rna:        transcript
            * protein:    amnio acid sequences
            * cds:        nucleotide coding sequences
            * gff3:       general feature file
            * gtf:        gene transfer format
            * gbff:       GenBank flat file
            * seq-report: sequence report file
            * none:       do not retrieve any sequence files
    dataset_fp: str, the filename for the output zip file, file must end with '.zip'
    """
    if not dataset_fp.endswith('.zip'):
        raise Exception(f'change {dataset_fp} to end with ".zip"')
    files_to_include = files_to_include.replace(' ','')
    download_dehydrated_ncbi_dataset(accessions_fp, files_to_include, dataset_fp)
    extract_path = dataset_fp.replace('.zip', '')
    unzip_ncbi_dir(dataset_fp, extract_path)
    rehydrate_ncbi_dir(extract_path)
    

def download_dehydrated_ncbi_dataset(accessions_fp, files_to_include, dataset_fp='ncbi_dataset.zip'):
    """
    Downloads an NCBI dataset using specified accessions and file types.
    
    accessions_fp: str, the path to the file containing ncbi assembly accessions
    files_to_include: str, the types of files to include (e.g., 'genome,protein,gff3').
            should not contain any spaces.
            * genome:     genomic sequence
            * rna:        transcript
            * protein:    amnio acid sequences
            * cds:        nucleotide coding sequences
            * gff3:       general feature file
            * gtf:        gene transfer format
            * gbff:       GenBank flat file
            * seq-report: sequence report file
            * none:       do not retrieve any sequence files
    dataset_fp: str, the filename for the output zip file
    """
    # Prepare the command to run
    command = [
        'datasets', 'download', 'genome', 'accession',
        '--inputfile', accessions_fp,
        '--include', files_to_include,
        '--dehydrated',
        '--filename', dataset_fp
    ]
    
    # Execute the command
    result = subprocess.run(command, capture_output=True, text=True)
    
    # Print the output and errors
    if result.stdout:
        print("Output:", result.stdout)
    if result.stderr:
        print("Error:", result.stderr)


def unzip_ncbi_dir(zip_path, extract_path):
    
    if not os.path.exists(extract_path):
        os.makedirs(extract_path)
        
    # Open the zip file
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        # Extract all the contents into the directory
        zip_ref.extractall(extract_path)
    print(f"Contents extracted to {extract_path}")


def rehydrate_ncbi_dir(dehydrated_dirname):
    command = ["datasets","rehydrate",
              "--directory", dehydrated_dirname
              ]
    subprocess.run(command)



def download_assembly_summaries(accessions_fp, output_file='./summaries.json'):
    """
    Runs the NCBI datasets summary command for genome accessions and outputs the result to a JSON file.
    
    accessions_fp: str, the path to the file containing accessions
    output_file: str, the path to the output JSON file
    """
    # Build the command
    command = [
        'datasets', 'summary', 'genome', 'accession',
        '--inputfile', accessions_fp
    ]
    
    # Open the output file in write mode
    with open(output_file, 'w') as file:
        # Execute the command and redirect stdout to the file
        result = subprocess.run(command, stdout=file, text=True)

    # Check for errors
    if result.returncode != 0:
        print("Error occurred during command execution")


def read_json(json_fp):
    with open(json_fp, 'r') as file:
        data_dict = json.load(file)
    return data_dict

def parse_ncbi_summary(json_fp):
    """"""
    ncbi_summary_d = read_json(json_fp)
    summary_d = {}   
    for report in ncbi_summary_d['reports']:
        accession = report['accession'] 
        summary_d[accession] = report
    return summary_d


def make_summary_df(json_fp):
    summary_d = parse_ncbi_summary(json_fp)
    cols = {'organism_name':[], 'tax_id':[]}
    accession_lst = []
    for acc, report in summary_d.items():
        accession_lst.append(acc)
        for col_name, lst in cols.items():
            lst.append(report['organism'][col_name])

    summary_df = pd.DataFrame(cols)
    summary_df.index = accession_lst
    return summary_df 

def write_summary_to_csv(summary_json_fp, out_fp):
    df = make_summary_df(summary_json_fp)
    df.to_csv(out_fp)


def download_genomes_from_acclist(accessions, 
                                  files_to_include='genome,protein,gff3',  
                                  dataset_fp="ncbi_dataset.zip"):
    
    """Download ncbi datasets genome assemblies from accession list or other iterable.

    goes through dehydration, unzipping and rehydration steps.

    accessions_fp: str, the path to the file containing ncbi assembly accessions
    files_to_include: str, the types of files to include (e.g., 'genome,protein,gff3').
            should not contain any spaces.
            * genome:     genomic sequence
            * rna:        transcript
            * protein:    amnio acid sequences
            * cds:        nucleotide coding sequences
            * gff3:       general feature file
            * gtf:        gene transfer format
            * gbff:       GenBank flat file
            * seq-report: sequence report file
            * none:       do not retrieve any sequence files
    dataset_fp: str, the filename for the output zip file, file must end with '.zip'
    """
    if not dataset_fp.endswith('.zip'):
        raise Exception(f'change {dataset_fp} to end with ".zip"')

    with tempfile.NamedTemporaryFile(mode='w+', delete=True) as tmp_file:
        print(f'temp_file created at {tmp_file.name}')
        for acc in accessions:
            tmp_file.write(acc+'\n')
        tmp_file.seek(0) #moves the file pointer back to start after writing
        download_genomes_from_accfile(tmp_file.name,files_to_include,dataset_fp)




def download_assembly_summaries_from_list(accessions, output_file='./summaries.json'):
    """
    Downloads ncbi datasets summaries from accessions and outputs the result to a JSON file.
    
    accessions_fp: str, the path to the file containing accessions
    output_file: str, the path to the output JSON file
    """
    with tempfile.NamedTemporaryFile(mode='w+', delete=True) as tmp_file:
        print(f'temp_file created at {tmp_file.name}')
        for acc in accessions:
            tmp_file.write(acc+'\n')
        tmp_file.seek(0) #moves the file pointer back to start after writing
        download_assembly_summaries(tmp_file.name, output_file)


def move_proteomes(data_dir, new_proteome_dir = './Proteomes'):
    """Moves proteomes from the nested ncbi data dirs to own dir and renames with accesion.faa"""
    
    os.makedirs(new_proteome_dir, exist_ok=True)
    no_proteome = []
    accs = [dir for dir in os.listdir(data_dir) if dir.startswith('GC')]
    for acc in accs:
        fp = os.path.join(data_dir, acc, 'protein.faa')
        if os.path.exists(fp):
            new_fp = os.path.join(new_proteome_dir, f'{acc}.faa')
            shutil.move(fp,new_fp )
        else:
            no_proteome.append(acc)
    return no_proteome

def move_genomes(data_dir, new_genome_dir = './Genomes'):
    """Moves genomes from the nested ncbi data dirs to own dir and renames with accesion.fna"""
    
    os.makedirs(new_genome_dir, exist_ok=True)
    no_genome = []
    accs = [dir for dir in os.listdir(data_dir) if dir.startswith('GC')]
    for acc in accs:
        for file in os.listdir(os.path.join(data_dir, acc)):
            if file.endswith('.fna'):
                print(file)
                fp = os.path.join(data_dir, acc, file)
                if os.path.exists(fp):
                    new_fp = os.path.join(new_genome_dir, f'{acc}.fna')
                    print(new_fp)
                    shutil.move(fp,new_fp )
                else:
                    no_genome.append(acc)
    return no_genome

def move_gffs(data_dir, new_gff_dir = './gff_files'):
    """Moves proteomes from the nested ncbi data dirs to own dir and renames with accesion.faa"""
    
    os.makedirs(new_gff_dir, exist_ok=True)
    no_gff = []
    accs = [dir for dir in os.listdir(data_dir) if dir.startswith('GC')]
    for acc in accs:
        fp = os.path.join(data_dir, acc, 'genomic.gff')
        if os.path.exists(fp):
            new_fp = os.path.join(new_gff_dir, f'{acc}.gff')
            shutil.move(fp,new_fp )
        else:
            no_gff.append(acc)
    return no_gff