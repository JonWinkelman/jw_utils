import subprocess
import os
import zipfile

def download_genomes_from_accfile(accessions_fp, files_to_include, dataset_fp):
    """Download ncbi datasets genomes from input file containing accession on each line

    
    
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
    files_to_include = files_to_include.replace(' ','')
    download_dehydrated_ncbi_dataset(accessions_fp, files_to_include, dataset_fp)
    extract_path = dataset_fp.replace('.zip', '')
    unzip_ncbi_dir(dataset_fp, extract_path)
    rehydrate_ncbi_dir(extract_path)
    

def download_dehydrated_ncbi_dataset(accessions_fp, files_to_include, dataset_fp):
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

