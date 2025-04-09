from concurrent.futures import ThreadPoolExecutor
import subprocess
import logging
from pathlib import Path
import os
import pandas as pd
from jw_utils import parse_fasta as pfa
from jw_utils import hmmer_utils as hu
from Bio import SeqIO
import os
from collections import defaultdict
from jw_utils import alignment_utils2 as au


def concatenate_hmm_files(directory, output_file, suffix='.HMM'):
    """
    Concatenates all files ending with '.HMM' in the given directory into one output file.

    Args:
        directory (str): Path to the directory containing .HMM files.
        output_file (str): Path to the output file where concatenated content will be saved.
    """
    with open(output_file, 'w') as outfile:
        
        for filename in os.listdir(directory):
            if filename.endswith(suffix):
                file_path = os.path.join(directory, filename)
                with open(file_path, 'r') as infile:
                    outfile.write(infile.read())
                    outfile.write('\n')  # Add a newline between files

    print(f"All .HMM files from '{directory}' concatenated into '{output_file}'.")




# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

def process_hmm_profiles(
    input_dir='bac120_hmm_profiles', 
    output_file='./results/bac120hmm_db/concat_bac120_profiles.hmm',
    run_concatenation=True
):
    """
    Processes HMM profiles by optionally concatenating files, ensuring the necessary 
    database directory exists, and running hmmpress.

    Parameters:
        input_dir (str): Directory containing the HMM profiles to concatenate.
        output_file (str): Path to save the concatenated HMM profile.
        run_concatenation (bool): Whether to run concatenate_hmm_files() before hmmpress.
    """
    
    # Use pathlib for better path handling
    output_path = Path(output_file).resolve()  # Get absolute path
    db_dir = output_path.parent  # Extract the directory where the output file is stored

    # Create database directory if it doesn't exist
    db_dir.mkdir(parents=True, exist_ok=True)
    
    if run_concatenation:
        logging.info(f"Concatenating HMM profiles from {input_dir} to {output_path}")
        concatenate_hmm_files(input_dir, str(output_path))  # Ensure correct type for function

    try:
        # Run hmmpress without changing directories
        hmmpress_cmd = ["hmmpress", str(output_path)]
        logging.info(f"Running command: {' '.join(hmmpress_cmd)}")
        
        result = subprocess.run(hmmpress_cmd, check=True, capture_output=True, text=True)
        logging.info(f"hmmpress output:\n{result.stdout}")

    except subprocess.CalledProcessError as e:
        logging.error(f"Error running hmmpress: {e.stderr if e.stderr else e}")



def run_hmmsearch(proteome_name,proteome_fp, hmm_db, output_dir):
    """Function to run hmmsearch with multiple hmms in an indexed hmm db on a single proteome"""
    output_file = os.path.join(output_dir, f'{proteome_name}_best_hits.txt')
    tblout_file = os.path.join(output_dir, f'{proteome_name}_table.txt')
    cmd = [
        'hmmsearch', '--tblout', tblout_file, '-E', '1e-5', 
        hmm_db, proteome_fp 
    ]
    with open(output_file, 'w') as outfile:
        subprocess.run(cmd, stdout=outfile)




def ThreadPool_hmm_search(hmm_db, output_dir, prot_name_fp_d, threads = 6):
    # Use ThreadPoolExecutor to run searches in parallel
    with ThreadPoolExecutor(max_workers=threads) as executor:  # Adjust max_workers based on CPU cores
        futures = [executor.submit(run_hmmsearch, 
                                   proteome_name,
                                   proteome_fp, 
                                   hmm_db, 
                                   output_dir) for  proteome_name, proteome_fp in prot_name_fp_d.items()]




def parse_tblout(tblout_file):
    hits = []
    with open(tblout_file) as f:
        for line in f:
            if not line.startswith('#'):  # Skip comments
                cols = line.split()
                target = cols[0]  # target sequence
                hmm_acc = cols[3]
                e_value = float(cols[4])  # e-value of the hit
                hits.append((target, hmm_acc, e_value))
    return hits



def aggregate_best_hits(prot_name_fp_d, hmm_output_dir):
    """
    Aggregates the best HMMER hits across all proteomes and saves them as CSV files.

    Parameters:
        prot_name_fp_d (dict): Dictionary mapping proteome names to file paths.
        hmm_output_dir (str): Directory containing the HMMER output files.
    """
    hits_dict = {}

    # Parse all HMMER output tables
    for name in prot_name_fp_d.keys():
        tblout_file = Path(hmm_output_dir) / f"{name}_table.txt"
        
        if tblout_file.exists():
            hits_dict[name] = parse_tblout(tblout_file)
        else:
            print(f"Warning: {tblout_file} not found. Skipping {name}.")

    # Create output directory for best hits
    best_hits_dir = Path(hmm_output_dir) / "best_hits"
    best_hits_dir.mkdir(exist_ok=True)

    # Set Pandas float display format
    pd.set_option("display.float_format", "{:.10e}".format)

    dfs = {}

    # Process hits for each proteome
    for name, hits in hits_dict.items():
        if not hits:
            print(f"Warning: No hits found for {name}. Skipping.")
            continue

        # Convert list of hits to DataFrame
        df = pd.DataFrame(hits, columns=["target_name", "hmm_acc", "e_value"])

        # Select best hit per HMM (lowest e-value)
        best_hits = df.sort_values("e_value").groupby("hmm_acc", as_index=False).first()

        # Store and save results
        dfs[name] = best_hits.set_index("target_name")
        best_hits.to_csv(best_hits_dir / f"{name}_best_hits.csv")

    return dfs  # Returning the dataframes in case they're needed later



def extract_bac120_proteins(prot_name_fp_d, dfs, output_dir='./fastTree/bac120_proteins'):
    """
    Extracts ~120 proteins from each proteome and writes them to FASTA files.

    Parameters:
        prot_name_fp_d (dict): Dictionary mapping proteome names to file paths.
        dfs (dict): Dictionary of best-hit dataframes, indexed by proteome name.
        output_dir (str): Directory where the extracted protein FASTA files will be stored.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)  # Ensure the output directory exists

    perc = 0.05  # Progress tracking
    total = len(prot_name_fp_d)

    for i, (name, proteome_fp) in enumerate(prot_name_fp_d.items(), start=1):
        # Print progress every 5%
        progress = i / total
        if progress > perc:
            print(f'{round(progress * 100, 2)}% finished.')
            perc += 0.05  # Increment progress threshold

        # Parse the proteome FASTA
        seq_d = pfa.get_seq_dict(proteome_fp)

        # Extract sequences for HMM hits
        hmm_seq_d = {str(row["hmm_acc"]): seq_d.get(protein_id, None) 
                     for protein_id, row in dfs[name].iterrows()}

        # Remove None values (proteins not found in the sequence dictionary)
        hmm_seq_d = {k: v for k, v in hmm_seq_d.items() if v is not None}

        if not hmm_seq_d:
            print(f"Warning: No matching proteins found for {name}. Skipping.")
            continue

        # Write to FASTA file
        fasta_path = output_dir / f"{name}_bac120hits.faa"
        pfa.write_to_fasta(hmm_seq_d, fasta_path)



def create_hmm_seq_files(prot_name_fp_d, hmm_bac120_ids, bac120_proteins, output_dir):
    """
    Extracts protein sequences for each HMM from each proteome and writes FASTA files.

    Parameters:
        prot_name_fp_d (dict): Dictionary mapping genome names to file paths.
        hmm_bac120_ids (list): List of HMM profile identifiers.
        bac120_proteins (str): Directory containing extracted protein FASTA files.
        output_dir (str): Directory where the aligned FASTA files will be stored.

    Returns:
        pres_abs_df (pd.DataFrame): Presence-absence matrix of HMM sequences across genomes.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)  # Ensure output directory exists

    names = list(prot_name_fp_d.keys())
    missing = {name: [] for name in names}

    for hmm in hmm_bac120_ids:
        seq_d = {}

        for name in names:
            fp = Path(bac120_proteins) / f"{name}_bac120hits.faa"
            seq = pfa.get_seq_dict(fp).get(hmm)

            if seq:
                seq_d[name] = seq
                missing[name].append(True)
            else:
                missing[name].append(False)

        fasta_output = output_dir / f"{hmm}.faa"
        pfa.write_to_fasta(seq_d, fasta_output)

    # Create presence-absence DataFrame
    pres_abs_df = pd.DataFrame.from_dict(missing, orient='index', columns=hmm_bac120_ids)
    
    return pres_abs_df

def run_hmm_alignments(hmm_bac120_ids, hmm_db, input_dir):
    """
    Runs HMMER alignment for each extracted HMM sequence file.

    Parameters:
        hmm_bac120_ids (list): List of HMM profile identifiers.
        hmm_db (str): Path to the HMM database.
        input_dir (str): Directory containing input FASTA files.

    Returns:
        list: File paths of generated Stockholm alignment files.
    """
    input_dir = Path(input_dir)
    sto_fps = []

    for hmm in hmm_bac120_ids:
        print(f'Aligning {hmm}...')

        seqs_fp = input_dir / f"{hmm}.faa"
        out_fp = seqs_fp.with_suffix('.sto')

        hu.run_hmm_alignment(hmm_db, str(seqs_fp), hmm, str(out_fp), output_format='Pfam', trim=True)
        sto_fps.append(out_fp)

    return sto_fps



def remove_pp_lines(hmm_alignment_fp):
    """
    Reads an HMM alignment file and removes lines starting with '#' or '/'.
    Reformats sequences by replacing '.' with '-' and ensures uppercase sequence output.

    Parameters:
        hmm_alignment_fp (str or Path): Path to the input HMM alignment file.

    Returns:
        dict: Reformatted sequence dictionary {sequence_id: sequence}.
    """
    new_aln_d = {}

    try:
        with open(hmm_alignment_fp, 'r') as f:
            for line in f:
                if line.startswith(('#', '/')) or len(line.strip()) <= 5:
                    continue
                line_lst = line.strip().split()
                if line_lst:
                    seq_id = line_lst[0]
                    seq = line_lst[-1].replace('.', '-').upper()
                    new_aln_d[seq_id] = seq

    except FileNotFoundError:
        print(f"Error: File not found - {hmm_alignment_fp}")
    
    return new_aln_d

def write_reformatted_aln_file(sto_alignment_fp, homogeneous_thresh=0.98, gap_threshold=0.8):
    """
    Reformats an HMM alignment by removing unnecessary columns and gaps.

    Parameters:
        sto_alignment_fp (str or Path): Path to the input Stockholm (.sto) alignment file.
        homogeneous_thresh (float): Threshold for removing homogeneous columns.
        gap_threshold (float): Threshold for removing gappy columns.

    Returns:
        dict: Trimmed sequence dictionary after formatting.
    """
    sto_alignment_fp = Path(sto_alignment_fp)
    reformatted_aligned_d = remove_pp_lines(sto_alignment_fp)

    if not reformatted_aligned_d:
        print(f"Warning: No valid sequences found in {sto_alignment_fp}. Skipping.")
        return {}

    # Apply column filtering
    trimmed_seq_d = au.remove_homogeneous_cols(reformatted_aligned_d, threshold=homogeneous_thresh)
    trimmed_seq_d = au.remove_gappy_columns(trimmed_seq_d, threshold=gap_threshold)

    # Write reformatted alignment to new file
    reformatted_align_fp = sto_alignment_fp.with_suffix('.sto.simple')
    pfa.write_to_fasta(trimmed_seq_d, reformatted_align_fp)

    return trimmed_seq_d

def process_hmm_alignments(sto_fps, homogeneous_thresh=0.98, gap_threshold=0.8):
    """
    Processes multiple HMM alignments by reformatting them and storing the results.

    Parameters:
        sto_fps (list): List of file paths to Stockholm (.sto) alignments.
        homogeneous_thresh (float): Threshold for removing homogeneous columns.
        gap_threshold (float): Threshold for removing gappy columns.

    Returns:
        list: File paths of reformatted alignment files.
    """
    simple_aln_fps = []

    for fp in sto_fps:
        trimmed_seq_d = write_reformatted_aln_file(fp, homogeneous_thresh, gap_threshold)

        if trimmed_seq_d:  # Ensure the file was processed
            
            simple_aln_fps.append(str(fp).replace('sto', 'sto.simple'))

    return simple_aln_fps



# def get_curated_accessions(names, remove):
#     """
#     Returns a set of curated genome accessions after removing unwanted ones.

#     Parameters:
#         names (set or list): A set of all genome names.
#         remove (set or list): A set of genome names to exclude.

#     Returns:
#         set: Curated accessions after exclusion.
#     """
#     return set(names).difference(remove)

def concatenate_sequences(fasta_dir, curated_accs, max_files=3000):
    """
    Concatenates sequences from multiple aligned FASTA files for each genome accession.

    Parameters:
        fasta_dir (str or Path): Directory containing aligned FASTA files.
        curated_accs (set): Set of curated genome accessions.
        max_files (int, optional): Maximum number of FASTA files to process (default: 30).

    Returns:
        dict: Dictionary mapping genome accessions to concatenated sequences.
    """
    fasta_dir = Path(fasta_dir)
    concatenated_sequences = defaultdict(str)

    # Get a sorted list of FASTA files
    fasta_files = sorted([f for f in fasta_dir.glob("*.simple")])[:max_files]

    for fasta_file in fasta_files:
        seq_d = pfa.get_seq_dict(fasta_file)

        if not seq_d:
            print(f"Warning: No valid sequences in {fasta_file}. Skipping.")
            continue

        aln_len = len(next(iter(seq_d.values())))  # Get alignment length from the first sequence

        for acc in curated_accs:
            seq = seq_d.get(acc, '-' * aln_len)
            concatenated_sequences[acc] += seq

    return concatenated_sequences

def write_concatenated_fasta(concatenated_sequences, output_fasta):
    """
    Writes concatenated sequences to a new FASTA file.

    Parameters:
        concatenated_sequences (dict): Dictionary of genome accessions and concatenated sequences.
        output_fasta (str or Path): Path to the output FASTA file.
    """
    output_fasta = Path(output_fasta)

    with output_fasta.open('w') as output_handle:
        for genome_accession, concatenated_seq in concatenated_sequences.items():
            output_handle.write(f'>{genome_accession}\n')
            output_handle.write(f'{concatenated_seq}\n')

def process_concatenated_alignment(names, remove, fasta_dir, output_fasta, max_files=3000):
    """
    Full pipeline to concatenate sequences from aligned FASTA files and save to a FASTA file.

    Parameters:
        names (set or list): Set of all genome names.
        remove (set or list): Set of genome names to exclude.
        fasta_dir (str or Path): Directory containing aligned FASTA files.
        output_fasta (str or Path): Path to the output concatenated FASTA file.
        max_files (int, optional): Maximum number of FASTA files to process (default: 30).
    """
    curated_accs = set(names).difference(remove)
    concatenated_sequences = concatenate_sequences(fasta_dir, curated_accs, max_files)
    write_concatenated_fasta(concatenated_sequences, output_fasta)



def run_fasttree(input_fasta, output_tree, **kwargs):
    """
    Runs FastTree on the given alignment file with user-specified flags.

    Parameters:
        input_fasta (str or Path): Path to the input concatenated alignment file.
        output_tree (str or Path): Path to save the output Newick tree file.
        kwargs (dict): Optional FastTree flags. Pass them as keyword arguments, e.g.:
                       run_fasttree("input.fasta", "output.nwk", quiet=True, log="logfile.txt")

    Returns:
        str: The standard output from FastTree execution.
    """
    input_fasta = Path(input_fasta).resolve()
    output_tree = Path(output_tree).resolve()

    # Base command
    cmd = ["fasttree"]

    # Add user-defined flags
    for flag, value in kwargs.items():
        if isinstance(value, bool):  # Flags that are just toggles
            if value:  
                cmd.append(f"-{flag}")
        else:  # Flags that require a value (e.g., -log logfile.txt)
            cmd.append(f"-{flag}")
            cmd.append(str(value))

    # Add input and output
    cmd.append(str(input_fasta))

    try:
        # Run FastTree and capture output
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)

        # Write output to the specified tree file
        with open(output_tree, "w") as tree_file:
            tree_file.write(result.stdout)

        print(f"FastTree successfully generated: {output_tree}")
        return result.stdout

    except subprocess.CalledProcessError as e:
        print(f"Error running FastTree: {e.stderr}")
        return None



def make_bac120_tree(hmm_profile_dir, proteome_dir_fp, proteome_suffix='.faa',
                     output_tree='results/bac120_tree/bac120_tree.nwk'):
    """
    Takes a list of HMM profiles, extracts the best hit from each proteome aligns, concatenates, and builds tree.
    
    hmm_profile_dir (path): This is the path to a directory containing individual text file of each hmm profile, e.g. PF00380.23.HMM
    proteome_dir_fp (path): This is the path to a directory containing individual fasta proteeomes, e.g. PF00380.23.HMM
    
    """
    
    # filepaths
    hmm_profile_dir = Path(hmm_profile_dir)
    proteome_dir_fp = Path(proteome_dir_fp)
    output_tree = Path(output_tree)
    parent_dir = output_tree.parent
    log = parent_dir / 'bac120tree_logfile.txt'
    concat_fasta = parent_dir / 'concatenated_alignment.fasta'
    hmm_db  = parent_dir / 'bac120hmm_db/concat_bac120_profiles.hmm'
    
    # dirs to make
    fasta_dir = parent_dir / 'hmm_alignments'
    hmm_output_dir = parent_dir / 'hmm_search_output/'
    bac_120_extacted_proteins_dir = parent_dir / 'bac_120_extacted_proteins'
    hmm_alignments_dir = parent_dir / 'hmm_alignments'

    if parent_dir.exists():
        raise Exception(f'"{parent_dir}" already exists, rename or delete if you want to run this function')
    
    for path in [parent_dir, fasta_dir, hmm_output_dir, hmm_db.parent, bac_120_extacted_proteins_dir, hmm_alignments_dir]:
        path.mkdir(exist_ok=True, parents=True)
    
    prot_name_fp_d = {f.strip(proteome_suffix).strip('.'):os.path.join(proteome_dir_fp, f ) for f in os.listdir(proteome_dir_fp) if f.endswith(proteome_suffix)}
    hmm_bac120_ids = [f.stem for f in hmm_profile_dir.glob('*.HMM')]
    if len(hmm_bac120_ids) < 1:
        raise Exception(f'are you sure files in {hmm_profile_dir} are present and formatted correctly? e.g. do the end with ".HMM"?')
    
    print("concatenating HMM profiles into one file and then pressing it into binary...")

    process_hmm_profiles(hmm_profile_dir, output_file= hmm_db, 
                         run_concatenation=True)

    
    print("searching each proteome with each HMM using hmmsearch...")
    ThreadPool_hmm_search(str(hmm_db), str(hmm_output_dir), prot_name_fp_d, threads=6)
    
    print("aggregating best hit from each proteome into dataframe...")
    dfs = aggregate_best_hits(prot_name_fp_d, hmm_output_dir)

    print("Extracting protein sequences for each HMM from each proteome and writing FASTA files...")
    
    extract_bac120_proteins(prot_name_fp_d, dfs, output_dir=bac_120_extacted_proteins_dir)
    
    

    
    pres_abs_df = create_hmm_seq_files(prot_name_fp_d, hmm_bac120_ids, bac_120_extacted_proteins_dir, hmm_alignments_dir)

    print(f'for each hmm, aligning all extracted proteins with that hmm')
    sto_fps =  run_hmm_alignments(hmm_bac120_ids, hmm_db, hmm_alignments_dir)
    simple_aln_fps = process_hmm_alignments(sto_fps)

    print(f'For each proteome, concatenating all hmm-aligned proteins into one long protein and saving as {concat_fasta}')
    process_concatenated_alignment(names=list(prot_name_fp_d.keys()), remove=[], fasta_dir=fasta_dir, output_fasta=concat_fasta)
    
    print(f'generating maximum likelihood tree using fasttree from concatenated protein sequences...')
    run_fasttree(input_fasta=concat_fasta, output_tree=output_tree, log=log)
    return pres_abs_df