import subprocess
import os
import zipfile
import json
import pandas as pd
import tempfile
import shutil
from pathlib import Path


def build_gtdb_style_taxonomy(
    selected_assemblies_fp,
    taxonomy_summary_fp,
    taxonomy_table_fp,
    gtdb_mapping_fp,
):
    """Join NCBI assembly metadata to fixed-rank, GTDB-style taxonomy.

    Parameters
    ----------
    selected_assemblies_fp : str or pathlib.Path
        Project assembly table containing ``Assembly Accession``,
        ``Organism Taxonomic ID``, ``Organism Name``, and optionally
        ``taxonomic_group``.
    taxonomy_summary_fp : str or pathlib.Path
        ``taxonomy_summary.tsv`` produced by NCBI Datasets.
    taxonomy_table_fp : str or pathlib.Path
        Destination for the wide rank-column table.
    gtdb_mapping_fp : str or pathlib.Path
        Destination for the two-column accession-to-taxonomy mapping.

    Returns
    -------
    pandas.DataFrame
        The wide taxonomy table. Missing formal ranks remain empty in the TSV
        and appear as an empty prefixed rank (for example, ``c__``) in the
        GTDB-style lineage.
    """
    selected_assemblies_fp = Path(selected_assemblies_fp)
    taxonomy_summary_fp = Path(taxonomy_summary_fp)
    taxonomy_table_fp = Path(taxonomy_table_fp)
    gtdb_mapping_fp = Path(gtdb_mapping_fp)

    assemblies = pd.read_csv(
        selected_assemblies_fp,
        sep="\t",
        dtype={"Organism Taxonomic ID": "string"},
    )
    ncbi_taxonomy = pd.read_csv(
        taxonomy_summary_fp,
        sep="\t",
        dtype={"Taxid": "string", "Query": "string"},
        keep_default_na=False,
    )

    required_assembly_columns = {
        "Assembly Accession",
        "Organism Taxonomic ID",
        "Organism Name",
    }
    missing_columns = required_assembly_columns.difference(assemblies.columns)
    if missing_columns:
        raise ValueError(
            "Assembly table is missing columns: "
            + ", ".join(sorted(missing_columns))
        )

    rank_sources = {
        "domain": "Domain/Realm name",
        "phylum": "Phylum name",
        "class": "Class name",
        "order": "Order name",
        "family": "Family name",
        "genus": "Genus name",
        "species": "Species name",
    }
    missing_taxonomy_columns = set(rank_sources.values()).difference(
        ncbi_taxonomy.columns
    )
    if missing_taxonomy_columns:
        raise ValueError(
            "NCBI taxonomy table is missing columns: "
            + ", ".join(sorted(missing_taxonomy_columns))
        )

    taxonomy_columns = ["Taxid", *rank_sources.values()]
    if ncbi_taxonomy["Taxid"].duplicated().any():
        duplicated = ncbi_taxonomy.loc[
            ncbi_taxonomy["Taxid"].duplicated(keep=False), "Taxid"
        ].unique()
        raise ValueError(f"Duplicate NCBI taxonomy records: {', '.join(duplicated)}")

    output = assemblies.merge(
        ncbi_taxonomy[taxonomy_columns],
        left_on="Organism Taxonomic ID",
        right_on="Taxid",
        how="left",
        validate="many_to_one",
        indicator=True,
    )
    unmatched = output.loc[output["_merge"] != "both", "Organism Taxonomic ID"]
    if not unmatched.empty:
        raise ValueError(
            "No NCBI taxonomy record for taxids: "
            + ", ".join(unmatched.astype(str))
        )

    output = output.rename(
        columns={
            "Assembly Accession": "assembly_accession",
            "Organism Taxonomic ID": "taxid",
            "Organism Name": "organism_name",
            **{source: rank for rank, source in rank_sources.items()},
        }
    )
    output = output.drop(columns=["Taxid", "_merge"])
    output[list(rank_sources)] = output[list(rank_sources)].fillna("")

    prefixes = {
        "domain": "d__",
        "phylum": "p__",
        "class": "c__",
        "order": "o__",
        "family": "f__",
        "genus": "g__",
        "species": "s__",
    }
    output["gtdb_taxonomy"] = output.apply(
        lambda row: ";".join(
            f"{prefixes[rank]}{row[rank]}" for rank in prefixes
        ),
        axis=1,
    )

    leading_columns = ["assembly_accession", "taxid", "organism_name"]
    if "taxonomic_group" in output.columns:
        leading_columns.append("taxonomic_group")
    output = output[
        leading_columns + list(rank_sources) + ["gtdb_taxonomy"]
    ].sort_values("assembly_accession").reset_index(drop=True)

    taxonomy_table_fp.parent.mkdir(parents=True, exist_ok=True)
    gtdb_mapping_fp.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(taxonomy_table_fp, sep="\t", index=False)
    output[["assembly_accession", "gtdb_taxonomy"]].to_csv(
        gtdb_mapping_fp,
        sep="\t",
        index=False,
        header=False,
    )
    return output


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


# def make_summary_df_full(json_fp):
#     summary_d = parse_ncbi_summary(json_fp)
#     cols = {'accession':[], 'organism_name':[], 'strain':[],'tax_id':[], 'completeness':[], 'contamination':[], 'contig_n50':[], 'assembly_level':[] }
#     for acc, report in summary_d.items():
#         cols['accession'].append(acc)
#         t = report.get('checkm_info')
#         if t:
#             cols['completeness'].append(t.get('completeness'))
#             cols['contamination'].append(t.get('contamination'))
#         else:
#             cols['completeness'].append(None)
#             cols['contamination'].append(None)
            
#         t = report.get('assembly_stats')
#         if t:
#             cols['contig_n50'].append(t.get('contig_n50'))
#         else:
#             cols['contig_n50'].append(None)
#         t = report['assembly_info']
#         if t:
#             cols['assembly_level'].append(t['assembly_level'])
#         else:
#             cols['assembly_level'].append(None)
            
#         cols['organism_name'].append(report['organism'].get('organism_name'))
#         cols['tax_id'].append(report['organism'].get('tax_id'))
#         t = report['organism'].get('infraspecific_names')
#         if t:
#             cols['strain'].append(t.get('strain', ''))
#         else:
#             cols['strain'].append('')
#     summary_df = pd.DataFrame(cols).set_index('accession')
#     summary_df['full_name'] = summary_df['organism_name'].str.cat(summary_df['strain'], sep='__')
#     return summary_df


def make_summary_df_full(json_fp):
    summary_d = parse_ncbi_summary(json_fp)

    rows = []

    for acc, report in summary_d.items():
        checkm_info = report.get("checkm_info") or {}
        assembly_stats = report.get("assembly_stats") or {}
        assembly_info = report.get("assembly_info") or {}
        organism = report.get("organism") or {}
        infraspecific_names = organism.get("infraspecific_names") or {}

        organism_name = organism.get("organism_name")
        strain = infraspecific_names.get("strain", "")

        rows.append(
            {
                "accession": acc,
                "organism_name": organism_name,
                "strain": strain,
                "tax_id": organism.get("tax_id"),
                "completeness": checkm_info.get("completeness"),
                "contamination": checkm_info.get("contamination"),
                "contig_n50": assembly_stats.get("contig_n50"),
                "assembly_level": assembly_info.get("assembly_level"),
            }
        )

    summary_df = pd.DataFrame(rows).set_index("accession")

    summary_df["full_name"] = summary_df.apply(
        lambda row: (
            f"{row['organism_name']}__{row['strain']}"
            if pd.notna(row["organism_name"]) and row["strain"]
            else row["organism_name"]
        ),
        axis=1,
    )

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
