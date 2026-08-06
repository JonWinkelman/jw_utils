import re
import pandas as pd
import os
import subprocess
import shlex
import tempfile
from collections.abc import Mapping
from jw_utils import parse_fasta as pfa
import requests
from time import sleep
from pathlib import Path
from tqdm import tqdm



from pathlib import Path
import subprocess

def run_hmmbuild(
    hmmfile,
    msafile,
    hmmbuild_exe="hmmbuild",
    **kwargs,
):
    """
    Run HMMER hmmbuild.

    Parameters
    ----------
    hmmfile : str or Path
        Output HMM file.
    msafile : str or Path
        Input multiple sequence alignment.
    hmmbuild_exe : str, default="hmmbuild"
        Path to hmmbuild executable.
    **kwargs
        Additional hmmbuild options.

        Boolean values become flags::

            amino=True
            fast=True
            hand=True

        Other values become option/value pairs::

            cpu=8
            n="AsnC_HOG1"
            symfrac=0.7
            informat="afa"
            seed=1

        Underscores are converted to hyphens::

            maxinsertlen=20
            mxfile="BLOSUM62.txt"

    Returns
    -------
    subprocess.CompletedProcess
    """
    cmd = [hmmbuild_exe]

    for key, value in kwargs.items():
        opt = "--" + key.replace("_", "-")

        if value is None or value is False:
            continue
        elif value is True:
            cmd.append(opt)
        else:
            cmd.extend([opt, str(value)])

    cmd.extend([str(hmmfile), str(msafile)])

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=True,
    )

    return result


def run_hmmscan(
    sequences,
    hmm_database,
    tblout_fp,
    domtblout_fp=None,
    raw_output_fp=None,
    hmmscan_exe="hmmscan",
    **kwargs,
):
    """Scan one or more protein sequences against an HMM profile database.

    Parameters
    ----------
    sequences : str, pathlib.Path, or mapping
        Input protein FASTA path, or a ``{sequence_name: sequence}`` mapping.
    hmm_database : str or pathlib.Path
        HMM database file. The database must have been indexed with
        ``hmmpress``; pass the base ``.hmm`` path, not an ``.h3*`` file.
    tblout_fp : str or pathlib.Path
        Destination for the per-sequence hit table (``--tblout``).
    domtblout_fp : str or pathlib.Path, optional
        Destination for the per-domain hit table (``--domtblout``). If not
        supplied, ``<tblout stem>.domtblout`` is used.
    raw_output_fp : str or pathlib.Path, optional
        Destination for HMMER's human-readable report. If not supplied,
        ``<tblout stem>.txt`` is used.
    hmmscan_exe : str or pathlib.Path, default="hmmscan"
        HMMER ``hmmscan`` executable name or path.
    **kwargs
        Additional HMMER options. Boolean values become flags; other values
        become option/value pairs. For example, ``cpu=4``, ``E=1e-5``, and
        ``noali=True`` become ``--cpu 4 -E 1e-05 --noali``.

    Returns
    -------
    subprocess.CompletedProcess
        The completed ``hmmscan`` process. Parsed tables can be obtained with
        :func:`parse_hmmsearch_output` and :func:`parse_domtblout_to_df`.

    Raises
    ------
    FileNotFoundError
        If the input FASTA, HMM database, or pressed database files are absent.
    TypeError
        If ``sequences`` is neither a FASTA path nor a mapping.
    subprocess.CalledProcessError
        If ``hmmscan`` exits unsuccessfully.
    """
    hmm_database = Path(hmm_database)
    tblout_fp = Path(tblout_fp)
    domtblout_fp = (
        Path(domtblout_fp)
        if domtblout_fp is not None
        else tblout_fp.with_suffix(".domtblout")
    )
    raw_output_fp = (
        Path(raw_output_fp)
        if raw_output_fp is not None
        else tblout_fp.with_suffix(".txt")
    )

    if not hmm_database.is_file():
        raise FileNotFoundError(f"HMM database not found: {hmm_database}")

    missing_indexes = [
        Path(f"{hmm_database}.{suffix}")
        for suffix in ("h3f", "h3i", "h3m", "h3p")
        if not Path(f"{hmm_database}.{suffix}").is_file()
    ]
    if missing_indexes:
        missing = ", ".join(path.name for path in missing_indexes)
        raise FileNotFoundError(
            f"HMM database is not pressed; missing {missing}. "
            f"Run: hmmpress {hmm_database}"
        )

    for output_path in (tblout_fp, domtblout_fp, raw_output_fp):
        output_path.parent.mkdir(parents=True, exist_ok=True)

    temporary_fasta = None
    if isinstance(sequences, Mapping):
        if not sequences:
            raise ValueError("sequences mapping cannot be empty")
        temporary_fasta = tempfile.NamedTemporaryFile(
            mode="w", suffix=".faa", delete=False
        )
        temporary_fasta.close()
        sequence_fp = Path(temporary_fasta.name)
        pfa.write_to_fasta(dict(sequences), sequence_fp)
    elif isinstance(sequences, (str, os.PathLike)):
        sequence_fp = Path(sequences)
        if not sequence_fp.is_file():
            raise FileNotFoundError(f"Input FASTA not found: {sequence_fp}")
    else:
        raise TypeError("sequences must be a FASTA path or a mapping of names to sequences")

    cmd = [
        str(hmmscan_exe),
        "--tblout", str(tblout_fp),
        "--domtblout", str(domtblout_fp),
        "-o", str(raw_output_fp),
    ]
    for key, value in kwargs.items():
        prefix = "-" if len(key) == 1 else "--"
        # HMMER spells its curated-cutoff flags with underscores; most other
        # long options use the key exactly as documented (for example domE).
        option_name = key if key in {"cut_ga", "cut_nc", "cut_tc"} else key.replace("_", "-")
        option = prefix + option_name
        if value is None or value is False:
            continue
        if value is True:
            cmd.append(option)
        else:
            cmd.extend([option, str(value)])
    cmd.extend([str(hmm_database), str(sequence_fp)])

    try:
        return subprocess.run(cmd, capture_output=True, text=True, check=True)
    finally:
        if temporary_fasta is not None:
            try:
                os.unlink(temporary_fasta.name)
            except OSError:
                pass



def run_simple_searches(hmm_profile_path, results_dir, proteome_dir, eval_thresh=1e-6):
    """
    Run HMMER searches against all proteomes in a directory and return
    filtered hits with their amino-acid sequences.

    Parameters
    ----------
    hmm_profile_path : str or pathlib.Path
        Path to the HMM profile file used as the query.

    results_dir : str or pathlib.Path
        Directory where raw hmmsearch output files will be written.

    proteome_dir : str or pathlib.Path
        Directory containing proteome FASTA files with the extension `.faa`.

    eval_thresh : float, default=1e-6
        Maximum E-value threshold for retaining HMM hits.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing filtered HMM hits from all proteomes.
        The index is `full_id`, constructed as:

            genome_acc__protein_id

        Columns include:

        - protein_id
        - HMM_name
        - E-value
        - description_of_target
        - AA_seq
        - genome_acc

    Notes
    -----
    This function:

    1. Finds all `.faa` files in `proteome_dir`.
    2. Runs a simple HMM search for each proteome using `run_simple_search`.
    3. Parses each HMM output using `parse_hmmsearch_output`.
    4. Filters hits by E-value.
    5. Adds amino-acid sequences from the corresponding proteome FASTA.
    6. Concatenates all per-genome results into one DataFrame.
    """

    hmm_profile_path = Path(hmm_profile_path)
    results_dir = Path(results_dir)
    proteome_dir = Path(proteome_dir)
    results_dir.mkdir(exist_ok=True, parents=True)
    proteome_fps = sorted(proteome_dir.glob("*.faa"))
    dfs = []
    for proteome_fp in tqdm(proteome_fps, desc="Running HMM searches"):
        hmm_out_fp = results_dir / f"{proteome_fp.stem}.txt"
        run_simple_search(
            str(hmm_profile_path),
            str(hmm_out_fp),
            proteome_fp=proteome_fp
        )
        df = parse_hmmsearch_output(hmm_out_fp)
        df = df[df["E-value"] < eval_thresh].copy()
        if df.empty:
            continue
        seq_df = (
            pd.Series(pfa.get_seq_dict(proteome_fp), name="AA_seq")
            .reset_index()
            .rename(columns={"index": "target_name"})
        )
        df = (
            pd.merge(df, seq_df, on="target_name")
            [[
                "target_name",
                "query_name",
                "E-value",
                "description_of_target",
                "AA_seq",
            ]]
        )
        df["genome_acc"] = proteome_fp.stem
        df["full_id"] = df["genome_acc"].str.cat(df["target_name"], sep="__")
        dfs.append(df)
    if not dfs:
        return pd.DataFrame(
            columns=[
                "protein_id",
                "HMM_name",
                "E-value",
                "description_of_target",
                "AA_seq",
                "genome_acc",
            ]
        )
    result_df = (
        pd.concat(dfs, ignore_index=True)
        .set_index("full_id")
        .rename(columns={
            "target_name": "protein_id",
            "query_name": "HMM_name",
        })
    )
    return result_df



def download_pfam_hmms(pfam_ids, out_dir, delay=0.2, overwrite=False):
    """
    Download Pfam HMMs from InterPro API.

    pfam_ids : list of strings (e.g., ['PF00210', 'PF00582'])
    out_dir  : directory to save .hmm files
    delay    : seconds between requests (avoid rate limiting)
    overwrite: overwrite existing files
    """
    os.makedirs(out_dir, exist_ok=True)

    base_url = "https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/{pfam}?annotation=hmm"

    for pfam in pfam_ids:
        out_path = os.path.join(out_dir, f"{pfam}.hmm.gz")

        if os.path.exists(out_path) and not overwrite:
            print(f"[SKIP] {pfam} already exists")
            continue

        url = base_url.format(pfam=pfam)

        try:
            print(f"[DOWNLOAD] {pfam}")
            r = requests.get(url)

            if r.status_code != 200:
                print(f"[ERROR] {pfam}: status {r.status_code}")
                continue

            with open(out_path, "wb") as f:
                f.write(r.content)

            sleep(delay)

        except Exception as e:
            print(f"[FAIL] {pfam}: {e}")


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


def parse_hmmscan_output(filepath):
    """Parse an HMMER ``hmmscan --tblout`` file into a typed DataFrame.

    Parameters
    ----------
    filepath : str or pathlib.Path
        Path to the sequence-level table written by ``hmmscan --tblout``.

    Returns
    -------
    pandas.DataFrame
        One row per profile hit. Column names reflect ``hmmscan`` semantics:
        the target is an HMM profile and the query is an input sequence.
        E-values, scores, and biases are floats; count columns are nullable
        integers. An output containing no hits returns an empty DataFrame with
        the same columns.
    """
    filepath = Path(filepath)
    if not filepath.is_file():
        raise FileNotFoundError(f"hmmscan table not found: {filepath}")

    columns = [
        "profile_name",
        "profile_accession",
        "sequence_name",
        "sequence_accession",
        "full_sequence_evalue",
        "full_sequence_score",
        "full_sequence_bias",
        "best_domain_evalue",
        "best_domain_score",
        "best_domain_bias",
        "expected_domains",
        "regions",
        "clusters",
        "overlaps",
        "envelopes",
        "domains",
        "reported_domains",
        "included_domains",
        "profile_description",
    ]
    records = []
    with filepath.open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split(maxsplit=len(columns) - 1)
            fields.extend([""] * (len(columns) - len(fields)))
            records.append(fields)

    dataframe = pd.DataFrame(records, columns=columns)
    float_columns = [
        "full_sequence_evalue",
        "full_sequence_score",
        "full_sequence_bias",
        "best_domain_evalue",
        "best_domain_score",
        "best_domain_bias",
        "expected_domains",
    ]
    integer_columns = [
        "regions",
        "clusters",
        "overlaps",
        "envelopes",
        "domains",
        "reported_domains",
        "included_domains",
    ]
    for column in float_columns:
        dataframe[column] = pd.to_numeric(dataframe[column], errors="coerce")
    for column in integer_columns:
        dataframe[column] = pd.to_numeric(
            dataframe[column], errors="coerce"
        ).astype("Int64")

    return dataframe


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
    float_columns = [
        "full_seq_Evalue", "full_seq_score", "full_seq_bias",
        "domain_cEvalue", "domain_iEvalue", "domain_score", "domain_bias",
        "acc",
    ]
    integer_columns = [
        "tlen", "qlen", "domain_number", "domains_in_target", "hmm_from",
        "hmm_to", "ali_from", "ali_to", "env_from", "env_to",
    ]
    for column in float_columns:
        df[column] = pd.to_numeric(df[column], errors="coerce")
    for column in integer_columns:
        df[column] = pd.to_numeric(df[column], errors="coerce").astype("Int64")
    return df


def extract_fasta_records(sequence_ids, fasta_fp, output_fp=None, strict=True):
    """Extract named protein sequences from a FASTA file.

    Parameters
    ----------
    sequence_ids : iterable of str
        Exact first-token FASTA identifiers to retrieve.
    fasta_fp : str or pathlib.Path
        Source FASTA file.
    output_fp : str or pathlib.Path, optional
        If supplied, write the extracted sequences to this FASTA path.
    strict : bool, default True
        Raise an error if any requested identifier is absent.

    Returns
    -------
    dict
        ``{sequence_id: amino_acid_sequence}`` in requested-ID order.
    """
    requested = list(dict.fromkeys(map(str, sequence_ids)))
    source = pfa.get_seq_dict(fasta_fp)
    missing = [sequence_id for sequence_id in requested if sequence_id not in source]
    if missing and strict:
        preview = ", ".join(missing[:5])
        raise KeyError(
            f"{len(missing)} requested sequence(s) absent from {fasta_fp}: {preview}"
        )
    extracted = {
        sequence_id: source[sequence_id]
        for sequence_id in requested
        if sequence_id in source
    }
    if output_fp is not None:
        output_fp = Path(output_fp)
        output_fp.parent.mkdir(parents=True, exist_ok=True)
        pfa.write_to_fasta(extracted, output_fp, line_size=80)
    return extracted


def summarize_pfam_architectures(domain_hits_df, sequence_ids=()):
    """Summarize ordered Pfam domain hits for each scanned protein.

    ``parse_domtblout_to_df`` uses hmmscan terminology: ``query_name`` is the
    protein and ``target_name`` is the Pfam profile. Coordinates are ordered
    along the protein sequence. Proteins without a reported domain can be
    retained by supplying ``sequence_ids``.
    """
    required = {"query_name", "target_name", "target_accession", "ali_from", "ali_to"}
    missing = required.difference(domain_hits_df.columns)
    if missing:
        raise ValueError(f"Domain table is missing columns: {sorted(missing)}")

    rows = []
    grouped = {name: group for name, group in domain_hits_df.groupby("query_name")}
    identifiers = list(dict.fromkeys([*map(str, sequence_ids), *map(str, grouped)]))
    for sequence_id in identifiers:
        group = grouped.get(sequence_id)
        if group is None or group.empty:
            rows.append({
                "sequence_id": sequence_id,
                "domain_count": 0,
                "pfam_architecture": "no_Pfam_domain",
                "pfam_accession_architecture": "",
            })
            continue
        ordered = group.sort_values(["ali_from", "ali_to", "domain_iEvalue"])
        rows.append({
            "sequence_id": sequence_id,
            "domain_count": len(ordered),
            "pfam_architecture": " | ".join(ordered["target_name"].astype(str)),
            "pfam_accession_architecture": " | ".join(
                ordered["target_accession"].astype(str).str.split(".").str[0]
            ),
        })
    return pd.DataFrame(rows)


def classify_hint_architectures(domain_hits_df, sequence_ids=()):
    """Assign conservative, architecture-based labels to Hint-hit proteins.

    These labels are candidate classes rather than functional assertions.
    Canonical Hedgehog candidates require an HH_signal domain N-terminal to a
    Hint domain. VWA-Hint proteins are separated from other non-Hedge Hint
    architectures; the latter include possible Hoglet/Hog-like proteins.
    """
    architectures = summarize_pfam_architectures(domain_hits_df, sequence_ids)
    hit_groups = {
        name: group.sort_values(["ali_from", "ali_to"])
        for name, group in domain_hits_df.groupby("query_name")
    }
    labels = []
    for sequence_id in architectures["sequence_id"]:
        group = hit_groups.get(sequence_id)
        if group is None or group.empty:
            labels.append("Hint hit; no additional Pfam domain")
            continue
        names = group["target_name"].astype(str)
        accessions = group["target_accession"].astype(str).str.split(".").str[0]
        hint_positions = group.loc[(names == "Hint") | (accessions == "PF01079"), "ali_from"]
        hedge_positions = group.loc[
            (names == "HH_signal") | (accessions == "PF01085"), "ali_from"
        ]
        vwa_positions = group.loc[names.str.contains("VWA|Vwa", case=False), "ali_from"]
        if not hint_positions.empty and not hedge_positions.empty and hedge_positions.min() < hint_positions.min():
            labels.append("canonical Hedgehog candidate")
        elif not hint_positions.empty and not vwa_positions.empty and vwa_positions.min() < hint_positions.min():
            labels.append("VWA-Hint candidate")
        elif not hint_positions.empty and len(group) == len(hint_positions):
            labels.append("Hint-only Hog/Hoglet-like candidate")
        elif not hint_positions.empty:
            labels.append("other Hog/Hoglet-like candidate")
        else:
            labels.append("Hint search hit without Pfam Hint call")
    architectures["architecture_class"] = labels
    return architectures


def select_hint_domain_representatives(
    domain_hits_df,
    protein_metadata_df,
    profile_accession="PF01079",
):
    """Select one best HINT-domain hit per annotated gene locus.

    Representatives are ranked by domain bit score, domain coverage, and then
    complete-protein length. This removes alternative isoforms without
    removing HINT domains encoded at distinct loci.
    """
    domains = domain_hits_df.copy()
    accessions = domains["target_accession"].astype(str).str.split(".").str[0]
    domains = domains.loc[accessions.eq(profile_accession)].copy()
    if domains.empty:
        raise ValueError(f"No {profile_accession} domain hits were found")
    metadata = protein_metadata_df.copy()
    sequence_column = "sequence_id" if "sequence_id" in metadata else "target_name"
    required = {sequence_column, "gene_locus_id"}
    missing = required.difference(metadata.columns)
    if missing:
        raise ValueError(f"Protein metadata is missing columns: {sorted(missing)}")
    keep_columns = [
        column for column in [
            sequence_column, "gene_locus_id", "genome_acc", "protein_id",
            "organism_name", "architecture_class", "pfam_architecture",
            "protein_length",
        ]
        if column in metadata.columns
    ]
    representatives = domains.merge(
        metadata[keep_columns].drop_duplicates(sequence_column),
        left_on="query_name",
        right_on=sequence_column,
        how="left",
        validate="many_to_one",
    )
    if representatives["gene_locus_id"].isna().any():
        missing_ids = representatives.loc[
            representatives["gene_locus_id"].isna(), "query_name"
        ].tolist()
        raise ValueError("No gene locus for: " + ", ".join(missing_ids[:5]))
    representatives["hmm_coverage"] = (
        representatives["hmm_to"] - representatives["hmm_from"] + 1
    ) / representatives["tlen"]
    if "protein_length" not in representatives:
        representatives["protein_length"] = representatives["qlen"]
    representatives = representatives.sort_values(
        ["gene_locus_id", "domain_score", "hmm_coverage", "protein_length"],
        ascending=[True, False, False, False],
    ).drop_duplicates("gene_locus_id")
    return representatives.reset_index(drop=True)


def extract_domain_sequences(
    representatives_df,
    protein_sequences,
    output_fp=None,
    coordinate_columns=("env_from", "env_to"),
):
    """Extract one-indexed inclusive domain coordinates from protein sequences."""
    if isinstance(protein_sequences, (str, os.PathLike)):
        protein_sequences = pfa.get_seq_dict(protein_sequences)
    start_column, end_column = coordinate_columns
    extracted = {}
    for row in representatives_df.itertuples(index=False):
        sequence_id = str(row.query_name)
        if sequence_id not in protein_sequences:
            raise KeyError(f"Protein sequence not found: {sequence_id}")
        start = int(getattr(row, start_column))
        end = int(getattr(row, end_column))
        sequence = protein_sequences[sequence_id]
        if start < 1 or end > len(sequence) or start > end:
            raise ValueError(
                f"Invalid {start_column}/{end_column} coordinates for "
                f"{sequence_id}: {start}-{end} of {len(sequence)}"
            )
        extracted[sequence_id] = sequence[start - 1:end]
    if output_fp is not None:
        output_fp = Path(output_fp)
        output_fp.parent.mkdir(parents=True, exist_ok=True)
        pfa.write_to_fasta(extracted, output_fp, line_size=80)
    return extracted


def fetch_hmm_profile(hmm_id, database_fp, output_fp, hmmfetch_exe="hmmfetch"):
    """Fetch one named HMM profile from an indexed HMM database."""
    output_fp = Path(output_fp)
    output_fp.parent.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(
        [str(hmmfetch_exe), str(database_fp), str(hmm_id)],
        capture_output=True,
        text=True,
        check=True,
    )
    if not output_fp.is_file() or output_fp.read_text() != result.stdout:
        output_fp.write_text(result.stdout)
    return output_fp


def run_profile_alignment(
    sequences,
    hmm_fp,
    output_fp,
    hmmalign_exe="hmmalign",
    trim=True,
):
    """Align protein domains to a profile HMM in A2M format."""
    output_fp = Path(output_fp)
    output_fp.parent.mkdir(parents=True, exist_ok=True)
    temporary_fasta = None
    if isinstance(sequences, Mapping):
        temporary_fasta = tempfile.NamedTemporaryFile(
            mode="w", suffix=".faa", delete=False
        )
        temporary_fasta.close()
        sequence_fp = Path(temporary_fasta.name)
        pfa.write_to_fasta(dict(sequences), sequence_fp, line_size=80)
    else:
        sequence_fp = Path(sequences)
    command = [
        str(hmmalign_exe), "--outformat", "A2M", "-o", str(output_fp)
    ]
    if trim:
        command.append("--trim")
    command.extend([str(hmm_fp), str(sequence_fp)])
    try:
        subprocess.run(command, capture_output=True, text=True, check=True)
    finally:
        if temporary_fasta is not None:
            os.unlink(temporary_fasta.name)
    return output_fp


def a2m_to_match_alignment(a2m_fp, output_fp):
    """Remove A2M insertion columns, retaining profile match-state columns."""
    # Do not use the project's general FASTA parser here: it normalizes
    # residues to uppercase, whereas A2M encodes insertions as lowercase.
    sequences = {}
    name = None
    sequence_parts = []
    with Path(a2m_fp).open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(sequence_parts)
                name = line[1:].split()[0]
                if name in sequences:
                    raise ValueError(f"Duplicate A2M sequence identifier: {name}")
                sequence_parts = []
            elif name is None:
                raise ValueError("A2M sequence data precedes the first header")
            else:
                sequence_parts.append(line)
    if name is not None:
        sequences[name] = "".join(sequence_parts)
    if not sequences:
        raise ValueError(f"No sequences found in {a2m_fp}")
    # A2M does not pad insertion residues across sequences. Remove each
    # sequence's lowercase insertion residues (and insertion-gap dots)
    # independently; the remaining uppercase/dash strings are the aligned
    # HMM match states and must consequently have equal lengths.
    match_alignment = {
        name: "".join(
            residue.upper()
            for residue in sequence
            if residue == "-" or residue.isupper()
        )
        for name, sequence in sequences.items()
    }
    match_lengths = {len(sequence) for sequence in match_alignment.values()}
    if len(match_lengths) != 1:
        raise ValueError(
            "Removing A2M insertions did not produce equal match-state lengths"
        )
    if not match_lengths or match_lengths == {0}:
        raise ValueError("No profile match states were found")
    output_fp = Path(output_fp)
    output_fp.parent.mkdir(parents=True, exist_ok=True)
    pfa.write_to_fasta(match_alignment, output_fp, line_size=80)
    return match_alignment


def parse_hmmsearch_output(
    hmm_results_fp,
    cut_ga=False,
    hmm_database=None,
    hmmfetch_exe="hmmfetch",
):
    """Parse an ``hmmsearch --tblout`` file into a DataFrame.

    Parameters
    ----------
    hmm_results_fp : str or pathlib.Path
        Table produced by ``hmmsearch --tblout``.
    cut_ga : bool, default=False
        If true, retain hits whose full-sequence score is at least the
        profile-specific Pfam gathering score. The score is retrieved from
        the profile's ``GA`` line rather than supplied manually.
    hmm_database : str or pathlib.Path, optional
        HMM database containing the query profile. Required when ``cut_ga``
        is true. For a Pfam search, pass the ``Pfam-A.hmm`` file.
    hmmfetch_exe : str or pathlib.Path, default="hmmfetch"
        HMMER ``hmmfetch`` executable used to retrieve the profile.

    Returns
    -------
    pandas.DataFrame
        One row per reported sequence hit, optionally filtered at the
        profile's full-sequence gathering threshold. The applied threshold is
        also available as ``dataframe.attrs["gathering_score"]``.

    Notes
    -----
    An HMM ``GA`` line contains full-sequence and per-domain scores. Because
    ``--tblout`` is a sequence-level table, this function uses the first
    (full-sequence) score. Use ``--domtblout`` for domain-level filtering.
    """
    hmm_results_fp = Path(hmm_results_fp)
    if not hmm_results_fp.is_file():
        raise FileNotFoundError(f"hmmsearch table not found: {hmm_results_fp}")

    
    lines = []
    with hmm_results_fp.open() as f:
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
        df = pd.DataFrame(columns=columns)
    else:
        df = pd.DataFrame(lines, columns=columns)

    float_columns = [
        'E-value', 'score', 'bias', 'E-value_2', 'score_2', 'bias_2', 'exp'
    ]
    integer_columns = ['reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc']
    df[float_columns] = df[float_columns].astype(float)
    df[integer_columns] = df[integer_columns].astype(int)

    if cut_ga:
        if hmm_database is None:
            raise ValueError("hmm_database is required when cut_ga=True")
        if df.empty:
            raise ValueError(
                "Cannot retrieve a gathering score because the hmmsearch "
                "table contains no reported hits"
            )

        query_names = df['query_name'].dropna().unique()
        if len(query_names) != 1:
            raise ValueError(
                "cut_ga=True requires exactly one query profile; found "
                f"{len(query_names)}: {', '.join(query_names)}"
            )

        fetch_result = subprocess.run(
            [str(hmmfetch_exe), str(hmm_database), str(query_names[0])],
            capture_output=True,
            text=True,
            check=True,
        )
        ga_match = re.search(
            r"^GA\s+([+-]?(?:\d+(?:\.\d*)?|\.\d+))\s+"
            r"([+-]?(?:\d+(?:\.\d*)?|\.\d+))\s*;?",
            fetch_result.stdout,
            flags=re.MULTILINE,
        )
        if ga_match is None:
            raise ValueError(
                f"Profile {query_names[0]!r} does not define a GA threshold"
            )

        gathering_score = float(ga_match.group(1))
        domain_gathering_score = float(ga_match.group(2))
        print(
            f"Gathering scores for {query_names[0]}: "
            f"sequence={gathering_score:g} bits, "
            f"domain={domain_gathering_score:g} bits"
        )
        df = df.loc[df['score'] >= gathering_score].reset_index(drop=True)
        df.attrs['gathering_score'] = gathering_score
        df.attrs['domain_gathering_score'] = domain_gathering_score

    return df


def run_hmm_fetch_search(
    hmm_id,
    db_fp,
    proteome_fp,
    output_fp,
    hmm_file_path=None,
    domtblout_fp=None,
    raw_output_fp=None,
    hmmfetch_exe="hmmfetch",
    hmmsearch_exe="hmmsearch",
    **kwargs,
):
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
    output_fp = Path(output_fp)
    domtblout_fp = (
        Path(domtblout_fp)
        if domtblout_fp is not None
        else output_fp.with_suffix(".domtblout")
    )
    raw_output_fp = (
        Path(raw_output_fp)
        if raw_output_fp is not None
        else output_fp.with_suffix(".txt")
    )
    for path in (output_fp, domtblout_fp, raw_output_fp):
        path.parent.mkdir(parents=True, exist_ok=True)

    search_command = [
        str(hmmsearch_exe),
        "--tblout", str(output_fp),
        "--domtblout", str(domtblout_fp),
        "-o", str(raw_output_fp),
    ]
    for key, value in kwargs.items():
        if value is None or value is False:
            continue
        prefix = "-" if len(key) == 1 else "--"
        option_name = key if key in {"cut_ga", "cut_nc", "cut_tc"} else key.replace("_", "-")
        search_command.append(prefix + option_name)
        if value is not True:
            search_command.append(str(value))

    if hmm_file_path is not None:
        search_command.extend([str(hmm_file_path), str(proteome_fp)])
        return subprocess.run(search_command, capture_output=True, text=True, check=True)

    search_command.extend(["-", str(proteome_fp)])
    fetch_command = [str(hmmfetch_exe), str(db_fp), str(hmm_id)]
    fetch_process = subprocess.Popen(
        fetch_command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    try:
        result = subprocess.run(
            search_command,
            stdin=fetch_process.stdout,
            capture_output=True,
            text=True,
            check=True,
        )
    finally:
        if fetch_process.stdout is not None:
            fetch_process.stdout.close()
    fetch_stderr = fetch_process.communicate()[1]
    if fetch_process.returncode:
        raise subprocess.CalledProcessError(
            fetch_process.returncode,
            fetch_command,
            stderr=fetch_stderr,
        )
    return result


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
    #convert paths to strings
    hmm_fp    = str(hmm_fp)
    output_fp = str(output_fp)

    
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
        seqs_file = str(seqs)

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
