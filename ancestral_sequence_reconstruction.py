import math
import pandas as pd
import matplotlib.pyplot as plt
import logomaker
import tempfile
import os
import subprocess
from jw_utils import parse_fasta as pfa
from jw_utils import alignment_utils2 as au 


def infer_ASRseq_w_gap(
    anc_fp,
    tree_fp,
    alignment_fp,
    node_ID,
    gap_freq_thresh: float = 0.8,
):
    """
    Infer an ancestral sequence for `node_ID` from an IQ-TREE .state file, then
    post-process it using ONLY descendant clade gap frequency:

      - If gap frequency at a site >= gap_freq_thresh → call it a gap ('-').
      - Else → keep the ML amino-acid state from the ASR file.

    Parameters
    ----------
    anc_fp : str | pathlib.Path
        Path to IQ-TREE ancestral states file (.state; tab-delimited, '#' comments).
        Must contain at least: ['Node', 'Site', 'State'].
    tree_fp : str | pathlib.Path
        Newick tree used to define the clade; passed to `make_alignment_subset_from_clade`.
    alignment_fp : str | pathlib.Path
        Alignment file corresponding to the tree; used to compute per-site gap frequency.
    node_ID : str
        Internal node identifier (e.g., "Node38") matching the 'Node' column in `anc_fp`.
    gap_freq_thresh : float, default 0.8
        Column is considered 'gap' if gap frequency >= this threshold.

    Returns
    -------
    str
        The post-processed ancestral sequence for `node_ID`.

    Notes
    -----
    - Requires helper functions:
        - `make_alignment_subset_from_clade(tree_fp, anc_fp, node_ID, alignment_fp, ...)`
        - `au.get_frequency_matrix(aln_dict, molecule_type='prot')`
          which should return a DataFrame with rows as characters (including '-' if present)
          and columns as alignment positions that can be aligned to `anc_df['Site']`.
    """
    # Build sub-alignment for the clade to compute gap frequencies
    subset_aln_d = make_alignment_subset_from_clade(tree_fp, anc_fp, node_ID, alignment_fp)

    # Load ASR table and subset to the node
    anc_df = pd.read_csv(anc_fp, sep="\t", comment="#")
    required_cols = {"Node", "Site", "State"}
    if not required_cols.issubset(anc_df.columns):
        missing = required_cols - set(anc_df.columns)
        raise ValueError(f"ASR file missing required columns: {missing}")

    if node_ID not in set(anc_df["Node"]):
        raise KeyError(f"{node_ID} not found in 'Node' column of {anc_fp}")

    node_df = anc_df.groupby("Node").get_group(node_ID).copy()

    # Compute per-column gap frequency for the clade
    freq_df = au.get_frequency_matrix(subset_aln_d, molecule_type="prot")
    if "-" in freq_df.index:
        gap_freq_series = freq_df.loc["-", :]
    else:
        # No gaps at all in clade — set zeros of appropriate length
        gap_freq_series = pd.Series(0.0, index=freq_df.columns)
    gap_freq_series.index = gap_freq_series.index+1
    
    # Align gap frequencies to the ASR site's positions using 'Site' column.
    if "Site" not in node_df.columns:
        raise ValueError("Expected 'Site' column in ASR file to align positions.")

    def _to_int_or_str(x):
        try:
            return int(x)
        except Exception:
            return str(x)

    site_labels = node_df["Site"].apply(_to_int_or_str).tolist()
    gap_freq_map = { _to_int_or_str(col): float(gap_freq_series[col]) for col in gap_freq_series.index }

    gap_freqs_aligned = []
    for pos in site_labels:
        if pos in gap_freq_map:
            gap_freqs_aligned.append(gap_freq_map[pos])
        else:
            # Graceful off-by-one fallback (common if labels differ by 1)
            if isinstance(pos, int):
                if (pos + 1) in gap_freq_map:
                    gap_freqs_aligned.append(gap_freq_map[pos + 1])
                    continue
                if (pos - 1) in gap_freq_map:
                    gap_freqs_aligned.append(gap_freq_map[pos - 1])
                    continue
            # Default: treat as not gappy
            gap_freqs_aligned.append(0.0)

    # Use the ML state per site as the base AA
    AA_states = node_df["State"].astype(str).tolist()
    if len(AA_states) != len(gap_freqs_aligned):
        raise ValueError("Length mismatch between gap frequencies and ASR sites.")

    # Build the sequence: gap if highly gappy, else keep AA
    out_seq = []
    for aa, gap_freq in zip(AA_states, gap_freqs_aligned):
        out_seq.append("-" if gap_freq >= gap_freq_thresh else aa)

    return "".join(out_seq)



def make_cluster_df(clstr_fp):
    clusters = []
    with open(clstr_fp) as f:
        cid = None
        members = []
        rep = None
        for line in f:
            if line.startswith(">Cluster"):
                # save previous
                if cid is not None:
                    clusters.append({
                        "cluster_id": cid,
                        "rep": rep,
                        "member_count": len(members),
                        "members": ";".join(members)
                    })
                # start new
                cid = int(line.split()[1])
                members = []
                rep = None
            else:
                seq = line.split(">")[1].split("...")[0]
                members.append(seq)
                if "*" in line:
                    rep = seq
        # last cluster
        clusters.append({
            "cluster_id": cid,
            "rep": rep,
            "member_count": len(members),
            "members": ";".join(members)
        })
    return pd.DataFrame(clusters).sort_values("cluster_id").set_index('cluster_id')




def open_alignment_in_jalview(alignment_dict, jalview_path='/Applications/Jalview.app/Contents/MacOS/Jalview', colour='clustal'):
    # Create a temporary FASTA file
    with tempfile.NamedTemporaryFile(suffix=".fasta", delete=False) as tmp:
        temp_file_path = tmp.name
    
    try:
        # Write alignment to FASTA
        pfa.write_to_fasta(alignment_dict, temp_file_path)
        print(f"Temporary FASTA saved at {temp_file_path}")
        aln_d = pfa.get_seq_dict(temp_file_path)
        
        # Launch Jalview
        subprocess.run([jalview_path, "-open", temp_file_path, "-colour", colour])
    
    finally:
        # Remove file when Jalview closes
        if os.path.exists(temp_file_path):
            os.remove(temp_file_path)
            print(f"Temporary file {temp_file_path} deleted.")
    return aln_d


import os
import subprocess

def run_iqtree_asr(alignment_file, tree_file=None, model='MFP',
                   threads=2, prefix='my_asr', output_dir=None,
                   sh_alrt_reps=1000, ufboot_reps=None, outgroups=None):
    """
    Runs IQ-TREE ancestral sequence reconstruction with SH-aLRT support.
    Optionally also runs ultrafast bootstrap (UFBoot) and supports outgroup rooting.

    Parameters:
    alignment_file: Path to alignment (FASTA, PHYLIP, etc.).
    tree_file: Optional Newick tree. If None, IQ-TREE infers one.
    model: Substitution model (default 'MFP').
    threads: Number of CPU threads.
    prefix: Output file prefix.
    output_dir: Directory for output files (optional).
    sh_alrt_reps: SH-aLRT replicates (default: 1000).
    ufboot_reps: UFBoot replicates (default: None → skip).
    outgroups: Outgroup taxa (comma-separated string or path to text file with one name per line).
    """

    # Ensure output directory exists
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        prefix_path = os.path.join(output_dir, prefix)
    else:
        prefix_path = prefix
    
    # Build base command
    command = [
        "iqtree2",
        "-s", alignment_file,
        "-m", model,
        "--ancestral",
        "-T", str(threads),
        "--alrt", str(sh_alrt_reps),
        "--prefix", prefix_path
    ]
    
    # Add bootstrap if requested
    if ufboot_reps:
        command.extend(["-B", str(ufboot_reps)])
    
    # Add tree file if provided
    if tree_file:
        command.extend(["-t", tree_file])
    
    # Handle outgroups
    if outgroups:
        if os.path.isfile(outgroups):  # If it's a file, read and join
            with open(outgroups) as f:
                outgroup_list = [line.strip() for line in f if line.strip()]
            outgroup_str = ",".join(outgroup_list)
        else:
            outgroup_str = outgroups  # Assume it's already comma-separated
        command.extend(["-o", outgroup_str])

    # Run IQ-TREE
    subprocess.run(command, check=True)
    print(f"IQ-TREE ASR completed. Results saved with prefix '{prefix_path}'")

import pandas as pd

def parse_blast_output(blast_out_fp):
    """
    Parses BLAST tabular output from:
    blastp -query $inp -db nr -remote -out $out_fp -max_target_seqs 500 \
           -entrez_query "Azotobacter[ORGN]" \
           -outfmt "6 sacc qseqid sseqid pident length evalue bitscore staxids sscinames"
    
    Parameters
    ----------
    blast_out_fp : str
        Path to BLAST output file
    
    Returns
    -------
    pd.DataFrame
        DataFrame containing parsed BLAST results with proper column names and types
    """
    
    # Define column names in the exact order as your -outfmt
    cols = [
        "sacc",        # Subject accession
        "qseqid",      # Query sequence ID
        "sseqid",      # Subject sequence ID
        "pident",      # Percent identity
        "length",      # Alignment length
        "evalue",      # Expect value
        "bitscore",    # Bit score
        "staxids",     # Subject tax IDs
        "sscinames"    # Subject scientific names
    ]
    
    # Read the tab-delimited file
    df = pd.read_csv(blast_out_fp, sep="\t", names=cols, comment="#", dtype=str)
    
    # Convert numeric columns to appropriate types
    numeric_cols = ["pident", "length", "evalue", "bitscore"]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    
    return df

import requests


def fetch_proteins(ids, out_fp, batch_size=100):
    """fetches protein sequences in batches to avoid error and writes a fasta file
    IDs (str): protein IDs, eg wp_#########.1
    """
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    
    with open(out_fp, "w") as out_f:
        for i in range(0, len(ids), batch_size):
            batch_ids = ids[i:i+batch_size]
            id_str = ",".join(batch_ids)
            
            params = {
                "db": "protein",
                "id": id_str,
                "rettype": "fasta",
                "retmode": "text"
            }
            r = requests.get(base_url, params=params)
            r.raise_for_status()
            out_f.write(r.text)



def _relabel_clades_in_alignment(subclade_d, aln_d, intNodes_simp_2_tree_d, tree):
    """
    subclade_d (dict): clade_name: new name to be givein across clade, e.g {'Node768':'all_psuedomonas', 'Node737':'all_aeruginosa', 'Node41':'Node41'}
    """
    renamed_aln_subset_d = {}
    for node, new_name in subclade_d.items():
        
        subcl = tree.find_any(intNodes_simp_2_tree_d[node])
        subcl_names_d = {cl.name:new_name for cl in subcl.get_terminals()}
        t_d = {f'{subcl_names_d.get(name)}_{i}':seq for i, (name, seq) in enumerate(aln_d.items()) if subcl_names_d.get(name) }
        renamed_aln_subset_d = renamed_aln_subset_d | t_d
    return renamed_aln_subset_d


import pandas as pd
from Bio import Phylo

def make_alignment_subset_from_clade(
    tree_fp,
    anc_fp,
    node_ID,
    alignment_fp,
    thresh=0.7,
    subclade_d=None,
):
    """
    Build a dict of sequences for all terminal taxa under a given internal node,
    and append the node's ASR sequence (with low-confidence sites masked as 'X').
    Useful for extracting a sub-alignment to view in Jalview.

    Parameters
    ----------
    tree_fp : str or pathlib.Path
        Path to the Newick tree used for ASR.
    anc_fp : str or pathlib.Path
        Path to IQ-TREE `.state` file (tab-delimited; lines starting with '#' ignored).
        Expected columns: ['Node', 'Site', 'State', <probability columns...>].
    node_ID : str
        Internal node ID of interest (e.g., "Node38"). Must match the 'Node' column in `anc_fp`.
        IQ-TREE often formats internal names like "Node38/..." in the tree; we handle that.
    alignment_fp : str or pathlib.Path
        Path to the alignment file that the ASR tree was generated from.
        Must be readable by `pfa.get_seq_dict` (user-defined helper).
    thresh : float, default 0.7
        Posterior probability threshold; sites with max probability <= thresh are replaced by 'X'.
    subclade_d : dict, optional
        Mapping of {child_internal_node_simple_name: replacement_label}.
        If provided, sequences belonging to those child clades are relabeled using the dict values.
        Requires helper `_relabel_clades_in_alignment(...)` (note: function name as implemented).

    Returns
    -------
    dict
        { sequence_name: sequence_string } for all leaves under `node_ID`,
        plus one entry for the reconstructed ancestor: f"{node_ID}_ASR".

    Raises
    ------
    KeyError
        If `node_ID` is not found in the 'Node' column of `anc_fp`.
    ValueError
        If the internal node cannot be found in the tree.

    Notes
    -----
    - ASR is derived by taking the maximum posterior probability across the amino acid
      probability columns per site; sites below/at `thresh` are masked to 'X'.
    - Internal node names in the tree are mapped by stripping any suffix after '/' so
      "Node38/..." matches `node_ID == "Node38"`.
    - Missing leaf names in the alignment are skipped (with a warning comment).

    Example
    -------
    >>> subset = make_alignment_subset_from_clade(
    ...     tree_fp="asr.treefile",
    ...     anc_fp="asr.state",
    ...     node_ID="Node38",
    ...     alignment_fp="aln.fasta",
    ...     thresh=0.8,
    ...     subclade_d={"Node40": "Clade_A", "Node42": "Clade_B"},
    ... )
    >>> list(subset)[:3]
    ['taxon1', 'taxon2', 'Node38_ASR']
    """
    # --- Load ASR state table
    anc_df = pd.read_csv(anc_fp, sep="\t", comment="#")

    if 'Node' not in anc_df.columns or 'State' not in anc_df.columns:
        raise ValueError("Expected columns 'Node' and 'State' not found in ASR file.")

    # Validate node presence in ASR table
    if node_ID not in set(anc_df['Node']):
        raise KeyError(f"{node_ID} not found in 'Node' column of {anc_fp}.")

    # (Optional) your graph function; kept as-is if you rely on side-effects/plotting
    # make_node_anc_graph(anc_df, node_name=node_ID)

    # --- Build ASR sequence for the node
    node_df = anc_df.groupby('Node').get_group(node_ID)
    # Probability columns start at col index 3 -> max per site over AA probs
    max_prob_per_site = node_df.iloc[:, 3:].max(axis=1).tolist()
    ASR_seq = ''.join(
        aa if prob > thresh else 'X'
        for aa, prob in zip(node_df['State'], max_prob_per_site)
    )

    # --- Read tree and locate the internal clade
    tree = Phylo.read(str(tree_fp), format="newick")

    # Map "NodeXX" -> full internal name in the tree (e.g., "NodeXX/...")
    int_map = {
        cl.name.split('/')[0]: cl.name
        for cl in tree.get_nonterminals() if cl.name
    }
    internal_name = int_map.get(node_ID, None)
    if internal_name is None:
        # Fallback: maybe the tree already stores the simple name
        internal_name = node_ID

    clade = tree.find_any(internal_name)
    if clade is None:
        raise ValueError(f"Internal node '{node_ID}' not found in tree.")

    # --- Collect terminal names under the clade
    clade_leaf_names = [t.name for t in clade.get_terminals() if t.name]

    # --- Load alignment dict (expects {name: sequence})
    aln_d = pfa.get_seq_dict(alignment_fp)

    # Keep only leaves present in the alignment
    present = [n for n in clade_leaf_names if n in aln_d]
    missing = [n for n in clade_leaf_names if n not in aln_d]
    if missing:
        # You might prefer logging.warning here
        # print(f"[make_alignment_subset_from_clade] Missing in alignment: {len(missing)} names")
        pass

    aln_subset_d = {k: aln_d[k] for k in present}

    # --- Optionally relabel child clades
    if subclade_d:
        aln_subset_d = _relabel_clades_in_alignment(subclade_d, aln_subset_d, int_map, tree)

    # --- Append reconstructed ancestor
    aln_subset_d[f"{node_ID}_ASR"] = ASR_seq

    return aln_subset_d




def make_ASR_fasta(node_list, thresh, anc_fp, rename_d=None):
    """
    Generate ancestral sequence reconstructions (ASR) from an IQ-TREE
    ancestral states file, replacing low-confidence sites with 'X'.

    Parameters
    ----------
    node_list : list of str
        List of node identifiers (matching 'Node' column in the ancestral states file)
        for which to generate sequences.
    thresh : float
        Posterior probability threshold. Sites with maximum state probability 
        below this threshold are assigned 'X'.
    anc_fp : str
        Path to the IQ-TREE ancestral states file (tab-delimited, '#'-prefixed lines are ignored).
    rename_d : dict, optional
        Dictionary mapping node names to alternative names in the output.
        If None, original node names are used.

    Returns
    -------
    dict
        Dictionary mapping node names (or renamed keys) to reconstructed sequences
        as strings of amino acids, with 'X' for low-confidence sites.

    Notes
    -----
    - Assumes the IQ-TREE output file has columns:
        ['Node', 'Site', 'State', <probability columns for each AA>]
    - The highest probability among amino acid columns is compared to `thresh`.
    - Use 'X' for low-confidence sites to indicate ambiguity.

    Examples
    --------
    >>> node_list = ["Node1", "Node2"]
    >>> thresh = 0.8
    >>> anc_fp = "asr.state"
    >>> rename_d = {"Node1": "AncA", "Node2": "AncB"}
    >>> make_ASR_fasta(node_list, thresh, anc_fp, rename_d)
    {'AncA': 'MKT...X', 'AncB': 'MKTA...'}
    """
    anc_df = pd.read_csv(anc_fp, sep='\t', comment='#')
    ASR_seq_d = {}

    for node in node_list:
        node_df = anc_df.groupby('Node').get_group(node)
        max_prob_per_site = node_df.iloc[:, 3:].max(axis=1).tolist()
        ASR_seq = ''.join(
            AA if prob > thresh else 'X'
            for AA, prob in zip(node_df['State'], max_prob_per_site)
        )
        name = rename_d.get(node, node) if rename_d else node
        ASR_seq_d[name] = ASR_seq

    return ASR_seq_d


def make_node_logo(
    anc_fp,
    node,
    thresh,
    n_panels=1,
    figsize_per_panel=(10, 2.5),
    share_y=True,
    alphabet_order=None,
    title=None,
):
    """
    Plot posterior-probability sequence logos for a given node, split into n_panels
    vertically stacked subplots.

    Parameters
    ----------
    anc_fp : str
        Path to TSV with columns including 'Node', 'Site', and p_* residue probs.
    node : str or int
        Node identifier to select from the table.
    thresh : float
        Zero out probabilities <= thresh (for visual sparsity).
    n_panels : int, optional (default=1)
        Number of vertical subplots to split the sequence across.
    figsize_per_panel : (w, h), optional
        Size of each panel. Total figure size scales by n_panels in height.
    share_y : bool, optional
        Share y-axis across panels (nice if all panels are comparable).
    alphabet_order : list[str], optional
        If provided, reorder columns to this residue/alphabet order.
    title : str, optional
        Overall figure title.

    Returns
    -------
    matplotlib.figure.Figure
        The Matplotlib Figure containing the stacked logos.
    """
    # Load
    anc_df = pd.read_csv(anc_fp, sep="\t", comment="#")

    # Select node and massage columns like original
    try:
        post_prob_mat = anc_df.groupby('Node').get_group(node).copy()
    except KeyError as e:
        raise ValueError(f"Node {node!r} not found in {anc_fp}") from e

    post_prob_mat.columns = [c.replace('p_', '') for c in post_prob_mat.columns]
    post_prob_mat = post_prob_mat.set_index('Site').iloc[:, 2:]  # keep only prob columns

    # Optional column/alphabet ordering (e.g., standard AA order)
    if alphabet_order is not None:
        missing = [a for a in alphabet_order if a not in post_prob_mat.columns]
        if missing:
            raise ValueError(f"alphabet_order has columns not in data: {missing}")
        post_prob_mat = post_prob_mat[alphabet_order]

    # Thresholding
    post_prob_mat = post_prob_mat.mask(post_prob_mat <= thresh, 0.0)

    # Split rows into panels
    n_panels = max(1, int(n_panels))
    n_sites = post_prob_mat.shape[0]
    if n_panels > n_sites:
        n_panels = n_sites  # cap so we don't create empty axes

    # Build figure
    fig_w, fig_h = figsize_per_panel
    fig, axes = plt.subplots(
        n_panels, 1,
        figsize=(fig_w, fig_h * n_panels),
        sharex=False,
        sharey=share_y
    )
    if n_panels == 1:
        axes = [axes]

    # Determine chunk sizes by site index order
    # (preserve existing order of 'Site' index)
    sites = post_prob_mat.index.to_list()
    chunk_size = math.ceil(len(sites) / n_panels)
    chunks = [sites[i:i+chunk_size] for i in range(0, len(sites), chunk_size)]

    # Plot each chunk
    for i, ax in enumerate(axes):
        if i >= len(chunks):  # safety
            ax.axis('off')
            continue
        rows = chunks[i]
        sub = post_prob_mat.loc[rows]

        # Logomaker uses columns = alphabet, index = positions
        _ = logomaker.Logo(sub, ax=ax)

        # Cosmetics
        start_site, end_site = rows[0], rows[-1]
        ax.set_xlabel(f"Sites {start_site}–{end_site}")
        ax.set_ylabel("Posterior prob" if (i == 0 or not share_y) else "")
        ax.set_title(f"{title or f'Node {node}'} — segment {i+1}/{n_panels}", fontsize=10)
        ax.margins(x=0.01)

    fig.tight_layout(h_pad=0.6)
    return fig