from pathlib import Path
import pandas as pd
import numpy as np
from ete4 import Tree
import re
from Bio import AlignIO
from Bio.Seq import Seq
import subprocess
from Bio.Align import MultipleSeqAlignment
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO, Phylo
from copy import deepcopy
import subprocess
from Bio.Phylo.BaseTree import Clade


def label_foreground_clade_for_codeml(tree, foreground_tips, out_tree_fp, label="#1"):
    """
    Write a Newick tree with the MRCA of `foreground_tips` labeled for codeml.

    The foreground tips must be monophyletic. If they are not, the function
    raises a ValueError rather than labeling an MRCA that contains non-foreground
    taxa.

    Parameters
    ----------
    tree : Bio.Phylo.BaseTree.Tree
        Tree whose tip labels match the codon alignment.
    foreground_tips : iterable of str
        Tip names defining the foreground clade.
    out_tree_fp : fp to output newick foreground labeled tree.
    label : str, default="#1"
        codeml foreground label.

    Returns
    -------
    None
    """
    tree = deepcopy(tree)
    foreground_tips = set(foreground_tips)

    tree_tips = {tip.name for tip in tree.get_terminals()}
    missing = foreground_tips - tree_tips

    if missing:
        raise ValueError(
            f"{len(missing)} foreground tips not found in tree. "
            f"First few: {sorted(missing)[:10]}"
        )

    mrca = tree.common_ancestor(list(foreground_tips))
    mrca_tip_names = {tip.name for tip in mrca.get_terminals()}

    if mrca_tip_names != foreground_tips:
        extra = mrca_tip_names - foreground_tips
        missing_from_mrca = foreground_tips - mrca_tip_names

        raise ValueError(
            "Foreground tips are not monophyletic. "
            f"MRCA contains {len(extra)} extra tips and is missing "
            f"{len(missing_from_mrca)} foreground tips. "
            f"First extra tips: {sorted(extra)[:10]}"
        )
    for clade in tree.get_nonterminals():
        clade.name = None
        
    mrca.name = label
    Phylo.write(tree, str(out_tree_fp), "newick")
    return Path(out_tree_fp)
 

def resolve_polytomies_binary(tree, branch_length=1e-6):
    """
    Return a copy of a Bio.Phylo tree with all polytomies resolved into
    arbitrary binary nodes.

    This is mainly for tools like codeml that require bifurcating trees.
    The added internal branches are assigned a tiny branch length.
    """
    tree = deepcopy(tree)

    for clade in tree.find_clades(order="postorder"):
        while len(clade.clades) > 2:
            child1 = clade.clades.pop()
            child2 = clade.clades.pop()

            new_internal = Clade(
                branch_length=branch_length,
                clades=[child1, child2],
            )

            clade.clades.append(new_internal)

    return tree



def parse_codeml_output(outfile):
    """
    Parse key results from a PAML/codeml output file.

    Returns
    -------
    dict
        Parsed values including lnL, np, kappa, and omega estimates.
    """
    text = Path(outfile).read_text()

    results = {}

    # Example:
    # lnL(ntime: 1725  np: 1728):  -12345.678
    lnL_match = re.search(
        r"lnL\(ntime:\s*\d+\s+np:\s*(\d+)\):\s*([-\d.]+)",
        text
    )

    if lnL_match:
        results["np"] = int(lnL_match.group(1))
        results["lnL"] = float(lnL_match.group(2))

    # kappa
    kappa_match = re.search(r"kappa \(ts/tv\)\s*=\s*([.\deE+-]+)", text)
    if kappa_match:
        results["kappa"] = float(kappa_match.group(1))

    # One-ratio model often has:
    # omega (dN/dS) = 0.12345
    omega_match = re.search(r"omega \(dN/dS\)\s*=\s*([.\deE+-]+)", text)
    if omega_match:
        results["omega"] = float(omega_match.group(1))

    # Branch model often has:
    # w (dN/dS) for branches:  0.0500  0.1200
    branch_omega_match = re.search(
        r"w \(dN/dS\) for branches:\s*([^\n]+)",
        text
    )

    if branch_omega_match:
        omegas = [
            float(x)
            for x in branch_omega_match.group(1).split()
        ]
        results["branch_omegas"] = omegas

        if len(omegas) >= 2:
            results["omega_background"] = omegas[0]
            results["omega_foreground"] = omegas[1]

    return results


def write_codeml_control_file(
    seqfile,
    treefile,
    outfile,
    ctlfile,
    model=0,
    codon_freq=2,
    omega=0.1,
    fix_omega=0,
    cleandata=0,
):
    """
    Write a PAML/codeml control file.

    Parameters
    ----------
    seqfile : str or Path
        Codon alignment file in PAML/PHYLIP format.
    treefile : str or Path
        Newick tree file. For branch models, foreground branches should be
        labeled with codeml labels such as #1.
    outfile : str or Path
        codeml output file.
    ctlfile : str or Path
        Path to write the control file.
    model : int, default=0
        codeml model. Use 0 for one-ratio model and 2 for branch model.
    codon_freq : int, default=2
        Codon frequency model. 2 = F3x4.
    omega : float, default=0.1
        Initial omega value.
    fix_omega : int, default=0
        0 estimates omega; 1 fixes omega.
    cleandata : int, default=0
        0 keeps ambiguous/gapped sites; 1 removes them.
    """
    seqfile = Path(seqfile)
    treefile = Path(treefile)
    outfile = Path(outfile)
    ctlfile = Path(ctlfile)

    ctl_text = f"""
seqfile = {seqfile}
treefile = {treefile}
outfile = {outfile}

noisy = 9
verbose = 1
runmode = 0

seqtype = 1
CodonFreq = {codon_freq}

model = {model}
NSsites = 0

icode = 0

fix_kappa = 0
kappa = 2

fix_omega = {fix_omega}
omega = {omega}

cleandata = {cleandata}

fix_alpha = 1
alpha = 0
Malpha = 0
ncatG = 8

getSE = 0
RateAncestor = 0

Small_Diff = .5e-6

method = 0
""".strip()

    ctlfile.parent.mkdir(parents=True, exist_ok=True)
    ctlfile.write_text(ctl_text + "\n")

    return ctlfile




def run_codeml(ctlfile, codeml_exe="codeml"):
    """
    Run PAML codeml using a control file.

    Parameters
    ----------
    ctlfile : str or Path
        Path to codeml .ctl file.
    codeml_exe : str
        codeml executable name or path.

    Returns
    -------
    subprocess.CompletedProcess
        stdout/stderr from codeml.
    """
    ctlfile = Path(ctlfile)

    result = subprocess.run(
        [codeml_exe, str(ctlfile)],
        capture_output=True,
        text=True,
        check=True,
    )

    return result
 







def make_paml_safe_alignment_and_tree(
    aln_fasta,
    tree,
    output_prefix,
):
    aln = AlignIO.read(aln_fasta, "fasta")
    tree = deepcopy(tree)
    for clade in tree.get_nonterminals():
        clade.name = None
    aln_ids = [rec.id for rec in aln]
    tree_ids = [tip.name for tip in tree.get_terminals()]

    if set(aln_ids) != set(tree_ids):
        raise ValueError(
            f"Tree/alignment mismatch: "
            f"{len(set(aln_ids) - set(tree_ids))} in alignment not tree; "
            f"{len(set(tree_ids) - set(aln_ids))} in tree not alignment."
        )

    id_map = {
        old_id: f"T{i:06d}"
        for i, old_id in enumerate(sorted(aln_ids), start=1)
    }

    for rec in aln:
        rec.id = id_map[rec.id]
        rec.name = id_map[rec.name] if rec.name in id_map else rec.id
        rec.description = rec.id

    for tip in tree.get_terminals():
        tip.name = id_map[tip.name]

    aln_out = Path(f"{output_prefix}.phy")
    tree_out = Path(f"{output_prefix}.nwk")
    map_out = Path(f"{output_prefix}.id_map.tsv")

    AlignIO.write(aln, aln_out, "phylip-sequential")
    Phylo.write(tree, tree_out, "newick")

    pd.DataFrame(
        [{"old_id": old, "paml_id": new} for old, new in id_map.items()]
    ).to_csv(map_out, sep="\t", index=False)

    return aln_out, tree_out, map_out


def get_alignment_ids(alignment_fp, fmt="fasta"):
    """Return sequence IDs from a multiple sequence alignment."""
    aln = AlignIO.read(alignment_fp, fmt)
    return [rec.id for rec in aln]


def get_tree_tip_names(tree):
    """Return terminal/tip names from a Bio.Phylo tree."""
    return [tip.name for tip in tree.get_terminals()]


def check_tree_alignment_match(tree, alignment_fp, fmt="fasta"):
    """
    Check whether tree tip labels exactly match alignment sequence IDs.

    Returns
    -------
    dict
        Summary of shared, missing, and extra labels.
    """
    aln_ids = set(get_alignment_ids(alignment_fp, fmt=fmt))
    tree_ids = set(get_tree_tip_names(tree))

    return {
        "n_alignment_ids": len(aln_ids),
        "n_tree_tips": len(tree_ids),
        "n_shared": len(aln_ids & tree_ids),
        "alignment_not_in_tree": sorted(aln_ids - tree_ids),
        "tree_not_in_alignment": sorted(tree_ids - aln_ids),
        "exact_match": aln_ids == tree_ids,
    }


def prune_tree_to_alignment(tree, alignment, fmt="fasta"):
    """
    Return a copy of `tree` pruned to the sequence IDs present in an alignment.

    tree: Bio.Phylo 
    alignment: dict or fp
    """

    
    tree = deepcopy(tree)
    if type(alignment) == dict:
        aln_ids = set([k for k in alignment])
    else:
        aln_ids = set(get_alignment_ids(alignment_fp, fmt=fmt))

    for tip in list(tree.get_terminals()):
        if tip.name not in aln_ids:
            tree.prune(tip)

    return tree



from pathlib import Path
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment


def prune_alignment_to_tree(
    alignment,
    tree,
    output_fp=None,
    fmt="fasta",
):
    """
    Return an alignment pruned to the terminal tips present in a tree.

    Removes sequences from an alignment whose IDs are not present as
    terminal nodes in the tree. This is useful before analyses such as
    codeml, which require identical sequence IDs between the tree and
    codon alignment.

    Parameters
    ----------
    alignment : dict, str, pathlib.Path, or Bio.Align.MultipleSeqAlignment
        Alignment to prune.

        Accepted inputs:
        - dict:
            Keys are sequence IDs and values are sequences.
            Returns a pruned dictionary.

        - str or Path:
            Path to an alignment file readable by Bio.AlignIO.
            Returns a MultipleSeqAlignment.

        - MultipleSeqAlignment:
            Returns a MultipleSeqAlignment.

    tree : Bio.Phylo.BaseTree.Tree
        Tree whose terminal node names define the IDs to retain.

    output_fp : str or pathlib.Path, optional
        If provided, write the pruned alignment to this file.
        Only applies to MultipleSeqAlignment output.

    fmt : str, default="fasta"
        Alignment format used for reading/writing files.

    Returns
    -------
    dict or Bio.Align.MultipleSeqAlignment
        Alignment containing only sequences present in the tree.
        Return type matches the input type.
    """

    tree_ids = {tip.name for tip in tree.get_terminals()}

    # -------------------------
    # Dictionary input
    # -------------------------
    if isinstance(alignment, dict):
        pruned = {
            seq_id: seq
            for seq_id, seq in alignment.items()
            if seq_id in tree_ids
        }

        return pruned

    # -------------------------
    # File input
    # -------------------------
    if isinstance(alignment, (str, Path)):
        aln = AlignIO.read(alignment, fmt)

    # -------------------------
    # Biopython alignment input
    # -------------------------
    elif isinstance(alignment, MultipleSeqAlignment):
        aln = alignment

    else:
        raise TypeError(
            "alignment must be dict, file path, or MultipleSeqAlignment"
        )

    pruned_records = [
        rec for rec in aln
        if rec.id in tree_ids
    ]

    pruned_aln = MultipleSeqAlignment(pruned_records)

    if output_fp is not None:
        output_fp = Path(output_fp)
        output_fp.parent.mkdir(parents=True, exist_ok=True)
        AlignIO.write(pruned_aln, output_fp, fmt)

    return pruned_aln



def run_pal2nal(
    protein_alignment: str | Path,
    cds_fasta: str | Path,
    output_alignment: str | Path,
    output_format: str = "paml",
    pal2nal: str = "./PAL2NAL/pal2nal.pl",
):
    """
    Run PAL2NAL to generate a codon alignment from a protein alignment.

    Parameters
    ----------
    protein_alignment
        Protein multiple sequence alignment.
    cds_fasta
        Nucleotide CDS sequences corresponding to the aligned proteins.
    output_alignment
        Output codon alignment.
    output_format
        PAL2NAL output format (e.g. 'paml', 'fasta', 'clustal').
    pal2nal
        Path to the PAL2NAL executable or script.
    """

    protein_alignment = Path(protein_alignment)
    cds_fasta = Path(cds_fasta)
    output_alignment = Path(output_alignment)
    output_alignment.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_alignment, "w") as fout:
        subprocess.run(
            [
                pal2nal,
                str(protein_alignment),
                str(cds_fasta),
                "-output",
                output_format,
            ],
            stdout=fout,
            check=True,
        )




def qc_codon_alignment(
    alignment_fp,
    fmt="fasta",
    genetic_code=11,
    gap_char="-",
    gap_rich_threshold=0.8,
):
    """
    Run QC checks on a codon alignment before codeml/PAML analysis.

    Fatal checks include:
    - alignment length not divisible by 3
    - per-sequence ungapped length not divisible by 3
    - internal stop codons
    - frameshift-like gap runs whose lengths are not multiples of 3

    Warning checks include:
    - ambiguous nucleotide bases
    - duplicate sequence IDs
    - high gap fraction
    - gap-rich columns
    - fully gapped columns

    The function also reports variable codons and parsimony-informative codons,
    which are useful for assessing whether the alignment has enough signal for
    dN/dS model fitting.

    Parameters
    ----------
    alignment_fp : str or pathlib.Path
        Path to codon alignment file.
    fmt : str, default="fasta"
        Alignment format readable by Bio.AlignIO.
    genetic_code : int, default=11
        NCBI translation table. Use 11 for bacterial coding sequences.
    gap_char : str, default="-"
        Character used for alignment gaps.
    gap_rich_threshold : float, default=0.8
        Fraction of sequences with a gap required to flag a column as gap-rich.

    Returns
    -------
    dict
        Dictionary with:
        - summary : dict
        - sequence_qc : pandas.DataFrame
        - gap_column_qc : pandas.DataFrame
        - codon_qc : pandas.DataFrame
        - summary_text : str
    """
    alignment_fp = Path(alignment_fp)
    aln = AlignIO.read(alignment_fp, fmt)

    seq_ids = [rec.id for rec in aln]
    aln_len = aln.get_alignment_length()
    n_seqs = len(aln)

    duplicate_ids = sorted({
        seq_id for seq_id in seq_ids
        if seq_ids.count(seq_id) > 1
    })

    seq_rows = []

    for rec in aln:
        aligned_seq = str(rec.seq).upper()
        ungapped_seq = aligned_seq.replace(gap_char, "")

        gap_runs = re.findall(f"{re.escape(gap_char)}+", aligned_seq)
        frameshift_gap_runs = [
            run for run in gap_runs
            if len(run) % 3 != 0
        ]

        length_divisible_by_3 = len(ungapped_seq) % 3 == 0

        internal_stops = []
        terminal_stop = False

        if length_divisible_by_3 and len(ungapped_seq) > 0:
            protein = str(Seq(ungapped_seq).translate(table=genetic_code))
            terminal_stop = protein.endswith("*")
            internal_stops = [
                i for i, aa in enumerate(protein[:-1], start=1)
                if aa == "*"
            ]

        ambiguous_count = sum(
            base not in {"A", "C", "G", "T", gap_char}
            for base in aligned_seq
        )

        seq_rows.append({
            "seq_id": rec.id,
            "aligned_length": len(aligned_seq),
            "ungapped_length": len(ungapped_seq),
            "ungapped_length_divisible_by_3": length_divisible_by_3,
            "n_internal_stops": len(internal_stops),
            "internal_stop_positions_aa": internal_stops,
            "terminal_stop": terminal_stop,
            "n_frameshift_gap_runs": len(frameshift_gap_runs),
            "frameshift_gap_lengths": [len(run) for run in frameshift_gap_runs],
            "n_ambiguous_bases": ambiguous_count,
            "gap_fraction": aligned_seq.count(gap_char) / len(aligned_seq),
        })

    sequence_qc = pd.DataFrame(seq_rows).set_index("seq_id")

    aln_array = np.array([list(str(rec.seq).upper()) for rec in aln])
    gap_fraction_by_col = (aln_array == gap_char).mean(axis=0)

    gap_column_qc = pd.DataFrame({
        "column": np.arange(1, aln_len + 1),
        "gap_fraction": gap_fraction_by_col,
        "fully_gapped": gap_fraction_by_col == 1.0,
        "gap_rich": gap_fraction_by_col >= gap_rich_threshold,
    })

    # Codon-level QC
    codon_rows = []

    if aln_len % 3 == 0:
        for codon_idx, start in enumerate(range(0, aln_len, 3), start=1):
            codons = [
                str(rec.seq[start:start + 3]).upper()
                for rec in aln
            ]

            non_gap_codons = [
                codon for codon in codons
                if codon != gap_char * 3 and gap_char not in codon
            ]

            codon_counts = pd.Series(non_gap_codons).value_counts()

            variable = codon_counts.shape[0] > 1
            parsimony_informative = (codon_counts >= 2).sum() >= 2

            codon_rows.append({
                "codon_position": codon_idx,
                "alignment_start_col": start + 1,
                "alignment_end_col": start + 3,
                "gap_fraction": sum(gap_char in codon for codon in codons) / n_seqs,
                "n_observed_codons": len(non_gap_codons),
                "n_unique_codons": codon_counts.shape[0],
                "variable": variable,
                "parsimony_informative": parsimony_informative,
            })

    codon_qc = pd.DataFrame(codon_rows)

    fatal_issues = {
        "alignment_length_not_divisible_by_3": aln_len % 3 != 0,
        "sequences_not_divisible_by_3": int((~sequence_qc["ungapped_length_divisible_by_3"]).sum()),
        "sequences_with_internal_stops": int((sequence_qc["n_internal_stops"] > 0).sum()),
        "sequences_with_frameshift_gaps": int((sequence_qc["n_frameshift_gap_runs"] > 0).sum()),
    }

    warning_issues = {
        "duplicate_ids": len(duplicate_ids),
        "sequences_with_ambiguous_bases": int((sequence_qc["n_ambiguous_bases"] > 0).sum()),
        "total_ambiguous_bases": int(sequence_qc["n_ambiguous_bases"].sum()),
        "gap_rich_columns": int(gap_column_qc["gap_rich"].sum()),
        "fully_gapped_columns": int(gap_column_qc["fully_gapped"].sum()),
        "mean_gap_fraction": float(sequence_qc["gap_fraction"].mean()),
    }

    summary = {
        "alignment_file": str(alignment_fp),
        "n_sequences": n_seqs,
        "alignment_length": aln_len,
        "n_codons": aln_len // 3 if aln_len % 3 == 0 else None,
        "alignment_length_divisible_by_3": aln_len % 3 == 0,
        "n_duplicate_ids": len(duplicate_ids),
        "duplicate_ids": duplicate_ids,
        "n_sequences_with_internal_stops": fatal_issues["sequences_with_internal_stops"],
        "n_sequences_with_frameshift_gaps": fatal_issues["sequences_with_frameshift_gaps"],
        "n_sequences_not_divisible_by_3": fatal_issues["sequences_not_divisible_by_3"],
        "n_sequences_with_ambiguous_bases": warning_issues["sequences_with_ambiguous_bases"],
        "total_ambiguous_bases": warning_issues["total_ambiguous_bases"],
        "mean_gap_fraction": warning_issues["mean_gap_fraction"],
        "n_gap_rich_columns": warning_issues["gap_rich_columns"],
        "n_fully_gapped_columns": warning_issues["fully_gapped_columns"],
        "n_variable_codons": int(codon_qc["variable"].sum()) if not codon_qc.empty else None,
        "n_parsimony_informative_codons": int(codon_qc["parsimony_informative"].sum()) if not codon_qc.empty else None,
        "passes_fatal_qc": not any(
            bool(v) if isinstance(v, bool) else v > 0
            for v in fatal_issues.values()
        ),
        "fatal_issues": fatal_issues,
        "warning_issues": warning_issues,
    }

    summary_text = _format_codon_qc_summary(summary)

    return {
        "summary": summary,
        "sequence_qc": sequence_qc,
        "gap_column_qc": gap_column_qc,
        "codon_qc": codon_qc,
        "summary_text": summary_text,
    }


def _format_codon_qc_summary(summary):
    """Format a codon-alignment QC summary as readable text."""
    status = "PASS" if summary["passes_fatal_qc"] else "FAIL"

    lines = [
        "Codon alignment QC",
        "------------------",
        f"Status:                         {status}",
        f"Alignment file:                 {summary['alignment_file']}",
        f"Sequences:                      {summary['n_sequences']}",
        f"Alignment length:               {summary['alignment_length']} bp",
        f"Codons:                         {summary['n_codons']}",
        f"Length divisible by 3:          {summary['alignment_length_divisible_by_3']}",
        f"Variable codons:                {summary['n_variable_codons']}",
        f"Parsimony-informative codons:   {summary['n_parsimony_informative_codons']}",
        "",
        "Fatal checks",
        "------------",
        f"Sequences not divisible by 3:   {summary['n_sequences_not_divisible_by_3']}",
        f"Sequences with internal stops:  {summary['n_sequences_with_internal_stops']}",
        f"Sequences with frameshift gaps: {summary['n_sequences_with_frameshift_gaps']}",
        "",
        "Warnings",
        "--------",
        f"Duplicate IDs:                  {summary['n_duplicate_ids']}",
        f"Sequences with ambiguous bases: {summary['n_sequences_with_ambiguous_bases']}",
        f"Total ambiguous bases:          {summary['total_ambiguous_bases']}",
        f"Mean gap fraction:              {summary['mean_gap_fraction']:.3f}",
        f"Gap-rich columns:               {summary['n_gap_rich_columns']}",
        f"Fully gapped columns:           {summary['n_fully_gapped_columns']}",
    ]

    return "\n".join(lines)





def trim_gap_rich_codons(alignment_fp, output_fp, fmt="fasta", max_gap_fraction=0.8):
    aln = AlignIO.read(alignment_fp, fmt)
    aln_len = aln.get_alignment_length()

    if aln_len % 3 != 0:
        raise ValueError("Alignment length is not divisible by 3.")

    keep_positions = []

    for codon_start in range(0, aln_len, 3):
        codon_cols = range(codon_start, codon_start + 3)

        # fraction of sequences with any gap in this codon
        gap_fraction = sum(
            "-" in str(rec.seq[codon_start:codon_start + 3])
            for rec in aln
        ) / len(aln)

        if gap_fraction <= max_gap_fraction:
            keep_positions.extend(codon_cols)

    trimmed_records = []
    for rec in aln:
        trimmed_seq = "".join(str(rec.seq)[i] for i in keep_positions)
        trimmed_records.append(
            SeqRecord(
                seq=rec.seq.__class__(trimmed_seq),
                id=rec.id,
                name=rec.name,
                description=rec.description,
            )
        )

    trimmed = MultipleSeqAlignment(trimmed_records)
    AlignIO.write(trimmed, output_fp, fmt)

    return {
        "original_length": aln_len,
        "trimmed_length": trimmed.get_alignment_length(),
        "codons_kept": trimmed.get_alignment_length() // 3,
        "codons_removed": (aln_len - trimmed.get_alignment_length()) // 3,
    }


def calculate_relative_protein_rates(gene_tree, species_tree):
    """
    Calculate the relative evolutionary rate of each protein by comparing its
    root-to-tip distance in the gene tree to the corresponding root-to-tip
    distance in the species tree.

    Protein identifiers are assumed to begin with a 15-character genome
    accession (e.g. 'GCF_000001405'), which is converted to the accession
    format used in the species tree (e.g. 'GCF.000001405').

    Parameters
    ----------
    gene_tree : Bio.Phylo.BaseTree.Tree
        Gene tree with branch lengths representing amino acid substitutions.
    species_tree : Bio.Phylo.BaseTree.Tree
        Species tree containing the corresponding genome accessions.

    Returns
    -------
    pandas.DataFrame
        DataFrame indexed by protein ID with columns:

        - gene_root_dist
        - species_root_dist
        - relative_rate
        - genome_acc
    """

    def protein_to_genome_acc(protein_id):
        """Convert a protein ID to its genome accession."""
        acc = protein_id[:15]
        return acc[::-1].replace("_", ".", 1)[::-1]

    # Root-to-tip distances in the gene tree
    gene_root_dist = {
        clade.name: gene_tree.distance(clade)
        for clade in gene_tree.get_terminals()
    }

    # Root-to-tip distances in the species tree
    species_root_dist = {}

    for protein_id in gene_root_dist:
        genome_acc = protein_to_genome_acc(protein_id)
        species_clade = species_tree.find_any(name=genome_acc)

        if species_clade is not None:
            species_root_dist[protein_id] = species_tree.distance(species_clade)

    dist_df = pd.DataFrame({
        "gene_root_dist": pd.Series(gene_root_dist),
        "species_root_dist": pd.Series(species_root_dist),
    })

    dist_df["relative_rate"] = (
        dist_df["gene_root_dist"] /
        dist_df["species_root_dist"]
    )

    dist_df["genome_acc"] = (
        dist_df.index.to_series()
        .apply(protein_to_genome_acc)
    )
    return dist_df