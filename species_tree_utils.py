"""Utilities for reproducible protein-based species-tree inference."""

from __future__ import annotations

import os
import shutil
import subprocess
import re
import sys
from itertools import combinations
from pathlib import Path
from urllib.parse import unquote

import dendropy
import pandas as pd
from Bio import SeqIO
from dendropy.calculate import treecompare

from jw_utils import parse_fasta as pfa
from jw_utils import alignment_utils2 as au


BUSCO_STATUSES = ("Complete", "Duplicated", "Fragmented", "Missing")
BUSCO_STATUS_PRIORITY = {
    "Duplicated": 4,
    "Complete": 3,
    "Fragmented": 2,
    "Missing": 1,
}


# Dependency discovery


def resolve_executable(name, additional_candidates=()):
    """Return an executable path from ``PATH`` or explicit candidates."""
    executable = shutil.which(str(name))
    if executable is not None:
        return Path(executable)
    for candidate in map(Path, additional_candidates):
        if candidate.is_file():
            return candidate
    return None


def check_dependencies(tool_candidates):
    """Resolve command-line tools and return paths plus availability status.

    Parameters
    ----------
    tool_candidates : mapping
        ``{tool_name: [candidate_path, ...]}``. An empty candidate collection
        searches only the current ``PATH``.
    """
    resolved = {
        name: resolve_executable(name, candidates)
        for name, candidates in tool_candidates.items()
    }
    status = {name: path is not None for name, path in resolved.items()}
    return resolved, status


# BUSCO execution


def parse_gff_attributes(attributes):
    """Parse a GFF3 attribute field into an unescaped key/value mapping."""
    parsed = {}
    for item in attributes.rstrip().split(";"):
        if not item:
            continue
        key, separator, value = item.partition("=")
        if separator:
            parsed[unquote(key)] = unquote(value)
    return parsed


def map_gff_proteins_to_genes(gff_fp):
    """Map GFF3 protein IDs to gene loci through CDS/transcript parents.

    The returned gene identifiers represent annotation loci, not gene symbols.
    A protein is therefore considered an isoform only when its CDS resolves to
    the same explicit GFF3 gene feature as another protein.
    """
    transcript_to_gene = {}
    cds_records = []
    with Path(gff_fp).open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            feature_type = fields[2]
            attributes = parse_gff_attributes(fields[8])
            if feature_type in {"mRNA", "transcript"}:
                transcript_id = attributes.get("ID")
                gene_id = attributes.get("Parent")
                if transcript_id and gene_id:
                    transcript_to_gene[transcript_id] = gene_id.split(",")[0]
            elif feature_type == "CDS":
                cds_records.append(attributes)

    protein_to_gene = {}
    conflicts = set()
    for attributes in cds_records:
        protein_id = attributes.get("protein_id") or attributes.get("Name")
        parent = attributes.get("Parent", "").split(",")[0]
        gene_id = transcript_to_gene.get(parent)
        if gene_id is None and parent.startswith("gene-"):
            gene_id = parent
        if not protein_id or not gene_id:
            continue
        previous = protein_to_gene.setdefault(protein_id, gene_id)
        if previous != gene_id:
            conflicts.add(protein_id)
    for protein_id in conflicts:
        protein_to_gene.pop(protein_id, None)
    return protein_to_gene


def select_longest_isoform_per_gene(proteome_fp, gff_fp):
    """Select the longest annotated protein isoform for every gene locus.

    Proteins that cannot be resolved to a gene are retained individually.
    Length ties are resolved deterministically by protein ID.
    """
    proteome_fp = Path(proteome_fp)
    protein_to_gene = map_gff_proteins_to_genes(gff_fp)
    records = list(SeqIO.parse(proteome_fp, "fasta"))
    groups = {}
    rows = []
    for record in records:
        gene_id = protein_to_gene.get(record.id)
        group_id = gene_id if gene_id is not None else f"unmapped:{record.id}"
        groups.setdefault(group_id, []).append(record)

    selected_ids = set()
    for group_id, candidates in groups.items():
        selected = sorted(candidates, key=lambda x: (-len(x.seq), x.id))[0]
        selected_ids.add(selected.id)
        for record in candidates:
            rows.append({
                "protein_id": record.id,
                "gene_id": None if group_id.startswith("unmapped:") else group_id,
                "length": len(record.seq),
                "isoforms_at_gene": len(candidates),
                "selected": record.id == selected.id,
                "selection_reason": (
                    "unmapped_retained"
                    if group_id.startswith("unmapped:")
                    else "longest_isoform" if record.id == selected.id
                    else "shorter_isoform"
                ),
            })
    return records, selected_ids, pd.DataFrame(rows)


def reduce_proteome_to_longest_isoforms(proteome_fp, gff_fp, output_fp):
    """Write one longest protein per annotated gene and return its audit table."""
    output_fp = Path(output_fp)
    output_fp.parent.mkdir(parents=True, exist_ok=True)
    records, selected_ids, audit = select_longest_isoform_per_gene(
        proteome_fp, gff_fp
    )
    SeqIO.write((r for r in records if r.id in selected_ids), output_fp, "fasta")
    audit.insert(0, "assembly", Path(proteome_fp).stem)
    return audit


def reduce_proteomes_to_longest_isoforms(
    proteome_dir, gff_dir, output_dir, audit_fp=None
):
    """Create accession-named, gene-level proteomes for a database.

    A matching ``<accession>.gff`` is required for every ``<accession>.faa``.
    The combined audit table records every retained and removed protein.
    """
    output_dir = Path(output_dir)
    audits = []
    outputs = {}
    for proteome_fp in discover_proteomes(proteome_dir):
        accession = proteome_fp.stem
        gff_fp = Path(gff_dir) / f"{accession}.gff"
        if not gff_fp.is_file():
            raise FileNotFoundError(f"No matching GFF3 for {accession}: {gff_fp}")
        output_fp = output_dir / proteome_fp.name
        audits.append(
            reduce_proteome_to_longest_isoforms(proteome_fp, gff_fp, output_fp)
        )
        outputs[accession] = output_fp
    audit = pd.concat(audits, ignore_index=True)
    if audit_fp is not None:
        audit_fp = Path(audit_fp)
        audit_fp.parent.mkdir(parents=True, exist_ok=True)
        audit.to_csv(audit_fp, sep="\t", index=False)
    return outputs, audit


def discover_proteomes(proteome_dir, pattern="*.faa"):
    """Return sorted proteome paths and require at least one match."""
    proteome_dir = Path(proteome_dir)
    proteomes = sorted(proteome_dir.glob(pattern))
    if not proteomes:
        raise FileNotFoundError(
            f"No proteomes matching {pattern!r} found in {proteome_dir}"
        )
    return proteomes


def busco_run_is_complete(busco_run_dir):
    """Return whether a BUSCO run contains its final JSON summary."""
    return any(Path(busco_run_dir).glob("short_summary*.json"))


def build_busco_command(
    proteome_fp,
    output_dir,
    lineage="eukaryota_odb12.2",
    busco_exe="busco",
    cpu=8,
    force=False,
):
    """Construct a BUSCO protein-mode command for one proteome."""
    proteome_fp = Path(proteome_fp)
    busco_exe = Path(busco_exe) if Path(str(busco_exe)).is_file() else busco_exe
    launcher = [str(busco_exe)]
    if isinstance(busco_exe, Path):
        with busco_exe.open(errors="ignore") as handle:
            if handle.readline().strip() == "#!/usr/bin/env python3":
                launcher = [sys.executable, str(busco_exe)]
    command = launcher + [
        "--in", str(proteome_fp),
        "--out", proteome_fp.stem,
        "--out_path", str(Path(output_dir)),
        "--mode", "proteins",
        "--lineage_dataset", str(lineage),
        "--cpu", str(cpu),
    ]
    if force:
        command.append("--force")
    return command


def run_busco_proteome(
    proteome_fp,
    output_dir,
    lineage="eukaryota_odb12.2",
    busco_exe="busco",
    cpu=8,
    force=False,
):
    """Run or safely resume BUSCO for one accession-named proteome."""
    proteome_fp = Path(proteome_fp)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    run_dir = output_dir / proteome_fp.stem
    if busco_run_is_complete(run_dir) and not force:
        print(f"BUSCO already complete; skipping {proteome_fp.stem}")
        return run_dir

    command = build_busco_command(
        proteome_fp=proteome_fp,
        output_dir=output_dir,
        lineage=lineage,
        busco_exe=busco_exe,
        cpu=cpu,
        force=force,
    )
    print(f"Running BUSCO: {proteome_fp.stem}")
    environment = os.environ.copy()
    if Path(str(busco_exe)).is_file():
        environment["PATH"] = (
            f"{Path(busco_exe).parent}{os.pathsep}{environment.get('PATH', '')}"
        )
    subprocess.run(command, check=True, env=environment)
    return run_dir


def discover_busco_runs(output_dir, require_complete=True):
    """Reconstruct an accession-to-run mapping from an existing BUSCO folder."""
    output_dir = Path(output_dir)
    runs = {
        run_dir.name: run_dir
        for run_dir in sorted(output_dir.iterdir())
        if run_dir.is_dir()
        and (not require_complete or busco_run_is_complete(run_dir))
    }
    if not runs:
        qualifier = "complete " if require_complete else ""
        raise FileNotFoundError(f"No {qualifier}BUSCO runs found in {output_dir}")
    return runs


def run_busco_proteomes(
    proteome_dir,
    output_dir,
    lineage="eukaryota_odb12.2",
    busco_exe="busco",
    cpu=8,
    force=False,
):
    """Run BUSCO protein mode separately on every accession-named proteome.

    Completed runs are skipped unless ``force=True``. Separate runs make each
    BUSCO result unambiguously traceable to an NCBI assembly accession.

    Returns
    -------
    dict
        ``{assembly_accession: output_directory}`` for all proteomes.
    """
    output_dir = Path(output_dir)
    return {
        proteome_fp.stem: run_busco_proteome(
            proteome_fp=proteome_fp,
            output_dir=output_dir,
            lineage=lineage,
            busco_exe=busco_exe,
            cpu=cpu,
            force=force,
        )
        for proteome_fp in discover_proteomes(proteome_dir)
    }


# BUSCO parsing and marker selection


def find_busco_full_table(busco_run_dir):
    """Locate the single BUSCO ``full_table.tsv`` under one run directory."""
    matches = list(Path(busco_run_dir).glob("run_*/full_table.tsv"))
    if len(matches) != 1:
        raise FileNotFoundError(
            f"Expected one BUSCO full_table.tsv under {busco_run_dir}; "
            f"found {len(matches)}"
        )
    return matches[0]


def parse_busco_full_table(full_table_fp):
    """Parse a BUSCO v6 ``full_table.tsv`` into a typed DataFrame."""
    full_table_fp = Path(full_table_fp)
    records = []
    with full_table_fp.open() as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            fields.extend([""] * (8 - len(fields)))
            records.append(fields[:8])
    columns = [
        "busco_id",
        "status",
        "sequence",
        "score",
        "length",
        "orthodb_url",
        "description",
        "extra",
    ]
    table = pd.DataFrame(records, columns=columns)
    table["score"] = pd.to_numeric(table["score"], errors="coerce")
    table["length"] = pd.to_numeric(table["length"], errors="coerce").astype(
        "Int64"
    )
    return table.drop(columns="extra")


def collapse_busco_records(table):
    """Collapse repeated BUSCO rows while preserving duplication evidence."""
    if not table["busco_id"].duplicated().any():
        return table.copy()
    return (
        table.assign(
            _status_priority=table["status"].map(BUSCO_STATUS_PRIORITY).fillna(0)
        )
        .sort_values(
            ["busco_id", "_status_priority", "score"],
            ascending=[True, False, False],
        )
        .drop_duplicates("busco_id")
        .drop(columns="_status_priority")
    )


def busco_table_to_series(table):
    """Return indexed status and sequence-ID series from one BUSCO table."""
    table = collapse_busco_records(table).set_index("busco_id")
    return table["status"], table["sequence"]


def build_busco_status_matrix(busco_runs):
    """Build marker-by-assembly BUSCO status and sequence-ID matrices.

    Parameters
    ----------
    busco_runs : mapping
        ``{assembly_accession: BUSCO_run_directory}``.

    Returns
    -------
    tuple of pandas.DataFrame
        Status matrix followed by the matched sequence-ID matrix.
    """
    status_by_assembly = {}
    sequence_by_assembly = {}
    for accession, run_dir in busco_runs.items():
        table = parse_busco_full_table(find_busco_full_table(run_dir))
        status, sequence = busco_table_to_series(table)
        status_by_assembly[accession] = status
        sequence_by_assembly[accession] = sequence

    status_matrix = pd.DataFrame(status_by_assembly).fillna("Missing")
    sequence_matrix = pd.DataFrame(sequence_by_assembly).fillna("")
    status_matrix.index.name = "busco_id"
    sequence_matrix.index.name = "busco_id"
    return status_matrix.sort_index(), sequence_matrix.reindex(status_matrix.index)


def summarize_busco_statuses(status_matrix):
    """Count each BUSCO status across assemblies for every marker."""
    summary = pd.DataFrame(index=status_matrix.index)
    for status in BUSCO_STATUSES:
        summary[status.lower()] = status_matrix.eq(status).sum(axis=1)
    summary.index.name = "busco_id"
    return summary


def summarize_busco_assemblies(status_matrix):
    """Summarize BUSCO completeness and status percentages per assembly.

    ``complete_percent`` follows BUSCO's definition of complete markers and
    therefore includes both single-copy and duplicated complete calls.
    """
    total_markers = status_matrix.shape[0]
    if total_markers == 0:
        raise ValueError("BUSCO status matrix contains no markers")

    records = []
    for accession in status_matrix.columns:
        counts = status_matrix[accession].value_counts()
        single_copy = int(counts.get("Complete", 0))
        duplicated = int(counts.get("Duplicated", 0))
        fragmented = int(counts.get("Fragmented", 0))
        missing = int(counts.get("Missing", 0))
        records.append(
            {
                "assembly_accession": accession,
                "total_markers": total_markers,
                "complete_single_copy": single_copy,
                "complete_duplicated": duplicated,
                "complete_total": single_copy + duplicated,
                "fragmented": fragmented,
                "missing": missing,
                "complete_percent": 100 * (single_copy + duplicated) / total_markers,
                "single_copy_percent": 100 * single_copy / total_markers,
                "duplicated_percent": 100 * duplicated / total_markers,
                "fragmented_percent": 100 * fragmented / total_markers,
                "missing_percent": 100 * missing / total_markers,
            }
        )
    return (
        pd.DataFrame.from_records(records)
        .set_index("assembly_accession")
        .sort_values("complete_percent", ascending=False)
    )


def filter_assemblies_by_busco_completeness(
    status_matrix,
    minimum_complete_percent=80.0,
):
    """Return retained accessions, excluded accessions, and assembly QC.

    Complete duplicated markers count as present for proteome-completeness QC,
    but duplicated markers are still rejected later during locus selection.
    Assemblies at or above the threshold are retained.
    """
    if not 0 <= minimum_complete_percent <= 100:
        raise ValueError("minimum_complete_percent must be between 0 and 100")
    summary = summarize_busco_assemblies(status_matrix)
    summary["passes_completeness_filter"] = (
        summary["complete_percent"] >= minimum_complete_percent
    )
    summary["qc_category"] = summary["complete_percent"].map(
        lambda percent: (
            "exclude"
            if percent < minimum_complete_percent
            else "preferred"
            if percent >= 90
            else "marginal"
        )
    )
    retained = summary.index[summary["passes_completeness_filter"]].tolist()
    excluded = summary.index[~summary["passes_completeness_filter"]].tolist()
    return retained, excluded, summary


def validate_marker_selection_inputs(
    status_matrix,
    minimum_occupancy,
    outgroup_accessions,
):
    """Validate occupancy and outgroup settings before marker selection."""
    if not 1 <= minimum_occupancy <= status_matrix.shape[1]:
        raise ValueError(
            f"minimum_occupancy must be between 1 and {status_matrix.shape[1]}"
        )
    unknown_outgroups = set(outgroup_accessions).difference(status_matrix.columns)
    if unknown_outgroups:
        raise ValueError(
            "Outgroups absent from status matrix: "
            + ", ".join(sorted(unknown_outgroups))
        )


def marker_selection_mask(
    status_matrix,
    summary,
    minimum_occupancy,
    outgroup_accessions=(),
    reject_any_duplication=True,
):
    """Return the boolean mask implementing the marker-selection policy."""
    keep = summary["complete"] >= minimum_occupancy
    if reject_any_duplication:
        keep &= summary["duplicated"].eq(0)
    if outgroup_accessions:
        outgroup_complete = status_matrix.loc[:, outgroup_accessions].eq("Complete")
        keep &= outgroup_complete.any(axis=1)
    return keep


def select_single_copy_buscos(
    status_matrix,
    minimum_occupancy=14,
    outgroup_accessions=(),
    reject_any_duplication=True,
):
    """Select ubiquitous single-copy BUSCO markers for species-tree analysis.

    A retained marker must be ``Complete`` in at least ``minimum_occupancy``
    assemblies. By default it cannot be ``Duplicated`` in any assembly and it
    must be complete in at least one requested outgroup.
    """
    outgroup_accessions = tuple(outgroup_accessions)
    validate_marker_selection_inputs(
        status_matrix, minimum_occupancy, outgroup_accessions
    )
    summary = summarize_busco_statuses(status_matrix)
    summary["selected"] = marker_selection_mask(
        status_matrix=status_matrix,
        summary=summary,
        minimum_occupancy=minimum_occupancy,
        outgroup_accessions=outgroup_accessions,
        reject_any_duplication=reject_any_duplication,
    )
    return summary.sort_values(
        ["selected", "complete", "duplicated"],
        ascending=[False, False, True],
    )


def write_busco_selection(
    status_matrix,
    sequence_matrix,
    selection_summary,
    output_dir,
):
    """Write BUSCO QC matrices and the selected marker list."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    status_matrix.to_csv(output_dir / "busco_status_matrix.tsv", sep="\t")
    sequence_matrix.to_csv(output_dir / "busco_sequence_matrix.tsv", sep="\t")
    selection_summary.to_csv(output_dir / "busco_marker_summary.tsv", sep="\t")
    selected = selection_summary.index[selection_summary["selected"]]
    (output_dir / "selected_buscos.txt").write_text(
        "".join(f"{busco_id}\n" for busco_id in selected)
    )
    return selected.tolist()


# Sequence extraction and alignment


def find_single_copy_busco_directory(busco_run_dir):
    """Locate one run's BUSCO single-copy protein directory."""
    matches = list(
        Path(busco_run_dir).glob(
            "run_*/busco_sequences/single_copy_busco_sequences"
        )
    )
    if len(matches) != 1:
        raise FileNotFoundError(
            f"Expected one single-copy BUSCO directory under {busco_run_dir}; "
            f"found {len(matches)}"
        )
    return matches[0]


def read_single_fasta(fasta_fp):
    """Read a FASTA that must contain exactly one sequence."""
    sequences = pfa.get_seq_dict(fasta_fp)
    if len(sequences) != 1:
        raise ValueError(f"Expected one sequence in {fasta_fp}; found {len(sequences)}")
    return next(iter(sequences.items()))


def extract_busco_locus(busco_id, busco_runs, output_fp):
    """Collect one BUSCO protein per available assembly into one FASTA."""
    output_fp = Path(output_fp)
    sequences = {}
    records = []
    for accession, run_dir in sorted(busco_runs.items()):
        source_fp = find_single_copy_busco_directory(run_dir) / f"{busco_id}.faa"
        if not source_fp.is_file():
            continue
        original_id, sequence = read_single_fasta(source_fp)
        sequences[accession] = sequence
        records.append(
            {
                "busco_id": busco_id,
                "assembly_accession": accession,
                "original_protein_id": original_id,
                "sequence_length": len(sequence),
                "source_fasta": str(source_fp),
            }
        )
    if not sequences:
        raise FileNotFoundError(f"No single-copy sequences found for {busco_id}")
    output_fp.parent.mkdir(parents=True, exist_ok=True)
    pfa.write_to_fasta(sequences, output_fp, line_size=80)
    return pd.DataFrame.from_records(records)


def extract_selected_buscos(busco_ids, busco_runs, output_dir):
    """Extract selected BUSCO loci and write a combined sequence manifest."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    manifests = [
        extract_busco_locus(
            busco_id,
            busco_runs,
            output_dir / f"{busco_id}.faa",
        )
        for busco_id in busco_ids
    ]
    manifest = pd.concat(manifests, ignore_index=True)
    manifest.to_csv(output_dir / "busco_sequence_manifest.tsv", sep="\t", index=False)
    return manifest


def run_mafft_alignment(
    input_fp,
    output_fp,
    mafft_exe="mafft",
    threads=1,
    force=False,
):
    """Align one protein locus with MAFFT L-INS-i."""
    input_fp = Path(input_fp)
    output_fp = Path(output_fp)
    if output_fp.is_file() and not force:
        return output_fp
    output_fp.parent.mkdir(parents=True, exist_ok=True)
    command = [
        str(mafft_exe),
        "--localpair",
        "--maxiterate", "1000",
        "--thread", str(threads),
        str(input_fp),
    ]
    with output_fp.open("w") as output_handle:
        subprocess.run(
            command,
            stdout=output_handle,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
        )
    return output_fp


def align_busco_loci(
    input_dir,
    output_dir,
    mafft_exe="mafft",
    threads=1,
    force=False,
):
    """Align every extracted BUSCO FASTA and return output paths."""
    input_fastas = discover_proteomes(input_dir)
    return [
        run_mafft_alignment(
            input_fp,
            Path(output_dir) / f"{input_fp.stem}.aligned.faa",
            mafft_exe=mafft_exe,
            threads=threads,
            force=force,
        )
        for input_fp in input_fastas
    ]


def run_clipkit_alignment(
    input_fp,
    output_fp,
    clipkit_exe="clipkit",
    mode="smart-gap",
    force=False,
):
    """Trim one alignment with ClipKIT."""
    input_fp = Path(input_fp)
    output_fp = Path(output_fp)
    if output_fp.is_file() and not force:
        return output_fp
    output_fp.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [str(clipkit_exe), str(input_fp), "-m", mode, "-o", str(output_fp)],
        capture_output=True,
        text=True,
        check=True,
    )
    return output_fp


def trim_busco_alignments(
    alignment_paths,
    output_dir,
    clipkit_exe="clipkit",
    mode="smart-gap",
    force=False,
):
    """Trim BUSCO alignments and return trimmed paths."""
    output_dir = Path(output_dir)
    return [
        run_clipkit_alignment(
            alignment_fp,
            output_dir / alignment_fp.name.replace(".aligned.faa", ".trimmed.faa"),
            clipkit_exe=clipkit_exe,
            mode=mode,
            force=force,
        )
        for alignment_fp in map(Path, alignment_paths)
    ]


def calculate_alignment_stats(alignment_fp):
    """Calculate basic sequence, length, gap, and variable-site alignment QC."""
    alignment_fp = Path(alignment_fp)
    sequences = pfa.get_seq_dict(alignment_fp)
    lengths = {len(sequence) for sequence in sequences.values()}
    if len(lengths) != 1:
        raise ValueError(f"Sequences are not aligned in {alignment_fp}")
    alignment_length = lengths.pop()
    columns = zip(*sequences.values())
    variable_sites = 0
    for column in columns:
        residues = {residue for residue in column if residue not in "-.?X"}
        variable_sites += len(residues) > 1
    total_cells = len(sequences) * alignment_length
    gap_cells = sum(sequence.count("-") for sequence in sequences.values())
    return {
        "locus": alignment_fp.name.split(".", 1)[0],
        "alignment_file": str(alignment_fp),
        "taxa": len(sequences),
        "alignment_length": alignment_length,
        "variable_sites": variable_sites,
        "gap_percent": 100 * gap_cells / total_cells if total_cells else 0,
    }


def summarize_alignments(alignment_paths, output_fp=None):
    """Summarize multiple alignments and optionally write a TSV report."""
    summary = pd.DataFrame(
        calculate_alignment_stats(path) for path in map(Path, alignment_paths)
    ).sort_values("locus").reset_index(drop=True)
    if output_fp is not None:
        output_fp = Path(output_fp)
        output_fp.parent.mkdir(parents=True, exist_ok=True)
        summary.to_csv(output_fp, sep="\t", index=False)
    return summary


def remove_gappy_alignment_columns(
    input_fp,
    output_fp,
    threshold=0.87,
    force=False,
):
    """Filter an alignment with ``au.remove_gappy_columns`` and write FASTA."""
    input_fp = Path(input_fp)
    output_fp = Path(output_fp)
    if output_fp.is_file() and not force:
        return output_fp
    sequences = pfa.get_seq_dict(input_fp)
    filtered = au.remove_gappy_columns(sequences, threshold=threshold)
    if not filtered or not next(iter(filtered.values())):
        raise ValueError(f"Gap filtering removed every column from {input_fp}")
    output_fp.parent.mkdir(parents=True, exist_ok=True)
    pfa.write_to_fasta(filtered, output_fp, line_size=80)
    return output_fp


def filter_gappy_busco_alignments(
    alignment_paths,
    output_dir,
    threshold=0.87,
    force=False,
):
    """Apply ``au.remove_gappy_columns`` to multiple BUSCO alignments."""
    output_dir = Path(output_dir)
    return [
        remove_gappy_alignment_columns(
            alignment_fp,
            output_dir / alignment_fp.name.replace(".trimmed.faa", ".filtered.faa"),
            threshold=threshold,
            force=force,
        )
        for alignment_fp in map(Path, alignment_paths)
    ]


# Gene-tree inference


def run_iqtree_gene_tree(
    alignment_fp,
    output_dir,
    iqtree_exe="iqtree3",
    threads=1,
    bootstrap_replicates=1000,
    alrt_replicates=1000,
    force=False,
):
    """Infer one maximum-likelihood protein gene tree with IQ-TREE."""
    alignment_fp = Path(alignment_fp)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    locus = alignment_fp.name.split(".", 1)[0]
    prefix = output_dir / locus
    tree_fp = Path(f"{prefix}.treefile")
    if tree_fp.is_file() and not force:
        return tree_fp
    command = [
        str(iqtree_exe),
        "-s", str(alignment_fp),
        "-m", "MFP",
        "-bb", str(bootstrap_replicates),
        "-alrt", str(alrt_replicates),
        "-nt", str(threads),
        "--prefix", str(prefix),
    ]
    if force:
        command.append("-redo")
    subprocess.run(command, capture_output=True, text=True, check=True)
    return tree_fp


def infer_busco_gene_trees(
    alignment_paths,
    output_dir,
    iqtree_exe="iqtree3",
    threads=1,
    bootstrap_replicates=1000,
    alrt_replicates=1000,
    force=False,
):
    """Infer restartable IQ-TREE gene trees for all BUSCO alignments."""
    trees = []
    for index, alignment_fp in enumerate(map(Path, alignment_paths), start=1):
        print(f"Gene tree {index}: {alignment_fp.stem}")
        trees.append(
            run_iqtree_gene_tree(
                alignment_fp=alignment_fp,
                output_dir=output_dir,
                iqtree_exe=iqtree_exe,
                threads=threads,
                bootstrap_replicates=bootstrap_replicates,
                alrt_replicates=alrt_replicates,
                force=force,
            )
        )
    return trees


# Concatenated species-tree inference


def concatenate_protein_alignments(
    alignment_paths,
    output_fasta_fp,
    partition_fp,
    taxa=None,
):
    """Concatenate aligned protein loci and write an IQ-TREE partition file.

    Taxa absent from a locus are represented by gaps. Alignment columns and
    partitions follow the sorted input-path order, making the supermatrix
    deterministic and auditable.

    Returns
    -------
    tuple
        ``(output_fasta_path, partition_path, partition_dataframe)``.
    """
    alignment_paths = sorted(map(Path, alignment_paths))
    if not alignment_paths:
        raise ValueError("No protein alignments were supplied")

    loci = []
    observed_taxa = set()
    for alignment_fp in alignment_paths:
        sequences = pfa.get_seq_dict(alignment_fp)
        lengths = {len(sequence) for sequence in sequences.values()}
        if len(lengths) != 1:
            raise ValueError(f"Sequences are not aligned in {alignment_fp}")
        length = lengths.pop()
        if length == 0:
            raise ValueError(f"Alignment contains zero columns: {alignment_fp}")
        loci.append((alignment_fp.stem.split(".filtered")[0], sequences, length))
        observed_taxa.update(sequences)

    taxa = sorted(observed_taxa if taxa is None else map(str, taxa))
    unexpected = observed_taxa.difference(taxa)
    if unexpected:
        raise ValueError("Alignments contain unexpected taxa: " + ", ".join(sorted(unexpected)))

    concatenated = {taxon: [] for taxon in taxa}
    partitions = []
    start = 1
    for locus, sequences, length in loci:
        end = start + length - 1
        partitions.append({
            "locus": locus,
            "start": start,
            "end": end,
            "sites": length,
            "taxa_present": len(sequences),
        })
        for taxon in taxa:
            concatenated[taxon].append(sequences.get(taxon, "-" * length))
        start = end + 1

    output_fasta_fp = Path(output_fasta_fp)
    partition_fp = Path(partition_fp)
    output_fasta_fp.parent.mkdir(parents=True, exist_ok=True)
    partition_fp.parent.mkdir(parents=True, exist_ok=True)
    fasta_sequences = {
        taxon: "".join(parts) for taxon, parts in concatenated.items()
    }
    fasta_content = "".join(
        f">{taxon}\n"
        + "\n".join(
            sequence[index:index + 80]
            for index in range(0, len(sequence), 80)
        )
        + "\n"
        for taxon, sequence in fasta_sequences.items()
    )
    if not output_fasta_fp.is_file() or output_fasta_fp.read_text() != fasta_content:
        output_fasta_fp.write_text(fasta_content)
    partition_content = "".join(
        f"LG, {row['locus']} = {row['start']}-{row['end']}\n"
        for row in partitions
    )
    if not partition_fp.is_file() or partition_fp.read_text() != partition_content:
        partition_fp.write_text(partition_content)
    return output_fasta_fp, partition_fp, pd.DataFrame(partitions)


def run_partitioned_iqtree_species_tree(
    alignment_fp,
    partition_fp,
    output_prefix,
    iqtree_exe="iqtree3",
    threads=8,
    bootstrap_replicates=1000,
    alrt_replicates=1000,
    model="MFP+MERGE",
    force=False,
):
    """Infer a partitioned concatenated protein tree with IQ-TREE.

    Existing output is reused only when it is newer than both the supermatrix
    and partition definition. The complete command is saved beside the tree.
    """
    alignment_fp = Path(alignment_fp)
    partition_fp = Path(partition_fp)
    output_prefix = Path(output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    tree_fp = Path(f"{output_prefix}.treefile")
    inputs_mtime = max(alignment_fp.stat().st_mtime, partition_fp.stat().st_mtime)
    if tree_fp.is_file() and not force and tree_fp.stat().st_mtime >= inputs_mtime:
        print(f"Partitioned IQ-TREE output is current; skipping {output_prefix.name}")
        return tree_fp

    command = [
        str(iqtree_exe),
        "-s", str(alignment_fp),
        "-p", str(partition_fp),
        "-m", str(model),
        "-B", str(bootstrap_replicates),
        "--alrt", str(alrt_replicates),
        "-T", str(threads),
        "--prefix", str(output_prefix),
    ]
    if force:
        command.append("-redo")
    Path(f"{output_prefix}.command.txt").write_text(
        " ".join(map(str, command)) + "\n"
    )
    subprocess.run(command, check=True)
    if not tree_fp.is_file():
        raise FileNotFoundError(f"IQ-TREE did not create {tree_fp}")
    return tree_fp


def prepare_iqtree_display_tree(
    input_tree_fp,
    taxonomy_fp,
    output_tree_fp,
):
    """Create a plotting copy with organism names and IQ-TREE support labels."""
    tree = dendropy.Tree.get(
        path=str(input_tree_fp),
        schema="newick",
        rooting="force-rooted",
        preserve_underscores=True,
    )
    taxonomy = pd.read_csv(taxonomy_fp, sep="\t", dtype=str, keep_default_na=False)
    name_by_accession = taxonomy.set_index("assembly_accession")["organism_name"].to_dict()
    for leaf in tree.leaf_node_iter():
        if leaf.taxon.label not in name_by_accession:
            raise ValueError(f"No organism name for {leaf.taxon.label}")
        leaf.taxon.label = name_by_accession[leaf.taxon.label]
    output_tree_fp = Path(output_tree_fp)
    output_tree_fp.parent.mkdir(parents=True, exist_ok=True)
    tree.write(
        path=str(output_tree_fp),
        schema="newick",
        suppress_rooting=True,
    )
    return output_tree_fp


# Gene-tree discordance


def read_unrooted_tree(tree_fp, taxon_namespace=None):
    """Read a Newick tree as unrooted using an optional shared namespace."""
    return dendropy.Tree.get(
        path=str(tree_fp),
        schema="newick",
        rooting="force-unrooted",
        taxon_namespace=taxon_namespace,
        preserve_underscores=True,
    )


def tree_tip_labels(tree_fp):
    """Return the set of terminal taxon labels in a Newick tree."""
    tree = read_unrooted_tree(tree_fp)
    return {taxon.label for taxon in tree.taxon_namespace}


def read_pruned_unrooted_tree(tree_fp, retained_taxa, taxon_namespace):
    """Read and prune a tree, then migrate it into a shared taxon namespace."""
    tree = read_unrooted_tree(tree_fp)
    tree.retain_taxa_with_labels(retained_taxa)
    return dendropy.Tree.get(
        data=tree.as_string(
            schema="newick",
            suppress_rooting=True,
            unquoted_underscores=True,
        ),
        schema="newick",
        rooting="force-unrooted",
        taxon_namespace=taxon_namespace,
        preserve_underscores=True,
    )


def normalized_robinson_foulds_distance(tree_a_fp, tree_b_fp):
    """Calculate normalized unrooted RF distance on pairwise shared taxa."""
    shared_taxa = tree_tip_labels(tree_a_fp).intersection(tree_tip_labels(tree_b_fp))
    if len(shared_taxa) < 4:
        return float("nan"), len(shared_taxa)
    namespace = dendropy.TaxonNamespace(sorted(shared_taxa))
    tree_a = read_pruned_unrooted_tree(tree_a_fp, shared_taxa, namespace)
    tree_b = read_pruned_unrooted_tree(tree_b_fp, shared_taxa, namespace)
    tree_a.encode_bipartitions()
    tree_b.encode_bipartitions()
    rf_distance = treecompare.symmetric_difference(tree_a, tree_b)
    maximum_rf = 2 * (len(shared_taxa) - 3)
    return rf_distance / maximum_rf, len(shared_taxa)


def gene_tree_rf_matrices(tree_paths):
    """Return pairwise normalized RF-distance and shared-taxon matrices."""
    tree_paths = [Path(path) for path in tree_paths]
    loci = [path.stem for path in tree_paths]
    rf_matrix = pd.DataFrame(0.0, index=loci, columns=loci)
    shared_taxa_matrix = pd.DataFrame(0, index=loci, columns=loci, dtype=int)
    for path, locus in zip(tree_paths, loci):
        shared_taxa_matrix.loc[locus, locus] = len(tree_tip_labels(path))
    for (path_a, locus_a), (path_b, locus_b) in combinations(
        zip(tree_paths, loci), 2
    ):
        distance, shared_taxa = normalized_robinson_foulds_distance(path_a, path_b)
        rf_matrix.loc[locus_a, locus_b] = distance
        rf_matrix.loc[locus_b, locus_a] = distance
        shared_taxa_matrix.loc[locus_a, locus_b] = shared_taxa
        shared_taxa_matrix.loc[locus_b, locus_a] = shared_taxa
    return rf_matrix, shared_taxa_matrix


def summarize_gene_tree_discordance(rf_matrix):
    """Rank loci by their mean and maximum normalized RF distance."""
    summary = pd.DataFrame(
        {
            "mean_normalized_rf": rf_matrix.apply(
                lambda row: row.drop(row.name).mean(), axis=1
            ),
            "median_normalized_rf": rf_matrix.apply(
                lambda row: row.drop(row.name).median(), axis=1
            ),
            "maximum_normalized_rf": rf_matrix.apply(
                lambda row: row.drop(row.name).max(), axis=1
            ),
        }
    )
    summary.index.name = "locus"
    return summary.sort_values("mean_normalized_rf", ascending=False)


def build_common_taxon_consensus(tree_paths, output_fp, minimum_frequency=0.5):
    """Build an unrooted consensus after pruning all trees to common taxa."""
    tree_paths = [Path(path) for path in tree_paths]
    common_taxa = set.intersection(*(tree_tip_labels(path) for path in tree_paths))
    if len(common_taxa) < 4:
        raise ValueError("Fewer than four taxa occur in every gene tree")
    namespace = dendropy.TaxonNamespace(sorted(common_taxa))
    trees = dendropy.TreeList(taxon_namespace=namespace)
    for path in tree_paths:
        trees.append(read_pruned_unrooted_tree(path, common_taxa, namespace))
    consensus = trees.consensus(min_freq=minimum_frequency)
    output_fp = Path(output_fp)
    output_fp.parent.mkdir(parents=True, exist_ok=True)
    consensus.write(
        path=str(output_fp),
        schema="newick",
        suppress_rooting=True,
        unquoted_underscores=True,
    )
    return consensus, sorted(common_taxa)


def parse_iqtree_branch_support(label):
    """Parse an IQ-TREE ``SH-aLRT/UFBoot`` internal-node label."""
    if label is None:
        return None
    match = re.fullmatch(
        r"\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+))/"
        r"([+-]?(?:\d+(?:\.\d*)?|\.\d+))\s*",
        str(label),
    )
    if match is None:
        return None
    return float(match.group(1)), float(match.group(2))


def collapse_weak_iqtree_branches(
    input_tree_fp,
    output_tree_fp,
    minimum_alrt=80.0,
    minimum_ufboot=95.0,
):
    """Collapse branches failing either IQ-TREE support threshold."""
    tree = read_unrooted_tree(input_tree_fp)
    internal_nodes = [
        node
        for node in tree.postorder_node_iter()
        if node is not tree.seed_node and not node.is_leaf()
    ]
    retained = 0
    collapsed = 0
    missing_support = 0
    for node in internal_nodes:
        support = parse_iqtree_branch_support(node.label)
        if support is None:
            missing_support += 1
            node.edge.collapse()
            collapsed += 1
            continue
        alrt, ufboot = support
        if alrt >= minimum_alrt and ufboot >= minimum_ufboot:
            retained += 1
        else:
            node.edge.collapse()
            collapsed += 1

    output_tree_fp = Path(output_tree_fp)
    output_tree_fp.parent.mkdir(parents=True, exist_ok=True)
    tree.write(
        path=str(output_tree_fp),
        schema="newick",
        suppress_rooting=True,
        unquoted_underscores=True,
    )
    return {
        "locus": Path(input_tree_fp).stem,
        "internal_branches": len(internal_nodes),
        "retained_branches": retained,
        "collapsed_branches": collapsed,
        "missing_support_labels": missing_support,
        "retained_percent": (
            100 * retained / len(internal_nodes) if internal_nodes else 0
        ),
        "collapsed_tree": str(output_tree_fp),
    }


def collapse_busco_gene_trees(
    tree_paths,
    output_dir,
    minimum_alrt=80.0,
    minimum_ufboot=95.0,
):
    """Collapse weak branches in multiple IQ-TREE trees and return QC."""
    output_dir = Path(output_dir)
    records = []
    for tree_fp in map(Path, tree_paths):
        records.append(
            collapse_weak_iqtree_branches(
                tree_fp,
                output_dir / f"{tree_fp.stem}.collapsed.tree",
                minimum_alrt=minimum_alrt,
                minimum_ufboot=minimum_ufboot,
            )
        )
    return pd.DataFrame.from_records(records).sort_values(
        "retained_percent", ascending=False
    ).reset_index(drop=True)


# ASTRAL species-tree inference


def write_newick_tree_list(tree_paths, output_fp):
    """Write one normalized Newick gene tree per line for ASTRAL.

    An unchanged file is left untouched so its modification time can serve as
    a reliable restart/provenance signal for downstream programs.
    """
    output_fp = Path(output_fp)
    output_fp.parent.mkdir(parents=True, exist_ok=True)
    newick_lines = []
    for tree_fp in map(Path, tree_paths):
        newick = tree_fp.read_text().strip()
        if not newick.endswith(";"):
            newick += ";"
        newick_lines.append(newick)
    content = "\n".join(newick_lines) + "\n"
    if not output_fp.is_file() or output_fp.read_text() != content:
        output_fp.write_text(content)
    return output_fp


def run_astral(
    tree_paths,
    output_dir,
    analysis_name,
    astral_exe="astral",
    annotation_level=2,
    force=False,
):
    """Infer an ASTRAL species tree from unrooted gene trees."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    input_fp = write_newick_tree_list(
        tree_paths, output_dir / f"{analysis_name}.gene_trees.tre"
    )
    output_tree_fp = output_dir / f"{analysis_name}.astral.tree"
    log_fp = output_dir / f"{analysis_name}.astral.log"
    if (
        output_tree_fp.is_file()
        and not force
        and output_tree_fp.stat().st_mtime >= input_fp.stat().st_mtime
    ):
        print(f"ASTRAL output is current; skipping {analysis_name}")
        return output_tree_fp
    if output_tree_fp.is_file() and not force:
        print(f"ASTRAL input changed; rerunning {analysis_name}")
    command = [
        str(astral_exe),
        "-i", str(input_fp),
        "-o", str(output_tree_fp),
        "-t", str(annotation_level),
    ]
    result = subprocess.run(command, capture_output=True, text=True, check=True)
    log_fp.write_text(result.stdout + result.stderr)
    return output_tree_fp


def root_tree_on_outgroup_clade(input_tree_fp, output_tree_fp, outgroup_taxa):
    """Root a Newick tree on an exclusive, monophyletic outgroup clade."""
    tree = read_unrooted_tree(input_tree_fp)
    outgroup_taxa = set(outgroup_taxa)
    observed_taxa = {taxon.label for taxon in tree.taxon_namespace}
    missing = outgroup_taxa.difference(observed_taxa)
    if missing:
        raise ValueError("Outgroup taxa absent from tree: " + ", ".join(sorted(missing)))

    outgroup_node = None
    for node in tree.postorder_node_iter():
        descendant_taxa = {leaf.taxon.label for leaf in node.leaf_iter()}
        if descendant_taxa == outgroup_taxa:
            outgroup_node = node
            break
    if outgroup_node is None or outgroup_node.edge.tail_node is None:
        raise ValueError(
            "Requested outgroups do not form an exclusive clade: "
            + ", ".join(sorted(outgroup_taxa))
        )

    edge_length = outgroup_node.edge.length
    half_length = edge_length / 2 if edge_length is not None else None
    tree.reroot_at_edge(
        outgroup_node.edge,
        length1=half_length,
        length2=half_length,
        update_bipartitions=False,
    )
    output_tree_fp = Path(output_tree_fp)
    output_tree_fp.parent.mkdir(parents=True, exist_ok=True)
    tree.write(
        path=str(output_tree_fp),
        schema="newick",
        suppress_rooting=False,
        unquoted_underscores=True,
    )
    return output_tree_fp


def duplication_loss_reconciliation_score(gene_tree, species_tree, gene_to_species):
    """Return LCA-reconciliation duplication and loss counts for rooted trees.

    The implementation uses the standard duplication-loss parsimony model.
    A gene-tree node is a duplication when its species-tree LCA mapping equals
    the mapping of at least one child. Losses along a gene edge equal the
    species-tree distance minus one after a speciation, or the full distance
    after a duplication.
    """
    species_leaves = {
        node.taxon.label: node for node in species_tree.leaf_node_iter()
    }
    missing_species = set(gene_to_species.values()).difference(species_leaves)
    if missing_species:
        raise ValueError(
            "Gene mappings reference species absent from the species tree: "
            + ", ".join(sorted(missing_species))
        )

    def species_lca(node_a, node_b):
        ancestors = set()
        node = node_a
        while node is not None:
            ancestors.add(node)
            node = node.parent_node
        node = node_b
        while node not in ancestors:
            node = node.parent_node
        return node

    def species_distance(ancestor, descendant):
        distance = 0
        node = descendant
        while node is not ancestor:
            if node is None:
                raise ValueError("Species mapping is not below its inferred LCA")
            node = node.parent_node
            distance += 1
        return distance

    mappings = {}
    duplications = 0
    losses = 0
    for node in gene_tree.postorder_node_iter():
        if node.is_leaf():
            gene_id = node.taxon.label
            if gene_id not in gene_to_species:
                raise ValueError(f"No species mapping for gene-tree tip: {gene_id}")
            mappings[node] = species_leaves[gene_to_species[gene_id]]
            continue
        children = list(node.child_node_iter())
        mapped_node = mappings[children[0]]
        for child in children[1:]:
            mapped_node = species_lca(mapped_node, mappings[child])
        mappings[node] = mapped_node
        is_duplication = any(mappings[child] is mapped_node for child in children)
        duplications += int(is_duplication)
        for child in children:
            distance = species_distance(mapped_node, mappings[child])
            losses += distance if is_duplication else max(distance - 1, 0)
    return duplications, losses


def root_gene_tree_by_duplication_loss(
    gene_tree_fp,
    rooted_species_tree_fp,
    gene_to_species,
    output_tree_fp,
    scores_fp=None,
    excluded_gene_ids=(),
    duplication_weight=1.0,
    loss_weight=1.0,
):
    """Root a gene tree on every edge and minimize reconciliation cost.

    All tied minimum-cost edges are retained in the returned score table. The
    first optimum is written as a deterministic representative rooted tree.
    Genes excluded because their species is absent from the species tree are
    pruned before scoring and do not appear in the rooted output.
    """
    gene_tree = read_unrooted_tree(gene_tree_fp)
    excluded_gene_ids = set(excluded_gene_ids)
    retained_gene_ids = set(gene_to_species).difference(excluded_gene_ids)
    gene_tree.retain_taxa_with_labels(retained_gene_ids)
    observed_gene_ids = {
        node.taxon.label for node in gene_tree.leaf_node_iter()
    }
    if observed_gene_ids != retained_gene_ids:
        missing = retained_gene_ids.difference(observed_gene_ids)
        raise ValueError("Mapped genes absent from gene tree: " + ", ".join(missing))
    species_tree = dendropy.Tree.get(
        path=str(rooted_species_tree_fp),
        schema="newick",
        rooting="force-rooted",
        preserve_underscores=True,
    )
    all_genes = frozenset(observed_gene_ids)
    candidates = []
    candidate_sides = []
    for node in gene_tree.preorder_node_iter():
        if node is gene_tree.seed_node:
            continue
        side = frozenset(leaf.taxon.label for leaf in node.leaf_iter())
        complement = all_genes.difference(side)
        canonical_side = min((side, complement), key=lambda x: (len(x), sorted(x)))
        if canonical_side in candidate_sides:
            continue
        candidate_sides.append(canonical_side)
        rooted = gene_tree.clone(depth=2)
        target = next(
            candidate
            for candidate in rooted.preorder_node_iter()
            if candidate is not rooted.seed_node
            and frozenset(leaf.taxon.label for leaf in candidate.leaf_iter()) == side
        )
        edge_length = target.edge.length
        half_length = None if edge_length is None else edge_length / 2
        rooted.reroot_at_edge(
            target.edge,
            length1=half_length,
            length2=half_length,
            update_bipartitions=False,
        )
        duplications, losses = duplication_loss_reconciliation_score(
            rooted, species_tree, gene_to_species
        )
        candidates.append(
            {
                "root_edge": len(candidates) + 1,
                "root_side_size": len(canonical_side),
                "root_side_genes": ";".join(sorted(canonical_side)),
                "duplications": duplications,
                "losses": losses,
                "reconciliation_cost": (
                    duplication_weight * duplications + loss_weight * losses
                ),
                "_tree": rooted,
            }
        )
    scores = pd.DataFrame(candidates).sort_values(
        ["reconciliation_cost", "duplications", "losses", "root_edge"]
    ).reset_index(drop=True)
    minimum_cost = scores["reconciliation_cost"].min()
    scores["is_optimal_root"] = scores["reconciliation_cost"].eq(minimum_cost)
    best_tree = scores.iloc[0]["_tree"]
    output_tree_fp = Path(output_tree_fp)
    output_tree_fp.parent.mkdir(parents=True, exist_ok=True)
    best_tree.write(
        path=str(output_tree_fp),
        schema="newick",
        suppress_rooting=False,
        unquoted_underscores=True,
    )
    scores = scores.drop(columns="_tree")
    if scores_fp is not None:
        scores_fp = Path(scores_fp)
        scores_fp.parent.mkdir(parents=True, exist_ok=True)
        scores.to_csv(scores_fp, sep="\t", index=False)
    return output_tree_fp, scores


def parse_astral_local_posterior(label):
    """Extract ASTRAL's main-resolution local posterior probability (pp1)."""
    if label is None:
        return None
    match = re.search(
        r"(?:^|[;\[])pp1=([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)",
        str(label),
    )
    return float(match.group(1)) if match else None


def prepare_astral_display_tree(
    input_tree_fp,
    taxonomy_fp,
    output_tree_fp,
    support_decimals=0,
):
    """Create a plotting tree with organism tips and pp1 support percentages."""
    tree = dendropy.Tree.get(
        path=str(input_tree_fp),
        schema="newick",
        rooting="force-rooted",
        preserve_underscores=True,
    )
    taxonomy = pd.read_csv(taxonomy_fp, sep="\t", dtype=str, keep_default_na=False)
    required = {"assembly_accession", "organism_name"}
    missing_columns = required.difference(taxonomy.columns)
    if missing_columns:
        raise ValueError(
            "Taxonomy table is missing columns: "
            + ", ".join(sorted(missing_columns))
        )
    name_by_accession = taxonomy.set_index("assembly_accession")[
        "organism_name"
    ].to_dict()
    missing_accessions = []
    for leaf in tree.leaf_node_iter():
        accession = leaf.taxon.label
        if accession not in name_by_accession:
            missing_accessions.append(accession)
            continue
        leaf.taxon.label = name_by_accession[accession]
    if missing_accessions:
        raise ValueError(
            "No organism name for accessions: "
            + ", ".join(sorted(missing_accessions))
        )

    for node in tree.postorder_node_iter():
        if node.is_leaf() or node is tree.seed_node:
            continue
        posterior = parse_astral_local_posterior(node.label)
        node.label = (
            f"{100 * posterior:.{support_decimals}f}"
            if posterior is not None
            else None
        )

    output_tree_fp = Path(output_tree_fp)
    output_tree_fp.parent.mkdir(parents=True, exist_ok=True)
    tree.write(
        path=str(output_tree_fp),
        schema="newick",
        suppress_rooting=True,
    )
    return output_tree_fp
