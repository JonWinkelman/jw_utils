from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import subprocess
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import to_tree
from Bio import AlignIO


def make_genome_sketch(fp, out_dir):
    fp = Path(fp)
    outfp = out_dir / fp.with_suffix(".msh").name

    if outfp.exists():
        return fp, "skipped"

    subprocess.run(
        ["mash", "sketch", str(fp), "-o", str(outfp)],
        check=True,
        capture_output=True,
        text=True,
    )

    return fp, "done"


def make_genome_sketches(input_genome_dir, out_dir, threads=8):
    input_genome_dir = Path(input_genome_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(exist_ok=True, parents=True)

    genome_paths = list(input_genome_dir.glob("*.fna"))

    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [
            executor.submit(make_genome_sketch, fp, out_dir)
            for fp in genome_paths
        ]

        for future in tqdm(as_completed(futures), total=len(futures), desc="Sketching genomes"):
            try:
                fp, status = future.result()
            except subprocess.CalledProcessError as e:
                print(f"\nFailed: {e.cmd}")
                print(e.stderr)


def mash_output_to_dataframe(dist_mat_path):
    """Return bottom triangle df from bottom triangle tab-delimeted text file"""
    with open(dist_mat_path, 'r') as f: 
        triangle = f.readlines()
    col_names = []

    num_cols_rows  = int(triangle[0][1:-1])
    triangle = triangle[1:]
    dist_matrix = np.full([len(triangle),len(triangle)], np.nan)
    index_col_names = []
    for i, row in enumerate(triangle):
        row_list = row.strip().split('\t')
        for col in range(num_cols_rows):
            if col == 0:
                name =  row_list[0].split('/')[-1][0:15]
                index_col_names.append(name)
            else:
                if col < len(row_list):
                    dist_matrix[i][col-1] = row_list[col]
    df = pd.DataFrame(dist_matrix, columns =  index_col_names)
    df.index = index_col_names
    return df            

            



def _write_textfile_for_mash(sketch_fps_lst, out_fp="temp_sketch_file_names.txt"):
    out_fp = Path(out_fp)
    with out_fp.open("w") as f:
        for fp in sketch_fps_lst:
            f.write(f"{Path(fp)}\n")
    return out_fp


from pathlib import Path
import subprocess


def make_dist_matrix_from_sketches(
    dist_mat_out_fp,
    sketch_fps_lst=None,
    sketch_dir=None,
    threads=8,
    overwrite=False,
):
    """
    Generate a pairwise Mash distance matrix from either a list of .msh files
    or all .msh files in a directory.
    """

    dist_mat_out_fp = Path(dist_mat_out_fp)
    sketch_file_names_fp = Path('./temp_sketch_fp_namse.txt')

    if dist_mat_out_fp.exists() and not overwrite:
        print(f"Using existing distance matrix: {dist_mat_out_fp}")
        return mash_output_to_dataframe(dist_mat_out_fp)

    if sketch_fps_lst is not None and sketch_dir is not None:
        raise ValueError("Provide either sketch_fps_lst or sketch_dir, not both.")

    if sketch_fps_lst is None and sketch_dir is None:
        raise ValueError("Provide either sketch_fps_lst or sketch_dir.")

    if sketch_fps_lst is not None:
        sketch_fps = [Path(fp) for fp in sketch_fps_lst]
    else:
        sketch_dir = Path(sketch_dir)
        sketch_fps = sorted(sketch_dir.glob("*.msh"))

    if not sketch_fps:
        raise ValueError("No .msh sketch files found.")

    _write_textfile_for_mash(sketch_fps, sketch_file_names_fp)

    cmd = [
        "mash",
        "triangle",
        "-p",
        str(threads),
        "-l",
        str(sketch_file_names_fp),
    ]

    with dist_mat_out_fp.open("w") as out_f:
        subprocess.run(
            cmd,
            stdout=out_f,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
        )
    sketch_file_names_fp.unlink(missing_ok=True)
    print(f'Parsing {dist_mat_out_fp} to a dataframe...')
    df = mash_output_to_dataframe(dist_mat_out_fp)
    return df





def lower_triangle_to_full(distance_df):
    # Force writable numeric copy
    arr = distance_df.to_numpy(dtype=float, copy=True)

    np.fill_diagonal(arr, 0)

    # Convert lower triangle to full symmetric matrix
    lower = np.tril(arr)
    full_arr = lower + lower.T
    np.fill_diagonal(full_arr, 0)

    return pd.DataFrame(
        full_arr,
        index=distance_df.index,
        columns=distance_df.columns
    )


def single_linkage_clusters_scipy(distance_df, threshold=0.98):
    genomes = distance_df.index.to_numpy()
    if threshold >= 1: 
        raise ValueError(f"Threshold {threshold} is >= 1. It needs to be less than 1 ")
    dist_thresh = 1 - threshold

    arr = distance_df.to_numpy()

    mask = (arr <= dist_thresh) & ~np.isnan(arr)

    # Remove self-connections
    np.fill_diagonal(mask, False)

    graph = csr_matrix(mask)

    n_components, labels = connected_components(
        csgraph=graph,
        directed=False,
        return_labels=True
    )

    return dict(zip(genomes, labels))

def save_upgma_to_newick(distance_df, newick_output_fp): 
    Z, labels  = build_upgma_tree(distance_df)
    with open(newick_output_fp, "w") as f:
        f.write(linkage_to_newick(Z, labels))

def build_upgma_tree(distance_df):
    """
    Construct a UPGMA (average-linkage) hierarchical clustering tree from
    a lower-triangular pairwise distance matrix.

    Parameters
    ----------
    distance_df : pandas.DataFrame
        Square lower-triangular distance matrix with genomes/samples as both
        the index and columns. Values should represent pairwise distances
        (e.g., Mash distances). The upper triangle may contain NaN values.

    Returns
    -------
    Z : numpy.ndarray
        SciPy linkage matrix of shape (n_samples - 1, 4). This encodes the
        complete UPGMA hierarchical clustering tree and can be used with
        scipy.cluster.hierarchy functions such as:

        - dendrogram()
        - fcluster()
        - to_tree()

    labels : list[str]
        Sample labels corresponding to the leaf order used in the linkage
        matrix. These labels can be used to map cluster assignments back
        to genome names or export the tree to Newick format.

    Notes
    -----
    This function:

    1. Converts the lower-triangular distance matrix into a full symmetric
       distance matrix using `lower_triangle_to_full()`.
    2. Converts the full distance matrix into SciPy's condensed distance
       format using `squareform()`.
    3. Performs average-linkage hierarchical clustering via
       `linkage(..., method='average')`, which is equivalent to UPGMA.

    Examples
    --------
    >>> Z, labels = build_upgma_tree(distance_df)

    Generate flat clusters:

    >>> from scipy.cluster.hierarchy import fcluster
    >>> clusters = fcluster(Z, t=0.05, criterion='distance')

    Plot a dendrogram:

    >>> from scipy.cluster.hierarchy import dendrogram
    >>> dendrogram(Z, labels=labels)
    """
    full_df = lower_triangle_to_full(distance_df)

    condensed = squareform(full_df.values)

    Z = linkage(condensed, method="average")  # UPGMA / average linkage

    return Z, full_df.index.tolist()



def linkage_to_newick(Z, labels):
    """
    Convert a SciPy linkage matrix into a Newick-formatted tree string.

    Parameters
    ----------
    Z : numpy.ndarray
        SciPy linkage matrix produced by
        `scipy.cluster.hierarchy.linkage()`. This matrix encodes the full
        hierarchical clustering tree (e.g., a UPGMA tree generated using
        `method='average'`).

    labels : list[str]
        Leaf labels corresponding to the original observations used to
        construct the linkage matrix. The order of labels must match the
        order of samples used to generate `Z`.

    Returns
    -------
    str
        Tree in Newick format. Branch lengths correspond to the linkage
        distances stored in the hierarchy and are formatted to six decimal
        places.

    Notes
    -----
    The linkage matrix stores cluster heights rather than explicit branch
    lengths. This function converts cluster heights into branch lengths by
    subtracting each node's height from its parent's height:

        branch_length = parent_height - node_height

    The resulting Newick string can be visualized using many phylogenetic
    tree tools, including:

    - FigTree
    - iTOL
    - ETE Toolkit
    - Biopython Phylo

    Examples
    --------
    Build a UPGMA tree and export to Newick:

    >>> Z, labels = build_upgma_tree(distance_df)
    >>> newick = linkage_to_newick(Z, labels)

    Save to file:

    >>> with open("tree.nwk", "w") as f:
    ...     f.write(newick)

    Example output:

    >>> print(newick)
    ((genome_A:0.012345,genome_B:0.012345):0.034567,
     genome_C:0.046912);
    """
    tree = to_tree(Z, rd=False)

    def build_newick(node, parent_dist):
        if node.is_leaf():
            name = labels[node.id]
            branch_length = parent_dist - node.dist
            return f"{name}:{branch_length:.6f}"

        left = build_newick(node.left, node.dist)
        right = build_newick(node.right, node.dist)

        branch_length = parent_dist - node.dist
        return f"({left},{right}):{branch_length:.6f}"

    return build_newick(tree, tree.dist) + ";"




def protein_alignment_distance_matrix(alignment_fp):
    """
    Compute pairwise amino-acid distances from a protein multiple sequence alignment.

    Distance is calculated as:

        1 - amino_acid_identity

    ignoring positions where either sequence has a gap.
    """

    alignment = AlignIO.read(alignment_fp, "fasta")

    names = [record.id for record in alignment]
    seqs = np.array([list(str(record.seq)) for record in alignment])

    n = len(names)
    dist = np.zeros((n, n), dtype=float)

    for i in range(n):
        for j in range(i):
            s1 = seqs[i]
            s2 = seqs[j]

            valid = (s1 != "-") & (s2 != "-")

            if valid.sum() == 0:
                d = np.nan
            else:
                identity = (s1[valid] == s2[valid]).mean()
                d = 1 - identity

            dist[i, j] = d
            dist[j, i] = d

    return pd.DataFrame(dist, index=names, columns=names)



def make_fclusters(Z, labels, cluster_distance):
    """
    Convert a hierarchical clustering tree into flat clusters using a
    distance threshold.

    Parameters
    ----------
    Z : numpy.ndarray
        SciPy linkage matrix produced by
        `scipy.cluster.hierarchy.linkage()`, typically generated from a
        pairwise distance matrix using UPGMA (`method='average'`) or
        another hierarchical clustering method.

    labels : list[str]
        Labels corresponding to the leaves of the linkage tree. The order
        must match the order of samples used to generate `Z`.

    cluster_distance : float
        Maximum linkage distance allowed within a cluster. The dendrogram
        is cut at this height using `criterion='distance'`.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing cluster assignments with columns:

        - genome : genome/sample identifier
        - cluster : integer cluster ID

    Notes
    -----
    This function uses:

        fcluster(Z, t=cluster_distance, criterion='distance')

    to convert a hierarchical clustering tree into flat clusters.

    Genomes are grouped together only if their clusters merge at a
    linkage distance less than or equal to `cluster_distance` in the
    dendrogram.

    For UPGMA clustering, the linkage distance corresponds to the
    average distance between the members of merged clusters.

    Examples
    --------
    Create clusters from a UPGMA tree:

    >>> cluster_df = make_fclusters(
    ...     Z,
    ...     labels,
    ...     cluster_distance=0.03
    ... )

    Count cluster sizes:

    >>> cluster_df['cluster'].value_counts()

    Merge cluster assignments back into metadata:

    >>> summary_df = summary_df.merge(
    ...     cluster_df,
    ...     left_index=True,
    ...     right_on='genome'
    ... )
    """

    fclusters = fcluster(
        Z,
        t=cluster_distance,
        criterion="distance"
    )

    cluster_df = pd.DataFrame({
        "genome": labels,
        "cluster": fclusters
    })

    return cluster_df


def pick_cluster_medoid(distance_df, members):
    """
    Return the sequence closest to the center of a cluster.

    Parameters
    ----------
    distance_df : pd.DataFrame
        Full pairwise distance matrix.

    members : list
        Sequence IDs belonging to a cluster.
    """
    sub = distance_df.loc[members, members]
    mean_dist = sub.mean(axis=1)
    return mean_dist.idxmin()
    
