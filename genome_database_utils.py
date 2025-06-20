from ete3 import Tree
import pandas as pd


def prune_tree_to_rank(db_tree_taxonomy_df, ete3_tree, rank='phylum'):
    """
    Prunes an ETE3 tree to one representative genome per specified taxonomic rank 
    (e.g., phylum) and renames tree leaf nodes to rank labels.

    Parameters:
        db_tree_taxonomy_df (pd.DataFrame): A DataFrame indexed by NCBI genome assembly 
            accessions (length-15 strings starting with "GC"), with taxonomy information.
        ete3_tree (ete3.Tree): An ETE3 Tree object whose leaves are genome accessions.
        rank (str): The taxonomic rank to prune by (e.g., "phylum", "class"). This rank name
            needs to match a db_tree_taxonomy_df column name. 

    Returns:
        pd.DataFrame: A filtered version of db_tree_taxonomy_df containing one genome per rank,
            matching leaves present in the tree.
    """
    # Validate index format
    accession_index = db_tree_taxonomy_df.index
    if not all(accession_index.str.len() == 15):
        raise ValueError("All DataFrame indices must be 15 characters long (NCBI accessions).")
    if not all(accession_index.str.startswith("GC")):
        raise ValueError("All DataFrame indices must start with 'GC' (NCBI accessions).")

    # Standardize tree leaf names to last 15 characters (genome accession)
    # GTDB often adds 'GB_' or 'RS_' to beginning of accession
    for clade in ete3_tree:
        clade.name = clade.name[-15:]

    # Get accessions present in both tree and taxonomy DataFrame
    tree_leaves = {leaf.name for leaf in ete3_tree}
    rank_df = db_tree_taxonomy_df.dropna(subset=[rank]).drop_duplicates(subset=[rank])
    accessions_to_keep = list(tree_leaves.intersection(rank_df.index))

    # Subset the taxonomy dataframe to only those kept in the tree
    rank_pruned_df = rank_df.loc[accessions_to_keep]

    # Prune tree to one genome per rank and rename leaves to rank name
    ete3_tree.prune(accessions_to_keep, preserve_branch_length=True)
    rank_map = rank_pruned_df[rank].to_dict()
    for clade in ete3_tree:
        clade.name = rank_map.get(clade.name, clade.name)

    return rank_pruned_df
