from ete3 import ncbi_taxonomy
ncbi_tax = ncbi_taxonomy.NCBITaxa()
from ete3 import Tree
from Bio import Phylo



def get_lineage_rank_dict(taxID):
    """Return a lineage dict {rank:taxid} from a single taxid
    
    ranks = 'no rank', 'superkingdom','family', 'genus', 'phylum', 'class', 'order', 'species', 
    
    taxid (int): ncbi taxonomic ID that can be parsed by the ete3 library
    """
    d = ncbi_tax.get_rank(ncbi_tax.get_lineage(taxID))
    return  {key:val for val, key in d.items()} 




def get_list_of_all_nodes(subtree, clades=None):
    """Get a list of all clade objects, including leaves, in a Bio.Phylo tree"""
    if not isinstance(subtree, (Phylo.Newick.Tree, Phylo.Newick.Clade)):
        raise TypeError(f'object entered needs to be of type {Phylo.Newick.Tree} or {Phylo.Newick.Clade}, you entered {type(subtree)}' )
        
    if clades is None:
        clades=[]
    for cl in subtree.root:
        clades.append(cl)
        if cl.is_terminal():
             clades.append(cl)
        else:
            get_list_of_all_nodes(cl, clades)
    return clades
