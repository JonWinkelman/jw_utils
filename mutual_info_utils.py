import pandas as pd
from sklearn.metrics import mutual_info_score


def check_freq_dict(freq_dict, char_dict, mol_type):
    """"""
    chars = char_dict[mol_type]
    freq_d_ed = {k.upper():val for k, val in freq_dict.items()}
    for char in chars:
        if char not in freq_d_ed.keys():
            freq_d_ed[char] = 0
    for char in freq_d_ed.keys():
        if char not in char_dict[mol_type]:
            raise KeyError(f'{char} not in set of {mol_type} characters: {char_dict[mol_type]}') 
    return freq_d_ed



def create_expectedfreq_matrix(freq_d1, freq_d2, moltype1='prot',moltype2='rna' ):
    """Return 20x4 matrix with expected frequencies p(x)p(y) of each pair, e.g. for a col in an alignment
    
     * dicts should contain the observed frequencies (marginal probabilities) of
     each character for a column in alignment.
     
    freq_d1 (dict): {character1:frequency1,...}
    freq_d2 (dict): {character1:frequency1, ...}
    """
    char_dict = {
        'prot': ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 
             'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V','-'],  # 20 standard amino acids
        'dna': ['A', 'T', 'C', 'G','-'],  # DNA nucleotides
        'rna': ['A', 'U', 'C', 'G','-']   # RNA nucleotides
    }
    ## Check dicts
    freq_d1 =  check_freq_dict(freq_d1, char_dict, moltype1)
    freq_d2 =  check_freq_dict(freq_d2, char_dict, moltype2)
    # Convert nucleotide and amino acid lists to sorted lists to ensure consistent order
    d1_chars = char_dict[moltype1]
    d2_chars = char_dict[moltype2]
    expected_matrix = np.zeros((len(char_dict[moltype1]), len(char_dict[moltype2])))

    # Fill the matrix with expected frequencies
    for i, d1_char in enumerate(d1_chars):
        for j, d2_char in enumerate(d2_chars):
            expected_matrix[i, j] = freq_d2[d2_char] * freq_d1[d1_char]
    df = pd.DataFrame(expected_matrix)
    df.index = d1_chars
    df.columns = d2_chars
    return df


def create_joint_probability_matrix(df, moltype1='prot',moltype2='rna' ):
    """
    
    df (pd.DataFrame): containing all observations, each row is a genome,  
                    column 1 id char state of molecule1, col 2 is character state
                    of moelcule2 in that genome
    """
    char_dict = {
        'prot': ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 
             'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V','-'],  # 20 standard amino acids
        'dna': ['A', 'T', 'C', 'G','-'],  # DNA nucleotides
        'rna': ['A', 'U', 'C', 'G','-']   # RNA nucleotides
    }
    all_pairs = pd.MultiIndex.from_product([char_dict[moltype1], char_dict[moltype2]], names=['molecule1', 'molecule2'])
    all_pairs_df = pd.DataFrame(index=all_pairs).reset_index()
    
    df.columns = ['molecule1', 'molecule2']
    pair_counts = df.groupby(['molecule1', 'molecule2']).size().reset_index(name='count')
    total_observations = len(df)
    pair_counts['joint_prob'] = pair_counts['count'] #/ total_observations
    joint_prob_df = pd.merge(all_pairs_df, pair_counts[['molecule1', 'molecule2', 'joint_prob']], how='left', on=['molecule1', 'molecule2'])
    joint_prob_df['joint_prob'] = joint_prob_df['joint_prob'].fillna(0)
    return  joint_prob_df.sort_values('joint_prob', ascending=False)