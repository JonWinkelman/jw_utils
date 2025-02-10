import plotly.graph_objects as go
from jw_utils import parse_gff as pgf
import os
import pandas as pd
from scipy import stats
import numpy as np

def create_merged_genecounts_df(fps, file_suffix='.txt', substring_to_remove='.bam'):
    "Create a merged df from all sample genecounts txt files"
    dfs = []
    for fp in fps:
        df = pd.read_csv(fp, skiprows=1, sep='\t').iloc[:,[0,6]]
        cols = list(df.columns)
        cols[1] = cols[1].replace(substring_to_remove, '')
        df.columns = cols
        dfs.append(df.set_index('Geneid'))
    merged_df = pd.concat(dfs, axis=1)
    return merged_df
    

def reads_per_mil_normalization(df):
    """divides each count by the (total reads in millions)

    e.g. if there are 2 million total reads in the sample, 
    then each count, e.g a count of 12, would be divided by 2, 6 reads per million total reads
    """
    df_norm=pd.DataFrame()
    for col in df:
        df_norm[col] = df[col]/(df.sum()/1000000)[col]
    return df_norm


def make_comparison_df_dict(samp_info_df, contrast_df, 
                            subread_dir='./results/subread', 
                            file_suffix = '_T1.featureCounts.txt',
                           ):
    """genereates a dict of dfs where each key contains the samples to be compared
    
    
    samp_info_df (df): two cols: ['sample', 'group']
                        "sample" values should match a featureCount file once the file_suffix has been removed
                        featureCount file name example: far_red_1_T1.featureCounts.txt in the subread dir
                        'group' variable groups replicates in the sample column. E.g. far_red contains
                        two replicate samples, far_red_1 and far_red_2
        e.g.
        sample	    group
        ___________________
        far_red_1	far_red
        far_red_2	far_red
        dark_1	    dark
        dark_2	    dark
        dark_3	    dark

    contrast_df: communicates which groups to compare


    id                 variable	     up	  reference
    _______________________________________________
    far_red-vs-dark     group	  far_red   dark
                        
    """
    comp_dfs = {}
    for index, row in contrast_df.iterrows():
        up_samples = list(samp_info_df.groupby('group').get_group(row['up'])['sample'])
        ref_samples =  list(samp_info_df.groupby('group').get_group(row['reference'])['sample'])
        up_sample_fps = [os.path.join(subread_dir, f'{f}{file_suffix}') for f in up_samples]
        ref_sample_fps = [os.path.join(subread_dir, f'{f}{file_suffix}') for f in ref_samples]
        up_df = reads_per_mil_normalization(create_merged_genecounts_df(up_sample_fps) )
        ref_df = reads_per_mil_normalization(create_merged_genecounts_df(ref_sample_fps)  )
        name=row['id']
        comp_dfs[name] = up_df,ref_df
    return comp_dfs
    

def log2_fold_change(up_vals, reference_vals, fraction=0.001):
    """
    Compute log2 fold change with pseudocount as a fraction of the minimum nonzero value.
    
    Parameters:
    - reference_vals (array-like): Baseline (control) expression values
    - up_vals (array-like): Experimental expression values
    - fraction (float): Fraction of the minimum nonzero value to use as pseudocount
    
    Returns:
    - log2 fold change values (numpy array)
    """
    reference_vals = np.array(reference_vals, dtype=np.float64)
    up_vals = np.array(up_vals, dtype=np.float64)
    
    # Find the minimum nonzero value across both arrays
    min_nonzero = np.min(np.concatenate([reference_vals[reference_vals > 0], up_vals[up_vals > 0]]))
    
    # Define the pseudocount as a fraction of the minimum nonzero value
    pseudocount = fraction * min_nonzero
    
    # Apply pseudocount to avoid division by zero
    adjusted_reference = reference_vals + pseudocount
    adjusted_up = up_vals + pseudocount
    
    # Compute log2 fold change
    log2_fc = np.log2(adjusted_up / adjusted_reference)

    return log2_fc


def make_stats_df(comp_dfs_d):
    """
    *Up df is first in comp_dfs_d.values() tuple, ref df is the second.
    """
    stat_dfs = {}
    for name, dfs in comp_dfs_d.items():
        up = name.split('-vs-')[0]
        ref = name.split('-vs-')[1]
        df_stat = pd.DataFrame()
        up_colname, ref_colname = f'{up}_mean_rpm', f'{ref}_mean_rpm'
        df_stat[up_colname] =  dfs[0].apply(lambda x: x.sum()/len(x), axis=1)
        df_stat[ref_colname] =  dfs[1].apply(lambda x: x.sum()/len(x), axis=1)
        up_means = dfs[0].apply(lambda x: list(x), axis=1).to_dict()
        ref_means = dfs[1].apply(lambda x: list(x), axis=1).to_dict()
        pval_d= {}
        for (up_gene, up_means), (ref_gene, ref_means) in zip(up_means.items(),ref_means.items() ):
            t_stat, p_value  = stats.ttest_ind(up_means, ref_means)
            pval_d[up_gene]= p_value
        pval_df = pd.DataFrame.from_dict(pval_d, orient='index', columns = ['pvals'])
        df_stat = df_stat.join(pval_df)
        df_stat['log2_FC'] = log2_fold_change(df_stat[up_colname], df_stat[ref_colname], fraction=0.0001)
        df_stat['-log10_pval'] = -np.log10(df_stat['pvals'])
        stat_dfs[name] = df_stat
    return stat_dfs
 


def make_text_annot(df, gff_fp, feature_type='gene'):
    """
    Return a list of plotly annotation strings derived from stats_df and gff annotation file.

    df (pd.DataFrame): with first 4 cols: up mean, ref mean, log2_FC, and -log10_pval. Index is gene name
    index of df should be a gene name witht the same format used in the gff file, e.g. gene-<geneid>
    
    """

    feature_names = list(df.index)
    if not feature_names[0].startswith(feature_type):
        df.index = [f'{feature_type.lower()}-' + i for i in df.index]
        feature_names = list(df.index)
    up_colname = df.columns[0]
    up_means = list(df[up_colname])
    ref_colname = df.columns[1]
    ref_means = list(df[ref_colname])
    
    annot_df = pgf.make_simple_annot_df(gff_fp)

    
    annot_dict = annot_df.apply(lambda x: list(x), axis=1).to_dict()
    text_annots = []
    for feat, up_mean, ref_mean in zip(feature_names, up_means, ref_means):
        s=""
        up_mean=round(up_mean, 2)
        ref_mean=round(ref_mean, 2)
        annots = annot_dict.get(feat, ['nan','nan','nan'])
        annot_str = f'<br>gene:{annots[1]}<br>protein:{annots[0]}<br>product:{annots[2]}<br>{up_colname}:{up_mean}<br> {ref_colname}: {ref_mean}'
        text_annots.append(annot_str)
    
    return text_annots



def get_simple_volc_trace(xvals, yvals, 
                          text=None, 
                          name=None, 
                          marker_size=2, 
                          marker_color='rgba(100,100,100,0.2)',
                         **kwargs):
    
    return go.Scatter(x=xvals, y=yvals, text=text, mode='markers',
                      marker={'size':marker_size, 'color':marker_color}, name=name,**kwargs)

# def _get_significant_gene_filt(df, pval_cutoff=0.05, FC_cutoff=5): 
#     fc_filt = (df.iloc[:, 2] > np.log2(FC_cutoff)) | (df.iloc[:, 2] < np.log2(1/FC_cutoff))
#     sig_vals_filt = fc_filt & (df.iloc[:, 3] > (-np.log10(pval_cutoff))  )
#     return sig_vals_filt


def _get_significant_gene_filt(df, pval_cutoff=0.05, FC_cutoff=5):
    """Return two filters, up and down regulated genes

    df (pd.DataFrame): with first 4 cols: up mean, ref mean, log2_FC, and -log10_pval. Index is gene name
    index of df should be a gene name witht the same format used in the gff file, e.g. gene-<geneid>
    """
    
    pval_filt = df.iloc[:, 3] > -np.log10(pval_cutoff)
    up_regulated = pval_filt & (df.iloc[:, 2] > np.log2(FC_cutoff))
    down_regulated = pval_filt & (df.iloc[:, 2] < -np.log2(FC_cutoff))
    return up_regulated, down_regulated



def get_significant_genes_df(df, pval_cutoff, FC_cutoff, gff_fp=None):
    """
    
    df (pd.DataFrame): with first 4 cols: up mean, ref mean, log2_FC, and -log10_pval. Index is gene name
    """
      
    up_f, down_f = _get_significant_gene_filt(df, pval_cutoff,FC_cutoff)
    sig_df = df[up_f | down_f]
    if gff_fp:
        return sig_df.join(pgf.make_simple_annot_df(gff_fp))
    else:
        return sig_df

def get_volc_traces(df, text_annots=None, gff_fp=None, pval_cutoff = 0.05, FC_cutoff=2,
                   up_reg_color='blue',down_reg_color='red', marker_color_nonsig='rgba(100,100,100,0.3)', marker_size_sig=4,):
    """
    df (pd.DataFrame): with first 4 cols: up mean, ref mean, log2_FC, and -log10_pval. Index is gene name 
    
    In 0 indexing for iloc, FC_col=2, -log10_pval=3
    """


    #get_significant_genes
    up_genes_filt, down_genes_filt = _get_significant_gene_filt(df, pval_cutoff, FC_cutoff)
    df_up = df[up_genes_filt]
    df_down = df[down_genes_filt]
    up_xvals, up_yvals = list(df_up.iloc[:, 2]), list(df_up.iloc[:, 3])
    down_xvals, down_yvals = list(df_down.iloc[:, 2]), list(df_down.iloc[:, 3])
    if text_annots:
        up_text_annots = make_text_annot(df_up, gff_fp)
        down_text_annots = make_text_annot(df_down, gff_fp)
        
    up_sig_volcano_trace = get_simple_volc_trace(up_xvals, up_yvals, text=up_text_annots, name='up_reg_genes', 
                                              marker_color=up_reg_color, marker_size=marker_size_sig)

    down_sig_volcano_trace = get_simple_volc_trace(down_xvals, down_yvals, text=down_text_annots, name='down_reg_genes', 
                                              marker_color=down_reg_color, marker_size=marker_size_sig)

    # get non-significant trace
    df_nonsignificant = df[~(up_genes_filt | down_genes_filt)]
    xvals, yvals = list(df_nonsignificant.iloc[:, 2]), list(df_nonsignificant.iloc[:, 3])
    nonsig_volcano_trace = get_simple_volc_trace(xvals, yvals, text=None, name='non significant genes',
                                                marker_color=marker_color_nonsig, hoverinfo="skip")

   

    return up_sig_volcano_trace, down_sig_volcano_trace,  nonsig_volcano_trace



def add_cutoff_traces(fig, pval_cutoff, FC_cutoff, row, col):
    """add cutoff traces to a plotly figure"""
    
    yline= [-np.log10(pval_cutoff), -np.log10(pval_cutoff)]
    VR_xline= [np.log2(FC_cutoff), np.log2(FC_cutoff)]
    VL_xline= [-np.log2(FC_cutoff), -np.log2(FC_cutoff)]
    
    htrace = go.Scatter(y=yline, x=fig.layout.xaxis.range, mode='lines', line={'color':'rgba(50,0,0,0.1)'}, showlegend=False)
    VR_trace = go.Scatter(y=fig.layout.yaxis.range, x=VR_xline, mode='lines', line={'color':'rgba(50,0,0,0.1)'}, showlegend=False)
    VL_trace = go.Scatter(y=fig.layout.yaxis.range, x=VL_xline, mode='lines', line={'color':'rgba(50,0,0,0.1)'}, showlegend=False)
    fig.add_trace(htrace, row, col)
    fig.add_trace(VR_trace, row, col)
    fig.add_trace(VL_trace, row, col)
    return fig

    
