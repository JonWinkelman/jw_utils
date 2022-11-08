#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 09:29:33 2022

@author: jonwinkelman
"""
import pandas as pd
import bisect
from plotly import offline as pyo
from plotly import graph_objects as go
import numpy as np
import os
from jw_utils import app_functions as afns
from jw_utils import parse_gff as pgf
from jw_utils import parse_fasta as pf
from jw_utils import genome_utils as gu
from jw_utils import parse_gbk as pgb
from jw_utils import plotly_utils as pu
from jw_utils import file_utils as fu
class GeneProfile:
    """Class for mapping reads to features in a contig using a genbank file."""
    
    def __init__(self, path_to_genbank, path_to_value_counts):
        """Contig object for mapping reads to features with a genbank file."""
        self.contig_name = path_to_value_counts.split('/')[-1].split('_')[0]
        self.geneObject_dict = pgb.build_genbank_dict(path_to_genbank)[self.contig_name]
        #self.genome_size = self.get_genome_descriptors()['genome_size']
        self.genes = list(self.geneObject_dict.keys())
        self.path_to_value_counts = path_to_value_counts
        self.value_counts = self.read_value_counts()
        self.plus_strand_genes = [gene for gene in self.genes if self.geneObject_dict[gene].strand == 1]
        self.minus_strand_genes = [gene for gene in self.genes if self.geneObject_dict[gene].strand == -1]
        self.total_mapped_reads = self.total_mapped_reads()



    def read_value_counts(self):
        """Return counts of 3 prime ends at each genomoic position."""
        df = pd.read_csv(self.path_to_value_counts)
        df.columns = ['genomic_position', "5'_(+)", "5''_(-)"]
        df = df.fillna(0)
        df = df.set_index('genomic_position')
        return df#df.reindex(list(range(1,self.genome_size+1)),fill_value=0)



    def total_mapped_reads(self, *args, **kwargs):
        """Return total number of reads present in value counts df."""
        df = self.value_counts
        return df[df.columns[0]].sum() + df[df.columns[1]].sum()


    def NormalizeData(self, data):
        """Scale data between 0 and 1, unless all values are 0, then return all zeros."""
        if (np.max(data) - np.min(data)) == 0:
            return np.full(len(data),0)
        else:
            return (data - np.min(data)) / (np.max(data) - np.min(data))


    def gene_hits(self, gene_name, start_offset=100, stop_offset=100, reads_per_mil=False, normalize=False):
        """Map reads from value counts file to a gene region."""
        df = self.value_counts
        gene_obj = self.geneObject_dict[gene_name]
        if gene_obj.strand == 1 :
            start_codon = gene_obj.location.start
            stop_codon = gene_obj.location.end
            start_index = bisect.bisect_left(df.index,(start_codon - start_offset))
            end_index = bisect.bisect_right(df.index, (stop_codon + stop_offset)) 
            df2 = df.iloc[start_index:end_index,0]

        elif gene_obj.strand == -1:
            start_codon = gene_obj.location.end
            stop_codon = gene_obj.location.start
            start_ind = bisect.bisect_left(df.index,(stop_codon - stop_offset))
            end_ind = bisect.bisect_right(df.index, (start_codon + start_offset)) 
            df2 = df.iloc[start_ind:end_ind,1]
    
            
        if reads_per_mil:
            mil_reads = self.total_mapped_reads/1000000
            df2 = df2.div(mil_reads)
        if normalize:
            if (np.max(df2) - np.min(df2)) != 0:
                df2 = self.NormalizeData(df2)

        return pd.DataFrame(df2), start_codon, stop_codon
    
    
    def plot_gene(self, gene, normalize = False, plot = True, reads_per_mil=False,
                  vert_shift=-1):
        """
        Return df, start and stop codon pos, plot gene hits with the gene arrow.
        
        Arguments:
        gene (str)
        normalize (bool)
        plot (bool)
        reads_per_mil (bool)
        """
        gene_obj = self.geneObject_dict[gene]
        if gene_obj.strand == -1: strand = '-'
        else: strand = '+'
        df2, start_codon, stop_codon = self.gene_hits(gene, normalize = normalize, reads_per_mil=reads_per_mil)
        if plot:
            max_val = df2.loc[:,df2.columns[0]].max()
            fig = go.Figure()
            fig.add_trace(go.Scatter(
                x=df2.index,
                y=df2.loc[:,df2.columns[0]],
                name = gene))
            fig.add_trace(go.Scatter(
                x=[start_codon,start_codon], 
                y=[0,max_val],
                line=dict(color="rgb(0,255,0)"),
                name = 'Start'))
            fig.add_trace(go.Scatter(
                x=[stop_codon,stop_codon], 
                y=[0,max_val],
                line=dict(color="rgb(255,0,0)"),
                name = 'Stop'))
            x,y = pu.make_gene_arrow_coords(start_codon, stop_codon, height = 1,vert_shift = vert_shift)
            fig.add_trace(go.Scatter(
                x = x, 
                y = y,
                line=dict(color="rgb(100,100,100)"),
                name = 'strand: '+strand))
            pyo.plot(fig)
        return df2, start_codon, stop_codon
    
    
    
    def graph_gene_region(self, feature_id, up_input, down_input,path_annot_file,
                          annot_type='gbk', contig_name=None):
        """Plot the reads in region of the input feature with upstream and downstream genes."""
        fig = go.Figure()
        #if up_input or down_input:     
        #gff_obj=pgf.make_seq_object_dict(path_to_gff, feature_type='CDS')
        
        start_offset, stop_offset, trimmed_df = afns.get_offsets(feature_id, up_input, 
                                            down_input, path_annot_file, annot_type=annot_type,
                                            contig_name=contig_name)
        l = []
            
        df, start, stop = self.gene_hits(feature_id,
                                 reads_per_mil=True,start_offset=start_offset,
                                 stop_offset=stop_offset)
        fig.add_trace(go.Bar(
            x=df.index,
            y=df.loc[:,df.columns[0]],
            name = f'{contig_name}\n{feature_id}'
            ))
        l.append(df.loc[:,df.columns[0]].max())
        max_val = max(l)
        fig.add_trace(go.Scatter(
            x=[stop,stop],
            y=[0,max_val*1.2],
            name = 'stop',
            mode = 'lines',
            line = {'width': 3,
                    'dash':'dash',
                    'color':'rgba(200,100,100,0.5)'}
            ))
        fig.add_trace(go.Scatter(
            x=[start,start],
            y=[0,max_val*1.2],
            name = 'start',
            mode = 'lines',
            line = {'width': 3,
                    'dash':'dash',
                    'color':'rgba(100,200,100,0.5)'}
            ))
        for ind in trimmed_df.index:
            if trimmed_df.loc[ind,'strand'] == '-':
                stop = trimmed_df.loc[ind,'start']
                start = trimmed_df.loc[ind,'end']
            else:
                start = trimmed_df.loc[ind,'start']
                stop = trimmed_df.loc[ind,'end']   
            fig.add_trace(go.Scatter(
                x = pu.make_gene_arrow_coords(start, stop, max_val/10, -max_val/8)[0],
                y = pu.make_gene_arrow_coords(start, stop, max_val/10, -max_val/8)[1],
                fill="toself",
                mode = 'lines',
                line = {'width': 1.5,
                        'color':'rgb(200,200,200)',
                        },
                name=feature_id,
                showlegend=False,
                #text = f'{gff_obj[ind].ID}<br>{gff_obj[ind].locus_tag}<br>{gff_obj[ind].Name}<br>{gff_obj[ind].product}'
                ))
        fig = afns.update_fig_style(fig,
                   title='Gene (selected from termination score graph)',
                   ylabel = "Reads (per million)",
                   xlabel = "Genomic position"
                       )
        return fig
