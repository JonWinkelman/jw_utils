#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 08:26:06 2022

@author: jonwinkelman
"""
import pandas as pd
import numpy as np
import logomaker


class LogoMatrix():
    ''' 
    generates LOGO matrices from aligned seqeunce data
    '''
    def __init__(self, aligned_sequences, aln_length = None, char_type = 'AA'):
        '''
        Parameters:
        aligned_sequences (Series or list or other iterable): contains aligned sequences. if
        all alignments are not the same lenght, the seq will be padded on the right with dashes
            
        aln_length (int): number of characters you want to include in the processing, 
        starting at the beginning of the sequence   
        '''
        self.char_type = char_type
        if aln_length:
            self.seq_length = aln_length
            
        else:
            self.seq_length = len(aligned_sequences[0])
        self.aligned_sequences = aligned_sequences
        self.padded_seqs = self.pad_seq()
        if char_type == 'AA':
            self.freq_matrix = self.get_freq_matrix_AA()
        if char_type == 'nt':
            self.freq_matrix = self.get_freq_matrix_nucleotide()
            
        
    def pad_seq(self):
        seq_add_dash = []
        for seq in self.aligned_sequences:
            seq_add_dash.append(seq.ljust(self.seq_length, '-'))
        return seq_add_dash

    
    #generate a df matrix for bits at each position
    def get_bit_matrix_nucleotide(self):   
        def observed_uncertainty():
            '''generate a matrix conatianing the uncertainty (entropy) associated with each position'''
            obs_uncertainty_matrix = -(self.freq_matrix['A']*np.log2(self.freq_matrix['A']) + 
                                        self.freq_matrix['T']*np.log2(self.freq_matrix['T']) + 
                                        self.freq_matrix['C']*np.log2(self.freq_matrix['C']) +
                                        self.freq_matrix['G']*np.log2(self.freq_matrix['G']) )
            return obs_uncertainty_matrix
        
        
        max_uncertainty = np.log2(4)  #4 is the number of symbols available, e.g protein would be 20
        change_in_uncertainty =  max_uncertainty - observed_uncertainty()
        bit_matrix_df = pd.DataFrame(columns = ['A','T','C','G'])
        bit_matrix_df['A'] = self.freq_matrix['A'] * change_in_uncertainty
        bit_matrix_df['T'] = self.freq_matrix['T'] * change_in_uncertainty
        bit_matrix_df['C'] = self.freq_matrix['C'] * change_in_uncertainty
        bit_matrix_df['G'] = self.freq_matrix['G'] * change_in_uncertainty
        
        return bit_matrix_df
    
    
    
    #generate a df matrix for bits at each position
    def get_bit_matrix_AA(self): 
        AAs = ['G','P','A','V','L','I','M','C','F','Y','W','H','K','R','Q','N','E','D','S','T']
        def observed_uncertainty():
            '''generate a matrix conatianing the uncertainty (entropy) associated with each position'''

            obs_uncertainty = 0
            for AA in AAs:
                t =  -(self.freq_matrix[AA]*np.log2(self.freq_matrix[AA]))
                obs_uncertainty = obs_uncertainty + t
            return obs_uncertainty
                
        
        max_uncertainty = np.log2(20)  #4 is the number of symbols available, e.g protein would be 20
        change_in_uncertainty =  max_uncertainty - observed_uncertainty()
        bit_matrix_df = pd.DataFrame(columns = AAs)
        for AA in AAs:
            bit_matrix_df[AA] = self.freq_matrix[AA] * change_in_uncertainty
        return bit_matrix_df
    

        

    def get_freq_matrix_nucleotide(self):
    
        ''' takes aligned sequences in pandas series and returns matrix for generating logo '''

        logo_matrix = pd.DataFrame()
        A_fraction = []
        T_fraction = []
        C_fraction = []
        G_fraction = []
        #dash_fraction = []

        # top loop corresponds to the length of a sequence.
        for i in range(0, self.seq_length):
            pos_x = []
            #nested loop appends the ith nucleotide from each sequence to the pos_x list
            #each loop
            for seq in self.padded_seqs:
                seq = seq.upper()
                pos_x.append(seq[i])   
            A_fraction.append(pos_x.count('A')/len(pos_x) )
            T_fraction.append(pos_x.count('T')/len(pos_x) )
            C_fraction.append(pos_x.count('C')/len(pos_x) )
            G_fraction.append(pos_x.count('G')/len(pos_x) )
            #dash_fraction.append(pos_x.count('-')/len(pos_x) )
        logo_matrix = pd.DataFrame()
        logo_matrix['A'] = A_fraction
        logo_matrix['T'] = T_fraction
        logo_matrix['C'] = C_fraction
        logo_matrix['G'] = G_fraction
        #logo_matrix['-'] = dash_fraction
        logo_matrix.index = range(1,self.seq_length+1)
        return logo_matrix
    
    
    def get_freq_matrix_AA(self):
    
        ''' takes aligned sequences in pandas series and returns matrix for generating logo '''
        AAs =['G','P','A','V','L','I','M','C','F','Y','W','H','K','R','Q','N','E','D','S','T']
        logo_matrix = pd.DataFrame()
        AA_fractions = {AA:[] for AA in AAs}
        #dash_fraction = []

        # top loop corresponds to the length of a sequence.
        for i in range(0, self.seq_length):
            pos_x = []
            #nested loop appends the ith nucleotide from each sequence to the pos_x list
            #each loop
            for seq in self.padded_seqs:
                seq = seq.upper()
                pos_x.append(seq[i])  
            for AA in AAs: 
                AA_fractions[AA].append( pos_x.count(AA)/len(pos_x) )
            #dash_fraction.append(pos_x.count('-')/len(pos_x) )
        logo_matrix = pd.DataFrame()   
        for AA in AAs:
            logo_matrix[AA] = AA_fractions[AA]
        logo_matrix.index = range(1,self.seq_length+1) #change to 1-based numbering
        return logo_matrix    
    

    
    def print_logo(self, logo_type = 'bit', title= '') :
        """
        Generate a logo from a matrix from the logos class attribute
        
        Parameters:
        logo_type (str): choose from 'bit' or 'freq' 
        char_type  (str): choose from 'AA' or 'nt'       
            
        """
        logo_type_options = ['bit', 'freq']
        char_type_options = ['AA', 'nt']
        if self.char_type not in char_type_options:
            raise Exception(f'{self.char_type} is not an option')
        if logo_type not in logo_type_options:
            raise Exception(f'{logo_type} is not an option')
            
        if self.char_type == 'AA' and logo_type == 'freq':
            freq_matrix = self.get_freq_matrix_AA()
            title = f'AA Logo {title}'
            ylabel = 'Freq'
        elif self.char_type == 'nt' and logo_type == 'freq':
            freq_matrix = self.get_freq_matrix_nucleotide()
            title = f'nt Logo {title}'
            ylabel = 'Freq'
        elif self.char_type == 'AA' and logo_type == 'bit':
            #freq_matrix = self.get_bit_matrix_AA()
            title = f'AA Logo {title}'
            ylabel = 'Bits'
            freq_matrix = logomaker.alignment_to_matrix(self.padded_seqs, to_type ='information' )
        elif self.char_type == 'nt' and logo_type == 'bit':
            freq_matrix = self.get_bit_matrix_nucleotide()
            freq_matrix = logomaker.alignment_to_matrix(self.padded_seqs, to_type ='information' )
            
            title = f'nt Logo {title}' 
            ylabel = 'Bits'
        
        
            
        logo = logomaker.Logo(
                        freq_matrix,
                        shade_below=.5,
                        fade_below=.5,
                        font_name='Arial Rounded MT Bold'
                        )
        logo.ax.set_ylabel(ylabel, fontsize = 20)
        logo.ax.set_title(title, fontsize = 20)
        logo.ax.set_xticks(range(1,self.seq_length+1))
            
            

            
