#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 15:04:36 2022

@author: jonwinkelman
"""
import os
    
    
def listdir(path_to_dir, 
                            begin_to_keep = None, 
                            end_to_keep = None,
                            end_to_exclude = None,
                            beg_to_exclude = None):
    '''
    return full filtered filenames in the given directory
    
    parameters:
        path_to_dir (str): path to the directory of interest
        
        ext_to_keep (str): files ending with this substring will be kept

    '''
    if end_to_keep:
        return [file for file in os.listdir(path_to_dir) if file.endswith(end_to_keep)]

    elif begin_to_keep:
        return [file for file in os.listdir(path_to_dir) if file.startswith(begin_to_keep)]
         
    elif beg_to_exclude:
        return [file for file in os.listdir(path_to_dir) if not file.startswith(beg_to_exclude)]
        
    elif end_to_exclude:
        return [file for file in os.listdir(path_to_dir) if not file.endswith(end_to_exclude)]
        
    else:
        return [file for file in os.listdir(path_to_dir)]
        
        

def get_filepaths_in_dir(path_to_dir, 
                            begin_to_keep = None, 
                            end_to_keep = None,
                            end_to_exclude = None,
                            beg_to_exclude = None):
    
    '''
    return path from working directory to all files in the given directory
    
    parameters:
        path_to_dir (str): path to the directory of interest
        
        ext_to_keep (str): files ending with this substring will be kept

    '''
    files = listdir(path_to_dir, 
                            begin_to_keep = begin_to_keep, 
                            end_to_keep = end_to_keep,
                            end_to_exclude = end_to_exclude,
                            beg_to_exclude = beg_to_exclude)
                            
    return [os.path.join(path_to_dir, file) for file in files]


                            
    
    
    
    



