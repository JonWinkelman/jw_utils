#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 16:13:54 2022

@author: jonwinkelman
"""
import numpy as np
import plotly.offline as pyo


def normalize_data(data):
    if (np.max(data) - np.min(data)) == 0:
        return np.full(len(data),0)
    else:
        return (data - np.min(data)) / (np.max(data) - np.min(data))
        
def linear_rescale(maximum, num_list):
    norm_data = normalize_data(num_list)
    norm_data = norm_data * maximum
    return norm_data


def linear_func(m,x,b):
    return m*x +b


def exp_func(x, a, b, c):
    return a * np.exp(-b * x) + c