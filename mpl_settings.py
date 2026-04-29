from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np
from matplotlib.colors import to_hex

def get_colors(num_colors, pallette='viridis', print_available_pallets=False, hexcolors=False):
    """
    Heatmaps / scores      → 'viridis' or 'cividis'
    Enrichment / depletion → 'RdBu' or 'coolwarm'
    Phylogeny categories   → 'tab10' or 'Set3'
    """
    if print_available_pallets:
        print(plt.colormaps())
    cmap = plt.get_cmap(pallette)
    colors = cmap(np.linspace(0, 1, num_colors))
    if hexcolors:
        return [to_hex(c, keep_alpha=True) for c in colors]
    return colors



def show_palette(name, n=10):
    
    cmap = plt.get_cmap(name)
    colors = cmap(np.linspace(0, 1, n))
    
    fig, ax = plt.subplots(figsize=(n, 1))
    
    for i, c in enumerate(colors):
        ax.add_patch(plt.Rectangle((i, 0), 1, 1, color=c))
    
    ax.set_xlim(0, n)
    ax.set_ylim(0, 1)
    ax.axis('off')
    plt.title(name)



def update_settings(): 
    mpl.rcParams['svg.fonttype'] = 'none'
    plt.rcParams.update({
        "font.size": 8,
        "axes.linewidth": 0.8,
        "axes.labelsize": 8,
        "axes.titlesize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "lines.linewidth": 1,
        "figure.dpi": 300,
        'svg.fonttype':'none'
    })

def format_subplots(axes):
    for ax in axes.flatten():
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(direction='out', length=3, width=0.8)

def format_ax(ax):   
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(direction='out', length=3, width=0.8)


def make_jitterplot(ax): 
    """"""
    
    return None