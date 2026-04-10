from matplotlib import pyplot as plt
import matplotlib as mpl






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