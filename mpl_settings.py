from matplotlib import pyplot as plt
import matplotlib as mpl
import numpy as np
from matplotlib.colors import to_hex
import numpy as np
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
from scipy.stats import norm

def make_gaussian_PDF(
    data,
    ax,
    histo_bins,
    n_components,
    random_state=1,
    linewidth=1,
    fitline_color='grey',
    histcolor='grey',
    fitline='-',
    plot_total=True, 
    rbarwidth=0.8,
):
    """
    Fit a Gaussian Mixture Model (GMM) to 1D data, plot the histogram and
    component probability density functions (PDFs), and return the total PDF.

    Parameters
    ----------
    data : array-like, shape (n_samples,) or (n_samples, 1)
        Input data. If 1D, it will be reshaped to (n_samples, 1) for fitting.
    ax : matplotlib.axes.Axes
        Axes object on which to plot the histogram and PDFs.
    n_components : int
        Number of Gaussian components in the mixture model.
    random_state : int, default=1
        Random seed for reproducibility of the GMM fit.
    linewidth : float, default=1
        Line width for plotted PDFs.
    fitline_color : str, default='grey'
        Color used for plotting individual component PDFs (and total PDF if enabled).
    histcolor : str, default='grey'
        Color of the histogram.
    fitline : str, default='-'
        Line style passed to `ax.plot` (e.g., '-', '--').
    plot_total : bool, default=True
        If True, plot the summed mixture PDF in addition to individual components.

    Returns
    -------
    total_pdf : np.ndarray, shape (1000,)
        The summed probability density function evaluated over `x`.
    ax : matplotlib.axes.Axes
        The axes object with the plotted histogram and PDFs.

    Notes
    -----
    - Uses `sklearn.mixture.GaussianMixture` for fitting.
    - Assumes a 1D distribution; multi-dimensional inputs are not supported.
    - The histogram is normalized (`density=True`) to match the scale of the PDFs.
    - Covariance handling assumes diagonal/spherical structure for 1D data.
    """

    # Ensure proper shape for sklearn (n_samples, n_features)
    data = np.asarray(data)
    if data.ndim == 1:
        data_reshaped = data.reshape(-1, 1)
    else:
        data_reshaped = data

    gmm = GaussianMixture(n_components=n_components, random_state=random_state)
    gmm.fit(data_reshaped)

    means = gmm.means_.flatten()

    # Handle covariance shapes robustly
    if gmm.covariance_type in ['full', 'tied']:
        covariances = np.array([c[0, 0] for c in gmm.covariances_])
    else:  # 'diag' or 'spherical'
        covariances = gmm.covariances_.flatten()

    sds = np.sqrt(covariances)
    weights = gmm.weights_.flatten()

    x = np.linspace(data.min(), data.max(), 1000)

    ax.hist(data, bins=histo_bins, density=True, alpha=0.5, rwidth=rbarwidth, color=histcolor)

    total_pdf = np.zeros_like(x)

    for i, (mean, sd, weight) in enumerate(zip(means, sds, weights), start=1):
        print(f"mean of gauss {i}: {mean}, weight of gauss {i}: {weight}")
        pdf = weight * norm.pdf(x, mean, sd)
        total_pdf += pdf
        ax.plot(x, pdf, fitline, linewidth=linewidth, color=fitline_color)

    if plot_total:
        ax.plot(x, total_pdf, '-', linewidth=linewidth + 1, color='black')

    return total_pdf, ax, gmm


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