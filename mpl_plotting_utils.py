import numpy as np
from sklearn.mixture import GaussianMixture
from scipy.stats import norm
from scipy.stats import gaussian_kde


def density_based_x_jitter(
    yvals: np.ndarray,
    center: float,
    max_spread: float = 0.03,
) -> np.ndarray:
    """
    Generate x-coordinates for a scatter plot by jittering points according to
    their local density along the y-axis.

    Points located in dense regions of the y-distribution receive larger
    horizontal jitter, while isolated points remain close to the center.
    This produces a swarm-like appearance without requiring iterative point
    placement.

    Parameters
    ----------
    yvals : array-like
        Y-values to plot.
    center : float
        Central x-coordinate around which points are jittered.
    max_spread : float, default=0.03
        Maximum standard deviation of the horizontal jitter. Points in the
        highest-density region receive approximately this amount of jitter.

    Returns
    -------
    tuple of xvals and yvals
    np.ndarray
        Jittered x-coordinates corresponding to `yvals` a
    """
    yvals = np.asarray(yvals, dtype=float)
    ilen=len(yvals)
    yvals = yvals[~np.isnan(yvals)]
    if len(yvals) < ilen:
        print('There were nan vals in yvals that were filered out')

    if len(yvals) == 0:
        return np.array([])

    if len(yvals) == 1:
        return np.array([center])

    density = gaussian_kde(yvals)(yvals)
    density /= density.max()

    jitter_sd = max_spread * density

    return np.random.normal(loc=center, scale=jitter_sd), yvals


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
