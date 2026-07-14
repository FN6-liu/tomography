"""
Utility functions for plotting and analysis.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm


def plot_phase_space(density, x_grid, xp_grid, title=None, ax=None, cmap='viridis'):
    """
    Plot 2D phase space density with contours.
    
    Parameters
    ----------
    density : np.ndarray
        2D density array
    x_grid : np.ndarray
        x coordinates [m]
    xp_grid : np.ndarray
        x' coordinates [rad]
    title : str, optional
        Plot title
    ax : matplotlib.axes.Axes, optional
        Axes to plot on
    cmap : str
        Colormap name
    
    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 5))
    
    # Calculate levels for contours
    max_val = np.max(density)
    levels = np.linspace(0.01 * max_val, max_val, 8)
    
    im = ax.imshow(density.T, origin='lower', 
                   extent=[x_grid[0], x_grid[-1], xp_grid[0], xp_grid[-1]],
                   aspect='auto', cmap=cmap, norm=PowerNorm(gamma=0.5))
    
    ax.contour(x_grid, xp_grid, density.T, levels=levels[1:], 
               colors='white', alpha=0.5, linewidths=0.5)
    
    ax.set_xlabel('x [m]')
    ax.set_ylabel("x' [rad]")
    if title:
        ax.set_title(title)
    
    plt.colorbar(im, ax=ax, label='Density')
    
    return ax


def plot_sinogram(projections, angles, title=None, ax=None):
    """
    Plot sinogram (projections vs angle).
    
    Parameters
    ----------
    projections : list of np.ndarray
        List of projection arrays
    angles : np.ndarray
        Array of angles [radians]
    title : str, optional
        Plot title
    ax : matplotlib.axes.Axes, optional
        Axes to plot on
    
    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))
    
    sinogram = np.array(projections)
    extent = [0, len(projections[0]), np.rad2deg(angles[0]), np.rad2deg(angles[-1])]
    
    im = ax.imshow(sinogram, origin='lower', extent=extent, aspect='auto', cmap='viridis')
    
    ax.set_xlabel('Projection bin')
    ax.set_ylabel('Angle [deg]')
    if title:
        ax.set_title(title)
    
    plt.colorbar(im, ax=ax, label='Intensity')
    
    return ax


def plot_comparison(original, fbp_recon, art_recon, x_grid, xp_grid,
                    titles=None, save_path=None):
    """
    Plot comparison of original and reconstructed distributions.
    
    Parameters
    ----------
    original : np.ndarray
        Original density
    fbp_recon : np.ndarray
        FBP reconstruction
    art_recon : np.ndarray
        ART reconstruction
    x_grid : np.ndarray
        x coordinates [m]
    xp_grid : np.ndarray
        x' coordinates [rad]
    titles : list, optional
        Titles for subplots
    save_path : str, optional
        Path to save figure
    """
    if titles is None:
        titles = ['Original', 'FBP', 'ART']
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    plot_phase_space(original, x_grid, xp_grid, titles[0], axes[0])
    plot_phase_space(fbp_recon, x_grid, xp_grid, titles[1], axes[1])
    plot_phase_space(art_recon, x_grid, xp_grid, titles[2], axes[2])
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    
    return fig, axes


def plot_projections_comparison(original, projections, x_grid, angles,
                                n_show=6, save_path=None):
    """
    Plot original phase space and selected projections.
    
    Parameters
    ----------
    original : np.ndarray
        Original phase space density
    projections : list of np.ndarray
        List of projections
    x_grid : np.ndarray
        x coordinates [m]
    angles : np.ndarray
        Array of angles [radians]
    n_show : int
        Number of projections to show
    save_path : str, optional
        Path to save figure
    """
    fig = plt.figure(figsize=(14, 6))
    
    # Phase space
    ax1 = fig.add_subplot(2, 4, 1)
    plot_phase_space(original, x_grid, np.linspace(-1, 1, 64), 
                     'Original Phase Space', ax1)
    
    # Projections
    n_proj = len(projections)
    idx = np.linspace(0, n_proj-1, min(n_show, 8), dtype=int)
    
    for i, idx_i in enumerate(idx):
        ax = fig.add_subplot(2, 4, i+5)
        ax.plot(x_grid, projections[idx_i])
        ax.set_xlabel('x [m]')
        ax.set_ylabel('Intensity')
        ax.set_title(f'θ = {np.rad2deg(angles[idx_i]):.0f}°')
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    
    return fig


def plot_residuals(original, reconstructed, x_grid, xp_grid, title=None, save_path=None):
    """
    Plot residuals between original and reconstructed distributions.
    
    Parameters
    ----------
    original : np.ndarray
        Original density
    reconstructed : np.ndarray
        Reconstructed density
    x_grid : np.ndarray
        x coordinates [m]
    xp_grid : np.ndarray
        x' coordinates [rad]
    title : str, optional
        Plot title
    save_path : str, optional
        Path to save figure
    """
    residuals = original - reconstructed
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Residual plot
    im = axes[0].imshow(residuals.T, origin='lower',
                        extent=[x_grid[0], x_grid[-1], xp_grid[0], xp_grid[-1]],
                        aspect='auto', cmap='RdBu_r')
    axes[0].set_xlabel('x [m]')
    axes[0].set_ylabel("x' [rad]")
    axes[0].set_title('Residuals')
    plt.colorbar(im, ax=axes[0])
    
    # Histogram of residuals
    axes[1].hist(residuals.flatten(), bins=50, edgecolor='black', alpha=0.7)
    axes[1].axvline(0, color='red', linestyle='--', label='Zero')
    axes[1].set_xlabel('Residual')
    axes[1].set_ylabel('Count')
    axes[1].set_title('Residual Distribution')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    if title:
        fig.suptitle(title)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
    
    return fig, axes


def compute_beam_parameters(density, x_grid, xp_grid):
    """
    Compute beam parameters (emittance, Twiss) from phase space density.
    
    Parameters
    ----------
    density : np.ndarray
        2D phase space density
    x_grid : np.ndarray
        x coordinates [m]
    xp_grid : np.ndarray
        x' coordinates [rad]
    
    Returns
    -------
    dict
        Dictionary with emittance, alpha, beta, gamma
    """
    from optics import cov_to_twiss, compute_emittance
    
    # Create temporary distribution
    from beam import BeamDistribution2D
    temp = BeamDistribution2D(density, x_grid, xp_grid)
    
    # Get moments
    cov = temp.get_moments()
    
    # Compute parameters
    alpha, beta, gamma, eps = cov_to_twiss(cov)
    
    return {
        'emittance': eps,
        'alpha': alpha,
        'beta': beta,
        'gamma': gamma,
        'covariance': cov
    }


def compute_reconstruction_metrics(original, reconstructed):
    """
    Compute quality metrics for reconstruction.
    
    Parameters
    ----------
    original : np.ndarray
        Original density
    reconstructed : np.ndarray
        Reconstructed density
    
    Returns
    -------
    dict
        Dictionary with RMSE, correlation, etc.
    """
    # Flatten
    orig_flat = original.flatten()
    recon_flat = reconstructed.flatten()
    
    # RMSE
    rmse = np.sqrt(np.mean((orig_flat - recon_flat)**2))
    
    # NRMSE (normalized by range)
    nrmse = rmse / (np.max(orig_flat) - np.min(orig_flat))
    
    # Correlation coefficient
    corr = np.corrcoef(orig_flat, recon_flat)[0, 1]
    
    # Relative error in total intensity
    total_error = abs(np.sum(recon_flat) - np.sum(orig_flat)) / np.sum(orig_flat)
    
    return {
        'rmse': rmse,
        'nrmse': nrmse,
        'correlation': corr,
        'total_error': total_error
    }


def print_beam_parameters(params, label=''):
    """
    Print beam parameters in a formatted way.
    
    Parameters
    ----------
    params : dict
        Dictionary from compute_beam_parameters
    label : str
        Label for the parameter set
    """
    print(f"\n=== Beam Parameters {label} ===")
    print(f"Emittance: {params['emittance']:.2e} m·rad")
    print(f"Beta:      {params['beta']:.4f} m")
    print(f"Alpha:     {params['alpha']:.4f}")
    print(f"Gamma:     {params['gamma']:.4f} 1/m")