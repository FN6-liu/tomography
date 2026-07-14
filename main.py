"""
Main script for beam tomography simulation and reconstruction.
"""

import numpy as np
import matplotlib.pyplot as plt

from optics import (
    compute_emittance, cov_to_twiss, 
    normalization_matrix, inverse_normalization_matrix
)
from beam import BeamDistribution2D
from tomography import fbp, art, art_fast
from utils import (
    plot_phase_space, plot_sinogram, plot_comparison,
    compute_beam_parameters, print_beam_parameters,
    compute_reconstruction_metrics, plot_residuals
)


def main():
    """Main function to run the tomography simulation."""
    
    print("=" * 60)
    print("Beam Phase Space Tomography")
    print("FBP vs ART Reconstruction")
    print("=" * 60)
    
    # ============================================================
    # 1. Simulation Parameters
    # ============================================================
    
    # Beam parameters (physical space)
    emittance = 1.0e-6          # Geometric emittance [m·rad]
    alpha0 = 0.0                # Twiss alpha at reconstruction point
    beta0 = 1.0                 # Twiss beta at reconstruction point [m]
    gamma0 = (1.0 + alpha0**2) / beta0
    
    # Phase space grid
    nx = 64                     # Number of grid points in x
    nxp = 64                    # Number of grid points in x'
    sigma_x = np.sqrt(emittance * beta0)
    sigma_xp = np.sqrt(emittance * gamma0)
    x_range = 4.0 * sigma_x     # Range of x [m]
    xp_range = 4.0 * sigma_xp   # Range of x' [rad]
    
    x_grid = np.linspace(-x_range, x_range, nx)
    xp_grid = np.linspace(-xp_range, xp_range, nxp)
    
    # Tomography parameters
    n_angles = 18               # Number of projection angles
    n_iter_art = 50             # ART iterations
    
    print(f"\n--- Simulation Parameters ---")
    print(f"Emittance: {emittance:.2e} m·rad")
    print(f"Alpha0: {alpha0:.3f}, Beta0: {beta0:.3f} m")
    print(f"Grid: {nx}x{nxp}")
    print(f"Number of angles: {n_angles}")
    print(f"ART iterations: {n_iter_art}")
    
    # ============================================================
    # 2. Generate Beam Distribution
    # ============================================================
    
    print("\n--- Generating Beam Distribution ---")
    
    # Create Gaussian beam in physical phase space
    beam_phys = BeamDistribution2D.from_twiss(
        emittance, alpha0, beta0, x_grid, xp_grid
    )
    
    # Transform to normalized phase space (where rotation is pure)
    T_norm = normalization_matrix(alpha0, beta0)
    
    # For normalized space, we create a new grid in normalized coordinates
    # The distribution in normalized space is just a rotated version
    # For a matched beam, it should be circular
    beam_norm = beam_phys.rotate(0)  # Same distribution, but in normalized coords
    
    # ============================================================
    # 3. Generate Projections (Sinogram)
    # ============================================================
    
    print("\n--- Generating Projections ---")
    
    angles = np.linspace(0, np.pi, n_angles, endpoint=False)
    projections = []
    
    for theta in angles:
        proj = beam_norm.project(theta)
        projections.append(proj)
    
    print(f"Generated {len(projections)} projections")
    
    # ============================================================
    # 4. Reconstruct using FBP
    # ============================================================
    
    print("\n--- FBP Reconstruction ---")
    
    recon_fbp = fbp(projections, angles, x_grid, xp_grid, window='hann')
    
    # ============================================================
    # 5. Reconstruct using ART
    # ============================================================
    
    print("\n--- ART Reconstruction ---")
    
    recon_art = art_fast(projections, angles, x_grid, xp_grid, 
                         n_iter=n_iter_art, relax=1.0)
    
    # ============================================================
    # 6. Compute Beam Parameters
    # ============================================================
    
    print("\n--- Computing Beam Parameters ---")
    
    # Original parameters
    params_orig = compute_beam_parameters(beam_norm.density, x_grid, xp_grid)
    
    # FBP parameters
    params_fbp = compute_beam_parameters(recon_fbp, x_grid, xp_grid)
    
    # ART parameters
    params_art = compute_beam_parameters(recon_art, x_grid, xp_grid)
    
    print_beam_parameters(params_orig, "(Original)")
    print_beam_parameters(params_fbp, "(FBP)")
    print_beam_parameters(params_art, "(ART)")
    
    # ============================================================
    # 7. Compute Reconstruction Metrics
    # ============================================================
    
    metrics_fbp = compute_reconstruction_metrics(
        beam_norm.density, recon_fbp
    )
    metrics_art = compute_reconstruction_metrics(
        beam_norm.density, recon_art
    )
    
    print(f"\n--- Reconstruction Metrics ---")
    print(f"FBP - RMSE: {metrics_fbp['rmse']:.4e}, Correlation: {metrics_fbp['correlation']:.4f}")
    print(f"ART - RMSE: {metrics_art['rmse']:.4e}, Correlation: {metrics_art['correlation']:.4f}")
    
    # ============================================================
    # 8. Plot Results
    # ============================================================
    
    print("\n--- Generating Plots ---")
    
    # Comparison plot
    fig1, axes = plot_comparison(
        beam_norm.density, recon_fbp, recon_art,
        x_grid, xp_grid,
        titles=['Original', 'FBP Reconstruction', 'ART Reconstruction']
    )
    fig1.savefig('comparison.png', dpi=150, bbox_inches='tight')
    print("Saved: comparison.png")
    
    # Sinogram
    fig2, ax = plt.subplots(figsize=(8, 5))
    plot_sinogram(projections, angles, 'Sinogram', ax)
    fig2.savefig('sinogram.png', dpi=150, bbox_inches='tight')
    print("Saved: sinogram.png")
    
    # Residuals
    fig3, _ = plot_residuals(
        beam_norm.density, recon_fbp,
        x_grid, xp_grid,
        title='FBP Residuals'
    )
    fig3.savefig('residuals_fbp.png', dpi=150, bbox_inches='tight')
    print("Saved: residuals_fbp.png")
    
    fig4, _ = plot_residuals(
        beam_norm.density, recon_art,
        x_grid, xp_grid,
        title='ART Residuals'
    )
    fig4.savefig('residuals_art.png', dpi=150, bbox_inches='tight')
    print("Saved: residuals_art.png")
    
    # ============================================================
    # 9. Show plots
    # ============================================================
    
    plt.show()
    
    print("\n--- Done! ---")


if __name__ == "__main__":
    main()