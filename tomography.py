"""
Tomography reconstruction algorithms: FBP and ART.
"""

import numpy as np
from scipy.fft import fft, ifft, fftfreq
from scipy.sparse import lil_matrix, csr_matrix
import warnings


def fbp(projections, angles, x_grid, xp_grid, window='hann'):
    """
    Filtered Back Projection reconstruction.
    
    Parameters
    ----------
    projections : list of np.ndarray
        List of 1D projection arrays
    angles : np.ndarray
        Array of projection angles [radians]
    x_grid : np.ndarray
        1D array of x coordinates [m]
    xp_grid : np.ndarray
        1D array of x' coordinates [rad]
    window : str or None
        Window function to apply to ramp filter ('hann', 'hamming', or None)
    
    Returns
    -------
    np.ndarray
        2D reconstructed phase space density
    """
    n_angles = len(angles)
    n_proj = len(projections[0])
    
    # Create ramp filter
    freq = fftfreq(n_proj)
    ramp = np.abs(freq)
    
    # Apply window if specified
    if window == 'hann':
        ramp *= np.hanning(n_proj)
    elif window == 'hamming':
        ramp *= np.hamming(n_proj)
    
    # Filter each projection
    filtered_projections = []
    for proj in projections:
        fproj = fft(proj)
        fproj_filtered = fproj * ramp
        filtered = np.real(ifft(fproj_filtered))
        filtered_projections.append(filtered)
    
    # Back-project
    nx = len(x_grid)
    nxp = len(xp_grid)
    recon = np.zeros((nx, nxp))
    
    dx = x_grid[1] - x_grid[0]
    
    for i, theta in enumerate(angles):
        for j, x in enumerate(x_grid):
            for k, xp in enumerate(xp_grid):
                # Projection position
                r = x * np.cos(theta) + xp * np.sin(theta)
                # Find index in projection
                idx = int(np.round((r - x_grid[0]) / dx))
                if 0 <= idx < n_proj:
                    recon[j, k] += filtered_projections[i][idx]
    
    # Normalize
    recon *= np.pi / n_angles
    
    # Ensure non-negative
    recon[recon < 0] = 0
    
    return recon


def art(projections, angles, x_grid, xp_grid, n_iter=50, relax=1.0):
    """
    Algebraic Reconstruction Technique using sparse matrix.
    
    Parameters
    ----------
    projections : list of np.ndarray
        List of 1D projection arrays
    angles : np.ndarray
        Array of projection angles [radians]
    x_grid : np.ndarray
        1D array of x coordinates [m]
    xp_grid : np.ndarray
        1D array of x' coordinates [rad]
    n_iter : int
        Number of iterations
    relax : float
        Relaxation parameter (0 < relax <= 2)
    
    Returns
    -------
    np.ndarray
        2D reconstructed phase space density
    """
    n_angles = len(angles)
    n_proj = len(projections[0])
    nx = len(x_grid)
    nxp = len(xp_grid)
    n_pixels = nx * nxp
    n_rays = n_angles * n_proj
    
    # Flatten projections
    b = np.concatenate(projections)
    
    # Build projection matrix P (sparse)
    print("Building projection matrix...")
    P = lil_matrix((n_rays, n_pixels))
    
    dx = x_grid[1] - x_grid[0]
    
    ray_idx = 0
    for i, theta in enumerate(angles):
        for m in range(n_proj):
            # x coordinate of this projection bin
            x_proj = x_grid[0] + m * dx
            for j, x in enumerate(x_grid):
                for k, xp in enumerate(xp_grid):
                    # Calculate projection position
                    r = x * np.cos(theta) + xp * np.sin(theta)
                    # Find bin index
                    bin_idx = int(np.round((r - x_grid[0]) / dx))
                    if bin_idx == m:
                        P[ray_idx, j * nxp + k] = 1.0
            ray_idx += 1
    
    P = P.tocsr()
    print(f"Projection matrix built: {P.shape}")
    
    # Initial guess (uniform distribution)
    rho = np.ones(n_pixels) / n_pixels
    
    # Iterative reconstruction (Kaczmarz method)
    print(f"Running ART with {n_iter} iterations...")
    for it in range(n_iter):
        for ray in range(n_rays):
            # Get non-zero elements
            row = P[ray]
            if row.nnz == 0:
                continue
            
            # Get indices and values
            indices = row.indices
            values = row.data
            
            # Compute dot product
            dot = np.sum(values * rho[indices])
            
            # Compute residual
            residual = b[ray] - dot
            
            # Compute norm squared
            norm2 = np.sum(values**2)
            
            if norm2 > 0:
                # Update
                update = relax * residual / norm2
                rho[indices] += update * values
        
        # Enforce non-negativity
        rho[rho < 0] = 0
        
        if (it + 1) % 10 == 0:
            print(f"Iteration {it + 1}/{n_iter}")
    
    # Reshape
    recon = rho.reshape((nx, nxp))
    
    # Normalize
    recon /= np.sum(recon)
    
    return recon


def art_fast(projections, angles, x_grid, xp_grid, n_iter=50, relax=1.0):
    """
    Faster ART implementation using projection/backprojection operations
    instead of explicit matrix construction.
    
    Parameters
    ----------
    projections : list of np.ndarray
        List of 1D projection arrays
    angles : np.ndarray
        Array of projection angles [radians]
    x_grid : np.ndarray
        1D array of x coordinates [m]
    xp_grid : np.ndarray
        1D array of x' coordinates [rad]
    n_iter : int
        Number of iterations
    relax : float
        Relaxation parameter (0 < relax <= 2)
    
    Returns
    -------
    np.ndarray
        2D reconstructed phase space density
    """
    nx = len(x_grid)
    nxp = len(xp_grid)
    n_proj = len(projections[0])
    n_angles = len(angles)
    
    dx = x_grid[1] - x_grid[0]
    dxp = xp_grid[1] - xp_grid[0]
    area_element = dx * dxp
    
    # Flatten projections
    b = np.concatenate(projections)
    
    # Initial guess (uniform distribution)
    rho = np.ones((nx, nxp)) / (nx * nxp)
    
    print(f"Running fast ART with {n_iter} iterations...")
    
    for it in range(n_iter):
        ray_idx = 0
        
        for i, theta in enumerate(angles):
            for m in range(n_proj):
                # x coordinate of this projection bin
                x_proj = x_grid[0] + m * dx
                
                # Forward project: compute projection value at bin m
                proj_val = 0.0
                
                # Find which phase space points project to this bin
                for j in range(nx):
                    for k in range(nxp):
                        # Calculate projection position
                        r = x_grid[j] * np.cos(theta) + xp_grid[k] * np.sin(theta)
                        bin_idx = int(np.round((r - x_grid[0]) / dx))
                        if bin_idx == m:
                            proj_val += rho[j, k] * area_element
                
                # Compute residual
                residual = b[ray_idx] - proj_val
                
                # Count number of points contributing to this bin
                count = 0
                if abs(residual) > 1e-12:
                    for j in range(nx):
                        for k in range(nxp):
                            r = x_grid[j] * np.cos(theta) + xp_grid[k] * np.sin(theta)
                            bin_idx = int(np.round((r - x_grid[0]) / dx))
                            if bin_idx == m:
                                count += 1
                    
                    # Update only points that project to this bin
                    if count > 0:
                        update = relax * residual / count
                        for j in range(nx):
                            for k in range(nxp):
                                r = x_grid[j] * np.cos(theta) + xp_grid[k] * np.sin(theta)
                                bin_idx = int(np.round((r - x_grid[0]) / dx))
                                if bin_idx == m:
                                    rho[j, k] += update
                
                ray_idx += 1
        
        # Enforce non-negativity
        rho[rho < 0] = 0
        
        # Normalize to preserve total intensity
        rho /= np.sum(rho)
        
        if (it + 1) % 10 == 0:
            print(f"Iteration {it + 1}/{n_iter}")
    
    return rho


def art_memory_efficient(projections, angles, x_grid, xp_grid, n_iter=50, relax=1.0):
    """
    Memory-efficient ART that avoids building the full projection matrix.
    Uses a more efficient forward/backward projection implementation.
    
    Parameters
    ----------
    projections : list of np.ndarray
        List of 1D projection arrays
    angles : np.ndarray
        Array of projection angles [radians]
    x_grid : np.ndarray
        1D array of x coordinates [m]
    xp_grid : np.ndarray
        1D array of x' coordinates [rad]
    n_iter : int
        Number of iterations
    relax : float
        Relaxation parameter (0 < relax <= 2)
    
    Returns
    -------
    np.ndarray
        2D reconstructed phase space density
    """
    nx = len(x_grid)
    nxp = len(xp_grid)
    n_proj = len(projections[0])
    n_angles = len(angles)
    
    dx = x_grid[1] - x_grid[0]
    dxp = xp_grid[1] - xp_grid[0]
    area_element = dx * dxp
    
    # Flatten projections
    b = np.concatenate(projections)
    
    # Initial guess
    rho = np.ones((nx, nxp)) / (nx * nxp)
    
    # Precompute cos and sin for each angle
    cos_theta = np.cos(angles)
    sin_theta = np.sin(angles)
    
    print(f"Running memory-efficient ART with {n_iter} iterations...")
    
    for it in range(n_iter):
        ray_idx = 0
        
        for i, theta in enumerate(angles):
            ct = cos_theta[i]
            st = sin_theta[i]
            
            for m in range(n_proj):
                x_proj = x_grid[0] + m * dx
                
                # Forward projection
                proj_val = 0.0
                for j in range(nx):
                    for k in range(nxp):
                        r = x_grid[j] * ct + xp_grid[k] * st
                        bin_idx = int(np.round((r - x_grid[0]) / dx))
                        if bin_idx == m:
                            proj_val += rho[j, k] * area_element
                
                residual = b[ray_idx] - proj_val
                
                if abs(residual) > 1e-12:
                    # Find all points that project to this bin
                    points = []
                    for j in range(nx):
                        for k in range(nxp):
                            r = x_grid[j] * ct + xp_grid[k] * st
                            bin_idx = int(np.round((r - x_grid[0]) / dx))
                            if bin_idx == m:
                                points.append((j, k))
                    
                    n_points = len(points)
                    if n_points > 0:
                        update = relax * residual / n_points
                        for j, k in points:
                            rho[j, k] += update
                
                ray_idx += 1
        
        # Enforce non-negativity
        rho[rho < 0] = 0
        rho /= np.sum(rho)
        
        if (it + 1) % 10 == 0:
            print(f"Iteration {it + 1}/{n_iter}")
    
    return rho