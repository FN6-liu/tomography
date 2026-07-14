"""
Tomography reconstruction algorithms: FBP and ART.
"""

import numpy as np
from scipy.fft import fft, ifft, fftfreq
from scipy.sparse import lil_matrix


def _normalize_density(density, dx, dxp):
    """Normalize a sampled phase-space density to unit integral."""
    total = np.sum(density) * dx * dxp
    if total > 0:
        return density / total
    return density


def _precompute_bin_indices(angles, x_grid, xp_grid):
    """Map each phase-space cell to its projection bin for every angle."""
    dx = x_grid[1] - x_grid[0]
    X, Xp = np.meshgrid(x_grid, xp_grid, indexing='ij')
    bin_maps = []

    for theta in angles:
        r = X * np.cos(theta) + Xp * np.sin(theta)
        bin_idx = np.rint((r - x_grid[0]) / dx).astype(int)
        valid = (bin_idx >= 0) & (bin_idx < len(x_grid))
        bin_maps.append((bin_idx, valid))

    return bin_maps


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

    freq = fftfreq(n_proj)
    ramp = np.abs(freq)

    if window == 'hann':
        ramp *= np.hanning(n_proj)
    elif window == 'hamming':
        ramp *= np.hamming(n_proj)

    filtered_projections = []
    for proj in projections:
        fproj = fft(proj)
        filtered = np.real(ifft(fproj * ramp))
        filtered_projections.append(filtered)

    nx = len(x_grid)
    nxp = len(xp_grid)
    recon = np.zeros((nx, nxp))
    dx = x_grid[1] - x_grid[0]
    dxp = xp_grid[1] - xp_grid[0]

    for i, theta in enumerate(angles):
        ct = np.cos(theta)
        st = np.sin(theta)
        for j, x in enumerate(x_grid):
            for k, xp in enumerate(xp_grid):
                r = x * ct + xp * st
                idx = int(np.rint((r - x_grid[0]) / dx))
                if 0 <= idx < n_proj:
                    recon[j, k] += filtered_projections[i][idx]

    recon *= np.pi / n_angles
    recon[recon < 0] = 0
    return _normalize_density(recon, dx, dxp)


def art(projections, angles, x_grid, xp_grid, n_iter=50, relax=1.0):
    """
    Algebraic Reconstruction Technique using an explicit sparse matrix.

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
    dx = x_grid[1] - x_grid[0]
    dxp = xp_grid[1] - xp_grid[0]
    area_element = dx * dxp

    b = np.concatenate(projections)
    print("Building projection matrix...")
    P = lil_matrix((n_rays, n_pixels))
    bin_maps = _precompute_bin_indices(angles, x_grid, xp_grid)

    ray_idx = 0
    for i in range(n_angles):
        bin_idx, valid = bin_maps[i]
        for m in range(n_proj):
            rows, cols = np.where(valid & (bin_idx == m))
            for j, k in zip(rows, cols):
                P[ray_idx, j * nxp + k] = area_element
            ray_idx += 1

    P = P.tocsr()
    print(f"Projection matrix built: {P.shape}")

    rho = np.ones(n_pixels) / (n_pixels * area_element)

    print(f"Running ART with {n_iter} iterations...")
    for it in range(n_iter):
        for ray in range(n_rays):
            row = P[ray]
            if row.nnz == 0:
                continue

            indices = row.indices
            values = row.data
            dot = np.sum(values * rho[indices])
            residual = b[ray] - dot
            norm2 = np.sum(values**2)

            if norm2 > 0:
                rho[indices] += relax * residual * values / norm2

        rho[rho < 0] = 0
        rho = _normalize_density(rho, dx, dxp)

        if (it + 1) % 10 == 0:
            print(f"Iteration {it + 1}/{n_iter}")

    recon = rho.reshape((nx, nxp))
    recon[recon < 0] = 0
    return _normalize_density(recon, dx, dxp)


def art_fast(projections, angles, x_grid, xp_grid, n_iter=50, relax=1.0):
    """
    ART implementation using precomputed projection-bin maps.

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
    dx = x_grid[1] - x_grid[0]
    dxp = xp_grid[1] - xp_grid[0]
    area_element = dx * dxp

    b = np.concatenate(projections)
    rho = np.ones((nx, nxp)) / (nx * nxp * area_element)
    bin_maps = _precompute_bin_indices(angles, x_grid, xp_grid)

    print(f"Running fast ART with {n_iter} iterations...")

    for it in range(n_iter):
        ray_idx = 0
        for i in range(len(angles)):
            bin_idx, valid = bin_maps[i]
            for m in range(n_proj):
                mask = valid & (bin_idx == m)
                count = np.count_nonzero(mask)
                proj_val = np.sum(rho[mask]) * area_element
                residual = b[ray_idx] - proj_val

                if count > 0 and abs(residual) > 1e-12:
                    rho[mask] += relax * residual / (count * area_element)

                ray_idx += 1

        rho[rho < 0] = 0
        rho = _normalize_density(rho, dx, dxp)

        if (it + 1) % 10 == 0:
            print(f"Iteration {it + 1}/{n_iter}")

    return rho


def art_memory_efficient(projections, angles, x_grid, xp_grid, n_iter=50, relax=1.0):
    """
    Memory-efficient ART using the same projection model as art_fast.

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
    return art_fast(projections, angles, x_grid, xp_grid, n_iter=n_iter, relax=relax)
