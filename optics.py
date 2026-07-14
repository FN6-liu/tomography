"""
Optics functions for beam transport and phase space normalization.
Based on linear accelerator optics (Courant-Snyder formalism).
"""

import numpy as np


def drift_matrix(L):
    """
    Transport matrix for a drift space of length L.
    
    Parameters
    ----------
    L : float
        Drift length [m]
    
    Returns
    -------
    np.ndarray
        2x2 transport matrix
    """
    return np.array([[1.0, L],
                     [0.0, 1.0]])


def quadrupole_matrix(f):
    """
    Transport matrix for a thin quadrupole with focal length f.
    f > 0 focusing in horizontal plane, f < 0 defocusing.
    
    Parameters
    ----------
    f : float
        Focal length [m]
    
    Returns
    -------
    np.ndarray
        2x2 transport matrix
    """
    return np.array([[1.0, 0.0],
                     [-1.0/f, 1.0]])


def thick_quadrupole_matrix(k, L):
    """
    Transport matrix for a thick quadrupole of length L with normalized
    gradient k (k > 0 focusing in horizontal plane).
    
    Parameters
    ----------
    k : float
        Normalized quadrupole strength [1/m^2]
    L : float
        Quadrupole length [m]
    
    Returns
    -------
    np.ndarray
        2x2 transport matrix
    """
    if k > 0:
        sqrt_k = np.sqrt(k)
        return np.array([[np.cos(sqrt_k * L), (1.0/sqrt_k) * np.sin(sqrt_k * L)],
                         [-sqrt_k * np.sin(sqrt_k * L), np.cos(sqrt_k * L)]])
    elif k < 0:
        sqrt_neg_k = np.sqrt(-k)
        return np.array([[np.cosh(sqrt_neg_k * L), (1.0/sqrt_neg_k) * np.sinh(sqrt_neg_k * L)],
                         [sqrt_neg_k * np.sinh(sqrt_neg_k * L), np.cosh(sqrt_neg_k * L)]])
    else:
        return drift_matrix(L)


def propagate_twiss(alpha, beta, M):
    """
    Propagate Twiss parameters through a transport matrix.
    
    Parameters
    ----------
    alpha : float
        Alpha parameter at start
    beta : float
        Beta parameter at start [m]
    M : np.ndarray
        2x2 transport matrix
    
    Returns
    -------
    tuple
        (alpha_end, beta_end)
    """
    C = M[0, 0]
    S = M[0, 1]
    Cp = M[1, 0]
    Sp = M[1, 1]
    
    gamma = (1.0 + alpha**2) / beta
    
    beta_end = C**2 * beta - 2.0 * C * S * alpha + S**2 * gamma
    alpha_end = -C * Cp * beta + (C * Sp + Cp * S) * alpha - S * Sp * gamma
    
    return alpha_end, beta_end


def normalization_matrix(alpha, beta):
    """
    Returns the transformation matrix from physical phase space
    to normalized phase space.
    
    Parameters
    ----------
    alpha : float
        Twiss alpha parameter
    beta : float
        Twiss beta parameter [m]
    
    Returns
    -------
    np.ndarray
        2x2 normalization matrix
    """
    return np.array([[1.0 / np.sqrt(beta), 0.0],
                     [alpha / np.sqrt(beta), np.sqrt(beta)]])


def inverse_normalization_matrix(alpha, beta):
    """
    Returns the transformation matrix from normalized phase space
    back to physical phase space.
    
    Parameters
    ----------
    alpha : float
        Twiss alpha parameter
    beta : float
        Twiss beta parameter [m]
    
    Returns
    -------
    np.ndarray
        2x2 inverse normalization matrix
    """
    return np.array([[np.sqrt(beta), 0.0],
                     [-alpha / np.sqrt(beta), 1.0 / np.sqrt(beta)]])


def normalized_rotation_angle(M_phys, alpha, beta):
    """
    Given a physical transport matrix M from point A to point B,
    and the Twiss parameters at point A, compute the rotation angle
    in normalized phase space.
    
    Parameters
    ----------
    M_phys : np.ndarray
        2x2 physical transport matrix
    alpha : float
        Twiss alpha at point A
    beta : float
        Twiss beta at point A [m]
    
    Returns
    -------
    float
        Rotation angle [radians]
    """
    T_A = normalization_matrix(alpha, beta)
    T_A_inv = np.linalg.inv(T_A)
    
    alpha_B, beta_B = propagate_twiss(alpha, beta, M_phys)
    T_B = normalization_matrix(alpha_B, beta_B)
    
    M_norm = T_B @ M_phys @ T_A_inv
    theta = np.arctan2(M_norm[0, 1], M_norm[0, 0])
    
    return theta


def compute_emittance(cov_matrix):
    """
    Compute geometric emittance from covariance matrix.
    
    Parameters
    ----------
    cov_matrix : np.ndarray
        2x2 covariance matrix [[<x^2>, <xx'>], [<xx'>, <x'^2>]]
    
    Returns
    -------
    float
        Geometric emittance [m·rad]
    """
    det = cov_matrix[0, 0] * cov_matrix[1, 1] - cov_matrix[0, 1]**2
    if det < 0:
        return np.nan
    return np.sqrt(det)


def cov_to_twiss(cov_matrix):
    """
    Compute Twiss parameters from covariance matrix.
    
    Parameters
    ----------
    cov_matrix : np.ndarray
        2x2 covariance matrix
    
    Returns
    -------
    tuple
        (alpha, beta, gamma, emittance)
    """
    eps = compute_emittance(cov_matrix)
    if np.isnan(eps) or eps <= 0:
        return np.nan, np.nan, np.nan, np.nan
    
    beta = cov_matrix[0, 0] / eps
    alpha = -cov_matrix[0, 1] / eps
    gamma = cov_matrix[1, 1] / eps
    
    return alpha, beta, gamma, eps