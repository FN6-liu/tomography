"""
Beam distribution models and projection generation for phase space tomography.
"""

import numpy as np
from scipy.stats import multivariate_normal
from scipy.interpolate import RegularGridInterpolator


class BeamDistribution2D:
    """
    2D transverse phase space distribution (x, x').
    """
    
    def __init__(self, density, x_grid, xp_grid):
        """
        Parameters
        ----------
        density : np.ndarray
            2D array of phase space density (shape: len(x_grid) x len(xp_grid))
        x_grid : np.ndarray
            1D array of x coordinates [m]
        xp_grid : np.ndarray
            1D array of x' coordinates [rad]
        """
        self.density = density
        self.x_grid = x_grid
        self.xp_grid = xp_grid
        self.dx = x_grid[1] - x_grid[0]
        self.dxp = xp_grid[1] - xp_grid[0]
        
        # Normalize density
        self.density /= np.sum(self.density)
    
    @classmethod
    def gaussian(cls, mean, cov, x_grid, xp_grid):
        """
        Create a Gaussian beam distribution.
        
        Parameters
        ----------
        mean : array-like
            Mean [x0, x0'] in physical phase space
        cov : array-like
            2x2 covariance matrix
        x_grid : np.ndarray
            1D array of x coordinates [m]
        xp_grid : np.ndarray
            1D array of x' coordinates [rad]
        
        Returns
        -------
        BeamDistribution2D
        """
        X, Xp = np.meshgrid(x_grid, xp_grid, indexing='ij')
        pos = np.stack([X, Xp], axis=-1)
        density = multivariate_normal(mean, cov).pdf(pos)
        return cls(density, x_grid, xp_grid)
    
    @classmethod
    def from_twiss(cls, emittance, alpha, beta, x_grid, xp_grid):
        """
        Create a Gaussian beam distribution from Twiss parameters.
        
        Parameters
        ----------
        emittance : float
            Geometric emittance [m·rad]
        alpha : float
            Twiss alpha parameter
        beta : float
            Twiss beta parameter [m]
        x_grid : np.ndarray
            1D array of x coordinates [m]
        xp_grid : np.ndarray
            1D array of x' coordinates [rad]
        
        Returns
        -------
        BeamDistribution2D
        """
        gamma = (1.0 + alpha**2) / beta
        mean = [0.0, 0.0]
        cov = [[emittance * beta, -emittance * alpha],
               [-emittance * alpha, emittance * gamma]]
        return cls.gaussian(mean, cov, x_grid, xp_grid)
    
    @classmethod
    def multi_gaussian(cls, components, x_grid, xp_grid):
        """
        Create a multi-Gaussian beam distribution.
        
        Parameters
        ----------
        components : list of dict
            List of dicts with 'weight', 'mean', 'cov' keys
        x_grid : np.ndarray
            1D array of x coordinates [m]
        xp_grid : np.ndarray
            1D array of x' coordinates [rad]
        
        Returns
        -------
        BeamDistribution2D
        """
        X, Xp = np.meshgrid(x_grid, xp_grid, indexing='ij')
        pos = np.stack([X, Xp], axis=-1)
        density = np.zeros_like(X, dtype=float)
        
        for comp in components:
            rv = multivariate_normal(comp['mean'], comp['cov'])
            density += comp['weight'] * rv.pdf(pos)
        
        return cls(density, x_grid, xp_grid)
    
    def rotate(self, theta):
        """
        Rotate the phase space distribution by angle theta [radians].
        
        Parameters
        ----------
        theta : float
            Rotation angle in normalized phase space [radians]
        
        Returns
        -------
        BeamDistribution2D
            New rotated distribution
        """
        # Rotate the grid
        X_rot = self.x_grid * np.cos(theta) - self.xp_grid * np.sin(theta)
        Xp_rot = self.x_grid * np.sin(theta) + self.xp_grid * np.cos(theta)
        
        # Interpolate density onto rotated grid
        interp = RegularGridInterpolator(
            (self.x_grid, self.xp_grid), self.density,
            bounds_error=False, fill_value=0.0
        )
        
        X_grid, Xp_grid = np.meshgrid(X_rot, Xp_rot, indexing='ij')
        points = np.stack([X_grid.ravel(), Xp_grid.ravel()], axis=-1)
        density_rot = interp(points).reshape(X_grid.shape)
        
        return BeamDistribution2D(density_rot, X_rot, Xp_rot)
    
    def project(self, theta=0.0):
        """
        Project the phase space density onto the x-axis after rotating by theta.
        
        Parameters
        ----------
        theta : float
            Rotation angle in normalized phase space [radians]
        
        Returns
        -------
        np.ndarray
            1D projection (profile) on the x-axis
        """
        if theta != 0.0:
            rotated = self.rotate(theta)
            projection = np.sum(rotated.density, axis=1) * rotated.dxp
            return projection
        else:
            projection = np.sum(self.density, axis=1) * self.dxp
            return projection
    
    def get_moments(self):
        """
        Compute the second-order moments of the distribution.
        
        Returns
        -------
        np.ndarray
            2x2 covariance matrix
        """
        total = np.sum(self.density)
        
        x_mesh, xp_mesh = np.meshgrid(self.x_grid, self.xp_grid, indexing='ij')
        mean_x = np.sum(x_mesh * self.density) / total
        mean_xp = np.sum(xp_mesh * self.density) / total
        
        x2 = np.sum((x_mesh - mean_x)**2 * self.density) / total
        xp2 = np.sum((xp_mesh - mean_xp)**2 * self.density) / total
        xxp = np.sum((x_mesh - mean_x) * (xp_mesh - mean_xp) * self.density) / total
        
        return np.array([[x2, xxp], [xxp, xp2]])
    
    def normalized_to_physical(self, alpha, beta):
        """
        Transform from normalized phase space to physical phase space.
        
        Parameters
        ----------
        alpha : float
            Twiss alpha parameter
        beta : float
            Twiss beta parameter [m]
        
        Returns
        -------
        BeamDistribution2D
        """
        from optics import inverse_normalization_matrix
        
        T = inverse_normalization_matrix(alpha, beta)
        
        X_phys = T[0, 0] * self.x_grid + T[0, 1] * self.xp_grid
        Xp_phys = T[1, 0] * self.x_grid + T[1, 1] * self.xp_grid
        
        interp = RegularGridInterpolator(
            (self.x_grid, self.xp_grid), self.density,
            bounds_error=False, fill_value=0.0
        )
        
        X_grid, Xp_grid = np.meshgrid(X_phys, Xp_phys, indexing='ij')
        points = np.stack([X_grid.ravel(), Xp_grid.ravel()], axis=-1)
        density_phys = interp(points).reshape(X_grid.shape)
        
        return BeamDistribution2D(density_phys, X_phys, Xp_phys)