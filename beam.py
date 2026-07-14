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
        self.density = np.array(density, dtype=float, copy=True)
        self.x_grid = np.array(x_grid, dtype=float, copy=True)
        self.xp_grid = np.array(xp_grid, dtype=float, copy=True)
        self.dx = self.x_grid[1] - self.x_grid[0]
        self.dxp = self.xp_grid[1] - self.xp_grid[0]

        # Normalize density as a continuous PDF on the phase-space grid.
        total = np.sum(self.density) * self.dx * self.dxp
        if total <= 0:
            raise ValueError("Density must have positive total integral.")
        self.density /= total
    
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
        ct = np.cos(theta)
        st = np.sin(theta)

        # rho_rot(x, x') = rho(R^-1 [x, x'])
        X_grid, Xp_grid = np.meshgrid(self.x_grid, self.xp_grid, indexing='ij')
        X_src = ct * X_grid + st * Xp_grid
        Xp_src = -st * X_grid + ct * Xp_grid

        interp = RegularGridInterpolator(
            (self.x_grid, self.xp_grid), self.density,
            bounds_error=False, fill_value=0.0
        )

        points = np.stack([X_src.ravel(), Xp_src.ravel()], axis=-1)
        density_rot = interp(points).reshape(X_grid.shape)

        return BeamDistribution2D(density_rot, self.x_grid, self.xp_grid)
    
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
        total = np.sum(self.density) * self.dx * self.dxp
        
        x_mesh, xp_mesh = np.meshgrid(self.x_grid, self.xp_grid, indexing='ij')
        mean_x = np.sum(x_mesh * self.density) * self.dx * self.dxp / total
        mean_xp = np.sum(xp_mesh * self.density) * self.dx * self.dxp / total
        
        x2 = np.sum((x_mesh - mean_x)**2 * self.density) * self.dx * self.dxp / total
        xp2 = np.sum((xp_mesh - mean_xp)**2 * self.density) * self.dx * self.dxp / total
        xxp = np.sum((x_mesh - mean_x) * (xp_mesh - mean_xp) * self.density) * self.dx * self.dxp / total
        
        return np.array([[x2, xxp], [xxp, xp2]])
    
    def transform(self, matrix, output_x_grid=None, output_xp_grid=None):
        """
        Apply a linear coordinate transform y = M x and resample onto
        a regular output grid.

        Parameters
        ----------
        matrix : np.ndarray
            2x2 transform matrix from source to target coordinates
        output_x_grid : np.ndarray, optional
            Target x grid
        output_xp_grid : np.ndarray, optional
            Target x' grid

        Returns
        -------
        BeamDistribution2D
        """
        matrix = np.asarray(matrix, dtype=float)
        matrix_inv = np.linalg.inv(matrix)

        if output_x_grid is None or output_xp_grid is None:
            corners = np.array([
                [self.x_grid[0], self.xp_grid[0]],
                [self.x_grid[0], self.xp_grid[-1]],
                [self.x_grid[-1], self.xp_grid[0]],
                [self.x_grid[-1], self.xp_grid[-1]],
            ])
            transformed = corners @ matrix.T
            if output_x_grid is None:
                output_x_grid = np.linspace(
                    transformed[:, 0].min(), transformed[:, 0].max(), len(self.x_grid)
                )
            if output_xp_grid is None:
                output_xp_grid = np.linspace(
                    transformed[:, 1].min(), transformed[:, 1].max(), len(self.xp_grid)
                )

        output_x_grid = np.array(output_x_grid, dtype=float, copy=True)
        output_xp_grid = np.array(output_xp_grid, dtype=float, copy=True)

        interp = RegularGridInterpolator(
            (self.x_grid, self.xp_grid), self.density,
            bounds_error=False, fill_value=0.0
        )

        X_out, Xp_out = np.meshgrid(output_x_grid, output_xp_grid, indexing='ij')
        target_points = np.stack([X_out.ravel(), Xp_out.ravel()], axis=-1)
        source_points = target_points @ matrix_inv.T
        density_out = interp(source_points).reshape(X_out.shape)

        return BeamDistribution2D(density_out, output_x_grid, output_xp_grid)

    def physical_to_normalized(self, alpha, beta, x_grid=None, xp_grid=None):
        """
        Transform from physical phase space to normalized phase space.

        Parameters
        ----------
        alpha : float
            Twiss alpha parameter
        beta : float
            Twiss beta parameter [m]
        x_grid : np.ndarray, optional
            Target normalized x grid
        xp_grid : np.ndarray, optional
            Target normalized x' grid

        Returns
        -------
        BeamDistribution2D
        """
        from optics import normalization_matrix

        T = normalization_matrix(alpha, beta)
        return self.transform(T, x_grid, xp_grid)

    def normalized_to_physical(self, alpha, beta, x_grid=None, xp_grid=None):
        """
        Transform from normalized phase space to physical phase space.
        
        Parameters
        ----------
        alpha : float
            Twiss alpha parameter
        beta : float
            Twiss beta parameter [m]
        x_grid : np.ndarray, optional
            Target physical x grid
        xp_grid : np.ndarray, optional
            Target physical x' grid
        
        Returns
        -------
        BeamDistribution2D
        """
        from optics import inverse_normalization_matrix
        
        T = inverse_normalization_matrix(alpha, beta)
        return self.transform(T, x_grid, xp_grid)
