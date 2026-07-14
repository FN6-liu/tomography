# Beam Tomography for Accelerator Physics

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

**A Python implementation of transverse phase space tomography for charged particle beams using Filtered Back Projection (FBP) and Algebraic Reconstruction Technique (ART).**

This project implements beam tomography in particle accelerators, reconstructing the 2D transverse phase space from beam profile projections at different rotation angles in normalized phase space.

---

## Table of Contents

- [Background](#background)
- [Theory Overview](#theory-overview)
- [Project Structure](#project-structure)
- [Installation](#installation)
- [Usage](#usage)
- [Code Walkthrough](#code-walkthrough)
- [Results](#results)
- [Extensions to 4D](#extensions-to-4d)
- [References](#references)
- [License](#license)

---

## Background

In particle accelerators, the transverse phase space distribution of a beam (x, x') determines its quality and performance. Profile monitors only record the projection onto coordinate space. By varying quadrupole magnet strengths, the phase space rotates, allowing collection of projections at multiple angles – analogous to medical CT scans. Reconstruction from these projections is **beam tomography**.

This project implements:

- **Filtered Back Projection (FBP)**: Fast reconstruction with uniform angles.
- **Algebraic Reconstruction Technique (ART)**: Iterative, robust with limited angles, enforces non-negativity.

Reconstruction is performed in **normalized phase space**, where transport is a pure rotation.

---

## Theory Overview

### Phase Space and Normalization
A particle is described by (x, x'), where x' = dx/ds. The beamline optics uses Twiss parameters (α, β, γ). Normalized coordinates (x_N, x'_N) are obtained via a linear transformation. In this frame, any drift or quadrupole becomes a pure rotation by the betatron phase advance μ. Changing quadrupole strengths changes μ, providing different projection angles.

*For detailed mathematical derivations, see docs/theory.md.*

### Tomography Setup
- Unknown beam distribution at reconstruction point.
- Measure transverse beam profiles on a downstream screen for many quadrupole settings.
- Recover 2D density from 1D projections at different angles.

### FBP
Applies a ramp filter to each projection, then back-projects onto the image plane.

### ART
Solves the linear system iteratively, updating pixel values to reduce residuals for each ray.

---

## Project Structure

```
beam_tomography/
│
├── README.md                 # Project documentation
├── requirements.txt          # Python package dependencies
│
├── optics.py                 # Transport matrices, Twiss propagation, normalization
├── beam.py                   # Beam distribution class (Gaussian, multi-Gaussian)
├── tomography.py             # FBP and ART reconstruction algorithms
├── utils.py                  # Plotting and helper functions
├── main.py                   # Main script: generate data, reconstruct, visualize
│
├── docs/
│   └── theory.md             # Detailed mathematical derivations with LaTeX
│
├── examples/
│   └── example_output.png    # Example reconstruction results
│
└── tests/
    ├── test_optics.py        # Unit tests for optics functions
    ├── test_beam.py          # Unit tests for beam distribution
    └── test_tomography.py    # Unit tests for reconstruction algorithms
```

---

## Installation

1. Clone or download the repository.
2. Install dependencies:

```bash
pip install -r requirements.txt
```

Dependencies:
- numpy>=1.20.0
- scipy>=1.7.0
- matplotlib>=3.4.0

---

## Usage

Run the main script:

```bash
python main.py
```

This generates:
- Original phase space distribution
- Sinogram (projections vs. angle)
- FBP and ART reconstructions
- Residual plots and emittance values

### Customization

Edit main.py to change:

| Variable | Description | Default |
|----------|-------------|---------|
| alpha0, beta0 | Twiss parameters at reconstruction point | 0.0, 1.0 m |
| emittance | Geometric emittance | 1e-6 m·rad |
| n_angles | Number of projection angles | 18 |
| n_iter_art | ART iterations | 50 |
| grid_size | Phase space grid resolution | 64x64 |

---

## Code Walkthrough

### optics.py - Beam Transport and Optics

| Function | Description |
|----------|-------------|
| drift_matrix(L) | 2x2 drift transport matrix |
| quadrupole_matrix(f) | Thin quadrupole matrix (f > 0 focusing) |
| thick_quadrupole_matrix(k, L) | Thick quadrupole matrix |
| propagate_twiss(alpha, beta, M) | Propagate Twiss through transport matrix |
| normalization_matrix(alpha, beta) | Physical to normalized phase space |
| inverse_normalization_matrix(alpha, beta) | Normalized to physical phase space |
| normalized_rotation_angle(M_phys, alpha, beta) | Compute rotation angle in normalized space |
| compute_emittance(cov_matrix) | Emittance from covariance matrix |
| cov_to_twiss(cov_matrix) | Twiss parameters from covariance matrix |

### beam.py - Beam Distribution

**Class: BeamDistribution2D**

| Method | Description |
|--------|-------------|
| gaussian(mean, cov, x_grid, xp_grid) | Create Gaussian distribution |
| from_twiss(emittance, alpha, beta, x_grid, xp_grid) | Create from Twiss parameters |
| multi_gaussian(components, x_grid, xp_grid) | Create multi-Gaussian distribution |
| rotate(theta) | Rotate phase space by theta (rad) |
| project(theta) | Project onto x-axis after rotation |
| get_moments() | Compute covariance matrix |
| normalized_to_physical(alpha, beta) | Transform to physical space |

### tomography.py - Reconstruction Algorithms

| Function | Description |
|----------|-------------|
| fbp(projections, angles, x_grid, xp_grid, window) | Filtered Back Projection |
| art(projections, angles, x_grid, xp_grid, n_iter) | Algebraic Reconstruction Technique |

**FBP Steps:**
1. Fourier transform each projection
2. Apply ramp filter (|ω|) with optional window
3. Inverse Fourier transform
4. Back-project filtered projections

**ART Steps:**
1. Build sparse projection matrix P
2. Iteratively solve P·ρ = p using Kaczmarz method
3. Enforce non-negativity after each iteration

### utils.py - Visualization and Analysis

| Function | Description |
|----------|-------------|
| plot_phase_space(density, x_grid, xp_grid, title) | Plot 2D phase space with contours |
| plot_sinogram(projections, angles, title) | Plot all projections as sinogram |
| plot_comparison(original, fbp_recon, art_recon, ...) | Compare original vs reconstructions |
| compute_beam_parameters(density, x_grid, xp_grid) | Calculate emittance, Twiss from reconstruction |
| plot_residuals(original, reconstructed, title) | Plot residual differences |

### main.py - Main Script

1. Set beam and simulation parameters
2. Generate beam distribution (Gaussian or multi-Gaussian)
3. Create projection angles (0 to 180 degrees)
4. Compute projections (sinogram)
5. Reconstruct using FBP
6. Reconstruct using ART
7. Compare results visually
8. Compute and display beam parameters (emittance, Twiss)
9. Calculate reconstruction quality metrics (RMSE, correlation)

---

## Extensions to 4D

Extend to two transverse degrees of freedom (x, y, x', y'):

- Projections become 2D images (x, y)
- Need quadrupole settings rotating both x and y spaces
- Larger data structures (N^4 voxels for 4D phase space)
- Computational challenges: memory scales as N^4
- Can use 4D ART or machine learning approaches

### Implementation Challenges

| Challenge | Solution |
|-----------|----------|
| Memory | Use sparse matrices, compression (DCT) |
| Computation | GPU acceleration, parallel processing |
| Data collection | Automated quadrupole scans |
| Reconstruction | ML algorithms for faster processing |

See references for experimental implementations.

---

## References

- K. M. Hock et al., "Beam tomography in transverse normalised phase space", Nucl. Instrum. Methods Phys. Res. A 642, 36 (2011).
- A. Wolski et al., "Transverse phase space characterization in an accelerator test facility", Phys. Rev. Accel. Beams 23, 032804 (2020).
- A. Wolski et al., "Transverse phase space tomography in an accelerator test facility using image compression and machine learning", Phys. Rev. Accel. Beams 25, 122803 (2022).
- V. Guo et al., "4D Beam Tomography at the UCLA Pegasus Laboratory", IBIC2021, TUPP15 (2021).
- E. Bravin, "Transverse emittance", CERN Accelerator School (2021).

---

## License

This project is licensed under the MIT License.

---

## Contributing

Contributions welcome! Please open an issue or submit a pull request.

---

*Happy reconstructing!* 🚀
