# Beam Tomography for Accelerator Physics

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

**A Python implementation of transverse phase space tomography for charged particle beams using Filtered Back Projection (FBP) and Algebraic Reconstruction Technique (ART).**

This project accompanies the study of beam tomography in particle accelerators, as discussed in the provided literature. It demonstrates how to reconstruct the 2D transverse phase space (one degree of freedom, e.g., horizontal) from a set of beam profile projections obtained at different rotation angles (betatron phase advances) in normalized phase space.

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

In particle accelerators, the transverse phase space distribution of a beam (e.g., \(x, x'\)) determines its quality and performance. Measuring the full distribution is challenging because profile monitors only record the projection onto coordinate space (the integral over angles). By varying the focusing strength of quadrupole magnets, the phase space rotates, allowing us to collect projections at multiple angles – a process analogous to medical CT scans. The reconstruction of the 2D (or 4D) phase space from these projections is known as **beam tomography**.

This project implements two standard reconstruction algorithms:

- **Filtered Back Projection (FBP)**: Fast and widely used when projections cover a full 180° range at uniform angular intervals.
- **Algebraic Reconstruction Technique (ART)**: Iterative and robust even with limited or non‑uniform angles; enforces non‑negativity.

Reconstruction is performed in **normalized phase space** (using Twiss parameters), where the transport between two points is a pure rotation, simplifying the geometry.

---

## Theory Overview

### Phase Space and Normalization
A particle in the horizontal plane is described by \((x, x')\), where \(x' = dx/ds\) is the angle relative to the reference trajectory. The linear optics of the beamline is characterized by the Twiss parameters \((\alpha, \beta, \gamma)\). The normalized coordinates are defined as:

\[
\begin{pmatrix} x_N \\ x'_N \end{pmatrix} =
\begin{pmatrix} 
1/\sqrt{\beta} & 0 \\
\alpha/\sqrt{\beta} & \sqrt{\beta}
\end{pmatrix}
\begin{pmatrix} x \\ x' \end{pmatrix}
\]

In this frame, a drift or a quadrupole corresponds to a pure rotation by the betatron phase advance \(\mu\). Changing the quadrupole strengths changes \(\mu\), effectively rotating the phase space and providing different projection angles.

### Tomography Setup
- The beam distribution at a **reconstruction point** (e.g., upstream of a set of quadrupoles) is unknown.
- By measuring the transverse beam profile on a downstream screen for many different quadrupole settings, we obtain projections at different angles.
- The goal is to recover the 2D density \(\rho(x_N, x'_N)\) from the set of 1D projections \(p_\theta(t) = \int \rho(t\cos\theta - u\sin\theta,\; t\sin\theta + u\cos\theta)\,du\).

### FBP
FBP applies a ramp filter in the frequency domain to each projection, then back‑projects the filtered projections onto the image plane. The reconstruction is:

\[
\rho(x, x') = \int_0^\pi q_\theta(x\cos\theta + x'\sin\theta)\,d\theta
\]

where \(q_\theta\) is the filtered projection.

### ART
ART solves the linear system \(\mathbf{P}\rho = \mathbf{p}\) iteratively. Each row of \(\mathbf{P}\) represents one ray (one projection bin). The Kaczmarz method updates the image pixel values to reduce the residual for each ray, one at a time, optionally with non‑negativity constraints.

---

## Project Structure

```
beam_tomography/
│
├── README.md                 # This file
├── requirements.txt          # Python dependencies
├── optics.py                 # Transport matrices, Twiss propagation, normalization
├── beam.py                   # Beam distribution class (Gaussian, multi‑peak, etc.)
├── tomography.py             # FBP and ART reconstruction functions
├── utils.py                  # Helper functions: plotting, emittance calculation, etc.
└── main.py                   # Main script: generate simulated data, reconstruct, plot
```

---

## Installation

1. **Clone the repository** (or download the source files).
2. Install the required packages (preferably in a virtual environment):

```bash
pip install -r requirements.txt
```

The dependencies are:
- `numpy`
- `scipy`
- `matplotlib`

---

## Usage

Run the main script to generate a synthetic beam, simulate projections, and compare FBP and ART reconstructions.

```bash
python main.py
```

You will see several figures:
- The original phase space distribution.
- The sinogram (all projections as a function of angle).
- Reconstructed distributions using FBP and ART.
- Residual plots and computed emittance values.

### Customizing the Simulation

You can modify the beam parameters, the number of projections, and the reconstruction settings inside `main.py`. Key variables include:

- `alpha0`, `beta0`: Twiss parameters at the reconstruction point.
- `emittance`: Geometric emittance (in m·rad).
- `n_angles`: Number of projection angles (default 18).
- `n_iter_art`: Iteration count for ART (default 50).

---

## Code Walkthrough

### `optics.py`
- `drift_matrix(L)`, `quad_matrix(f)`: 2×2 transport matrices.
- `propagate_twiss(alpha, beta, M)`: compute Twiss at a new location given the transfer matrix.
- `normalization_matrix(alpha, beta)`: returns the transform from physical to normalized coordinates.
- `rotation_angle(M, alpha0, beta0)`: given a physical transport matrix and initial Twiss, returns the rotation angle in normalized space.

### `beam.py`
- `BeamDistribution2D`: a class that holds a discrete phase space density on a 2D grid.
- Methods:
  - `project(theta)`: projects the distribution onto the x‑axis after rotating by `theta` (radians) in normalized space. Returns a 1D array.

### `tomography.py`
- `fbp(projections, angles, grid)`: performs FBP reconstruction. It applies a ramp filter (with optional window) and back‑projects.
- `art(projections, angles, grid, n_iter)`: performs ART by constructing a sparse projection matrix and iterating with the Kaczmarz algorithm.

### `main.py`
- Generates a test beam (Gaussian or user‑defined).
- Computes projections for the chosen angles.
- Reconstructs using both methods.
- Plots and compares results.

---

## Extensions to 4D

The same principles can be extended to **two transverse degrees of freedom** (x, y, x', y'), i.e., 4D phase space. This requires:

- More complex projection data: on a screen we observe 2D images (x,y), which are projections of the 4D distribution.
- A set of quadrupole settings that rotates both the x and y phase spaces.
- Larger data structures and reconstruction techniques (e.g., 4D ART or machine learning).

The framework presented here can serve as a starting point. For experimental results and advanced methods (including machine learning), see the referenced literature.

---

## References

- K. M. Hock et al., “Beam tomography in transverse normalised phase space”, *Nucl. Instrum. Methods Phys. Res. A* **642**, 36 (2011).
- A. Wolski et al., “Transverse phase space characterization in an accelerator test facility”, *Phys. Rev. Accel. Beams* **23**, 032804 (2020).
- A. Wolski et al., “Transverse phase space tomography in an accelerator test facility using image compression and machine learning”, *Phys. Rev. Accel. Beams* **25**, 122803 (2022).

---

## License

This project is licensed under the MIT License – see the [LICENSE](LICENSE) file for details.

---

## Contributing

Contributions and suggestions are welcome! Please open an issue or submit a pull request.

---

*Happy reconstructing!* 🚀
