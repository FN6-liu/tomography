markdown
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

In particle accelerators, the transverse phase space distribution of a beam (e.g., horizontal coordinate x and angle x') determines its quality and performance. Measuring the full distribution is challenging because profile monitors only record the projection onto coordinate space (the integral over angles). By varying the focusing strength of quadrupole magnets, the phase space rotates, allowing us to collect projections at multiple angles – a process analogous to medical CT scans. The reconstruction of the 2D (or 4D) phase space from these projections is known as **beam tomography**.

This project implements two standard reconstruction algorithms:

- **Filtered Back Projection (FBP)**: Fast and widely used when projections cover a full 180° range at uniform angular intervals.
- **Algebraic Reconstruction Technique (ART)**: Iterative and robust even with limited or non‑uniform angles; enforces non‑negativity.

Reconstruction is performed in **normalized phase space** (using Twiss parameters), where the transport between two points is a pure rotation, simplifying the geometry.

---

## Theory Overview

### Phase Space and Normalization
A particle in the horizontal plane is described by (x, x'), where x' = dx/ds is the angle relative to the reference trajectory. The linear optics of the beamline is characterized by the Twiss parameters (α, β, γ). The normalized coordinates (x_N, x'_N) are obtained by a linear transformation that depends on α and β. In this normalized frame, any drift or quadrupole corresponds to a pure rotation by the betatron phase advance μ. Changing the quadrupole strengths changes μ, thus rotating the phase space and yielding different projection angles.

*(For the exact mathematical definitions in LaTeX, please refer to the [detailed theory document](docs/theory.md).)*

### Tomography Setup
- The beam distribution at a **reconstruction point** (e.g., upstream of a set of quadrupoles) is unknown.
- By measuring the transverse beam profile on a downstream screen for many different quadrupole settings, we obtain projections at different angles.
- The goal is to recover the 2D density ρ(x_N, x'_N) from the set of 1D projections p_θ(t), where t is the coordinate along the projection axis.

### FBP (Filtered Back Projection)
FBP applies a ramp filter (in the frequency domain) to each projection, then back‑projects the filtered projections onto the image plane. The reconstruction is the integral over all angles of the filtered projection evaluated at the corresponding ray position.

### ART (Algebraic Reconstruction Technique)
ART solves the linear system **P·ρ = p** iteratively. Each row of P represents one ray (one projection bin). The Kaczmarz method updates the image pixel values to reduce the residual for each ray, one at a time, optionally enforcing non‑negativity.

---
markdown
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

In particle accelerators, the transverse phase space distribution of a beam (e.g., horizontal coordinate x and angle x') determines its quality and performance. Measuring the full distribution is challenging because profile monitors only record the projection onto coordinate space (the integral over angles). By varying the focusing strength of quadrupole magnets, the phase space rotates, allowing us to collect projections at multiple angles – a process analogous to medical CT scans. The reconstruction of the 2D (or 4D) phase space from these projections is known as **beam tomography**.

This project implements two standard reconstruction algorithms:

- **Filtered Back Projection (FBP)**: Fast and widely used when projections cover a full 180° range at uniform angular intervals.
- **Algebraic Reconstruction Technique (ART)**: Iterative and robust even with limited or non‑uniform angles; enforces non‑negativity.

Reconstruction is performed in **normalized phase space** (using Twiss parameters), where the transport between two points is a pure rotation, simplifying the geometry.

---

## Theory Overview

### Phase Space and Normalization
A particle in the horizontal plane is described by (x, x'), where x' = dx/ds is the angle relative to the reference trajectory. The linear optics of the beamline is characterized by the Twiss parameters (α, β, γ). The normalized coordinates (x_N, x'_N) are obtained by a linear transformation that depends on α and β. In this normalized frame, any drift or quadrupole corresponds to a pure rotation by the betatron phase advance μ. Changing the quadrupole strengths changes μ, thus rotating the phase space and yielding different projection angles.

*(For the exact mathematical definitions in LaTeX, please refer to the [detailed theory document](docs/theory.md).)*

### Tomography Setup
- The beam distribution at a **reconstruction point** (e.g., upstream of a set of quadrupoles) is unknown.
- By measuring the transverse beam profile on a downstream screen for many different quadrupole settings, we obtain projections at different angles.
- The goal is to recover the 2D density ρ(x_N, x'_N) from the set of 1D projections p_θ(t), where t is the coordinate along the projection axis.

### FBP (Filtered Back Projection)
FBP applies a ramp filter (in the frequency domain) to each projection, then back‑projects the filtered projections onto the image plane. The reconstruction is the integral over all angles of the filtered projection evaluated at the corresponding ray position.

### ART (Algebraic Reconstruction Technique)
ART solves the linear system **P·ρ = p** iteratively. Each row of P represents one ray (one projection bin). The Kaczmarz method updates the image pixel values to reduce the residual for each ray, one at a time, optionally enforcing non‑negativity.

---
beam_tomography/
│
├── README.md # This file
├── requirements.txt # Python dependencies
├── optics.py # Transport matrices, Twiss propagation, normalization
├── beam.py # Beam distribution class (Gaussian, multi‑peak, etc.)
├── tomography.py # FBP and ART reconstruction functions
├── utils.py # Helper functions: plotting, emittance calculation, etc.
├── main.py # Main script: generate simulated data, reconstruct, plot
└── docs/
└── theory.md # Detailed mathematical derivations with LaTeX
---

## Installation

1. **Clone the repository** (or download the source files).
2. Install the required packages (preferably in a virtual environment):

```bash
pip install -r requirements.txt
Usage
Run the main script to generate a synthetic beam, simulate projections, and compare FBP and ART reconstructions.

bash
python main.py
You will see several figures:

The original phase space distribution.

The sinogram (all projections as a function of angle).

Reconstructed distributions using FBP and ART.

Residual plots and computed emittance values.

Customizing the Simulation
You can modify the beam parameters, the number of projections, and the reconstruction settings inside main.py. Key variables include:

alpha0, beta0: Twiss parameters at the reconstruction point.

emittance: Geometric emittance (in m·rad).

n_angles: Number of projection angles (default 18).

n_iter_art: Iteration count for ART (default 50).

Code Walkthrough
optics.py
drift_matrix(L), quad_matrix(f): 2×2 transport matrices.

propagate_twiss(alpha, beta, M): compute Twiss at a new location given the transfer matrix.

normalization_matrix(alpha, beta): returns the transform from physical to normalized coordinates.

rotation_angle(M, alpha0, beta0): given a physical transport matrix and initial Twiss, returns the rotation angle in normalized space.

beam.py
BeamDistribution2D: a class that holds a discrete phase space density on a 2D grid.

Methods:

project(theta): projects the distribution onto the x‑axis after rotating by theta (radians) in normalized space. Returns a 1D array.

tomography.py
fbp(projections, angles, grid): performs FBP reconstruction. It applies a ramp filter (with optional window) and back‑projects.

art(projections, angles, grid, n_iter): performs ART by constructing a sparse projection matrix and iterating with the Kaczmarz algorithm.

main.py
Generates a test beam (Gaussian or user‑defined).

Computes projections for the chosen angles.

Reconstructs using both methods.

Plots and compares results.

Results
An example output for a Gaussian beam with 18 projections (10° step) is shown below:

Original	FBP Reconstruction	ART Reconstruction
https://docs/original.png	https://docs/fbp.png	https://docs/art.png
FBP produces a smooth reconstruction with typical streaking artefacts when angles are limited.

ART often preserves edges and details better, especially with few projections, but may be slower.

The emittance (computed from the covariance matrix of the reconstructed distribution) typically agrees within a few percent of the true value when using ≥18 angles.

Extensions to 4D
The same principles can be extended to two transverse degrees of freedom (x, y, x', y'), i.e., 4D phase space. This requires:

More complex projection data: on a screen we observe 2D images (x,y), which are projections of the 4D distribution.

A set of quadrupole settings that rotates both the x and y phase spaces.

Larger data structures and reconstruction techniques (e.g., 4D ART or machine learning).

The framework presented here can serve as a starting point. For experimental results and advanced methods (including machine learning), see the referenced literature.

References
K. M. Hock et al., “Beam tomography in transverse normalised phase space”, Nucl. Instrum. Methods Phys. Res. A 642, 36 (2011).

A. Wolski et al., “Transverse phase space characterization in an accelerator test facility”, Phys. Rev. Accel. Beams 23, 032804 (2020).

A. Wolski et al., “Transverse phase space tomography in an accelerator test facility using image compression and machine learning”, Phys. Rev. Accel. Beams 25, 122803 (2022).

License
This project is licensed under the MIT License – see the LICENSE file for details.

Contributing
Contributions and suggestions are welcome! Please open an issue or submit a pull request.

Happy reconstructing! 🚀

text

---

## 文件 2：`docs/theory.md`（包含完整的 LaTeX 公式）

```markdown
# Mathematical Details of Beam Tomography

This document provides the full mathematical formulation used in the beam tomography simulation and reconstruction. It is intended for readers who wish to understand the exact equations behind the code.

## 1. Phase Space and Normalization

A single particle in the transverse horizontal plane is described by the coordinates:

\[
(x, x') \quad \text{where } x' = \frac{dx}{ds}
\]

with \(s\) the longitudinal coordinate along the reference trajectory.

The linear optics of the beamline is characterized by the Courant–Snyder (Twiss) parameters \((\alpha, \beta, \gamma)\), which satisfy:

\[
\gamma = \frac{1 + \alpha^2}{\beta}
\]

The **normalized phase space** coordinates \((x_N, x'_N)\) are defined by:

\[
\begin{pmatrix} x_N \\ x'_N \end{pmatrix}
=
\begin{pmatrix}
1/\sqrt{\beta} & 0 \\
\alpha/\sqrt{\beta} & \sqrt{\beta}
\end{pmatrix}
\begin{pmatrix} x \\ x' \end{pmatrix}
\]

In this frame, the transport matrix for any linear element (drift, quadrupole) becomes a pure rotation:

\[
\mathbf{M}_N = \begin{pmatrix}
\cos \mu & \sin \mu \\
-\sin \mu & \cos \mu
\end{pmatrix}
\]

where \(\mu\) is the betatron phase advance.

## 2. The Tomography Problem

Let \(\rho(x_N, x'_N)\) be the unknown beam density in normalized phase space. A projection at angle \(\theta\) is the integral along lines perpendicular to the projection direction:

\[
p_\theta(t) = \int_{-\infty}^{\infty} \rho\left(t\cos\theta - u\sin\theta,\; t\sin\theta + u\cos\theta\right) du
\]

Here \(t\) is the coordinate on the detector (the projection axis), and \(\theta\) is the rotation angle of the phase space relative to the detector.

Given a set of projections \(p_\theta(t)\) for multiple \(\theta\), the goal is to reconstruct \(\rho\).

## 3. Filtered Back Projection (FBP)

The FBP reconstruction formula is:

\[
\rho(x_N, x'_N) = \int_0^\pi q_\theta\left(x_N\cos\theta + x'_N\sin\theta\right) d\theta
\]

where \(q_\theta(t)\) is the filtered projection:

\[
q_\theta(t) = \mathcal{F}^{-1}\left\{ |\omega| \, \mathcal{F}\{p_\theta\}(\omega) \right\}(t)
\]

Here \(\mathcal{F}\) denotes the Fourier transform, and \(|\omega|\) is the ramp filter. In practice, the filter is often windowed (e.g., with a Hann window) to reduce noise.

## 4. Algebraic Reconstruction Technique (ART)

We discretize the phase space into an \(N \times N\) grid, with the density represented as a vector \(\boldsymbol{\rho} \in \mathbb{R}^{N^2}\). Each projection bin provides a linear equation:

\[
\sum_{j} P_{ij} \rho_j = p_i
\]

where \(P_{ij}\) is the contribution of pixel \(j\) to ray \(i\) (usually the length of intersection or a weighting factor). The full system is:

\[
\mathbf{P} \boldsymbol{\rho} = \mathbf{p}
\]

ART solves this iteratively. The Kaczmarz method updates the solution by cycling through the rows of \(\mathbf{P}\):

\[
\boldsymbol{\rho}^{(k+1)} = \boldsymbol{\rho}^{(k)} + \frac{p_i - \mathbf{P}_{i:} \cdot \boldsymbol{\rho}^{(k)}}{\|\mathbf{P}_{i:}\|_2^2} \mathbf{P}_{i:}^\mathrm{T}
\]

Optionally, a non‑negativity constraint \(\rho_j \ge 0\) is applied after each full iteration.

## 5. Emittance and Twiss Parameters from Reconstructed Distribution

Given the reconstructed density \(\rho(x, x')\) in physical space (after inverse normalization), the second moments are computed as:

\[
\langle x^2 \rangle = \iint x^2 \rho(x, x') dx dx'
\]
\[
\langle x'^2 \rangle = \iint x'^2 \rho(x, x') dx dx'
\]
\[
\langle x x' \rangle = \iint x x' \rho(x, x') dx dx'
\]

The geometric emittance is:

\[
\epsilon = \sqrt{ \langle x^2 \rangle \langle x'^2 \rangle - \langle x x' \rangle^2 }
\]

The Twiss parameters at the reconstruction point are then:

\[
\beta = \frac{\langle x^2 \rangle}{\epsilon}, \quad \alpha = -\frac{\langle x x' \rangle}{\epsilon}, \quad \gamma = \frac{\langle x'^2 \rangle}{\epsilon}
\]

---

*This document supplements the main README and provides the mathematical backbone of the implementation.*
将这两个文件放入您的项目仓库即可。README.md 中已经链接到 docs/theory.md，方便需要详细公式的用户。如果您的平台支持 LaTeX 渲染（如 GitLab 或带有插件的 Markdown 阅读器），theory.md 中的公式将正确显示。

