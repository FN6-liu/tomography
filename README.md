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
