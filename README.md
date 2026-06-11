# CavitySpectroscopy.jl

[![CI](https://github.com/garrekstemo/CavitySpectroscopy.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/garrekstemo/CavitySpectroscopy.jl/actions/workflows/CI.yml)
[![docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://garrekstemo.github.io/CavitySpectroscopy.jl/dev/)
[![codecov](https://codecov.io/gh/garrekstemo/CavitySpectroscopy.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/garrekstemo/CavitySpectroscopy.jl)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![license](https://img.shields.io/github/license/garrekstemo/CavitySpectroscopy.jl)](LICENSE)

CavitySpectroscopy.jl analyzes Fabry-Pérot cavity transmission spectra and
polariton dispersion for light-matter strong coupling experiments
(vibrational strong coupling, exciton polaritons). It covers the full chain
from raw transmittance data to Rabi splitting and Hopfield coefficients:

1. Multi-oscillator Lorentz dielectric function (via CurveFitModels.jl)
2. Complex refractive index `n`, `k` from the dielectric function
3. Fabry-Pérot Airy transmittance with an absorbing intracavity medium
4. Coupled oscillator polariton model: branches, eigenvalues, mixing fractions
5. Nonlinear least-squares fitting (via CurveFit.jl) of spectra and dispersion

## Installation

```julia
using Pkg
Pkg.add("CavitySpectroscopy")
```

## Quick start

### Polariton physics

```julia
using CavitySpectroscopy

# Cavity photon energy vs incidence angle (radians)
E_cav = cavity_mode_energy([2040.0, 1.5], deg2rad.(0:5:30))

# Upper and lower polariton branches (2-level coupled oscillator model)
LP, UP = polariton_branches(E_cav, 2055.0, 25.0)

# Light-matter mixing fractions
h = hopfield_coefficients(E_cav, 2055.0, 25.0)
h.photon_LP   # photon fraction of the lower polariton at each angle

# N molecular modes coupled to one cavity mode
eigs = polariton_eigenvalues(2050.0, [2030.0, 2060.0], [15.0, 20.0])
```

### Fitting a cavity transmission spectrum

```julia
result = fit_cavity_spectrum(nu, T;
    oscillators = [(nu0 = 2055.0, Gamma = 23.0)],
    L = 12.0e-4,      # cavity length (cm)
    n_bg = 1.4)       # background refractive index

result.R                  # fitted mirror reflectance
result.polariton_peaks    # auto-extracted peak positions
predict(result)           # fitted curve on the data grid
```

### Fitting polariton dispersion

```julia
result = fit_dispersion(angles, lp_positions, up_positions;
    molecular_modes = 2055.0)

result.rabi_splitting     # Rabi splitting (with result.rabi_err)
result.hopfield_zero      # mixing fractions at zero detuning
```

Data can also be carried in a `CavitySpectrum(x, y; metadata...)`, which
`fit_cavity_spectrum` accepts directly (cavity length read from metadata,
percent transmittance normalized automatically).

## Conventions

- Wavenumber units (cm⁻¹) throughout; angles in radians.
- Model functions follow the `fn(p, x)` signature and are
  ForwardDiff-compatible, so they can be used directly with CurveFit.jl.
- Hopfield convention: with detuning `δ = E_cav − E_vib` and
  `θ = ½·atan(Ω, δ)`, the photon fraction of the lower polariton is
  `sin²θ = ½(1 − δ/√(δ² + Ω²))` — at far positive detuning the LP converges
  to the bare vibration (matter-like), verified against direct Hamiltonian
  diagonalization in the test suite.

## Related packages

- [TransferMatrix.jl](https://github.com/garrekstemo/TransferMatrix.jl) —
  general 4×4 transfer-matrix optics for multilayer structures
- [CurveFitModels.jl](https://github.com/garrekstemo/CurveFitModels.jl) —
  pure model functions (Lorentz oscillator, lineshapes)
- [JASCOFiles.jl](https://github.com/garrekstemo/JASCOFiles.jl) —
  read JASCO FTIR/Raman/UV-Vis files to get `x`/`y` data for fitting
