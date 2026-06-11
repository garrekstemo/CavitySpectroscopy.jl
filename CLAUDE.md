# CavitySpectroscopy.jl

Cavity polariton spectroscopy analysis: coupled-oscillator dispersion,
Hopfield coefficients, Fabry-Pérot cavity transmittance, and fitting for
cavity spectra and polariton dispersion. Public; registration pending
(initial version 0.1.0). Extracted from QPSTools.jl 2026-06-11.

## Design rules

- Depends on **CurveFit** (SciML) and **CurveFitModels** directly — never on
  OpticalSpectroscopy or JASCOFiles. This package must stay independently
  registerable.
- Model functions follow `fn(p, x)` and are ForwardDiff-compatible
  (`similar(p, N)`, never `zeros(N)`; promote Hamiltonian eltypes).
- `CavitySpectrum` is generic over plain data `(x, y, metadata)`.
  Instrument-specific wrappers (e.g. QPSTools' JASCO-backed `CavitySpectrum`
  and its metadata-aware `fit_cavity_spectrum` dispatch) live downstream.
- Plotting stays out — QPSTools owns the Makie layer. (A neutral Makie
  extension is a natural 0.2 feature.)

## Physics conventions

- θ = ½·atan2(Ω, δ) with detuning δ = E_cav − E_vib.
- Hopfield: photon fraction of LP = sin²θ = ½(1 − δ/√(δ² + Ω²));
  LP becomes matter-like as δ → +∞. Pinned by eigenvector regression tests
  against direct 2×2 diagonalization (test/test_physics.jl) — do not weaken
  these.
- `polariton_eigenvalues` handles N coupled molecular modes via a Symmetric
  Hamiltonian.

## Structure

```
src/
  CavitySpectroscopy.jl  # module, exports
  physics.jl             # cavity_transmittance, dispersion, Hopfield
  types.jl               # CavitySpectrum, CavityFitResult, DispersionFitResult
  fitting.jl             # fit_cavity_spectrum, fit_dispersion
test/                    # physics regressions + synthetic fit round-trips + Aqua
docs/                    # Documenter site
```

## Downstream

QPSTools re-exports the physics vocabulary as a documented exception and
bridges `OpticalSpectroscopy.format_results` to the result types here. The
exports `CavitySpectrum`, `wavenumber`, `transmittance`, `format_results`
clash with sibling packages by design — consumers qualify or pick one.
