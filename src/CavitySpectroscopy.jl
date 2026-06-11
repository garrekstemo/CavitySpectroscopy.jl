"""
    CavitySpectroscopy

Fabry-Pérot cavity and polariton spectroscopy analysis.

Provides the physics and fitting chain for vibrational strong coupling (VSC)
and polariton experiments:

- Fabry-Pérot transmittance with an absorbing intracavity medium
  ([`cavity_transmittance`](@ref), [`compute_cavity_transmittance`](@ref))
- Angle-dependent cavity mode dispersion ([`cavity_mode_energy`](@ref))
- Coupled oscillator polariton model ([`polariton_branches`](@ref),
  [`polariton_eigenvalues`](@ref), [`hopfield_coefficients`](@ref))
- Spectrum and dispersion fitting ([`fit_cavity_spectrum`](@ref),
  [`fit_dispersion`](@ref)) built on CurveFit.jl

Data enters as plain vectors or as a [`CavitySpectrum`](@ref)
(x, y, metadata). Model functions follow the `fn(p, x)` convention and are
ForwardDiff-compatible.
"""
module CavitySpectroscopy

using LinearAlgebra: Symmetric, eigvals
using Statistics: mean
using CurveFit: NonlinearCurveFitProblem, solve, coef, stderror
import CurveFit: predict, residuals
using CurveFitModels: dielectric_real, dielectric_imag

include("physics.jl")

# Physics
export cavity_transmittance, compute_cavity_transmittance
export refractive_index, extinction_coeff
export cavity_mode_energy, polariton_branches, polariton_eigenvalues
export hopfield_coefficients

end # module
