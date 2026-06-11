# Cavity and polariton physics.
#
# Physics chain for the transmittance model:
# 1. Multi-oscillator dielectric function (CurveFitModels: `dielectric_real`,
#    `dielectric_imag`)
# 2. Complex refractive index from the dielectric function
#    (`refractive_index`, `extinction_coeff`)
# 3. Absorption coefficient from the extinction coefficient
# 4. Fabry-Perot Airy function (`cavity_transmittance`)
#
# All functions are pure and ForwardDiff-compatible.

"""
    refractive_index(eps1, eps2)

Compute refractive index n from real (eps1) and imaginary (eps2) parts of
the dielectric function.

n = sqrt((sqrt(eps1^2 + eps2^2) + eps1) / 2)
"""
function refractive_index(eps1, eps2)
    @. sqrt(0.5 * (eps1 + sqrt(eps1^2 + eps2^2)))
end

"""
    extinction_coeff(eps1, eps2)

Compute extinction coefficient k from real (eps1) and imaginary (eps2) parts of
the dielectric function.

k = sqrt((sqrt(eps1^2 + eps2^2) - eps1) / 2)
"""
function extinction_coeff(eps1, eps2)
    @. sqrt(0.5 * (-eps1 + sqrt(eps1^2 + eps2^2)))
end

"""
    cavity_transmittance(p, ν)

Fabry-Perot cavity transmittance with an absorbing medium as a function of frequency.

# Arguments
- `p`: Parameters [n, α, L, R, ϕ]
  - `n`: Refractive index
  - `α`: Absorption coefficient
  - `L`: Cavity length
  - `R`: Mirror reflectance (T = 1 - R assumed)
  - `ϕ`: Phase shift upon reflection
- `ν`: Frequency (independent variable)

```math
\\begin{aligned}
    T(\\nu) = \\frac{(1-R)^2 e^{-\\alpha L}}{1 + R^2 e^{-2\\alpha L} - 2R e^{-\\alpha L} \\cos(4\\pi n L \\nu + 2\\phi)}
\\end{aligned}
```

[https://en.wikipedia.org/wiki/Fabry%E2%80%93P%C3%A9rot_interferometer](https://en.wikipedia.org/wiki/Fabry%E2%80%93P%C3%A9rot_interferometer)
"""
function cavity_transmittance(p, ν)
    n, α, L, R, ϕ = p[1], p[2], p[3], p[4], p[5]
    T = 1 - R
    e = exp(-α * L)
    @. T^2 * e / (1 + R^2 * e^2 - 2 * R * e * cos(4π * n * L * ν + 2 * ϕ))
end

"""
    compute_cavity_transmittance(nu::Number, nu0s, Gammas, As, R, L, n_bg, phi)

Compute cavity transmittance at a single frequency for multiple Lorentzian oscillators.

Builds the full physics chain:
1. Sum Lorentz oscillator dielectric contributions
2. Compute complex refractive index (n, k) from dielectric function
   ([`refractive_index`](@ref), [`extinction_coeff`](@ref))
3. Compute absorption coefficient alpha = 4pi * k * nu
4. Feed into the Fabry-Perot Airy function ([`cavity_transmittance`](@ref))

# Arguments
- `nu`: Frequency (cm^-1)
- `nu0s`: Vector of oscillator center frequencies (cm^-1)
- `Gammas`: Vector of oscillator linewidths (cm^-1)
- `As`: Vector of oscillator amplitudes/strengths (cm^-2)
- `R`: Mirror reflectivity
- `L`: Cavity length (cm)
- `n_bg`: Background refractive index
- `phi`: Phase shift upon reflection

ForwardDiff-compatible.
"""
function compute_cavity_transmittance(nu::Number, nu0s, Gammas, As, R, L, n_bg, phi)
    # Build dielectric function from all oscillators
    eps1 = n_bg^2
    eps2 = zero(eps1)

    for i in eachindex(nu0s, Gammas, As)
        eps1 += dielectric_real((As[i], nu0s[i], Gammas[i]), nu)
        eps2 += dielectric_imag((As[i], nu0s[i], Gammas[i]), nu)
    end

    # Complex refractive index
    n = refractive_index(eps1, eps2)
    k = extinction_coeff(eps1, eps2)

    # Absorption coefficient
    alpha = 4π * k * nu

    # Fabry-Perot transmittance
    return cavity_transmittance((n, alpha, L, R, phi), nu)
end

"""
    compute_cavity_transmittance(nus::AbstractArray, nu0s, Gammas, As, R, L, n_bg, phi)

Array dispatch: compute cavity transmittance for multiple frequencies.
"""
function compute_cavity_transmittance(nus::AbstractArray, nu0s, Gammas, As, R, L, n_bg, phi)
    return [compute_cavity_transmittance(nu, nu0s, Gammas, As, R, L, n_bg, phi) for nu in nus]
end

"""
    cavity_mode_energy(p, thetas)

Compute cavity photon energy as a function of incidence angle.

E_cav(theta) = E0 / sqrt(1 - (sin(theta) / n_eff)^2)

# Arguments
- `p`: Parameters [E0, n_eff] where E0 is normal-incidence cavity energy
  and n_eff is effective refractive index
- `thetas`: Incidence angles (radians)
"""
function cavity_mode_energy(p, thetas)
    E0, n_eff = p[1], p[2]
    @. E0 / sqrt(1 - (sin(thetas) / n_eff)^2)
end

"""
    polariton_branches(E_cav, E_vib, Omega)

Compute upper and lower polariton energies from the 2-level coupled oscillator model.

E_pm = (E_cav + E_vib) / 2 +/- sqrt(Omega^2 + (E_cav - E_vib)^2) / 2

# Arguments
- `E_cav`: Cavity photon energy (scalar or vector)
- `E_vib`: Vibrational mode energy (scalar)
- `Omega`: Rabi splitting (scalar)

# Returns
`(E_LP, E_UP)` — lower and upper polariton energies, same shape as `E_cav`.
"""
function polariton_branches(E_cav, E_vib, Omega)
    delta = @. sqrt(Omega^2 + (E_cav - E_vib)^2)
    E_LP = @. 0.5 * (E_vib + E_cav - delta)
    E_UP = @. 0.5 * (E_vib + E_cav + delta)
    return E_LP, E_UP
end

"""
    polariton_eigenvalues(E_cav, E_vibs, Omegas)

Compute polariton energies for N vibrational modes coupled to one cavity mode.

Builds the (N+1) x (N+1) Hamiltonian and returns sorted eigenvalues.
For a single vibrational mode, this reduces to [`polariton_branches`](@ref).

The Hamiltonian eltype is promoted from the inputs, so dual numbers
propagate through (ForwardDiff provides `eigvals` for `Symmetric` dual
matrices), which makes the multi-mode dispersion fit differentiable.

# Arguments
- `E_cav`: Cavity photon energy (scalar)
- `E_vibs`: Vector of vibrational mode energies
- `Omegas`: Vector of Rabi splittings (one per mode)

# Returns
Sorted vector of N+1 eigenvalues (polariton energies).
"""
function polariton_eigenvalues(E_cav, E_vibs, Omegas)
    N = length(E_vibs)
    @assert length(Omegas) == N "Need one Rabi splitting per vibrational mode"

    T = float(promote_type(typeof(E_cav), eltype(E_vibs), eltype(Omegas)))
    H = zeros(T, N + 1, N + 1)
    H[1, 1] = E_cav
    for i in eachindex(E_vibs)
        H[i + 1, i + 1] = E_vibs[i]
        H[1, i + 1] = Omegas[i] / 2
        H[i + 1, 1] = Omegas[i] / 2
    end

    return sort!(eigvals(Symmetric(H)))
end

"""
    hopfield_coefficients(E_cav, E_vib, Omega)

Compute Hopfield coefficients (light-matter mixing fractions) for the
2-level coupled oscillator model.

At a given detuning (E_cav - E_vib), returns the photonic and matter
fractions for the lower and upper polariton branches.

# Arguments
- `E_cav`: Cavity photon energy (scalar or vector)
- `E_vib`: Vibrational mode energy (scalar)
- `Omega`: Rabi splitting (scalar)

# Returns
Named tuple `(photon_LP, matter_LP, photon_UP, matter_UP)`.
Each element has the same shape as `E_cav`.
Values satisfy: photon + matter = 1 for each branch.

# Convention
With `theta = 0.5 * atan(Omega, delta)` and `delta = E_cav - E_vib`, the
photon fraction of the LP is `sin^2(theta) = (1 - delta/sqrt(delta^2 + Omega^2))/2`:
at far positive detuning (E_cav >> E_vib) the LP converges to the bare
vibration, so its photon fraction must vanish (theta -> 0, sin^2 -> 0).
Verified against direct Hamiltonian diagonalization in the test suite.
"""
function hopfield_coefficients(E_cav, E_vib, Omega)
    delta = @. E_cav - E_vib
    theta = @. 0.5 * atan(Omega, delta)

    photon_LP = @. sin(theta)^2
    matter_LP = @. cos(theta)^2
    photon_UP = matter_LP
    matter_UP = photon_LP

    return (photon_LP=photon_LP, matter_LP=matter_LP,
            photon_UP=photon_UP, matter_UP=matter_UP)
end
