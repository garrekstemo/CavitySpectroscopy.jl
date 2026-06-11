# Fitting: cavity transmission spectra and polariton dispersion.

"""
    _find_local_maxima(x, y; min_prominence=0.0)

Find local maxima in y(x). Returns vector of x positions sorted by
prominence (highest first).

Prominence is topographic: the height of the peak above the higher of the
two valleys separating it from higher terrain (or the data boundary), so
it is independent of the sampling density.
"""
function _find_local_maxima(x, y; min_prominence::Real=0.0)
    peaks = Float64[]
    prominences = Float64[]
    for i in 2:(length(y) - 1)
        if y[i] > y[i-1] && y[i] > y[i+1]
            # Walk left until terrain higher than the peak (or boundary),
            # tracking the lowest valley; same to the right.
            left_min = y[i]
            j = i - 1
            while j >= 1 && y[j] <= y[i]
                left_min = min(left_min, y[j])
                j -= 1
            end
            right_min = y[i]
            j = i + 1
            while j <= length(y) && y[j] <= y[i]
                right_min = min(right_min, y[j])
                j += 1
            end
            prom = y[i] - max(left_min, right_min)
            if prom > min_prominence
                push!(peaks, x[i])
                push!(prominences, prom)
            end
        end
    end
    # Sort by prominence (highest first)
    order = sortperm(prominences, rev=true)
    return peaks[order]
end

"""
    fit_cavity_spectrum(nu, T_data; oscillators, L, n_bg,
                        R_init=0.92, phi_init=0.3, A_init=3000.0,
                        scale_init=1.0, offset_init=0.0,
                        region=nothing, fit_nu0=false, fit_Gamma=false)

Fit a cavity transmission spectrum with a multi-oscillator Fabry-Perot model.

# Arguments
- `nu`: Wavenumber array (cm^-1)
- `T_data`: Transmittance data (fractional, 0-1)
- `oscillators`: Vector of named tuples `(nu0=..., Gamma=...)` defining oscillator
  center frequencies and linewidths. These are fixed by default.
- `L`: Cavity length (cm)
- `n_bg`: Background refractive index
- `R_init`: Initial guess for mirror reflectivity (default: 0.92)
- `phi_init`: Initial guess for phase shift (default: 0.3)
- `A_init`: Initial guess for oscillator amplitude (default: 3000.0)
- `scale_init`: Initial guess for scale factor (default: 1.0)
- `offset_init`: Initial guess for baseline offset (default: 0.0)
- `region`: Optional `(lo, hi)` tuple to restrict fitting range
- `fit_nu0`: If true, also fit oscillator center frequencies (default: false)
- `fit_Gamma`: If true, also fit oscillator linewidths (default: false)

# Returns
[`CavityFitResult`](@ref) with fitted parameters and auto-extracted polariton peaks.
"""
function fit_cavity_spectrum(nu::AbstractVector, T_data::AbstractVector;
    oscillators,
    L::Real,
    n_bg::Real,
    R_init::Real=0.92,
    phi_init::Real=0.3,
    A_init::Real=3000.0,
    scale_init::Real=1.0,
    offset_init::Real=0.0,
    region=nothing,
    fit_nu0::Bool=false,
    fit_Gamma::Bool=false)

    # Apply region mask
    if !isnothing(region)
        mask = region[1] .<= nu .<= region[2]
        nu = Float64.(nu[mask])
        T_data = Float64.(T_data[mask])
    else
        nu = Float64.(nu)
        T_data = Float64.(T_data)
    end

    n_osc = length(oscillators)

    # Build parameter vector: [R, phi, scale, offset, A1, A2, ..., (nu0_1, ...), (Gamma_1, ...)]
    p0 = Float64[R_init, phi_init, scale_init, offset_init]
    for _ in 1:n_osc
        push!(p0, A_init)
    end
    if fit_nu0
        for osc in oscillators
            push!(p0, Float64(osc.nu0))
        end
    end
    if fit_Gamma
        for osc in oscillators
            push!(p0, Float64(osc.Gamma))
        end
    end

    # Fixed values for nu0 and Gamma when not fitting
    fixed_nu0s = Float64[osc.nu0 for osc in oscillators]
    fixed_Gammas = Float64[osc.Gamma for osc in oscillators]

    function model(p, x)
        R_val = p[1]
        phi_val = p[2]
        scale_val = p[3]
        offset_val = p[4]

        As = p[5:5 + n_osc - 1]

        idx = 5 + n_osc
        if fit_nu0
            nu0s = p[idx:idx + n_osc - 1]
            idx += n_osc
        else
            nu0s = fixed_nu0s
        end

        if fit_Gamma
            Gammas = p[idx:idx + n_osc - 1]
        else
            Gammas = fixed_Gammas
        end

        T = compute_cavity_transmittance(x, nu0s, Gammas, As, R_val, L, n_bg, phi_val)
        return T .* scale_val .+ offset_val
    end

    prob = NonlinearCurveFitProblem(model, p0, nu, T_data)
    sol = solve(prob)
    c = coef(sol)

    # Extract fitted parameters
    R_fit = c[1]
    phi_fit = c[2]
    scale_fit = c[3]
    offset_fit = c[4]
    A_fits = c[5:5 + n_osc - 1]

    idx = 5 + n_osc
    if fit_nu0
        nu0_fits = c[idx:idx + n_osc - 1]
        idx += n_osc
    else
        nu0_fits = fixed_nu0s
    end

    if fit_Gamma
        Gamma_fits = c[idx:idx + n_osc - 1]
    else
        Gamma_fits = fixed_Gammas
    end

    # Build oscillator results
    osc_results = [(nu0=nu0_fits[i], Gamma=Gamma_fits[i], A=A_fits[i]) for i in 1:n_osc]

    # Compute R^2
    y_fit = model(c, nu)
    ss_res = sum((T_data .- y_fit).^2)
    ss_tot = sum((T_data .- mean(T_data)).^2)
    rsq = 1 - ss_res / ss_tot

    # Auto-extract polariton peaks from fitted curve
    polariton_peaks = _find_local_maxima(nu, y_fit; min_prominence=0.005 * maximum(y_fit))

    return CavityFitResult(R_fit, L, n_bg, phi_fit, scale_fit, offset_fit,
                           osc_results, polariton_peaks, rsq, nu, T_data, sol)
end

"""
    fit_cavity_spectrum(spec::CavitySpectrum; kwargs...)

Fit a [`CavitySpectrum`](@ref). Uses the spectrum's wavenumber/transmittance
data; the cavity length `L` defaults to `spec.metadata["cavity_length"]`
when present and not given as a keyword.

Transmittance is normalized from percent (0-100) to fractional (0-1)
automatically if the maximum value exceeds 1.5.
"""
function fit_cavity_spectrum(spec::CavitySpectrum; kwargs...)
    nu = spec.x
    T = spec.y

    # Auto-normalize percent transmittance to fractional
    if maximum(T) > 1.5
        T = T ./ 100.0
    end

    # Pull defaults from metadata if not provided in kwargs
    kw = Dict{Symbol, Any}(kwargs)
    if !haskey(kw, :L) && haskey(spec.metadata, "cavity_length")
        kw[:L] = spec.metadata["cavity_length"]
    end

    return fit_cavity_spectrum(nu, T; kw...)
end

"""
    fit_dispersion(lp_angles, lp_positions, up_angles, up_positions;
                   molecular_modes, E0_init=nothing, n_eff_init=1.5, Omega_init=20.0)

Fit the coupled oscillator model to polariton dispersion data.

LP and UP data can be measured at different angles (common in experiments where
only the photon-like branch is visible at large detuning).

For a single molecular mode, fits the analytic 2-level model:
E_LP, E_UP = (E_cav + E_vib)/2 +/- sqrt(Omega^2 + (E_cav - E_vib)^2)/2

where E_cav(theta) = E0 / sqrt(1 - (sin(theta)/n_eff)^2). For multiple
molecular modes, LP/UP are the lowest/highest eigenvalues of the
(N+1)-level Hamiltonian ([`polariton_eigenvalues`](@ref)).

# Arguments
- `lp_angles`: Incidence angles for LP data (radians)
- `lp_positions`: Lower polariton energies (cm^-1)
- `up_angles`: Incidence angles for UP data (radians)
- `up_positions`: Upper polariton energies (cm^-1)
- `molecular_modes`: Scalar or vector of molecular mode energies (cm^-1), fixed
- `E0_init`: Initial guess for normal-incidence cavity energy (default: min of LP - 10)
- `n_eff_init`: Initial guess for effective refractive index (default: 1.5)
- `Omega_init`: Initial guess for Rabi splitting (default: 20.0)

# Returns
[`DispersionFitResult`](@ref)
"""
function fit_dispersion(lp_angles::AbstractVector, lp_positions::AbstractVector,
                        up_angles::AbstractVector, up_positions::AbstractVector;
                        molecular_modes,
                        E0_init=nothing,
                        n_eff_init::Real=1.5,
                        Omega_init::Real=20.0)

    mol_modes = molecular_modes isa Number ? [Float64(molecular_modes)] : Float64.(molecular_modes)

    # Default E0 guess: slightly below the lowest LP position
    if isnothing(E0_init)
        E0_init = minimum(lp_positions) - 10.0
    end

    n_lp = length(lp_angles)

    # Stack LP and UP data: fit both branches simultaneously
    y_data = Float64.([lp_positions; up_positions])

    # p = [E0, n_eff, Omega]
    p0 = Float64[E0_init, n_eff_init, Omega_init]

    # x-axis: LP angles then UP angles
    x = Float64.([lp_angles; up_angles])

    function model(p, x)
        E0, n_eff, Omega = p[1], p[2], p[3]
        a_lp = x[1:n_lp]
        a_up = x[n_lp+1:end]

        E_cav_lp = cavity_mode_energy((E0, n_eff), a_lp)
        E_cav_up = cavity_mode_energy((E0, n_eff), a_up)

        if length(mol_modes) == 1
            lp, _ = polariton_branches(E_cav_lp, mol_modes[1], Omega)
            _, up = polariton_branches(E_cav_up, mol_modes[1], Omega)
        else
            # Comprehensions (not similar/setindex!) so ForwardDiff duals
            # propagate through the eigenvalue solver.
            Omegas = fill(Omega, length(mol_modes))
            lp = [polariton_eigenvalues(e, mol_modes, Omegas)[1] for e in E_cav_lp]
            up = [polariton_eigenvalues(e, mol_modes, Omegas)[end] for e in E_cav_up]
        end

        return [lp; up]
    end

    prob = NonlinearCurveFitProblem(model, p0, x, y_data)
    sol = solve(prob)
    c = coef(sol)
    errs = stderror(sol)

    E0_fit, n_eff_fit, Omega_fit = c[1], c[2], c[3]
    E0_err, n_eff_err, Omega_err = errs[1], errs[2], errs[3]

    # Compute R^2
    y_pred = model(c, x)
    ss_res = sum((y_data .- y_pred).^2)
    ss_tot = sum((y_data .- mean(y_data)).^2)
    rsq = 1 - ss_res / ss_tot

    # Hopfield coefficients at zero detuning (E_cav = E_vib).
    # For multi-mode, computed at E_cav = mean of the molecular modes.
    E_zero = mean(mol_modes)
    h = hopfield_coefficients(E_zero, E_zero, Omega_fit)
    hopfield_zero = (photon_LP=h.photon_LP, matter_LP=h.matter_LP,
                     photon_UP=h.photon_UP, matter_UP=h.matter_UP)

    return DispersionFitResult(Omega_fit, mol_modes, n_eff_fit, E0_fit,
                               Omega_err, n_eff_err, E0_err,
                               Float64.(lp_angles), Float64.(lp_positions),
                               Float64.(up_angles), Float64.(up_positions),
                               hopfield_zero, rsq, sol)
end

"""
    fit_dispersion(angles, lp_positions, up_positions; kwargs...)

Convenience method when LP and UP are measured at the same angles.
"""
function fit_dispersion(angles::AbstractVector, lp_positions::AbstractVector,
                        up_positions::AbstractVector; kwargs...)
    return fit_dispersion(angles, lp_positions, angles, up_positions; kwargs...)
end

"""
    fit_dispersion(results::Vector{CavityFitResult}; molecular_modes, angles)

Extract LP/UP peak positions from a vector of [`CavityFitResult`](@ref)s and fit
the coupled oscillator model.

# Arguments
- `results`: Vector of cavity fit results (one per detuning/angle)
- `molecular_modes`: Molecular mode energy or vector of energies (cm^-1)
- `angles`: Vector of incidence angles (radians). Must match length of `results`.
"""
function fit_dispersion(results::Vector{CavityFitResult};
                        molecular_modes,
                        angles::AbstractVector)
    @assert length(angles) == length(results) "Need one angle per CavityFitResult"

    mol = molecular_modes isa Number ? Float64(molecular_modes) : Float64.(molecular_modes)
    mol_center = mol isa Number ? mol : mean(mol)

    lp = Float64[]
    up = Float64[]
    valid_angles = Float64[]

    for (i, r) in enumerate(results)
        if length(r.polariton_peaks) >= 2
            sorted_peaks = sort(r.polariton_peaks)
            # LP = peak below molecular mode, UP = peak above
            below = filter(p -> p < mol_center, sorted_peaks)
            above = filter(p -> p >= mol_center, sorted_peaks)

            if !isempty(below) && !isempty(above)
                push!(lp, last(below))    # Highest peak below molecular mode
                push!(up, first(above))   # Lowest peak above molecular mode
                push!(valid_angles, angles[i])
            end
        end
    end

    if length(lp) < 3
        error("Need at least 3 valid LP/UP pairs for dispersion fitting, got $(length(lp))")
    end

    return fit_dispersion(valid_angles, lp, up; molecular_modes=molecular_modes)
end
