# Data and result types.

# =============================================================================
# CavitySpectrum
# =============================================================================

"""
    CavitySpectrum(x, y, metadata::Dict{String,Any})
    CavitySpectrum(x, y; metadata...)

Cavity transmission spectrum: plain data plus sample metadata.

# Fields
- `x::Vector{Float64}` - Wavenumber axis (cm⁻¹)
- `y::Vector{Float64}` - Transmittance
- `metadata::Dict{String,Any}` - Sample metadata (mirror, cavity_length, angle, ...)

Metadata keys used elsewhere in the package:
- `"cavity_length"` - picked up as the default `L` by [`fit_cavity_spectrum`](@ref)

# Examples
```julia
spec = CavitySpectrum(nu, T)
spec = CavitySpectrum(nu, T; mirror="Au", angle=10, cavity_length=12e-4)
spec.metadata["mirror"]
```

Loaders that wrap instrument-specific formats (e.g. JASCO FTIR files)
live in downstream lab packages; this type only carries data.
"""
struct CavitySpectrum
    x::Vector{Float64}
    y::Vector{Float64}
    metadata::Dict{String, Any}

    function CavitySpectrum(x::AbstractVector{<:Real}, y::AbstractVector{<:Real},
                            metadata::AbstractDict)
        length(x) == length(y) ||
            throw(ArgumentError("x and y must have the same length, got $(length(x)) and $(length(y))"))
        return new(Float64.(x), Float64.(y), Dict{String, Any}(metadata))
    end
end

function CavitySpectrum(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}; metadata...)
    md = Dict{String, Any}(String(k) => v for (k, v) in metadata)
    return CavitySpectrum(x, y, md)
end

"""
    wavenumber(s::CavitySpectrum) -> Vector{Float64}

Return the wavenumber axis (cm⁻¹).
"""
wavenumber(s::CavitySpectrum) = s.x

"""
    transmittance(s::CavitySpectrum) -> Vector{Float64}

Return the transmittance signal.
"""
transmittance(s::CavitySpectrum) = s.y

function Base.show(io::IO, spec::CavitySpectrum)
    label = get(spec.metadata, "_id", nothing)
    n = length(spec.x)
    if isnothing(label)
        print(io, "CavitySpectrum($n points)")
    else
        print(io, "CavitySpectrum(\"$label\", $n points)")
    end
end

function Base.show(io::IO, ::MIME"text/plain", spec::CavitySpectrum)
    println(io, "CavitySpectrum:")

    id = get(spec.metadata, "_id", nothing)
    !isnothing(id) && println(io, "  id: $id")

    for key in ["sample", "mirror", "cavity_length", "angle", "solute", "concentration", "solvent"]
        val = get(spec.metadata, key, nothing)
        !isnothing(val) && println(io, "  $key: $val")
    end

    x = spec.x
    if !isempty(x)
        println(io, "  range: $(round(minimum(x), digits=1)) - $(round(maximum(x), digits=1)) cm⁻¹")
    end
    print(io, "  points: $(length(x))")
end

# =============================================================================
# CavityFitResult
# =============================================================================

"""
    CavityFitResult

Result from fitting a single cavity transmission spectrum.

# Fields
- `R`: Mirror reflectivity
- `L`: Cavity length (cm)
- `n_bg`: Background refractive index
- `phi`: Phase shift
- `scale`: Scale factor applied to transmittance
- `offset`: Baseline offset
- `oscillators`: Vector of NamedTuples `(nu0, Gamma, A)` for each oscillator
- `polariton_peaks`: Vector of peak positions (cm^-1) extracted from fit
- `rsquared`: R^2 goodness of fit
- `_nu`: Wavenumber array used in fit (internal)
- `_T_data`: Transmittance data used in fit (internal)
- `_sol`: CurveFit solution object (internal)
"""
struct CavityFitResult
    R::Float64
    L::Float64
    n_bg::Float64
    phi::Float64
    scale::Float64
    offset::Float64
    oscillators::Vector{NamedTuple{(:nu0, :Gamma, :A), Tuple{Float64, Float64, Float64}}}
    polariton_peaks::Vector{Float64}
    rsquared::Float64
    _nu::Vector{Float64}
    _T_data::Vector{Float64}
    _sol::Any
end

"""
    wavenumber(r::CavityFitResult) -> Vector{Float64}

Return the wavenumber array used in the fit.
"""
wavenumber(r::CavityFitResult) = r._nu

"""
    transmittance(r::CavityFitResult) -> Vector{Float64}

Return the transmittance data used in the fit.
"""
transmittance(r::CavityFitResult) = r._T_data

"""
    predict(result::CavityFitResult)
    predict(result::CavityFitResult, nu)

Return fitted transmittance on the original wavenumber grid, or on a
custom wavenumber array `nu`. Extends `CurveFit.predict`.
"""
predict(result::CavityFitResult) = predict(result, result._nu)

function predict(result::CavityFitResult, nu::AbstractVector)
    nu0s = [osc.nu0 for osc in result.oscillators]
    Gammas = [osc.Gamma for osc in result.oscillators]
    As = [osc.A for osc in result.oscillators]
    T = compute_cavity_transmittance(nu, nu0s, Gammas, As,
                                     result.R, result.L, result.n_bg, result.phi)
    return T .* result.scale .+ result.offset
end

"""
    residuals(result::CavityFitResult)

Return residuals (data - fit) on the original wavenumber grid.
Extends `CurveFit.residuals`.
"""
residuals(result::CavityFitResult) = result._T_data .- predict(result)

function Base.show(io::IO, r::CavityFitResult)
    n_osc = length(r.oscillators)
    n_peaks = length(r.polariton_peaks)
    print(io, "CavityFitResult($n_osc oscillator$(n_osc == 1 ? "" : "s"), ",
          "$n_peaks polariton peak$(n_peaks == 1 ? "" : "s"), ",
          "R^2=$(round(r.rsquared, digits=4)))")
end

function Base.show(io::IO, ::MIME"text/plain", r::CavityFitResult)
    println(io, "Cavity Spectrum Fit")
    println(io, "=" ^ 50)

    println(io, "\nCavity parameters:")
    println(io, "  R       = $(round(r.R, digits=4))")
    println(io, "  L       = $(r.L) cm")
    println(io, "  n_bg    = $(round(r.n_bg, digits=3))")
    println(io, "  phi     = $(round(r.phi, digits=4))")
    println(io, "  scale   = $(round(r.scale, digits=4))")
    println(io, "  offset  = $(round(r.offset, digits=4))")

    if !isempty(r.oscillators)
        println(io, "\nOscillators:")
        for (i, osc) in enumerate(r.oscillators)
            println(io, "  [$i] nu0 = $(round(osc.nu0, digits=1)) cm^-1, Gamma = $(round(osc.Gamma, digits=1)) cm^-1, A = $(round(osc.A, digits=1))")
        end
    end

    if !isempty(r.polariton_peaks)
        println(io, "\nPolariton peaks:")
        for (i, pk) in enumerate(r.polariton_peaks)
            println(io, "  [$i] $(round(pk, digits=1)) cm^-1")
        end
    end

    print(io, "\nR^2 = $(round(r.rsquared, digits=6))")
end

"""
    format_results(result) -> String

Format a fit result ([`CavityFitResult`](@ref) or [`DispersionFitResult`](@ref))
as a markdown string, suitable for lab notebooks and reports.
"""
function format_results end

function format_results(r::CavityFitResult)
    lines = String[]
    push!(lines, "## Cavity Spectrum Fit\n")

    push!(lines, "| Parameter | Value |")
    push!(lines, "|-----------|-------|")
    push!(lines, "| R | $(round(r.R, digits=4)) |")
    push!(lines, "| L | $(r.L) cm |")
    push!(lines, "| n_bg | $(round(r.n_bg, digits=3)) |")
    push!(lines, "| phi | $(round(r.phi, digits=4)) |")
    push!(lines, "| scale | $(round(r.scale, digits=4)) |")
    push!(lines, "| offset | $(round(r.offset, digits=4)) |")
    push!(lines, "| R^2 | $(round(r.rsquared, digits=6)) |")

    if !isempty(r.oscillators)
        push!(lines, "\n### Oscillators\n")
        push!(lines, "| # | nu0 (cm^-1) | Gamma (cm^-1) | A |")
        push!(lines, "|---|-------------|---------------|---|")
        for (i, osc) in enumerate(r.oscillators)
            push!(lines, "| $i | $(round(osc.nu0, digits=1)) | $(round(osc.Gamma, digits=1)) | $(round(osc.A, digits=1)) |")
        end
    end

    if !isempty(r.polariton_peaks)
        push!(lines, "\n### Polariton Peaks\n")
        push!(lines, "| # | Position (cm^-1) |")
        push!(lines, "|---|-----------------|")
        for (i, pk) in enumerate(r.polariton_peaks)
            push!(lines, "| $i | $(round(pk, digits=1)) |")
        end
    end

    return join(lines, "\n")
end

# =============================================================================
# DispersionFitResult
# =============================================================================

"""
    DispersionFitResult

Result from fitting the coupled oscillator model to polariton dispersion data.

# Fields
- `rabi_splitting`: Rabi splitting Omega (cm^-1)
- `molecular_modes`: Vector of molecular mode energies (cm^-1)
- `n_eff`: Effective refractive index
- `E0`: Normal-incidence cavity energy (cm^-1)
- `rabi_err`: Uncertainty in Rabi splitting
- `n_eff_err`: Uncertainty in n_eff
- `E0_err`: Uncertainty in E0
- `lp_angles`: Incidence angles for LP data (radians)
- `lp_positions`: Lower polariton positions at each LP angle
- `up_angles`: Incidence angles for UP data (radians)
- `up_positions`: Upper polariton positions at each UP angle
- `hopfield_zero`: Hopfield coefficients at zero detuning
- `rsquared`: R^2 goodness of fit
- `_sol`: CurveFit solution object (internal)
"""
struct DispersionFitResult
    rabi_splitting::Float64
    molecular_modes::Vector{Float64}
    n_eff::Float64
    E0::Float64
    rabi_err::Float64
    n_eff_err::Float64
    E0_err::Float64
    lp_angles::Vector{Float64}
    lp_positions::Vector{Float64}
    up_angles::Vector{Float64}
    up_positions::Vector{Float64}
    hopfield_zero::NamedTuple{(:photon_LP, :matter_LP, :photon_UP, :matter_UP),
                              NTuple{4, Float64}}
    rsquared::Float64
    _sol::Any
end

function Base.show(io::IO, r::DispersionFitResult)
    print(io, "DispersionFitResult(Omega=$(round(r.rabi_splitting, digits=1)) cm^-1, ",
          "R^2=$(round(r.rsquared, digits=4)))")
end

function Base.show(io::IO, ::MIME"text/plain", r::DispersionFitResult)
    println(io, "Dispersion Fit (Coupled Oscillator Model)")
    println(io, "=" ^ 50)

    println(io, "\nFitted parameters:")
    println(io, "  Rabi splitting = $(round(r.rabi_splitting, digits=1)) +/- $(round(r.rabi_err, digits=1)) cm^-1")
    println(io, "  E0 (normal)    = $(round(r.E0, digits=1)) +/- $(round(r.E0_err, digits=1)) cm^-1")
    println(io, "  n_eff          = $(round(r.n_eff, digits=3)) +/- $(round(r.n_eff_err, digits=3))")

    println(io, "\nMolecular modes:")
    for (i, m) in enumerate(r.molecular_modes)
        println(io, "  [$i] $(round(m, digits=1)) cm^-1")
    end

    h = r.hopfield_zero
    println(io, "\nHopfield coefficients (zero detuning):")
    println(io, "  LP: photon = $(round(h.photon_LP, digits=3)), matter = $(round(h.matter_LP, digits=3))")
    println(io, "  UP: photon = $(round(h.photon_UP, digits=3)), matter = $(round(h.matter_UP, digits=3))")

    println(io, "\nR^2 = $(round(r.rsquared, digits=6))")
    print(io, "Data points: $(length(r.lp_angles)) LP, $(length(r.up_angles)) UP")
end

function format_results(r::DispersionFitResult)
    lines = String[]
    push!(lines, "## Dispersion Fit (Coupled Oscillator)\n")

    push!(lines, "| Parameter | Value | Uncertainty |")
    push!(lines, "|-----------|-------|-------------|")
    push!(lines, "| Rabi splitting | $(round(r.rabi_splitting, digits=1)) cm^-1 | $(round(r.rabi_err, digits=1)) |")
    push!(lines, "| E0 | $(round(r.E0, digits=1)) cm^-1 | $(round(r.E0_err, digits=1)) |")
    push!(lines, "| n_eff | $(round(r.n_eff, digits=3)) | $(round(r.n_eff_err, digits=3)) |")
    push!(lines, "| R^2 | $(round(r.rsquared, digits=6)) | |")

    push!(lines, "\n### Molecular Modes\n")
    for (i, m) in enumerate(r.molecular_modes)
        push!(lines, "- Mode $i: $(round(m, digits=1)) cm^-1")
    end

    h = r.hopfield_zero
    push!(lines, "\n### Hopfield Coefficients (zero detuning)\n")
    push!(lines, "| Branch | Photon | Matter |")
    push!(lines, "|--------|--------|--------|")
    push!(lines, "| LP | $(round(h.photon_LP, digits=3)) | $(round(h.matter_LP, digits=3)) |")
    push!(lines, "| UP | $(round(h.photon_UP, digits=3)) | $(round(h.matter_UP, digits=3)) |")

    return join(lines, "\n")
end
