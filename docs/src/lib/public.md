# Public API

```@docs
CavitySpectroscopy
```

## Types

```@docs
CavitySpectrum
CavityFitResult
DispersionFitResult
```

## Accessors

```@docs
wavenumber
transmittance
```

## Physics

```@docs
cavity_transmittance
compute_cavity_transmittance
refractive_index
extinction_coeff
cavity_mode_energy
polariton_branches
polariton_eigenvalues
hopfield_coefficients
```

## Fitting

```@docs
fit_cavity_spectrum
fit_dispersion
predict(::CavityFitResult)
residuals(::CavityFitResult)
format_results
```
