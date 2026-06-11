@testset "Types" begin

    @testset "CavitySpectrum construction" begin
        x = collect(1900.0:1.0:2200.0)
        y = rand(length(x))

        # Plain positional construction
        spec = CavitySpectrum(x, y)
        @test spec.x == x
        @test spec.y == y
        @test isempty(spec.metadata)

        # Keyword metadata convenience
        spec2 = CavitySpectrum(x, y; mirror="Au", angle=10, cavity_length=12e-4)
        @test spec2.metadata["mirror"] == "Au"
        @test spec2.metadata["angle"] == 10
        @test spec2.metadata["cavity_length"] == 12e-4

        # Explicit metadata dict
        md = Dict{String, Any}("sample" => "NH4SCN/DMF")
        spec3 = CavitySpectrum(x, y, md)
        @test spec3.metadata["sample"] == "NH4SCN/DMF"

        # Generic numeric input converts to Float64
        spec4 = CavitySpectrum([1, 2, 3], [4, 5, 6])
        @test spec4.x isa Vector{Float64}
        @test spec4.y isa Vector{Float64}

        # Length mismatch errors
        @test_throws ArgumentError CavitySpectrum([1.0, 2.0], [1.0])
    end

    @testset "CavitySpectrum accessors" begin
        x = collect(1900.0:1.0:2200.0)
        y = collect(range(0.0, 1.0, length=length(x)))
        spec = CavitySpectrum(x, y)

        @test wavenumber(spec) == x
        @test transmittance(spec) == y
    end

    @testset "CavitySpectrum show" begin
        x = collect(1900.0:1.0:2200.0)
        y = rand(length(x))
        spec = CavitySpectrum(x, y; sample="1M NH4SCN", mirror="Au", angle=10)

        @test occursin("CavitySpectrum", repr(spec))
        @test occursin("301", repr(spec))  # point count

        buf = IOBuffer()
        show(buf, MIME("text/plain"), spec)
        output = String(take!(buf))
        @test occursin("CavitySpectrum", output)
        @test occursin("1M NH4SCN", output)
        @test occursin("Au", output)
        @test occursin("1900", output)
        @test occursin("2200", output)
    end

    @testset "CavityFitResult predict and residuals" begin
        nu = collect(1900.0:1.0:2200.0)
        oscillators = [(nu0=2055.0, Gamma=23.0, A=3000.0)]
        R, L, n_bg, phi = 0.92, 12.0e-4, 1.4, 0.3
        T = compute_cavity_transmittance(nu, [2055.0], [23.0], [3000.0],
                                         R, L, n_bg, phi)

        result = CavityFitResult(R, L, n_bg, phi, 1.0, 0.0,
                                 oscillators, Float64[], 1.0, nu, T, nothing)

        # predict on the original grid reproduces the model exactly
        # (scale = 1, offset = 0)
        @test predict(result) ≈ T atol=1e-12

        # predict on a custom grid
        y_custom = predict(result, nu[1:10])
        @test length(y_custom) == 10
        @test y_custom ≈ T[1:10] atol=1e-12

        # residuals are zero for self-consistent data
        @test residuals(result) ≈ zeros(length(nu)) atol=1e-12

        # scale and offset are applied
        result2 = CavityFitResult(R, L, n_bg, phi, 2.0, 0.1,
                                  oscillators, Float64[], 1.0, nu, T, nothing)
        @test predict(result2) ≈ 2.0 .* T .+ 0.1 atol=1e-12

        # accessors
        @test wavenumber(result) == nu
        @test transmittance(result) == T
    end

    @testset "CavityFitResult show and format_results" begin
        nu = collect(1900.0:1.0:2200.0)
        oscillators = [(nu0=2055.0, Gamma=23.0, A=3000.0)]
        T = compute_cavity_transmittance(nu, [2055.0], [23.0], [3000.0],
                                         0.92, 12.0e-4, 1.4, 0.3)
        result = CavityFitResult(0.92, 12.0e-4, 1.4, 0.3, 1.0, 0.0,
                                 oscillators, [2030.0, 2080.0], 0.999, nu, T, nothing)

        @test occursin("CavityFitResult", repr(result))

        buf = IOBuffer()
        show(buf, MIME("text/plain"), result)
        output = String(take!(buf))
        @test occursin("Cavity Spectrum Fit", output)
        @test occursin("R^2", output)
        @test occursin("2055", output)

        md = format_results(result)
        @test md isa String
        @test occursin("## Cavity Spectrum Fit", md)
        @test occursin("| R |", md)
        @test occursin("Polariton Peaks", md)
    end

    @testset "DispersionFitResult show and format_results" begin
        h = hopfield_coefficients(2055.0, 2055.0, 25.0)
        result = DispersionFitResult(25.0, [2055.0], 1.5, 2040.0,
                                     0.5, 0.01, 1.0,
                                     [0.0, 0.1], [2030.0, 2032.0],
                                     [0.0, 0.1], [2080.0, 2085.0],
                                     h, 0.998, nothing)

        @test occursin("DispersionFitResult", repr(result))
        @test occursin("25", repr(result))

        buf = IOBuffer()
        show(buf, MIME("text/plain"), result)
        output = String(take!(buf))
        @test occursin("Dispersion Fit", output)
        @test occursin("Rabi splitting", output)
        @test occursin("Hopfield", output)

        md = format_results(result)
        @test md isa String
        @test occursin("## Dispersion Fit", md)
        @test occursin("Rabi splitting", md)
    end

end
