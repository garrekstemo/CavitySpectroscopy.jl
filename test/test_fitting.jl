@testset "Fitting" begin

    @testset "Cavity transmittance reflectance fit" begin
        nu = collect(0.0:0.005:2.0)
        p_true = [1.0, 0.0, 1.0, 0.9, 0.0]
        T_true = cavity_transmittance(p_true, nu)
        noise = 0.005 * randn(length(nu))
        T_data = clamp.(T_true .+ noise, 0.0, 1.0)

        # Fit reflectance with the raw model (keeping other params fixed)
        model(p, x) = cavity_transmittance([1.0, 0.0, 1.0, p[1], 0.0], x)
        p0 = [0.85]
        prob = CavitySpectroscopy.NonlinearCurveFitProblem(model, p0, nu, T_data)
        sol = CavitySpectroscopy.solve(prob)
        R_fit = CavitySpectroscopy.coef(sol)[1]

        @test isapprox(R_fit, p_true[4], atol=0.05)
    end

    @testset "Synthetic cavity spectrum round-trip" begin
        # Generate synthetic cavity spectrum with known parameters
        nu = collect(1900.0:0.5:2200.0)
        R_true = 0.92
        L_true = 12.0e-4
        n_bg_true = 1.4
        phi_true = 0.3
        A_true = 3000.0
        nu0_true = 2055.0
        Gamma_true = 23.0

        T_true = compute_cavity_transmittance(nu, [nu0_true], [Gamma_true], [A_true],
                                              R_true, L_true, n_bg_true, phi_true)

        # Add small noise
        T_noisy = clamp.(T_true .+ 0.002 .* randn(length(nu)), 0.0, 1.0)

        result = fit_cavity_spectrum(nu, T_noisy;
            oscillators=[(nu0=nu0_true, Gamma=Gamma_true)],
            L=L_true,
            n_bg=n_bg_true,
            R_init=0.9,
            phi_init=0.2,
            A_init=2500.0)

        @test result isa CavityFitResult
        @test result.rsquared > 0.95
        @test isapprox(result.R, R_true, atol=0.05)
        @test length(result.oscillators) == 1
        @test result.oscillators[1].nu0 ≈ nu0_true
        @test result.oscillators[1].Gamma ≈ Gamma_true

        # Strong coupling: fitted curve should expose two polariton peaks
        @test length(result.polariton_peaks) >= 2

        # predict should work
        y_fit = predict(result)
        @test length(y_fit) == length(nu)
        @test all(isfinite, y_fit)

        y_fit_custom = predict(result, nu[1:10])
        @test length(y_fit_custom) == 10

        # residuals
        res = residuals(result)
        @test length(res) == length(nu)
    end

    @testset "Region parameter" begin
        nu = collect(1900.0:0.5:2200.0)
        T = compute_cavity_transmittance(nu, [2055.0], [23.0], [3000.0],
                                         0.92, 12.0e-4, 1.4, 0.3)

        result = fit_cavity_spectrum(nu, T;
            oscillators=[(nu0=2055.0, Gamma=23.0)],
            L=12.0e-4, n_bg=1.4,
            region=(1950, 2150))

        @test result.rsquared > 0.99
        @test length(result._nu) < length(nu)
        @test all(1950 .<= result._nu .<= 2150)
    end

    @testset "fit_nu0 and fit_Gamma flags" begin
        nu = collect(1900.0:0.5:2200.0)
        T = compute_cavity_transmittance(nu, [2055.0], [23.0], [3000.0],
                                         0.92, 12.0e-4, 1.4, 0.3)

        # Start the oscillator slightly off and let the fit recover it
        result = fit_cavity_spectrum(nu, T;
            oscillators=[(nu0=2050.0, Gamma=20.0)],
            L=12.0e-4, n_bg=1.4,
            fit_nu0=true, fit_Gamma=true)

        @test result.rsquared > 0.99
        @test isapprox(result.oscillators[1].nu0, 2055.0, atol=1.0)
        @test isapprox(result.oscillators[1].Gamma, 23.0, atol=2.0)
    end

    @testset "fit_cavity_spectrum on CavitySpectrum" begin
        nu = collect(1900.0:0.5:2200.0)
        T = compute_cavity_transmittance(nu, [2055.0], [23.0], [3000.0],
                                         0.92, 12.0e-4, 1.4, 0.3)

        # Percent transmittance + cavity length from metadata
        spec = CavitySpectrum(nu, 100.0 .* T; cavity_length=12.0e-4)

        result = fit_cavity_spectrum(spec;
            oscillators=[(nu0=2055.0, Gamma=23.0)],
            n_bg=1.4)

        @test result isa CavityFitResult
        @test result.rsquared > 0.99
        # Transmittance was auto-normalized from percent to fractional
        @test maximum(result._T_data) <= 1.5
        # L came from metadata
        @test result.L == 12.0e-4

        # Explicit kwarg overrides metadata
        result2 = fit_cavity_spectrum(spec;
            oscillators=[(nu0=2055.0, Gamma=23.0)],
            L=12.0e-4, n_bg=1.4)
        @test result2.L == 12.0e-4
    end

    @testset "Dispersion round-trip" begin
        # Generate synthetic dispersion data
        E_vib = 2055.0
        Omega_true = 25.0
        E0_true = 2040.0
        n_eff_true = 1.5

        angles = collect(0.0:2.0:30.0) .* (pi / 180)

        E_cav = cavity_mode_energy([E0_true, n_eff_true], angles)
        lp_true, up_true = polariton_branches(E_cav, E_vib, Omega_true)

        # Add small noise
        lp_noisy = lp_true .+ 0.5 .* randn(length(angles))
        up_noisy = up_true .+ 0.5 .* randn(length(angles))

        result = fit_dispersion(angles, lp_noisy, up_noisy;
            molecular_modes=E_vib,
            E0_init=2035.0,
            n_eff_init=1.4,
            Omega_init=20.0)

        @test result isa DispersionFitResult
        @test result.rsquared > 0.95
        @test isapprox(result.rabi_splitting, Omega_true, atol=3.0)
        @test isapprox(result.E0, E0_true, atol=5.0)
        @test isapprox(result.n_eff, n_eff_true, atol=0.2)
        @test length(result.molecular_modes) == 1

        # Hopfield at zero detuning should be ~50/50
        @test result.hopfield_zero.photon_LP ≈ 0.5 atol=0.05

        # Stored data should match input
        @test length(result.lp_angles) == length(angles)
        @test length(result.up_angles) == length(angles)
        @test result.lp_positions ≈ lp_noisy atol=1e-10
        @test result.up_positions ≈ up_noisy atol=1e-10

        # Uncertainties should be finite and positive
        @test result.rabi_err > 0 && isfinite(result.rabi_err)
        @test result.E0_err > 0 && isfinite(result.E0_err)
        @test result.n_eff_err > 0 && isfinite(result.n_eff_err)
    end

    @testset "Dispersion textbook values (Omega=100)" begin
        # Textbook: cavity at 2000 cm^-1, Rabi splitting 100 cm^-1
        E_vib = 2000.0
        Omega_true = 100.0
        E0_true = 1950.0
        n_eff_true = 1.5

        angles = collect(0.0:2.0:40.0) .* (pi / 180)

        E_cav = cavity_mode_energy([E0_true, n_eff_true], angles)
        lp_true, up_true = polariton_branches(E_cav, E_vib, Omega_true)

        # Noiseless round-trip: should recover exact parameters
        result = fit_dispersion(angles, lp_true, up_true;
            molecular_modes=E_vib,
            E0_init=1940.0,
            n_eff_init=1.4,
            Omega_init=80.0)

        @test result isa DispersionFitResult
        @test result.rsquared > 0.999
        @test isapprox(result.rabi_splitting, Omega_true, atol=0.5)
        @test isapprox(result.E0, E0_true, atol=1.0)
        @test isapprox(result.n_eff, n_eff_true, atol=0.01)
    end

    @testset "Dispersion with different LP/UP angles" begin
        # In experiments, LP and UP may be measured at different angles
        E_vib = 2000.0
        Omega_true = 80.0
        E0_true = 1970.0
        n_eff_true = 1.5

        lp_angles = collect(0.0:3.0:25.0) .* (pi / 180)
        up_angles = collect(5.0:3.0:35.0) .* (pi / 180)

        E_cav_lp = cavity_mode_energy([E0_true, n_eff_true], lp_angles)
        E_cav_up = cavity_mode_energy([E0_true, n_eff_true], up_angles)
        lp_true, _ = polariton_branches(E_cav_lp, E_vib, Omega_true)
        _, up_true = polariton_branches(E_cav_up, E_vib, Omega_true)

        result = fit_dispersion(lp_angles, lp_true, up_angles, up_true;
            molecular_modes=E_vib,
            Omega_init=60.0)

        @test result isa DispersionFitResult
        @test result.rsquared > 0.999
        @test isapprox(result.rabi_splitting, Omega_true, atol=1.0)
        @test length(result.lp_angles) == length(lp_angles)
        @test length(result.up_angles) == length(up_angles)
    end

    @testset "Multi-mode dispersion fit" begin
        # Two molecular modes coupled to one cavity mode; LP/UP positions
        # generated from the (N+1)-level Hamiltonian eigenvalues.
        mol_modes = [2000.0, 2060.0]
        Omega_true = 60.0
        E0_true = 1960.0
        n_eff_true = 1.5

        angles = collect(0.0:2.0:35.0) .* (pi / 180)
        E_cav = cavity_mode_energy([E0_true, n_eff_true], angles)

        lp = [polariton_eigenvalues(e, mol_modes, [Omega_true, Omega_true])[1] for e in E_cav]
        up = [polariton_eigenvalues(e, mol_modes, [Omega_true, Omega_true])[end] for e in E_cav]

        result = fit_dispersion(angles, lp, up;
            molecular_modes=mol_modes,
            E0_init=1950.0,
            n_eff_init=1.4,
            Omega_init=50.0)

        @test result isa DispersionFitResult
        @test result.rsquared > 0.999
        @test isapprox(result.rabi_splitting, Omega_true, atol=1.0)
        @test isapprox(result.E0, E0_true, atol=2.0)
        @test length(result.molecular_modes) == 2
        # Uncertainties require ForwardDiff through the eigenvalue solver
        @test isfinite(result.rabi_err)
    end

    @testset "fit_dispersion from CavityFitResult vector" begin
        # Build synthetic CavityFitResults whose polariton_peaks follow a
        # known dispersion; the angle-resolved peak extraction should
        # recover the Rabi splitting.
        E_vib = 2055.0
        Omega_true = 25.0
        E0_true = 2040.0
        n_eff_true = 1.5

        angles = collect(0.0:5.0:30.0) .* (pi / 180)
        E_cav = cavity_mode_energy([E0_true, n_eff_true], angles)
        lp, up = polariton_branches(E_cav, E_vib, Omega_true)

        nu = collect(1900.0:1.0:2200.0)
        T = zeros(length(nu))
        results = [CavityFitResult(0.92, 12.0e-4, 1.4, 0.3, 1.0, 0.0,
                                   [(nu0=E_vib, Gamma=23.0, A=3000.0)],
                                   [lp[i], up[i]], 0.999, nu, T, nothing)
                   for i in eachindex(angles)]

        result = fit_dispersion(results; molecular_modes=E_vib, angles=angles)

        @test result isa DispersionFitResult
        @test isapprox(result.rabi_splitting, Omega_true, atol=1.0)
        @test isapprox(result.E0, E0_true, atol=2.0)

        # Mismatched lengths should error
        @test_throws AssertionError fit_dispersion(results;
            molecular_modes=E_vib, angles=angles[1:2])

        # Too few valid LP/UP pairs should error
        few = results[1:2]
        @test_throws ErrorException fit_dispersion(few;
            molecular_modes=E_vib, angles=angles[1:2])
    end

end
