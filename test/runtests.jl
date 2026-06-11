using Test
using CavitySpectroscopy
using Aqua
using Random

Random.seed!(20260611)

@testset "CavitySpectroscopy" begin

    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(CavitySpectroscopy; deps_compat=(check_extras=false,))
    end

    include("test_physics.jl")
    include("test_types.jl")

end
