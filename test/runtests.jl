using ClassicalSpinMonteCarlo
using Test

@testset "ClassicalSpinMonteCarlo.jl" begin
    include("test.jl")
    include("test_Lattice.jl")
    include("test_spin.jl")
    include("test_UnitCell.jl")
end
