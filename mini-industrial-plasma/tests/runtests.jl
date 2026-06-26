using Test
include(joinpath(@__DIR__, "..", "src", "IndustrialPlasma.jl"))
using .IndustrialPlasma
@testset "IndustrialPlasma" begin @test true end
