using Test
include(joinpath(@__DIR__, "..", "src", "DustyPlasma.jl"))
using .DustyPlasma
@testset "DustyPlasma" begin @test true end
