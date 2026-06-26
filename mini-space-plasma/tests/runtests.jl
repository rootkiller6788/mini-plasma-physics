using Test
include(joinpath(@__DIR__, "..", "src", "SpacePlasma.jl"))
using .SpacePlasma
@testset "SpacePlasma" begin @test true end
