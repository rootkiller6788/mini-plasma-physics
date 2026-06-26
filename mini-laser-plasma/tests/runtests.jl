using Test
include(joinpath(@__DIR__, "..", "src", "LaserPlasma.jl"))
using .LaserPlasma
@testset "LaserPlasma" begin @test true end
