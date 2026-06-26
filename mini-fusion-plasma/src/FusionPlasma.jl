# FusionPlasma.jl -- 聚变等离子体
# 参考: MIT 22.611 / Wesson

module FusionPlasma
export lawson_criterion
export triple_product
export fusion_power_density
export bremsstrahlung_power
export energy_confinement_IPB98
export greenwald_density_limit
export ignition_condition
include("core.jl")
end
