# IndustrialPlasma.jl -- 工业等离子体
# 参考: MIT 22.611 / Lieberman

module IndustrialPlasma
export child_langmuir_current
export debye_sheath_thickness
export bohm_criterion
export floating_potential
export electron_energy_dist_maxwellian
include("core.jl")
end
