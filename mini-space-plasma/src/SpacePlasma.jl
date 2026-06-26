# SpacePlasma.jl -- 空间等离子体
# 参考: MIT 22.611 / Kivelson & Russell

module SpacePlasma
export parker_spiral_B
export alfven_mach_number
export sonic_mach_number
export magnetosphere_convection
include("core.jl")
end
