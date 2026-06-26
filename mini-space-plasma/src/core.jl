# core.jl -- SpacePlasma
# MIT 22.611 / Kivelson & Russell

const K_B=1.380649e-23
const HBAR=1.054571817e-34
const C_L=299792458.0
const M_E=9.1093837e-31
const E_C=1.602176634e-19
const G=6.67430e-11
const R=8.314462618
const EPS0=8.8541878128e-12
const MU0=1.25663706212e-6
const N_A=6.02214076e23
const H=6.62607015e-34
const SIGMA_SB=5.670374419e-8

function parker_spiral_B(r,B0,r0,omega,vsw);Br=B0*(r0/r)^2;Bphi=-Br*omega*r/vsw;(Br,Bphi);end
alfven_mach_number(v,vA)=v/vA
sonic_mach_number(v,cs)=v/cs
function magnetosphere_convection(E,B);cross(E,B)/dot(B,B);end
