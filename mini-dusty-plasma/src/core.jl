# core.jl -- DustyPlasma
# MIT 22.611 / Shukla & Mamun

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

dust_charge(Zd)=Zd*E_C
function dust_plasma_freq(nd,md,Zd);sqrt(nd*Zd^2*E_C^2/(EPS0*md));end
function dust_acoustic_speed(Te,Zd,md);sqrt(Zd*K_B*Te/md);end
havnes_parameter(P,ne,Zd,nd)=Zd*nd/ne
function dust_coulomb_coupling(Qd,a,Td);Qd^2/(4*pi*EPS0*a*K_B*Td);end
dust_crystal_condition(Gamma)=Gamma>170
