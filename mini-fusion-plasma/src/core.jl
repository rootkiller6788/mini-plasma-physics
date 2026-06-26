# core.jl -- FusionPlasma
# MIT 22.611 / Wesson

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

lawson_criterion(n,T,tauE)=n*T*tauE
triple_product(n,T,tauE)=n*T*tauE
fusion_power_density(nD,nT,sigmav,Efus=17.6e6*E_C)=nD*nT*sigmav*Efus
function bremsstrahlung_power(ne,Zeff,Te);5.35e-37*ne^2*Zeff*sqrt(Te);end
function energy_confinement_IPB98(I,B,n,Ploss,R,a,kappa);0.0562*I^0.93*B^0.15*n^0.41*R^1.97*kappa^0.78;end
greenwald_density_limit(Ip,a)=Ip/(pi*a^2)
ignition_condition(Q)=Q
