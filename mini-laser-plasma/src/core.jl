# core.jl -- LaserPlasma
# MIT 22.611 / Kruer

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

function critical_density(lam);omega=2*pi*C_L/lam;EPS0*M_E*omega^2/E_C^2;end
function ponderomotive_potential(I,lam);omega=2*pi*C_L/lam;E0=sqrt(2*I/(EPS0*C_L));E_C^2*E0^2/(4*M_E*omega^2);end
function relativistic_intensity(lam);omega=2*pi*C_L/lam;EPS0*C_L*M_E^2*omega^2*C_L^2/E_C^2;end
wakefield_gradient(np)=sqrt(np/1e24)*100e9
function inverse_bremsstrahlung(ne,nc,nuei);ne^2*nuei/(nc*sqrt(1-ne/nc));end
