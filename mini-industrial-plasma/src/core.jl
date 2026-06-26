# core.jl -- IndustrialPlasma
# MIT 22.611 / Lieberman

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

function child_langmuir_current(V,d,A,mi=1.67e-27);4*EPS0/9*sqrt(2*E_C/mi)*A*V^1.5/d^2;end
function debye_sheath_thickness(Vbias,Te,lamD);sqrt(2/3)*lamD*(2*Vbias/Te)^0.75;end
bohm_criterion(vi,cs)=vi>=cs
function floating_potential(Te,mi);-K_B*Te/(2*E_C)*log(mi/(2*pi*M_E));end
function electron_energy_dist_maxwellian(E,Te);2*sqrt(E/pi)*(K_B*Te)^(-1.5)*exp(-E/(K_B*Te));end
