# 工业等离子体 Cheatsheet
> MIT 22.611 / Lieberman
```julia
function child_langmuir_current(V,d,A,mi=1.67e-27);4*EPS0/9*sqrt(2*E_C/mi)*A*V^1.5/d^2;end
function debye_sheath_thickness(Vbias,Te,lamD);sqrt(2/3)*lamD*(2*Vbias/Te)^0.75;end
bohm_criterion(vi,cs)=vi>=cs
function floating_potential(Te,mi);-K_B*Te/(2*E_C)*log(mi/(2*pi*M_E));end
function electron_energy_dist_maxwellian(E,Te);2*sqrt(E/pi)*(K_B*Te)^(-1.5)*exp(-E/(K_B*Te));end
```
