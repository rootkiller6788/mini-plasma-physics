# 激光等离子体 Cheatsheet
> MIT 22.611 / Kruer
```julia
function critical_density(lam);omega=2*pi*C_L/lam;EPS0*M_E*omega^2/E_C^2;end
function ponderomotive_potential(I,lam);omega=2*pi*C_L/lam;E0=sqrt(2*I/(EPS0*C_L));E_C^2*E0^2/(4*M_E*omega^2);end
function relativistic_intensity(lam);omega=2*pi*C_L/lam;EPS0*C_L*M_E^2*omega^2*C_L^2/E_C^2;end
wakefield_gradient(np)=sqrt(np/1e24)*100e9
function inverse_bremsstrahlung(ne,nc,nuei);ne^2*nuei/(nc*sqrt(1-ne/nc));end
```
