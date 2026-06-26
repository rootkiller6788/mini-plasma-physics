# 尘埃等离子体 Cheatsheet
> MIT 22.611 / Shukla & Mamun
```julia
dust_charge(Zd)=Zd*E_C
function dust_plasma_freq(nd,md,Zd);sqrt(nd*Zd^2*E_C^2/(EPS0*md));end
function dust_acoustic_speed(Te,Zd,md);sqrt(Zd*K_B*Te/md);end
havnes_parameter(P,ne,Zd,nd)=Zd*nd/ne
function dust_coulomb_coupling(Qd,a,Td);Qd^2/(4*pi*EPS0*a*K_B*Td);end
dust_crystal_condition(Gamma)=Gamma>170
```
