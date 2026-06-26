# 聚变等离子体 Cheatsheet
> MIT 22.611 / Wesson
```julia
lawson_criterion(n,T,tauE)=n*T*tauE
triple_product(n,T,tauE)=n*T*tauE
fusion_power_density(nD,nT,sigmav,Efus=17.6e6*E_C)=nD*nT*sigmav*Efus
function bremsstrahlung_power(ne,Zeff,Te);5.35e-37*ne^2*Zeff*sqrt(Te);end
function energy_confinement_IPB98(I,B,n,Ploss,R,a,kappa);0.0562*I^0.93*B^0.15*n^0.41*R^1.97*kappa^0.78;end
greenwald_density_limit(Ip,a)=Ip/(pi*a^2)
ignition_condition(Q)=Q
```
