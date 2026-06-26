# 空间等离子体 Cheatsheet
> MIT 22.611 / Kivelson & Russell
```julia
function parker_spiral_B(r,B0,r0,omega,vsw);Br=B0*(r0/r)^2;Bphi=-Br*omega*r/vsw;(Br,Bphi);end
alfven_mach_number(v,vA)=v/vA
sonic_mach_number(v,cs)=v/cs
function magnetosphere_convection(E,B);cross(E,B)/dot(B,B);end
```
