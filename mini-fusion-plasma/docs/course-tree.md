# Course Tree — mini-fusion-plasma

## Prerequisites
```
Classical Mechanics (L4: Newton, Lagrangian, Hamiltonian)
  └─ Electromagnetism (Maxwell equations, Lorentz force)
      └─ Statistical Mechanics (Maxwell-Boltzmann, kinetic theory)
          └─ Fluid Mechanics (Navier-Stokes, MHD)
              └─ Quantum Mechanics (atomic physics, cross-sections)
                  └─ mini-fusion-plasma (THIS MODULE)
```

## Internal Dependencies
```
fusion_plasma.h (types + L2 fundamentals)
  ├─ fusion_core.c (L2-L5 core physics)
  │   ├─ fusion_confinement.c (L4-L8 transport + stability)
  │   ├─ fusion_mhd.c (L3-L8 equilibrium + MHD)
  │   ├─ fusion_heating.c (L2-L7 heating + CD)
  │   └─ fusion_reactor.c (L7-L9 power plant design)
  └─ fusion_equilibrium.c (L6-L7 operating points)
```

## Downstream Modules
- mini-nuclear-physics (cross-sections, neutronics)
- mini-particle-physics (fusion as particle source)
- mini-astrophysics (stellar fusion, solar interior)
- mini-computational-physics (PIC, MHD codes)