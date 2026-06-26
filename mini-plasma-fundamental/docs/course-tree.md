# Course Tree — mini-plasma-fundamental

## Prerequisites
- Classical mechanics (Newton's laws, phase space)
- Electromagnetism (Maxwell equations, Lorentz force)
- Statistical mechanics (Maxwell-Boltzmann distribution)
- Fluid dynamics (Navier-Stokes equations)
- Numerical methods (finite differences, ODE integration)

## Dependency Graph

```
L1: Definitions
├── L2: Core Concepts
│   ├── L3: Math Structures
│   │   ├── L4: Fundamental Laws
│   │   │   ├── L5: Computational Methods
│   │   │   │   ├── L6: Canonical Systems
│   │   │   │   │   ├── L7: Applications
│   │   │   │   │   │   └── L8: Advanced Topics
│   │   │   │   │   │       └── L9: Research Frontiers
```

## Module Dependencies
- mini-classical-mechanics → Newton, phase space
- mini-electromagnetism → Maxwell, Lorentz
- mini-thermodynamics-statphys → distributions, equilibrium

## Internal Module Flow
1. Start with `plasma_params.c` (L1-L2) for core parameters
2. Understand `plasma_kinetic.c` for kinetic theory (L3-L4)
3. Study `plasma_mhd.c` for fluid description (L3-L4)
4. Explore `plasma_waves.c` for wave physics (L4-L6)
5. Learn `plasma_particle.c` for simulation methods (L5-L6)
6. Apply with `plasma_diagnostics.c` (L7)
7. Design with `plasma_confinement.c` (L7-L8)

## Research Frontiers (L9)
- Burning plasma (ITER, SPARC) → plasma_confinement.c
- Stellarator optimization → plasma_confinement.c
- Warm dense matter → plasma_params.c (regime classification)
- Plasma accelerators → documented only
