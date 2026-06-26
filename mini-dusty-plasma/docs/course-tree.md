# Course Tree -- mini-dusty-plasma

## Prerequisites

Plasma Fundamentals (mini-plasma-fundamental)
  |-- Debye length
  |-- Plasma frequency
  |-- Quasineutrality
  |-- Fluid equations
      |
      v
Dusty Plasma Physics (this module)
  |-- Dust charging (OML theory)
  |-- Dust-acoustic waves (DAW)
  |-- Dust-ion-acoustic waves (DIAW)
  |-- Dust lattice waves (DLW)
  |-- Coulomb crystallization
  |-- Yukawa potential & forces
  |-- Transport coefficients
  |-- Sheath & void physics
      |
      v
Advanced Plasma Modules
  |-- mini-fusion-plasma (dust in tokamaks)
  |-- mini-space-plasma (dust in rings/comets)
  |-- mini-industrial-plasma (dust contamination)
  |-- mini-waves-instabilities (DAW/DIAW instabilities)

## Internal Dependency Graph

dusty_constants.h (L0)
  |
  v
dusty_plasma.h (L1-L2) --> dusty_plasma_core.c
  |
  |-- dusty_charging.h --> dusty_charging.c
  |-- dusty_waves.h ---- > dusty_waves.c
  |-- dusty_crystal.h --> dusty_crystal.c
  |-- dusty_transport.h -> dusty_transport.c
  |-- dusty_forces.h ----> dusty_forces.c
  |-- dusty_potential.c (standalone)
  |-- dusty_dynamics.c (standalone)

dusty_plasma.lean (L1-L9 formalization, independent of C)

## Learning Path

1. Start with dusty_constants.h -- physical constants
2. Read dusty_plasma.h -- core definitions (L1-L2)
3. Study dusty_charging.c -- OML theory (L3-L4)
4. Explore dusty_waves.c -- wave modes (L4-L6)
5. Understand dusty_crystal.c -- phase transitions (L6-L8)
6. Learn dusty_transport.c -- collisions & diffusion (L3-L5)
7. Apply dusty_forces.c -- force balance (L2-L5)
8. Run dusty_dynamics.c -- trajectories (L5)
9. Simulate dusty_potential.c -- sheath & void (L4-L7)
10. Verify with dusty_plasma.lean -- formal statements
11. Run tests/test_dusty.c -- comprehensive testing
12. Study examples/ -- canonical problems (L6-L7)