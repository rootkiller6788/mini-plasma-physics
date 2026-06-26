# Course Tree -- mini-laser-plasma

## Prerequisite Dependency Tree

mini-classical-mechanics: Lagrangian/Hamiltonian, phase space
  -> mini-electromagnetism: Maxwell eqns, EM waves, Poynting
  -> mini-quantum-mechanics: atomic levels, tunneling (Keldysh)
  -> mini-optics-photonics: Gaussian beams, paraxial wave eq.
  -> mini-thermodynamics-statphys: Maxwell-Boltzmann, EOS
  -> mini-laser-plasma [THIS MODULE]

Downstream:
  -> mini-solid-state-physics: warm dense matter
  -> mini-nuclear-physics: ICF fusion, neutron yield
  -> mini-fluid-dynamics: plasma fluid models
  -> mini-particle-physics: QED at extreme a0
  -> mini-quantum-field-theory: strong-field QED

## Internal Dependency Graph

plasma_constants.h -> plasma_params.h -> laser_plasma.h
  laser_plasma.h -> laser_propagation.h, instabilities.h
  plasma_params.h -> wakefield.h, absorption.h, ionization.h
  particle_motion.h (independent)

## Knowledge Flow

L1 Definitions -> All layers: fundamental plasma parameters
L2 Concepts -> L5-L6: physical mechanisms -> numerical methods
L3 Math Structures -> L5: structures enable algorithms
L4 Fundamental Laws -> L2: laws validated through tests
L5 Computational Methods -> L6-L7: solve canonical problems
L6 Canonical Systems -> L7-L8: foundation for applications
