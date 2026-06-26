# Knowledge Graph — mini-plasma-fundamental

## L1: Definitions — COMPLETE

| # | Item | File | Status |
|---|------|------|--------|
| 1 | Debye length λ_D | plasma_params.h/c | Complete |
| 2 | Plasma frequency ω_p | plasma_params.h/c | Complete |
| 3 | Electron/ion plasma frequency | plasma_params.h/c | Complete |
| 4 | Thermal velocity v_th | plasma_params.h/c | Complete |
| 5 | Ion sound speed c_s (Bohm) | plasma_params.h/c | Complete |
| 6 | Alfven speed v_A | plasma_params.h/c | Complete |
| 7 | Plasma beta β | plasma_params.h/c | Complete |
| 8 | Gyrofrequency ω_c | plasma_params.h/c | Complete |
| 9 | Gyroradius (Larmor radius) r_L | plasma_params.h/c | Complete |
| 10 | Coulomb logarithm ln Λ | plasma_params.h/c | Complete |
| 11 | Spitzer resistivity η_∥ | plasma_params.h/c | Complete |
| 12 | Saha ionization equation | plasma_params.h/c | Complete |
| 13 | Maxwellian distribution | plasma_params.h/c | Complete |
| 14 | Plasma parameter g | plasma_params.h/c | Complete |
| 15 | Debye sphere count N_D | plasma_params.h/c | Complete |

## L2: Core Concepts — COMPLETE

| # | Item | File | Status |
|---|------|------|--------|
| 1 | Quasi-neutrality condition | plasma_params.h/c | Complete |
| 2 | Debye shielding (Yukawa potential) | plasma_params.h/c | Complete |
| 3 | Collective behavior | plasma_params.h/c | Complete |
| 4 | Plasma regime classification | plasma_params.h/c | Complete |
| 5 | ExB drift | plasma_particle.h/c | Complete |
| 6 | Grad-B drift | plasma_particle.h/c | Complete |
| 7 | Curvature drift | plasma_particle.h/c | Complete |
| 8 | Polarization drift | plasma_particle.h/c | Complete |
| 9 | MHD wave modes (Alfven, fast, slow) | plasma_mhd.h/c | Complete |
| 10 | Magnetic pressure and tension | plasma_mhd.h/c | Complete |

## L3: Mathematical Structures — COMPLETE

| # | Item | File | Status |
|---|------|------|--------|
| 1 | Vlasov equation (1D1V) | plasma_kinetic.h/c | Complete |
| 2 | Fokker-Planck collision operator | plasma_kinetic.h/c | Complete |
| 3 | Rosenbluth potentials | plasma_kinetic.h/c | Complete |
| 4 | Ideal MHD flux functions | plasma_mhd.h/c | Complete |
| 5 | Div B cleaning (Powell 8-wave) | plasma_mhd.h/c | Complete |
| 6 | Cold plasma dielectric tensor (Stix) | plasma_waves.h/c | Complete |
| 7 | Pitch-angle scattering operator | plasma_kinetic.h/c | Complete |
| 8 | Phase-space structures (1D1V, 2D) | plasma_kinetic.h | Complete |

## L4: Fundamental Laws — COMPLETE

| # | Item | File | Status |
|---|------|------|--------|
| 1 | Maxwell equations (plasma context) | plasma_mhd.h/c | Complete |
| 2 | Boltzmann equation with collisions | plasma_kinetic.h/c | Complete |
| 3 | Landau damping (collisionless) | plasma_kinetic.h/c | Complete |
| 4 | Ideal MHD equations (conservation form) | plasma_mhd.h/c | Complete |
| 5 | Grad-Shafranov equation (MHD equilibrium) | plasma_mhd.h/c | Complete |
| 6 | SOR solver for Grad-Shafranov | plasma_mhd.c | Complete |
| 7 | Resistive MHD (Lundquist number) | plasma_mhd.h/c | Complete |
| 8 | Sweet-Parker / Petschek reconnection | plasma_mhd.h/c | Complete |
| 9 | Kruskal-Shafranov kink stability | plasma_mhd.h/c | Complete |
| 10 | Spitzer resistivity (formal derivation) | plasma_formal.lean | Complete |
| 11 | Lawson criterion for fusion | plasma_diagnostics.c | Complete |
| 12 | Quasi-neutrality theorem | plasma_formal.lean | Complete |

## L5: Computational Methods — COMPLETE

| # | Item | File | Status |
|---|------|------|--------|
| 1 | Boris particle pusher | plasma_particle.h/c | Complete |
| 2 | PIC charge deposition (NGP, CIC) | plasma_particle.h/c | Complete |
| 3 | PIC field interpolation (CIC) | plasma_particle.h/c | Complete |
| 4 | PIC Poisson solver (tridiagonal) | plasma_particle.h/c | Complete |
| 5 | PIC full cycle | plasma_particle.h/c | Complete |
| 6 | Vlasov advection (operator splitting) | plasma_kinetic.h/c | Complete |
| 7 | Plasma dispersion function Z(ζ) | plasma_kinetic.h/c | Complete |
| 8 | Dispersion root finder (bisection) | plasma_waves.h/c | Complete |
| 9 | Langmuir probe I-V analysis | plasma_diagnostics.c | Complete |
| 10 | FDTD for Maxwell equations | plasma_kinetic.h/c | Complete |

## L6: Canonical Systems — COMPLETE

| # | Item | File | Status |
|---|------|------|--------|
| 1 | Langmuir waves (Bohm-Gross) | plasma_waves.h/c | Complete |
| 2 | Ion acoustic waves | plasma_waves.h/c | Complete |
| 3 | Alfven waves (shear) | plasma_waves.h/c | Complete |
| 4 | Fast/slow magnetosonic waves | plasma_waves.h/c | Complete |
| 5 | Whistler (R-wave, helicon) | plasma_waves.h/c | Complete |
| 6 | Lower hybrid waves | plasma_waves.h/c | Complete |
| 7 | Two-stream instability | plasma_waves.h/c | Complete |
| 8 | Weibel instability | plasma_waves.h/c | Complete |
| 9 | Firehose/mirror instabilities | plasma_waves.h/c | Complete |
| 10 | Kelvin-Helmholtz instability | plasma_waves.h/c | Complete |
| 11 | Rayleigh-Taylor instability | plasma_waves.h/c | Complete |
| 12 | Banana orbits (tokamak) | plasma_particle.h/c | Complete |
| 13 | Trapped particle fraction | plasma_particle.h/c | Complete |
| 14 | Tokamak safety factor | plasma_mhd.h/c | Complete |

## L7: Applications — COMPLETE (4 apps)

| # | Item | File | Status |
|---|------|------|--------|
| 1 | ITER parameters | plasma_confinement.c | Complete |
| 2 | SPARC (CFS/MIT) parameters | plasma_confinement.c | Complete |
| 3 | Solar wind diagnostics | plasma_diagnostics.c | Complete |
| 4 | Magnetopause distance (space plasma) | plasma_diagnostics.c | Complete |
| 5 | Wendelstein 7-X parameters | plasma_confinement.c | Complete |
| 6 | Fusion power/reactivity (Bosch-Hale) | plasma_diagnostics.c | Complete |
| 7 | Langmuir probe analysis | plasma_diagnostics.c | Complete |
| 8 | Microwave interferometry | plasma_diagnostics.c | Complete |
| 9 | Thomson scattering | plasma_diagnostics.c | Complete |
| 10 | Plasma parameter space explorer | demos/ | Complete |

## L8: Advanced Topics — PARTIAL (4 topics)

| # | Item | File | Status |
|---|------|------|--------|
| 1 | Three-wave resonant coupling | plasma_waves.h/c | Complete |
| 2 | Parametric decay instability | plasma_waves.h/c | Complete |
| 3 | H-mode pedestal physics | plasma_confinement.c | Complete |
| 4 | Mirror machine confinement (Pastukhov) | plasma_confinement.c | Complete |
| 5 | Nonlinear wave coupling (full theory) | — | Missing |
| 6 | Neoclassical transport (full matrix) | — | Missing |

## L9: Research Frontiers — DOCUMENTED ONLY

| # | Item | Status |
|---|------|--------|
| 1 | Burning plasma physics (ITER) | Documented in confinement |
| 2 | Stellarator optimization | Documented in confinement |
| 3 | Advanced divertor concepts | Not implemented |
| 4 | Plasma-based accelerators | Not implemented |
| 5 | Quantum plasma effects | Not implemented |
