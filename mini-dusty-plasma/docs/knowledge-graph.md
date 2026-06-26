# Knowledge Graph -- mini-dusty-plasma

## L1: Definitions Complete

| # | Concept | C Implementation | Lean Definition |
|---|---------|-----------------|-----------------|
| 1 | Dust grain (mass, radius, charge) | DustGrain struct in dusty_plasma.h | structure DustGrain |
| 2 | Plasma state (n_e, n_i, n_d, T_e, T_i, T_d) | DustPlasmaState struct | structure DustPlasmaState |
| 3 | Debye length (electron, ion, total) | dust_debye_electron/ion/total() | debyeLength, totalDebyeLength |
| 4 | Plasma frequency (electron, ion, dust) | dust_plasma_freq_electron/ion/dust() | dustPlasmaFreq |
| 5 | Havnes parameter P | dust_havnes_parameter() | havnesParameter |
| 6 | Coulomb coupling Gamma | dust_coulomb_coupling() | coulombCoupling |
| 7 | Yukawa coupling Gamma* | dust_yukawa_coupling() | -- |
| 8 | Dust acoustic speed c_da | dust_acoustic_speed() | dustAcousticSpeed |
| 9 | Ion acoustic speed c_s | dust_ion_acoustic_speed() | ionAcousticSpeed |
| 10 | Wave mode (omega, k, damping) | WaveMode struct | -- |
| 11 | Yukawa pair potential/force | YukawaPair struct, yukawa_*() | yukawaPotential, yukawaForceMagnitude |
| 12 | Dust charging current components | DustChargingCurrent struct | -- |
| 13 | Force components on a grain | DustForceResult struct | -- |
| 14 | Pair correlation function | PairCorrelation struct | -- |
| 15 | Dust crystal lattice | DustCrystal struct | -- |
| 16 | Grain trajectory state | DustTrajectory struct | -- |

## L2: Core Concepts Complete

| # | Concept | Implementation |
|---|---------|---------------|
| 1 | Quasineutrality (n_i = n_e + Z_d n_d) | dust_plasma_state_init(); Lean: quasineutrality |
| 2 | Dust charging (OML theory) | dust_oml_electron/ion_current() |
| 3 | Floating potential equilibrium | dust_floating_potential_solve() |
| 4 | Dust-acoustic wave (DAW) | dust_acoustic_wave_dispersion() |
| 5 | Dust-ion-acoustic wave (DIAW) | dust_ion_acoustic_wave_dispersion() |
| 6 | Dust lattice wave (DLW) | dust_lattice_wave_1d/2d_*() |
| 7 | Coulomb crystallization | dust_crystal_condition(); Lean: crystallizationCondition |
| 8 | Ion drag force | dust_ion_drag_collection/orbit/total() |
| 9 | Neutral drag (Epstein) | dust_neutral_drag_force(); Lean: epsteinCollisionFreq |
| 10 | Dust diffusion and mobility | dust_diffusion_coefficient(), dust_mobility(); Lean: einsteinDiffusion |
| 11 | Thermophoretic force | dust_thermophoretic_force() |
| 12 | Dust in plasma sheath | dust_sheath_electric_field/potential() |
| 13 | Dust void formation | dust_void_radius_estimate(); Lean: voidFormationCondition |
| 14 | Dust cyclotron motion | dust_cyclotron_frequency() |
| 15 | Debye-Huckel screening | dust_debye_huckel_potential/field() |

## L3: Mathematical Structures Complete

| # | Structure | Implementation |
|---|-----------|---------------|
| 1 | OML charging integral | Closed-form expressions in dusty_charging.c |
| 2 | Yukawa potential (exponential screening) | yukawa_potential() |
| 3 | Lattice sums (Madelung) | dust_yukawa_madelung() |
| 4 | Yukawa force vector calculus | yukawa_force_vector() |
| 5 | Dispersion relation (omega^2 vs k^2) | WaveMode struct, analytical dispersion |
| 6 | Coulomb logarithm | dust_ion_coulomb_logarithm() |
| 7 | Cunningham slip correction | dust_general_drag_force() |
| 8 | Richardson-Dushman thermionic emission | dust_thermionic_current() |
| 9 | Box-Muller Gaussian sampling | dust_box_muller() |
| 10 | Numerical derivative (central diff) | dust_net_current_derivative() |
| 11 | HNC closure (pair correlation) | dust_pair_correlation_hnc() |

## L4: Fundamental Laws Complete

| # | Law/Theorem | C Implementation | Lean Statement |
|---|------------|-----------------|----------------|
| 1 | Gauss law (spherical grain) | Q_d = 4*pi*eps0*a*phi_f | equilibriumCharge |
| 2 | Current balance (charging eq.) | I_e + I_i = 0 at phi_f | floatingPotentialEquilibrium |
| 3 | Dust momentum equation | dust_langevin_deterministic_accel() | -- |
| 4 | DAW dispersion relation | daw_dispersion() | dawDispersion |
| 5 | DIAW dispersion relation | dust_ion_acoustic_wave_dispersion() | -- |
| 6 | DLW dispersion (1D long./trans.) | dust_lattice_wave_1d_*() | -- |
| 7 | Yukawa force law | F = (Q^2/(4*pi*eps0*r^2))*(1+r/lambda)*exp(-r/lambda) | yukawaForceMagnitude |
| 8 | Epstein drag law | nu_dn = (8/3)*sqrt(2/pi)*a^2*n_n*v_thn*m_n/m_d | epsteinCollisionFreq |
| 9 | Einstein relation | D_d = k_B*T_d/(m_d*nu_dn) | einsteinRelation |
| 10 | Stefan-Boltzmann (radiative cooling) | dust_radiative_cooling_power() | -- |
| 11 | Ion drag (Barnes model) | dust_ion_drag_collection/orbit() | ionDragCollection |
| 12 | Modified Bohm criterion | dust_modified_bohm_velocity() | modifiedBohmVelocity |
| 13 | Landau damping (DAW kinetic) | dust_acoustic_wave_landau_damping() | -- |
| 14 | Read-Shockley grain boundary | dust_grain_boundary_energy() | -- |

## L5: Computational Methods Complete

| # | Method | Implementation |
|---|--------|---------------|
| 1 | Newton-Raphson root finding | dust_floating_potential_solve() |
| 2 | Bisection method | dust_steady_state_temperature() |
| 3 | Forward Euler integration (ODE) | dust_integrate_charge_dynamics() |
| 4 | RK4 integration (6D phase space) | dust_rk4_step() |
| 5 | Box-Muller (Gaussian RNG) | dust_box_muller() |
| 6 | Lattice sum computation | dust_yukawa_madelung() |
| 7 | Numerical integration (trapezoidal) | dust_structure_factor() |
| 8 | Central difference derivative | dust_net_current_derivative() |
| 9 | Minimum image convention | dust_crystal_total_energy() |

## L6: Canonical Systems Complete

| # | System | Implementation |
|---|--------|---------------|
| 1 | 1D dust-acoustic wave | example_daw.c |
| 2 | Dust charging equilibrium | example_charging.c |
| 3 | 2D dust crystal phase diagram | example_crystal.c |
| 4 | Dust levitation in sheath | example_sheath.c |
| 5 | 1D longitudinal DLW | dust_lattice_wave_1d_longitudinal() |
| 6 | 1D transverse DLW | dust_lattice_wave_1d_transverse() |
| 7 | 2D out-of-plane DLW | dust_lattice_wave_2d_out_of_plane() |
| 8 | 3D bcc Madelung constant | dust_yukawa_madelung(lattice_type=2) |
| 9 | Dust void in RF discharge | dust_void_radius_estimate() |

## L7: Applications Complete (4 apps)

| # | Application | Implementation |
|---|-------------|---------------|
| 1 | Semiconductor manufacturing | example_sheath.c -- levitation, ISO cleanroom |
| 2 | Fusion reactor dust (ITER) | example_sheath.c -- thermophoretic force |
| 3 | Space dusty plasmas (Saturn) | example_sheath.c -- Havnes parameter |
| 4 | Nanoparticle synthesis | dust_nucleation_rate() |

## L8: Advanced Topics Complete (7 topics)

| # | Topic | Implementation |
|---|-------|---------------|
| 1 | Strongly coupled Yukawa systems | dust_phase_determine(); Lean: stronglyCoupled |
| 2 | 2D dust lattice waves | dust_lattice_wave_2d_out_of_plane() |
| 3 | Phase transition dynamics | dust_nucleation_rate(), dust_critical_coupling_yukawa() |
| 4 | Pair correlation function (HNC) | dust_pair_correlation_hnc() |
| 5 | Static structure factor S(k) | dust_structure_factor() |
| 6 | Mach cone / ion focusing | dust_oml_ion_current_drift() |
| 7 | Grain boundary energy | dust_grain_boundary_energy() |

## L9: Research Frontiers Partial (3 frontiers)

| # | Frontier | Documentation |
|---|----------|--------------|
| 1 | Dust in magnetic confinement fusion (ITER) | tokamakDustForceBalance in Lean |
| 2 | Quantum dusty plasmas (nano-dust) | quantumDustyCondition in Lean |
| 3 | Active (smart) dust particles | Documented in knowledge-graph |