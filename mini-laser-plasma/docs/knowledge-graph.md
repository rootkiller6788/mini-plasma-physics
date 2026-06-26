# Knowledge Graph ！ mini-laser-plasma

## L1: Definitions (Complete)

| # | Concept | C Struct/Function | Lean Definition |
|---|---------|-------------------|-----------------|
| 1 | Plasma frequency omega_p | plasma_frequency() | PlasmaDerived.wp_over_w0 |
| 2 | Critical density n_c | critical_density() | critical_density_cutoff |
| 3 | Debye length lambda_D | debye_length() | debye_length_pos |
| 4 | Plasma parameter N_D | debye_sphere_particles() | plasma_parameter_collective |
| 5 | Electron thermal velocity v_the | electron_thermal_velocity() | ！ |
| 6 | Ion sound speed c_s | ion_sound_speed() | ！ |
| 7 | Normalized vector potential a0 | normalized_vector_potential() | LaserPulse.a0 |
| 8 | Ponderomotive potential Phi_p | ponderomotive_potential() | ponderomotive_quadratic |
| 9 | Relativistic intensity I_rel | relativistic_intensity() | ！ |
| 10 | Keldysh parameter gamma_K | keldysh_parameter() | ！ |
| 11 | ADK tunneling rate | adk_tunneling_rate() | ！ |
| 12 | Plasma beta | plasma_beta() | ！ |
| 13 | Cold wave-breaking field E_wb | cold_wavebreaking_field() | wakefield_gradient_bound |
| 14 | Plasma wavelength lambda_p | plasma_wavelength() | WakefieldRegime |

## L2: Core Concepts (Complete)

| # | Concept | Implementation |
|---|---------|---------------|
| 1 | Ponderomotive force | ponderomotive_force() |
| 2 | Relativistic self-focusing | relativistic_self_focusing_threshold() |
| 3 | Ponderomotive self-focusing | ponderomotive_self_focusing_threshold() |
| 4 | Relativistic transparency | relativistic_transparency_density() |
| 5 | Hole boring | hole_boring_velocity() |
| 6 | SRS (Stimulated Raman Scattering) | srs_growth_rate_backscatter/forward() |
| 7 | SBS (Stimulated Brillouin Scattering) | sbs_growth_rate() |
| 8 | TPD (Two-Plasmon Decay) | tpd_growth_rate() |
| 9 | Filamentation instability | filamentation_growth_rate() |
| 10 | Inverse bremsstrahlung absorption | inverse_bremsstrahlung_coefficient() |
| 11 | Resonance absorption | resonance_absorption_fraction() |
| 12 | Brunel effect (vacuum heating) | brunel_absorption_fraction() |
| 13 | JxB heating | jxb_heating_fraction() |
| 14 | Bubble/blowout regime | bubble_accelerating_field() |
| 15 | BSI threshold | barrier_suppression_intensity() |
| 16 | Avalanche ionization | avalanche_ionization_rate() |
| 17 | Ponderomotive scattering | ponderomotive_scattering_angle() |
| 18 | Convective vs absolute instability | srs_absolute_threshold_a0() |

## L3: Mathematical Structures (Complete)

| # | Structure | Implementation |
|---|-----------|---------------|
| 1 | Cold plasma dielectric epsilon(omega) | cold_plasma_permittivity() |
| 2 | Drude collisional dielectric | collisional_permittivity() |
| 3 | Bohm-Gross warm dielectric | warm_plasma_permittivity() |
| 4 | Plasma refractive index N | plasma_refractive_index() |
| 5 | Group velocity dispersion | group_velocity() |
| 6 | Paraxial wave equation state | ParaxialState struct |
| 7 | Particle 6D phase space | Particle3D struct |
| 8 | EM field tensor (E,B) | EMField struct |
| 9 | Gaussian beam optics | GaussianBeam struct |
| 10 | Lean dielectric structure | DielectricFunction |

## L4: Fundamental Laws (Complete)

| # | Law/Theorem | Verification |
|---|-------------|-------------|
| 1 | EM dispersion: omega^2 = wp^2 + c^2 k^2 | critical_density() + tests |
| 2 | Debye screening (Poisson-Boltzmann) | debye_length() + debye_length_pos |
| 3 | Spitzer collision frequency | electron_ion_collision_frequency() |
| 4 | Lorentz force (relativistic) | lorentz_force() |
| 5 | Ponderomotive potential from Lorentz | ponderomotive_potential() |
| 6 | SRS three-wave matching | srs_matching() + test_srs_matching |
| 7 | SBS frequency matching | sbs_matching() |
| 8 | Volkov analytic solution | plane_wave_orbit_analytic() |
| 9 | Wakefield: L_deph = (lambda_p/2)*(nc/ne) | dephasing_length_1D() |
| 10 | Wave-breaking: E_wb = m_e c omega_p / e | cold_wavebreaking_field() |
| 11 | ADK tunneling: Keldysh-Faisal-Reiss | adk_tunneling_rate() |
| 12 | Ginzburg function for resonance absorption | ginzburg_function() |

## L5: Computational Methods (Complete)

| # | Method | Implementation |
|---|--------|---------------|
| 1 | Boris relativistic particle pusher | boris_push() |
| 2 | RK4 particle integrator | rk4_push() |
| 3 | Split-step paraxial propagation | paraxial_step_vacuum/plasma() |
| 4 | Ray tracing eikonal step | ray_trace_step() |
| 5 | Snell law for stratified plasma | snells_law_plasma() |
| 6 | Rate equation integration (ionization) | integrate_ionization_rate() |
| 7 | Convective gain (Rosenbluth) | srs_convective_gain() |
| 8 | Combined absorption calculation | total_absorption_fraction() |
| 9 | Trajectory recording/diagnostics | particle_track_record() |
| 10 | PlasmaDerived one-shot computation | compute_all_derived() |

## L6: Canonical Systems (Complete)

| # | System | Example/Demo |
|---|--------|-------------|
| 1 | Critical density for key laser wavelengths | example_critical_density.c |
| 2 | LWFA in three operating regimes | example_wakefield.c |
| 3 | SRS/SBS/TPD analysis for ICF plasma | example_instabilities.c |
| 4 | Debye sphere in laser-produced plasma | test_debye_sphere |

## L7: Applications (Partial ！ 2+ items)

| # | Application | Reference |
|---|-------------|-----------|
| 1 | GeV-scale LWFA electron acceleration | example_wakefield.c (W_max calculation) |
| 2 | ICF hohlraum plasma instability analysis | example_instabilities.c (SRS/SBS/TPD) |

## L8: Advanced Topics (Partial ！ 1+ item)

| # | Topic | Implementation |
|---|-------|---------------|
| 1 | Rosenbluth convective gain theory | srs_convective_gain() |
| 2 | Brunel vacuum heating (collisionless) | brunel_absorption_fraction() |

## L9: Research Frontiers (Partial ！ documented)

| # | Topic | Status |
|---|-------|--------|
| 1 | Laser-driven ion acceleration (RPA/LS) | hole_boring_velocity() as building block | 
| 2 | Fast ignition fusion | Absorption models applicable |
