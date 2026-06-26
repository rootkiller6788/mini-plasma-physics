# Knowledge Graph — mini-mhd

## L1 — Definitions

| Item | Type | Location |
|------|------|----------|
| Physical constants (mu_0, epsilon_0, c, k_B, e, m_e, m_p) | #define | mhd_defs.h |
| MHDState (8 primitive fields) | typedef struct | mhd_defs.h |
| MHDConserved (8 conserved fields) | typedef struct | mhd_defs.h |
| MHDFlux (8 flux components) | typedef struct | mhd_defs.h |
| MHDWaveSpeeds (4 characteristic speeds) | typedef struct | mhd_defs.h |
| MHDResistiveSource (ohmic heating + diffusion) | typedef struct | mhd_defs.h |
| MHDCoordinateSystem (4 coordinate types) | typedef enum | mhd_defs.h |
| MHDGeometry (domain descriptor) | typedef struct | mhd_defs.h |
| PlasmaParameters (13 dimensionless numbers) | typedef struct | mhd_defs.h |
| MHDBCType (5 boundary condition types) | typedef enum | mhd_numerics.h |
| MHDState (Lean 4) | structure | mhd.lean |
| Vec3 (Lean 4) | structure | mhd.lean |

## L2 — Core Concepts

| Concept | Function | Location |
|---------|----------|----------|
| Sound speed c_s = sqrt(gamma*p/rho) | mhd_sound_speed | mhd_defs.h |
| Alfven speed v_A = B/sqrt(mu_0*rho) | mhd_alfven_speed | mhd_defs.h |
| Magnetic pressure p_mag = B^2/(2*mu_0) | mhd_magnetic_pressure | mhd_defs.h |
| Magnetic tension B^2/mu_0 | mhd_magnetic_tension | mhd_defs.h |
| Plasma beta beta = 2*mu_0*p/B^2 | mhd_plasma_beta | mhd_defs.h |
| Alfven Mach M_A = v/v_A | mhd_mach_alfven | mhd_defs.h |
| Sonic Mach M_s = v/c_s | mhd_mach_sonic | mhd_defs.h |
| Magnetic Reynolds R_m = v*L/eta | mhd_magnetic_reynolds | mhd_defs.h |
| Lundquist S = v_A*L/eta | mhd_lundquist | mhd_defs.h |
| Magnetic Prandtl P_m = nu/eta | mhd_magnetic_prandtl | mhd_defs.h |
| Hartmann Ha = B*L/sqrt(mu_0*rho*nu*eta) | mhd_hartmann | mhd_defs.h |
| Ion inertial length d_i = c/omega_pi | mhd_ion_inertial_length | mhd_defs.h |
| Ion Larmor radius rho_i | mhd_ion_larmor_radius | mhd_defs.h |
| Debye length lambda_D | mhd_debye_length | mhd_defs.h |
| Ion cyclotron frequency omega_ci | mhd_ion_cyclotron_frequency | mhd_defs.h |
| Plasma frequency omega_pe | mhd_plasma_frequency | mhd_defs.h |
| Elsasser number Lambda | mhd_elsasser | mhd_defs.h |
| Frozen-in flux theorem | mhd_alfven_theorem_check | mhd_eqns.c |
| Jeans criterion (MHD) | mhd_jeans_length | mhd_eqns.c |

## L3 — Mathematical Structures

| Structure | Implementation | Location |
|-----------|---------------|----------|
| Vector magnitude, dot, cross product | inline functions | mhd_defs.h |
| Energy densities (kinetic, magnetic, internal, total) | inline functions | mhd_defs.h |
| J = curl(B)/mu_0 (Ampere) | mhd_current_from_curlB | mhd_defs.h |
| Lorentz force J x B | mhd_lorentz_force | mhd_defs.h |
| Joule heating eta*J^2 | mhd_joule_heating | mhd_defs.h |
| Primitive <-> Conserved transform | conversion functions | mhd_defs.c |
| MHD flux vector computation | mhd_flux_compute | mhd_defs.c |
| Maximum wave speed (CFL) | mhd_max_wavespeed | mhd_defs.c |
| Characteristic decomposition | mhd_wave_decomposition | mhd_waves.c |
| Phase/group velocity | mhd_wave_phase_velocity / group_velocity | mhd_waves.c |

## L4 — Fundamental Laws

| Law | Implementation | Verification |
|-----|---------------|-------------|
| Continuity eq: d_t(rho) + div(rho*v) = 0 | mhd_continuity_residual | test |
| Momentum eq: d_t(rho*v) + div(...) = 0 | mhd_momentum_residual | test |
| Induction (ideal): d_t(B) = curl(v x B) | mhd_induction_residual | test |
| Induction (resistive): d_t(B) = curl(v x B) + eta*Lap(B) | mhd_induction_resistive_residual | test |
| Energy (ideal): d_t(E) + div(F_E) = 0 | mhd_energy_residual | test |
| Energy (resistive): d_t(E) + div(F_E) = eta*J^2 | mhd_energy_resistive_residual | test |
| Ohm's law: E + v x B = eta*J | mhd_ohm_law | test |
| Poynting vector S = E x B/mu_0 | mhd_poynting_vector | test |
| div(B) = 0 constraint | mhd_divB | test |
| Energy principle (Bernstein) | mhd_energy_principle_fluid | mhd_instability.c |
| Lorentz_no_work (Lean) | theorem lorentz_no_work | mhd.lean |
| div_curl_zero (Lean) | theorem div_curl_zero | mhd.lean |
| Helicity conservation (Lean) | theorem helicity_conservation_ideal | mhd.lean |

## L5 — Computational Methods

| Method | Function | Location |
|--------|----------|----------|
| FD gradient 2nd order | mhd_fd_gradient_2nd | mhd_numerics.c |
| FD gradient 4th order | mhd_fd_gradient_4th | mhd_numerics.c |
| FD Laplacian 2nd order | mhd_fd_laplacian_2nd | mhd_numerics.c |
| FD divergence 3D | mhd_fd_divergence_2nd | mhd_numerics.c |
| FD curl 3D | mhd_fd_curl_2nd | mhd_numerics.c |
| Minmod limiter | mhd_limiter_minmod | mhd_numerics.c |
| Superbee limiter | mhd_limiter_superbee | mhd_numerics.c |
| Van Leer limiter | mhd_limiter_van_leer | mhd_numerics.c |
| MC limiter | mhd_limiter_mc | mhd_numerics.c |
| MUSCL-Hancock | mhd_muscl_extrapolate | mhd_numerics.c |
| HLL simplified flux | mhd_hllc_simple | mhd_numerics.c |
| HLLD flux (Miyoshi-Kusano) | mhd_hlld_flux | mhd_numerics.c |
| Roe flux | mhd_roe_flux | mhd_numerics.c |
| RK2 (Heun) | mhd_rk2_step | mhd_numerics.c |
| RK3 (TVD Shu-Osher) | mhd_rk3_step | mhd_numerics.c |
| RK4 (classical) | mhd_rk4_step | mhd_numerics.c |
| Constrained Transport | mhd_ct_emf_2d, mhd_ct_update_B_2d | mhd_numerics.c |
| Powell 8-wave source | mhd_powell_source | mhd_numerics.c |
| Boundary conditions (5 types) | mhd_apply_boundary_1d/2d | mhd_numerics.c |
| CFL timestep | mhd_cfl_timestep | mhd_numerics.c |

## L6 — Canonical Systems

| System | Implementation | Example |
|--------|---------------|---------|
| Alfven wave dispersion | mhd_alfven_dispersion | alfven_wave.c |
| Magnetosonic dispersion | mhd_magnetosonic_dispersion | waves test |
| Friedrichs diagram | mhd_friedrichs_diagram | waves test |
| Z-pinch (Bennett) | mhd_zpinch_equilibrium | zpinch.c |
| Theta-pinch | mhd_thetapinch_equilibrium | — |
| Screw pinch | mhd_screwpinch_equilibrium | — |
| Bennett relation | mhd_bennett_temperature | zpinch.c |
| Soloviev tokamak | mhd_soloviev_psi, mhd_soloviev_field | tokamak.c |
| Safety factor q | mhd_safety_factor_q | tokamak.c |
| Kruskal-Shafranov limit | mhd_kruskal_shafranov_condition | tokamak.c |

## L7 — Applications

| Application | Location |
|-------------|----------|
| PlasmaParameter computation | mhd_compute_all_parameters (mhd_defs.c) |
| Bennett relation verification | zpinch.c example |
| Safety factor profile analysis | tokamak.c example |

## L8 — Advanced Topics

| Topic | Implementation | Status |
|-------|---------------|--------|
| Grad-Shafranov equation | mhd_grad_shafranov_operator | Complete |
| Taylor relaxation (Bessel function) | mhd_taylor_relaxation | Complete |
| Beltrami fields | mhd_beltrami_check | Complete |
| Ballooning parameter | mhd_ballooning_parameter | Complete |
| Mercier criterion | mhd_mercier_criterion | Complete |
| Tearing mode growth rate | mhd_tearing_growth_rate | Complete |
| MHD Rayleigh-Taylor | mhd_rayleigh_taylor_growth | Complete |
| MHD Kelvin-Helmholtz | mhd_kelvin_helmholtz_growth | Complete |
| Linear stability scan | mhd_linear_stability_scan | Complete |
| Beltrami theorem (Lean) | theorem beltrami_iterative | Complete |

## L9 — Research Frontiers

| Topic | Status |
|-------|--------|
| Magnetic helicity conservation | Documented (Lean theorem) |
| Fusion MHD (tokamak stability) | Documented (KS, ballooning, tearing) |
| Coronal heating problem | Documented |
| Dynamo theory | Documented (Elsasser number) |
| Magnetic reconnection | Documented (Lundquist, Sweet-Parker scaling) |
