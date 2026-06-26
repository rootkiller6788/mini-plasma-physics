# Knowledge Graph — mini-space-plasma

## L1: Definitions
- [x] Physical constants (CODATA 2018): kB, hbar, c, me, mp, e, eps0, mu0, G, NA, sigma_SB
- [x] Rydberg energy, AU, Earth radius/mass, Solar radius/mass
- [x] Earth equatorial magnetic field B0 = 3.12e-5 T, dipole moment mu_E = 7.78e22 A.m^2
- [x] plasma_state_t: n_e, n_i, T_e, T_i, B[3], v[3], p, Z, is_magnetized
- [x] particle_orbit_t: r_L, omega_c, v_perp, v_par, mu_m, pitch_angle, guiding_center
- [x] mhd_wave_t: v_A, c_s, v_fast, v_slow, beta, theta
- [x] debye_params_t: lambda_D, N_D, omega_pe, omega_pi, tau
- [x] magnetosphere_t: R_mp, R_ss, B_dipole, L_shell, pressure_sw, Dst
- [x] solar_wind_t: v_sw, n_sw, T_sw, B_imf, Ma, Ms, M_ms, beta_sw
- [x] reconnection_t: B_in, v_in, v_out, delta, L, eta, S, E_rec
- [x] distribution_t: f[], v_min, v_max, dv, n_v, n0, T, v_th
- [x] mhd_primitive_t / mhd_conserved_t / mhd_flux_t / mhd_eigensystem_t
- [x] stix_tensor_t: S, D, P, R, L
- [x] wave_mode_t enum (12 wave types)
- [x] parker_profile_t / parker_solution_type_t
- [x] plasma_region_t enum (12 magnetospheric regions)

## L2: Core Concepts
- [x] Debye shielding (lambda_D, N_D, Coulomb coupling Gamma)
- [x] Plasma frequency (electron, ion, upper/lower hybrid)
- [x] Gyrofrequency and Larmor radius (omega_c, r_L, T_c)
- [x] Alfven velocity, sound speed, magnetosonic speed
- [x] Plasma beta (thermal/magnetic pressure ratio)
- [x] Collision frequencies (electron-ion, Coulomb logarithm, Spitzer resistivity)
- [x] MHD conserved/primitive variable conversion
- [x] MHD flux computation (8-wave Powell formulation)
- [x] MHD wave eigenspeeds (7 eigenvalues)
- [x] Ideal MHD frozen-in flux condition (E + v x B = 0)
- [x] Magnetospheric E x B convection
- [x] Ring current (Dessler-Parker-Sckopke, energy density)
- [x] Plasmasphere/Plasmapause (Carpenter-Anderson model)
- [x] Solar wind Mach numbers (Alfven, sonic, magnetosonic, fast)

## L3: Mathematical Structures
- [x] MHD equilibrium solutions (force-free, Harris sheet, Z-pinch)
- [x] Grad-Shafranov elliptic operator (finite difference discretization)
- [x] Cold plasma dielectric tensor (Stix parameters S,D,P,R,L)
- [x] Appleton-Hartree dispersion relation (A n^4 - B n^2 + C = 0)
- [x] Wave polarization (E_x/E_y ratio, helicity, ellipticity)
- [x] CMA diagram classifier (Clemmow-Mullaly-Allis)
- [x] MHD stress tensor T_ij

## L4: Fundamental Laws
- [x] Parker isothermal solar wind equation (transcendental, Newton-Raphson solve)
- [x] Mass continuity for solar wind (n v r^2 = const)
- [x] Parker spiral magnetic field (B_r ~ r^-2, B_phi ~ r^-1)
- [x] Chapman-Ferraro magnetopause (ram-magnetic pressure balance)
- [x] Dipole magnetic field (B ~ r^-3)
- [x] Rankine-Hugoniot MHD shock jump conditions
- [x] Warm plasma dielectric function (Fried-Conte, Bohm-Gross)
- [x] Bennett relation (Z-pinch equilibrium)
- [x] Dessler-Parker-Sckopke relation (ring current Dst)
- [x] Cold plasma dispersion relation (Appleton-Hartree)

## L5: Computational Methods
- [x] Lax-Friedrichs MHD solver (1D, first-order TVD)
- [x] MHD CFL timestep condition
- [x] Newton-Raphson dispersion root solver (solve_dispersion_omega)
- [x] Group velocity via central finite differences
- [x] Grad-Shafranov operator (2D finite differences)
- [x] div(B) error monitoring for MHD
- [x] Parker equation Newton-Raphson solver

## L6: Canonical Systems
- [x] Parker solar wind transonic profile (class I-V)
- [x] Parker spiral at 1 AU (B field, spiral angle, sector polarity)
- [x] CIR formation (fast/slow wind interaction)
- [x] Earth dipole field (Cartesian, McIlwain L, invariant latitude)
- [x] Chapman-Ferraro magnetopause (standoff distance, CF current)
- [x] Bow shock (standoff distance, MHD shock jumps)
- [x] Harris current sheet (1D tearing mode equilibrium)
- [x] Z-pinch equilibrium
- [x] Force-free magnetic fields (Lundquist solution)
- [x] Particle drifts (grad-B, curvature, dipole total drift, bounce period)
- [x] Magnetosheath flow deflection
- [x] Tsyganenko-like tail model
- [x] Plasmasphere density profile
- [x] Langmuir waves (Bohm-Gross)
- [x] Ion acoustic waves (L-L damping)
- [x] Alfven waves (shear, kinetic)
- [x] Whistler waves (R-mode dispersion)
- [x] Fast/slow magnetosonic waves
- [x] Lower/upper hybrid waves

## L7: Applications
- [x] Parker solar wind profile computation (corona → 1 AU)
- [x] Magnetopause standoff for quiet/fast/CME conditions
- [x] Cross-polar cap potential (Boyle 1997 empirical formula)
- [x] Dungey cycle reconnection voltage
- [x] Ring current Dst prediction
- [x] CMA propagation band classification

## L8: Advanced Topics
- [x] Kinetic Alfven waves (finite Larmor radius corrections)
- [x] Warm plasma dielectric (Landau damping regime)
- [x] Wave polarization analysis (helicity, ellipticity)
- [x] Tsyganenko stretched field model (tail current sheet)

## L9: Research Frontiers
- [x] CME magnetic cloud parameters (documented)
- [x] Heliospheric current sheet tilt evolution (documented)
- [x] Bi-Maxwellian ring current energy (anisotropy effects)
