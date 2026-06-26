# Knowledge Graph — mini-fusion-plasma

## L1: Definitions (Complete)
- PlasmaParameters, DebyeSphere, LawsonCriteria, FusionReactionRate
- EnergyBalance, TokamakGeometry, ConfinementScaling
- GradShafranovSolution, HeatingSystem, TritiumBreedingRatio
- DivertorParameters, FusionFuelType, FluxSurface
- Physical constants (CODATA 2018): kB, e, me, mp, epsilon0, mu0, c
- Fusion-specific: E_DT, E_DD, E_DHe3, E_alpha

## L2: Core Concepts (Complete)
- Debye shielding: debye_length, plasma_parameter
- Plasma oscillation: plasma_frequency, ion_plasma_frequency
- Coulomb collisions: coulomb_logarithm, collision frequencies
- Thermal velocities: electron_thermal_velocity, ion_thermal_velocity
- Gyro-motion: ion/electron gyroradius, cyclotron frequencies
- MHD waves: alfven_speed, ion_sound_speed, fast/slow magnetosonic
- Plasma beta: plasma_beta, plasma_beta_poloidal
- Safety factor: safety_factor_cylindrical, safety_factor_from_current
- Tokamak geometry: aspect_ratio, plasma_volume, plasma_surface_area
- Density limit: greenwald_density_limit, greenwald_fraction

## L3: Mathematical Structures (Complete)
- 3D vectors (Vec3) and 3x3 matrices (Mat3x3)
- Flux coordinates: FluxSurface (psi, theta, zeta, q, shear)
- Collisionality: nu_star, trapped_particle_fraction
- Spitzer resistivity: spitzer_resistivity
- Energy equilibration: energy_equilibration_time
- MHD dimensionless numbers: magnetic_reynolds_number, lundquist_number
- Grad-Shafranov operator: delta_star_psi

## L4: Fundamental Laws (Complete)
- Lawson criterion: lawson_criterion, triple_product
- Ignition condition: ignition_condition, breakeven_condition
- Fusion power density: fusion_power_density_dt/dd/dhe3/pb11
- Radiation: bremsstrahlung, cyclotron, recombination, line radiation
- Heating: ohmic_heating_power, nbi_heating_power_deposited
- Alpha physics: alpha_power_density, alpha_slowing_down_time, alpha_critical_energy
- Fusion gain: fusion_gain_Q, engineering_fusion_gain, alpha_heating_fraction
- Energy confinement: energy_confinement_time, plasma_stored_energy
- IPB98(y,2) and ITER-89P scaling laws
- MHD equilibrium: force_balance_residual, virial_integrand

## L5: Computational Methods (Complete)
- Bosch-Hale D-T cross-section parameterization (5-parameter fit, <0.25% error)
- Bosch-Hale D-D and D-He3 reactivity
- Maxwellian reactivity integration (4-point Gauss-Laguerre quadrature)
- Neoclassical transport: banana, plateau, Pfirsch-Schluter regimes
- Anomalous transport: Bohm and Gyro-Bohm diffusion
- Equilibrium reconstruction: finite-difference GS operator

## L6: Canonical Systems (Complete)
- ITER baseline: 500 MW, Q=10, R=6.2m, B=5.3T
- DEMO: 2 GW, Q=25-50, R=9.0m, B=5.7T
- SPARC: compact high-field, R=1.85m, B=12.2T (HTS)
- JET D-T record: 16.1 MW, Q=0.67
- H-mode pedestal: pedestal_pressure, pedestal_width
- Sawtooth oscillations: sawtooth_period, sawtooth_mixing_radius

## L7: Applications (Complete)
- ITER power balance analysis (iter_power_balance)
- Fusion power plant net power: thermal_to_electric_efficiency
- Tritium breeding: tritium_breeding_ratio_needed, lithium_enrichment_for_tbr
- Neutron wall loading and energy multiplication
- Divertor heat flux: divertor_heat_flux
- Levelized cost of electricity: fusion_lcoe
- Capital cost scaling: fusion_capital_cost_busd
- Plant availability and waste disposal rating
- Tritium inventory: tritium_inventory_required, tritium_doubling_time

## L8: Advanced Topics (Complete)
- Burning plasma physics: alpha_power_to_electrons/ions, alpha_ash_confinement
- Neoclassical tearing modes: ntm_threshold_island_width
- Tearing mode stability: tearing_mode_delta_prime
- Ballooning stability: ballooning_alpha
- Mercier criterion: mercier_criterion
- Vertical stability: vertical_stability_growth_rate
- Resistive wall modes: resistive_wall_mode_growth_rate
- Halo currents: halo_current_fraction
- Liquid metal blankets: liquid_metal_tbr_estimate, beryllium_neutron_multiplier
- Structural materials: first_wall_dpa, helium_production
- Zonal flows: zonal_flow_shearing_rate, geodesic_acoustic_mode_frequency
- ITG/TEM/ETG turbulent transport models
- Critical gradient threshold for transport barrier formation

## L9: Research Frontiers (Partial)
- ARC-class compact fusion reactors: arc_magnet_cost
- Stellarator optimization: stellarator_coil_complexity, stellarator_alpha_confinement
- Inertial confinement fusion: icf_gain_from_rhoR, icf_driver_efficiency
- Fusion safety and lifecycle: fusion_safety_factor, fusion_carbon_intensity
- Economic competitiveness: fusion_competitiveness_index