# Coverage Report — mini-space-plasma

## Assessment Date: 2026-06-20

| Level | Name | Rating | Score | Notes |
|-------|------|--------|-------|-------|
| L1 | Definitions | **Complete** | 2 | 17+ struct/typedef/enum definitions, all CODATA constants |
| L2 | Core Concepts | **Complete** | 2 | 14+ core concepts implemented (Debye, Alfven, beta, MHD waves, etc.) |
| L3 | Math Structures | **Complete** | 2 | Grad-Shafranov operator, Stix tensor, CMA diagram, MHD stress tensor |
| L4 | Fundamental Laws | **Complete** | 2 | Parker equation, Chapman-Ferraro, dipole field, R-H shocks, Bennett relation |
| L5 | Compute Methods | **Complete** | 2 | Lax-Friedrichs MHD, Newton-Raphson dispersion, CFL, FD operators |
| L6 | Canonical Systems | **Complete** | 2 | Solar wind, magnetosphere, MHD equilibria, wave modes, particle drifts |
| L7 | Applications | **Complete** | 2 | 6 applications (Parker profile, magnetopause, CPCP, Dungey, Dst, CMA) |
| L8 | Advanced Topics | **Complete** | 2 | Kinetic Alfven, warm plasma Z-function, wave polarization, Tsyganenko model |
| L9 | Research Frontiers | **Partial** | 1 | CME parameters, HCS tilt, bi-Maxwellian energy (documented, no dedicated sim) |

**Total Score: 17/18 → COMPLETE**

## Detailed Coverage

### L1 Complete: 17+ struct/enum/constant groups
- Fundamental constants: 17 (SP_KB through SP_MU_EARTH)
- Struct types: plasma_state_t, particle_orbit_t, mhd_wave_t, debye_params_t, magnetosphere_t, solar_wind_t, reconnection_t, distribution_t, velocity_coord_t, dispersion_t
- MHD types: mhd_primitive_t, mhd_conserved_t, mhd_flux_t, mhd_eigensystem_t, mhd_grid_t
- Stix tensor: stix_tensor_t, wave_polarization_t
- Enums: wave_mode_t (12 values), resistivity_model_t, plasma_region_t (12 values), parker_solution_type_t (5 values)

### L2 Complete: 14+ concepts
Debye screening, plasma frequency, gyrofrequency, Alfven speed, sound speed, plasma beta, collision frequency, MHD prim/cons conversion, MHD flux, MHD eigenspeeds, frozen-in flux, ExB convection, ring current, plasmasphere, Mach numbers

### L3 Complete: 5 mathematical structures
MHD equilibria (force-free, Harris, Z-pinch), Grad-Shafranov operator, Stix cold plasma tensor, Appleton-Hartree dispersion, CMA diagram

### L4 Complete: 10 fundamental relations
Parker wind equation, mass continuity, Parker spiral, Chapman-Ferraro, dipole B-field, R-H shock jumps, warm plasma dielectric (Fried-Conte/Bohm-Gross), Bennett relation, Dessler-Parker-Sckopke, cold plasma dispersion

### L5 Complete: 7 numerical methods
Lax-Friedrichs MHD, CFL condition, Newton-Raphson dispersion solver, group velocity FD, Grad-Shafranov FD, div(B) monitoring, Parker equation Newton solve

### L6 Complete: 19 canonical problems
Parker wind profile (classes I-V), Parker spiral, CIR formation, Earth dipole, magnetopause, bow shock, Harris sheet, Z-pinch, force-free fields, particle drifts (grad-B, curvature, total), bounce period, magnetosheath flow, Tsyganenko model, plasmasphere profile, Langmuir/Bohm-Gross, ion acoustic, Alfven (shear+kinetic), whistler, fast/slow magnetosonic, lower/upper hybrid

### L7 Complete: 6 applications
Parker profile computation (corona→1AU), magnetopause standoff (quiet/fast/CME), cross-polar cap potential (Boyle 1997), Dungey cycle voltage, ring current Dst prediction, CMA propagation bands

### L8 Complete: 4 advanced topics
Kinetic Alfven waves (FLR), warm plasma dielectric (Landau damping), wave polarization (helicity), Tsyganenko stretched field

### L9 Partial: 3 research areas documented
CME magnetic cloud parameters, HCS tilt evolution, bi-Maxwellian ring current anisotropy

## Statistics
- include/*.h: 1920 lines (6 files)
- src/*.c: 2103 lines (5 files)
- Total (include+src): 4023 lines >= 3000 OK
- tests/test_plasma.c: 41 tests, 100% pass rate
- examples/: 3 end-to-end examples
