# Course Alignment — mini-space-plasma

## Reference Textbooks
- **Goldston & Rutherford** "Introduction to Plasma Physics" (1995) — Primary
- **Kivelson & Russell** "Introduction to Space Physics" (1995) — Primary
- **Chen** "Introduction to Plasma Physics and Controlled Fusion" (1984)
- **Priest** "Solar Magnetohydrodynamics" (1982)
- **Stix** "Waves in Plasmas" (1992)
- **Swanson** "Plasma Waves" (2003)

## Nine-School Course Mapping

| School | Course | Chapters Covered | Our Implementation |
|--------|--------|-----------------|-------------------|
| **MIT** | 22.611 Space Plasma Physics | L1-L22 (full) | plasma_parameters.c, solar_wind.c, magnetosphere.c, plasma_waves.c |
| **MIT** | 22.615 MHD Theory | L1-L20 (full) | mhd_core.c (equilibria, waves, Lax-Friedrichs) |
| **Stanford** | PHYSICS 370 Plasma Physics | Debye, MHD, waves | All modules |
| **Berkeley** | PHYS 242 Plasma Physics | Kinetic theory, MHD | plasma_waves.c, mhd_core.c |
| **Caltech** | Ph 106 Plasma Physics | Fundamentals | plasma_parameters.c |
| **Princeton** | AST 553 Plasma Astrophysics | MHD, solar wind | mhd_core.c, solar_wind.c |
| **Cambridge** | Part III Plasma Physics | Waves, MHD, space | plasma_waves.c, magnetosphere.c |
| **Oxford** | CMT Plasma | MHD theory | mhd_core.c (Grad-Shafranov, equilibria) |
| **ETH** | 402-0891 Plasma Physics | Space plasma | all modules |
| **Tokyo** | 宇宙プラズマ物理学 | Space plasma | solar_wind.c, magnetosphere.c |

## Chapter-Level Mapping

### Goldston & Rutherford
- Ch.2: Single Particle Motion → plasma_parameters.c (gyrofreq, Larmor radius)
- Ch.3: Plasma as a Fluid → mhd_core.c (prim/cons, flux)
- Ch.4: Waves in Cold Plasma → plasma_waves.c (Stix, dispersion)
- Ch.6: Collisions → plasma_parameters.c (Coulomb, Spitzer)
- Ch.9: MHD Equilibrium → mhd_core.c (Harris, Z-pinch, GS)

### Kivelson & Russell
- Ch.2: Single Particle Motion → magnetosphere.c (gradB/curvature drift)
- Ch.3: Trapped Particles → magnetosphere.c (bounce, drift, L-shell)
- Ch.4: MHD → mhd_core.c
- Ch.5: Plasma Waves → plasma_waves.c
- Ch.6: Solar Wind → solar_wind.c (Parker, spiral, CIR)
- Ch.7-8: Magnetosphere → magnetosphere.c (dipole, CF, plasmasphere)
- Ch.9: Magnetospheric Dynamics → magnetosphere.c (ExB, convection)
