# Course Dependency Tree — mini-space-plasma

## Prerequisites
```
General Physics (mechanics, EM, thermodynamics)
  └─ Classical Electrodynamics (Maxwell equations, vector calculus)
       └─ Plasma Physics Fundamentals
            ├── Single Particle Motion (Larmor, drifts)
            ├── Debye Shielding & Plasma Parameter
            ├── MHD Theory (fluid plasma)
            │    ├── Ideal MHD (flux freezing)
            │    ├── MHD Waves (Alfven, magnetosonic)
            │    └── MHD Equilibria (force-free, Grad-Shafranov)
            ├── Kinetic Theory (Vlasov, Landau damping)
            ├── Cold Plasma Waves (Stix, Appleton-Hartree)
            └── Collisions & Transport (Spitzer resistivity)
                 │
                 └─ Space Plasma Physics
                      ├── Solar Wind (Parker model, spiral)
                      │    └── Heliosphere (CIR, HCS, termination shock)
                      ├── Planetary Magnetospheres
                      │    ├── Dipole Field & Field Lines
                      │    ├── Magnetopause (Chapman-Ferraro)
                      │    ├── Bow Shock (Rankine-Hugoniot)
                      │    ├── Magnetotail (current sheet)
                      │    ├── Plasmasphere & Plasmapause
                      │    ├── Ring Current (Dessler-Parker-Sckopke)
                      │    └── Radiation Belts (particle trapping)
                      └── Reconnection & Energy Transfer
```

## Module Internal Dependencies
```
space_plasma.h (base constants + structs)
  ├─ plasma_parameters.h/.c (no deps beyond space_plasma.h)
  ├─ mhd_core.h/.c (depends on space_plasma.h)
  ├─ plasma_waves.h/.c (depends on space_plasma.h + plasma_parameters.h)
  ├─ solar_wind.h/.c (depends on space_plasma.h + plasma_parameters.h)
  └─ magnetosphere.h/.c (depends on space_plasma.h + plasma_parameters.h)
```

## Learning Path
1. Start with `plasma_parameters.c` — fundamental quantities (Debye, beta, collisions)
2. Move to `mhd_core.c` — fluid description of plasma
3. Study `plasma_waves.c` — wave phenomena in magnetized plasma
4. Apply to `solar_wind.c` — Parker's model and heliosphere
5. Apply to `magnetosphere.c` — Earth's magnetic environment
