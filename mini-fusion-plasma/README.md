# mini-fusion-plasma — Fusion Plasma Physics

> MIT 22.611/22.615 | Princeton PHY 563 | Cambridge Part III | Oxford CMT | ETH 402-0851

## Module Status: COMPLETE ✅

- **L1 Definitions**: Complete (14 struct/enum types)
- **L2 Core Concepts**: Complete (28 fundamental plasma functions)
- **L3 Math Structures**: Complete (15+ mathematical operations, flux coordinates)
- **L4 Fundamental Laws**: Complete (24 physics laws + Lean theorems)
- **L5 Computational Methods**: Complete (Bosch-Hale, Gauss-Laguerre, transport models)
- **L6 Canonical Systems**: Complete (ITER, DEMO, SPARC, JET, sawtooth)
- **L7 Applications**: Complete (power balance, TBR, NWL, LCOE, divertor)
- **L8 Advanced Topics**: Complete (NTM, RWM, ballooning, ITG/TEM/ETG, zonal flows)
- **L9 Research Frontiers**: Partial (ARC, stellarator, ICF, lifecycle — 6 functions)

**Score: 17/18 — COMPLETE** (threshold: 16)

**Line count**: include/ + src/ = **3354 lines** (threshold: 3000 ✅)

## Core Definitions (L1)

| Type | Description |
|------|-------------|
| `PlasmaParameters` | Plasma state: ne, Te, B, Ip, R, a, kappa, ... |
| `DebyeSphere` | Debye screening: lambda_D, N_D, ln(Lambda) |
| `LawsonCriteria` | Lawson criterion: n, T, tau_E, triple_product |
| `FusionReactionRate` | Fusion reaction rate: n1, n2, sigma_v, P_fus |
| `EnergyBalance` | Power balance: P_fus, P_alpha, P_aux, Q |
| `TokamakGeometry` | Geometry: R, a, kappa, delta, volume |
| `FluxSurface` | Flux coordinates: psi, theta, zeta, q, shear |
| `FusionFuelType` | Fuel enum: DT, DD, DHe3, pB11, CAT_DD, TT |
| `HeatingSystem` | Heating: P_nbi, P_icrf, P_ecrh, P_lh, P_alpha |
| `TritiumBreedingRatio` | TBR, required TBR, doubling time |

## Core Theorems (L4)

### Lawson Criterion (1957)
```
n * T * tau_E >= 3 x 10^21 m^-3 keV s  (D-T ignition)
```
- `lawson_criterion(n, T, tau_E)` — compute triple product
- `ignition_condition(n, T, tau_E)` — evaluate ignition (returns 0/1)
- Lean: `isIgnited`, `ignition_from_triple_product`, `ignition_monotonic_tau`

### Fusion Power Density
```
P_DT = nD * nT * <sigma v>_DT * E_fus
P_alpha = nD * nT * <sigma v>_DT * E_alpha
```

### Bremsstrahlung Radiation (Karzas & Latter, 1961)
```
P_brem = 5.35e-37 * Zeff * ne^2 * sqrt(Te)  [W/m^3]
```

### IPB98(y,2) Confinement Scaling (ITER Physics Basis, 1999)
```
tau_E = 0.0562 * Ip^0.93 * B^0.15 * n^0.41 * P^-0.69 * R^1.97 * eps^0.58 * kappa^0.78
```

### Bosch-Hale D-T Reactivity (1992)
```
<sigma v>_DT = C1 * theta * sqrt(xi/(mR*c^2*Ti^3)) * exp(-3*xi)
xi = (B_G^2/(4*theta))^(1/3), theta = Ti/(1 - Ti*(C2+...)/(1+Ti*(C3+...)))
```

## Core Algorithms (L5)

| Algorithm | Description | Source |
|-----------|-------------|--------|
| Bosch-Hale 5-param fit | D-T/D-D/D-He3 reactivity (0.5-200 keV, <0.25% error) | Nucl. Fusion 32, 611 |
| Gauss-Laguerre quadrature | Maxwellian reactivity integration (4-point) | Numerical Recipes |
| Neoclassical transport | Banana/plateau/Pfirsch-Schluter unified formula | Chang & Hinton (1982) |
| IPB98(y,2) scaling | ELMy H-mode tau_E prediction (RMS ~14%) | NF 39, 2137 |
| Finite-difference GS | Grad-Shafranov operator Delta* on 2D R-Z grid | Numerical |
| ITG/TEM/ETG models | Mixing-length turbulent diffusivity estimates | Various |

## Canonical Systems (L6)

| System | Key Parameters | Fusion Power |
|--------|---------------|-------------|
| **ITER** | R=6.2m, B=5.3T, Ip=15MA | 500 MW, Q=10 |
| **DEMO** | R=9.0m, B=5.7T, Ip=19.6MA | 2000 MW, Q=25-50 |
| **SPARC** | R=1.85m, B=12.2T, Ip=8.7MA | 50-140 MW, Q>=2 |
| **JET** | R=2.96m, B=3.8T, Ip=4MA | 16.1 MW, Q=0.67 |

## Nine-School Course Mapping

| University | Course | Key Topics |
|------------|--------|-----------|
| MIT | 22.611/22.615 | Debye, collisions, MHD, GS equation |
| Princeton | PHY 563 | Transport, neoclassical, tokamak physics |
| Cambridge | Part III Plasma/Fusion | Lawson, ignition, confinement |
| Oxford | CMT Plasma Theory | Grad-Shafranov, bootstrap, shear |
| ETH | 402-0851 | Waves, heating, current drive |
| Stanford | PHYSICS 370 CM | Kinetic theory, Fokker-Planck |
| Berkeley | PHYS 242 CM | MHD stability, transport |
| Caltech | Ph 205 GR / Ph 230 | MHD as fluid/field theory |
| Tokyo | Plasma Physics | Tokamak experiments, ITER |

## Build

```bash
cd mini-fusion-plasma
make          # build static library libfusionplasma.a
make test     # run all tests
make examples # run all examples
make count    # line count statistics
make audit    # safety check (filler/stub/TODO detection)
make clean    # clean build artifacts
```

## File Structure

```
mini-fusion-plasma/
  Makefile
  README.md                          # This file
  include/
    fusion_plasma.h                  # Main header (types, constants, API)
    fusion_confinement.h             # Confinement and transport
    fusion_transport.h               # Transport theory
    fusion_mhd.h                     # MHD equilibrium and stability
    fusion_heating.h                 # Heating and current drive
    fusion_equilibrium.h             # Operating points
  src/
    fusion_core.c                    # Core plasma physics (635 lines)
    fusion_confinement.c             # Confinement/transport (493 lines)
    fusion_mhd.c                     # MHD equilibrium (410 lines)
    fusion_heating.c                 # Heating systems (269 lines)
    fusion_equilibrium.c             # Operating points (344 lines)
    fusion_reactor.c                 # Power plant design (421 lines)
    fusion_cheatsheet.lean           # Lean 4 formalization (179 lines)
  tests/
    test_core.c                      # 25+ assert-based tests
  examples/
    example1_iter_power_balance.c     # ITER power balance
    example2_reactivity_scan.c        # Cross-section temperature scan
    example3_reactor_comparison.c     # ITER/DEMO/SPARC/JET comparison
  docs/
    knowledge-graph.md               # Nine-level knowledge map
    coverage-report.md               # Per-level completeness assessment
    gap-report.md                    # Missing knowledge items + priority
    course-alignment.md              # Nine-school course mapping
    course-tree.md                   # Prerequisite dependency tree
```

## References

- Goldston & Rutherford, *Introduction to Plasma Physics* (1995)
- Wesson, *Tokamaks* 4th ed. (2011)
- Freidberg, *Plasma Physics and Fusion Energy* (2007)
- Bosch & Hale, *Nucl. Fusion* 32, 611 (1992)
- ITER Physics Basis, *Nucl. Fusion* 39, 2137 (1999)
- Sorbom et al., *FED* 100, 378 (2015) — SPARC/ARC
- Maisonnier et al., *FED* 75-79, 1173 (2005) — DEMO
- Stacey, *Fusion: An Introduction* 2nd ed. (2010)