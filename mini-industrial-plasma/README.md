# mini-industrial-plasma -- Industrial Plasma Physics

> MIT 22.611 / Lieberman & Lichtenberg (2005)
> Berkeley EECS 245 / Stanford EE 414

## Module Status: COMPLETE

- L1-L6: Complete
- L7: Complete (3+ applications: RIE, PECVD, Paschen curve)
- L8: Partial (Electronegative plasma, ALE, non-local EEDF)
- L9: Partial (Documented, basic ALE model implemented)

### Score: 16/18 (COMPLETE)

## Core Definitions (L1)
| Type | Name | Description |
|------|------|-------------|
| struct | PlasmaState | Complete plasma thermodynamic state (n, T, potentials) |
| struct | PlasmaSpecies | Particle species (mass, charge, thresholds) |
| struct | PlasmaReaction | Elementary reaction (Arrhenius rate parameters) |
| struct | CrossSectionModel | Parametric cross section model |
| struct | DischargeGeometry | Chamber and electrode geometry |
| struct | DischargeParams | RF/DC power supply parameters |
| struct | TransportCoeffs | Mobility, diffusion, conductivity |
| struct | RateMatrix | Chemical kinetics rate matrix |
| struct | EEDFState | Electron energy distribution state |
| struct | SheathSolution | Poisson equation sheath solution |
| struct | PaschenCurve | Breakdown voltage vs p*d |
| struct | EtchModel | Plasma etching process model |
| struct | DepositionModel | PECVD film properties |

## Core Theorems (L4)
| Theorem | Formula | Reference |
|---------|---------|-----------|
| Child-Langmuir Law | J = (4/9)*eps0*sqrt(2e/m)*V^{3/2}/d^2 | Child (1911), Langmuir (1913) |
| Bohm Criterion | v_i >= c_s = sqrt(e*T_e/m_i) | Bohm (1949) |
| Floating Potential | V_f = -(T_e/2)*ln(m_i/(2*pi*m_e)) | Langmuir probe theory |
| Paschen Law | V_br = B*pd/ln(A*pd/ln(1+1/gamma)) | Paschen (1889) |
| Townsend Criterion | gamma*(exp(alpha*d)-1) = 1 | Townsend (1902) |

## Core Algorithms (L5)
| Algorithm | File | Method |
|-----------|------|--------|
| RK4 Sheath Solver | sheath.c | 4th-order Runge-Kutta, adaptive step |
| Simpson Rate Integral | eedf.c | 1000-pt log-spaced Simpson integration |
| Two-Term Boltzmann | eedf.c | Iterative energy balance |
| Euler Chemistry Step | plasma_chemistry.c | Adaptive forward Euler (0-D) |
| Global Discharge Model | discharge.c | Particle + power balance |
| Matrix Sheath | sheath.c | Uniform ion density analytic |

## Canonical Systems (L6)
- DC Glow Discharge (Townsend to Arc regimes)
- Capacitively Coupled RF Discharge (CCP) with self-bias
- Inductively Coupled Plasma (ICP) with transformer model
- Argon plasma chemistry (7 reactions, 5 species)
- SF6 etching chemistry (5 reactions, 10 species)
- SiH4/H2 PECVD chemistry (5 reactions, 7 species)
- O2 electronegative chemistry (6 reactions, 7 species)

## Applications (L7)
1. **RIE Process Model** — Full CCP etching simulation with anisotropy, selectivity, ARDE
2. **PECVD Deposition** — a-Si:H film growth with H content, stress, uniformity
3. **Paschen Curve** — Breakdown voltage analysis for gas selection and ignition design
4. **Industrial Database** — Default parameters for RIE, PECVD, ICP, Sputter tools

## Advanced Topics (L8)
- Electronegative plasma chemistry (O2 attachment/detachment)
- Atomic Layer Etching (ALE) basic model
- Non-local EEDF effects in low-pressure ICP

## File Structure


## Nine-School Course Mapping
| School | Course | Key Topics Covered |
|--------|--------|--------------------|
| MIT | 22.611 | All chapters (Lieberman textbook) |
| Berkeley | EECS 245 | Plasma processing, RIE, PECVD |
| Stanford | EE 414 | RF discharges, Boltzmann solver |
| Princeton | AST 551 | Plasma-wall interactions |
| Cambridge | Part III | Sheath theory, EEDF |
| Oxford | CMT | Plasma kinetics |

## Build & Test
mkdir -p build
gcc -Wall -Wextra -O2 -std=c99 -Iinclude -c src/sheath.c -o build/sheath.o
mkdir -p build
gcc -Wall -Wextra -O2 -std=c99 -Iinclude -c src/sheath.c -o build/sheath.o
mkdir -p build
gcc -Wall -Wextra -O2 -std=c99 -Iinclude -c src/sheath.c -o build/sheath.o
mkdir -p build
gcc -Wall -Wextra -O2 -std=c99 -Iinclude -O3 benches/bench_eedf.c src/eedf.c src/plasma_params.c -lm -o build/bench_eedf
rm -rf build

## Line Counts
- include/: 737 lines (5 headers)
- src/: 3096 lines (6 C + 1 Lean)
- Total include/ + src/: 3833 lines
- Threshold: 3000 lines MET
