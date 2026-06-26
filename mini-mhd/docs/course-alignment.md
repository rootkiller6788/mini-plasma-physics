# Course Alignment — mini-mhd

## MIT 22.611 — MHD Theory (Freidberg)

| Chapter | Topic | Implementation |
|---------|-------|---------------|
| 1 | MHD equations | mhd_eqns.c (ideal + resistive residuals) |
| 2 | MHD equilibrium | mhd_equilibrium.c (pinches, Soloviev) |
| 3 | MHD energy principle | mhd_instability.c (energy_principle_fluid) |
| 4 | MHD waves | mhd_waves.c (dispersion, Friedrichs) |
| 5 | Kink/sausage/interchange | mhd_instability.c (growth rates) |
| 6 | Ballooning modes | mhd_instability.c (ballooning_parameter) |
| 7 | Resistive instabilities | mhd_instability.c (tearing_growth_rate) |

## Princeton AST 554 — MHD in Astrophysics

| Topic | Implementation |
|-------|---------------|
| MHD waves | mhd_waves.c |
| Virial theorem | (implicit in equilibrium) |
| Magnetic helicity | mhd_eqns.c, mhd.lean |
| Dynamo action | mhd_defs.h (Elsasser number) |

## Cambridge Part III — MHD

| Topic | Implementation |
|-------|---------------|
| Alfven's theorem | mhd_eqns.c (alfven_theorem_check) |
| Force-free fields | mhd_equilibrium.c (beltrami_check) |
| Grad-Shafranov | mhd_equilibrium.c (grad_shafranov_operator) |
| Magnetic reconnection | (Sweet-Parker scaling in tearing mode) |

## Berkeley PHYS 242 — Classical Mechanics (MHD as continuum)

| Topic | Implementation |
|-------|---------------|
| Lagrangian for MHD | (documented in knowledge-graph) |
| Noether invariants | mhd_eqns.c (helicities, energy) |
| Hamiltonian structure | (documented) |

## Caltech Ph 135 — Electromagnetism

| Topic | Implementation |
|-------|---------------|
| Maxwell stress tensor | mhd_defs.h (magnetic pressure/tension) |
| Ohm's law for moving media | mhd_eqns.c (mhd_ohm_law) |
| Quasistatic approximation | mhd_defs.h (Ampere without displacement) |

## ETH 402-0891 — General Relativity

| Topic | Implementation |
|-------|---------------|
| Relativistic MHD | (documented in L9) |

## Oxford CMT — Condensed Matter / Plasma Theory

| Topic | Implementation |
|-------|---------------|
| Taylor relaxation | mhd_equilibrium.c (taylor_relaxation) |
| Magnetic helicity | mhd.lean (helicity_conservation_ideal) |
| Beltrami fields | mhd_equilibrium.c (bessel_J0/J1) |

## 东京大学 — Plasma Physics

| Topic | Implementation |
|-------|---------------|
| Tokamak equilibrium | mhd_equilibrium.c (Soloviev), tokamak.c |
| Kruskal-Shafranov | mhd_instability.c, tokamak.c |
| Bennett relation | mhd_equilibrium.c, zpinch.c |

## Cross-Cutting Topics (All Schools)

| Topic | Appears In |
|-------|-----------|
| Plasma beta | All plasma courses |
| Alfven waves | MIT, Princeton, Cambridge, Tokyo |
| Grad-Shafranov | MIT, Princeton, Tokyo |
| Magnetic reconnection | Princeton, Cambridge, Tokyo |
| Kink instability | MIT, Tokyo |
| Taylor relaxation | Oxford, Princeton |
