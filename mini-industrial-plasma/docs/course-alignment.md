# Course Alignment - mini-industrial-plasma
## MIT 22.611 (Plasma Science and Fusion Engineering)
| Section | Topic | Implementation |
|---------|-------|----------------|
| Ch.1 | Plasma parameters | plasma_state_init(), debye_length |
| Ch.2 | Single particle motion | electron_larmor_radius(), hall_parameter() |
| Ch.4 | Plasma as fluid | TransportCoeffs, plasma_beta() |
| Ch.6 | Sheath theory | sheath.c (full implementation) |
| Ch.7-8 | Plasma chemistry | plasma_chemistry.c (Ar,SF6,SiH4,O2) |
| Ch.9 | Etching | surface.c (RIE model) |
| Ch.10 | Global model | global_model_electron_temp() |
| Ch.11 | Capacitive discharges | discharge.c (CCP model) |
| Ch.12 | Inductive discharges | icp_skin_depth(), ICP model |
| Ch.14-16 | PECVD/Sputtering | surface.c (deposition, sputtering) |

## Berkeley EECS 245 (Plasma Processing)
| Topic | Implementation |
|-------|----------------|
| Sheath dynamics | solve_planar_sheath(), sheath_impedance() |
| RIE lag / ARDE | rie_lag_etch_rate() |
| PECVD a-Si:H | pecvd_deposition_rate(), pecvd_hydrogen_content() |
| Wafer uniformity | wafer_etch_uniformity() |

## Stanford EE 414 (Plasma Processing)
| Topic | Implementation |
|-------|----------------|
| RF discharges | CCP self-bias, global model |
| Boltzmann solver | solve_two_term_boltzmann() |
| Surface kinetics | surface_coverage_langmuir() |

## Corresponding Textbooks
- Lieberman & Lichtenberg, Principles of Plasma Discharges (2005)
- Chabert & Braithwaite, Physics of RF Plasmas (2011)
- Goldston & Rutherford, Plasma Physics (1995)
- von Engel, Ionized Gases (1965)
