/**
 * fusion_heating.h — Plasma Heating and Current Drive Header
 *
 * Refs: Stix "Waves in Plasmas" (1992),
 *       ITER Physics Basis (1999), Fisch RMP (1987)
 *
 * L2: NBI, ICRF, ECRH, LH heating
 * L4: Alpha heating, ohmic heating
 * L5: NBI deposition models, RF propagation
 * L7: ITER/DEMO heating mix, cost models
 */

#ifndef MINI_FUSION_HEATING_H
#define MINI_FUSION_HEATING_H

#include "fusion_plasma.h"

/* L2: NBI */
double nbi_stopping_cross_section(double E_beam_keV, double ne, double lnLambda);
double nbi_shine_through_fraction(double a, double ne, double sigma_stop);
double nbi_orbit_loss_fraction(double a, double rho_fast, double tangency_radius);
double nbi_current_drive_efficiency(double ne_20, double R, double Z_eff);

/* L2: RF */
double icrf_resonance_frequency(double B, double mi_kg, double Z);
double ecrh_resonance_frequency(double B);
double lh_resonance_frequency(double B, double mi_kg, double Z, double ne);
double rf_power_deposition(double r, double r_res, double delta_r,
                             double P_total, double a);

/* L4: Alpha */
double alpha_power_to_electrons(double P_alpha, double Te_eV);
double alpha_power_to_ions(double P_alpha, double Te_eV);
double alpha_ash_confinement_time(double tau_E, double S_He);

/* L7: System Integration */
double heating_system_cost(double P_nbi, double P_ecrh, double P_icrf, double P_lh);
void   optimum_heating_mix(HeatingSystem *hs);

#endif /* MINI_FUSION_HEATING_H */