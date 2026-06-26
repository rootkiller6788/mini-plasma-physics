#ifndef EEDF_H
#define EEDF_H

#include "plasma_types.h"
#include <stddef.h>

/*
 * eedf.h — Electron Energy Distribution Functions
 *
 * Reference: Lieberman §5.3, Hagelaar & Pitchford Plasma Sources
 *   Sci. Technol. 14, 722 (2005) (BOLSIG+)
 * Course: MIT 22.611, Stanford EE 414
 */

/* L1: Maxwellian EEDF — f(epsilon) = 2*sqrt(epsilon/pi)*(T_e^{-3/2})*exp(-epsilon/T_e) */
double maxwellian_eedf(double energy_eV, double T_e_eV);

double maxwellian_eedf_cumulative(double E_max, double T_e_eV);

/* Effective electron temperature from mean energy:
 * <epsilon> = 3/2 * T_e for Maxwellian */
double mean_energy_to_Te(double mean_energy_eV);

/* L1: Druyvesteyn EEDF — for high E/N, elastic collisions dominant
 * f(epsilon) = 0.565 * <epsilon>^{-3/2} * epsilon^{1/2} * exp(-0.243*(epsilon/<epsilon>)^2) */
double druyvesteyn_eedf(double energy_eV, double mean_energy_eV);

/* L2: Bi-Maxwellian EEDF — bulk + tail electrons */
double bimaxwellian_eedf(double energy_eV, double T_bulk, double T_tail,
                         double n_tail_frac);

/* L3: Generalized EEDF with shape parameter x
 * f(eps) = c1 * eps^{1/2} * exp(-c2 * (eps/T_eff)^x)
 * x=1 -> Maxwellian, x=2 -> Druyvesteyn */
double generalized_eedf(double energy_eV, double T_eff, double x);

/* L3: Normalization constant for generalized EEDF */
double generalized_eedf_norm(double T_eff, double x);

/*
 * L4: Boltzmann equation — rate coefficient calculation
 * k(T_e) = sqrt(2*e/m_e) * integral_0^inf sigma(eps) * f(eps) * eps * d(eps)
 */
double rate_coefficient_maxwellian(CrossSectionModel *cs, double T_e_eV);

double rate_coefficient_druyvesteyn(CrossSectionModel *cs, double mean_E);

/*
 * L5: Two-term Boltzmann solver — EEDF from cross sections
 * Solves the steady-state, spatially homogeneous Boltzmann equation
 * using the two-term spherical harmonic expansion.
 * Reference: Hagelaar & Pitchford (2005), Morgan & Penetrante (1990)
 */
int solve_two_term_boltzmann(EEDFState *state, CrossSectionModel **cs_list,
                              int n_cs, double *gas_fractions,
                              double E_over_N_Td, int max_iter);

/*
 * L5: Calculate transport coefficients from EEDF
 * Mobility: mu = -(e/3)*sqrt(2/m_e)*integral (eps/sigma_m)*(df/deps)*d(eps)
 * Diffusion: D = (1/3)*sqrt(2/m_e)*integral (eps/sigma_m)*f*d(eps)
 */
void compute_transport_from_eedf(EEDFState *state,
                                  CrossSectionModel *sigma_m,
                                  TransportCoeffs *coeffs);

/* L5: Ionization rate coefficient from tabulated cross sections */
double ionization_rate_bolsig(double E_over_N, double ionization_energy,
                              double *table, int table_len);

/*
 * L6: EEDF for specific gas mixtures relevant to semiconductor processing
 */
/* Argon EEDF (widely used in sputtering and ICP) */
void argon_default_cross_sections(CrossSectionModel *cs, int *n_cs);

/* CF4 EEDF (fluorocarbon etching gas) */
void cf4_default_cross_sections(CrossSectionModel *cs, int *n_cs);

/* O2 EEDF (ashing and surface activation) */
void o2_default_cross_sections(CrossSectionModel *cs, int *n_cs);

/* L8: Non-local EEDF — kinetic effects in low-pressure ICP */
double nonlocal_eedf_estimate(double epsilon, double T_e_bulk,
                              double potential_drop, double lambda_e);

#endif /* EEDF_H */
