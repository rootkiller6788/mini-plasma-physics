/**
 * fusion_confinement.h — Magnetic Confinement Physics
 *
 * Refs: Wesson (2011), Freidberg (2007),
 *       ITER Physics Basis (1999)
 *
 * L2: Bootstrap current, Pfirsch-Schluter factor
 * L4: Troyon/Greenwald limits, Kruskal-Shafranov
 * L5: Neoclassical transport regimes, Bohm/Gyro-Bohm
 * L6: H-mode pedestal, beta limits
 * L8: Tearing mode stability, ballooning, Mercier, disruptions
 */

#ifndef MINI_FUSION_CONFINEMENT_H
#define MINI_FUSION_CONFINEMENT_H

#include "fusion_plasma.h"

/* L4: Beta Limits */
double troyon_beta_limit(double beta_percent, double a, double B, double Ip_MA);
double normalized_beta(double beta_percent, double a, double B, double Ip_MA);
int    kruskal_shafranov_limit(double q_edge);

/* L2: Bootstrap Current */
double bootstrap_current_fraction(double epsilon, double beta_p);
double bootstrap_current_density(double epsilon, double B_pol, double dp_dr);
double pfirsch_schluter_factor(double q);

/* L5: Neoclassical Transport */
double banana_regime_diffusivity(double epsilon, double q, double rho_i, double nu_ii);
double plateau_regime_diffusivity(double rho_i, double v_thi, double R, double q);
double pfirsch_schluter_diffusivity(double q, double rho_i, double nu_ii);
double neoclassical_ion_heat_diffusivity(double epsilon, double q, double rho_i,
                                          double nu_ii, double v_thi, double R);

/* L5: Anomalous Transport */
double bohm_diffusion_coefficient(double Te_eV, double B);
double gyrobohm_diffusion_coefficient(double Te_eV, double B, double rho_i, double a);

/* L6: H-mode Physics */
double h_mode_power_threshold(double n_20, double B, double S);
double pedestal_pressure(double Ip_MA, double a, double B);
double pedestal_width(double beta_p_ped, double a);

/* L8: MHD Stability */
double tearing_mode_delta_prime(double psi_prime_plus, double psi_prime_minus, double psi_resonant);
double ballooning_alpha(double q, double R, double dp_dr, double B);
double mercier_criterion(double q, double shear, double magnetic_well, double pressure_gradient);
double magnetic_shear(double r, double q, double dq_dr);
double magnetic_well_depth(double V_prime_0, double V_prime_psi);
double sawtooth_period(double a, double eta, double v_A);
double sawtooth_mixing_radius(double r_q1);

/* L8: Disruptions */
double disruption_density_limit(double Ip, double a, double margin);
double ntm_threshold_island_width(double rho_i, double beta_p);
double halo_current_fraction(double I_halo, double I_p0);

/* L7: Divertor */
double divertor_heat_flux(double P_sol, double alpha_deg, double R, double lambda_q, double f_exp);
double power_exhaust_fraction(double P_rad, double P_heat);

#endif /* MINI_FUSION_CONFINEMENT_H */