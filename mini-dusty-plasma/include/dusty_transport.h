#ifndef DUSTY_TRANSPORT_H
#define DUSTY_TRANSPORT_H
/**
 * @file  dusty_transport.h
 * @brief Transport coefficients in dusty plasmas.
 *
 * Transport processes govern the macroscopic behavior: collisions with
 * neutrals and ions determine dust mobility, diffusion, and thermal
 * relaxation. In strongly coupled systems, non-Newtonian viscosity
 * and anomalous diffusion can arise.
 *
 * Refs:
 *   Bouchoule (1999), "Dusty Plasmas" (Wiley)
 *   Fortov et al. (2005), Phys. Rep. 421, 1 — review
 *   Khrapak & Morfill (2009), Contrib. Plasma Phys. 49, 148
 */
#include "dusty_plasma.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================
 * L2 — Collision Frequencies
 * ================================================================ */

/**
 * @brief Dust-neutral collision frequency (Epstein drag).
 *
 * For a spherical grain in the free-molecular regime
 * (Kn = lambda_mfp / a >> 1):
 *
 * nu_dn = (8/3) * sqrt(2/pi) * a^2 * n_n * v_thn * (m_n / m_d)
 *
 * where v_thn = sqrt(k_B * T_n / m_n).
 *
 * For the specular reflection case, the prefactor is (8/3)*sqrt(2/pi).
 * For diffuse reflection with accommodation coefficient alpha:
 * nu_dn = nu_dn_specular * (1 + pi * alpha / 8)
 *
 * Ref: Epstein (1924), Phys. Rev. 23, 710
 */
double dust_neutral_collision_freq(
    double a, double n_n, double T_n, double m_n, double m_d,
    double alpha_accommodation);

/**
 * @brief Dust-ion collision frequency (Coulomb scattering).
 *
 * nu_di = (4/3) * sqrt(2*pi) * n_i * a^2 * v_thi
 *         * (Z_d * e^2 / (4*pi*eps0 * a * k_B * T_i))^2
 *         * ln(Lambda)
 *
 * where ln(Lambda) is the Coulomb logarithm for dust-ion collisions:
 * ln(Lambda) ≈ ln(lambda_D / r_min) with r_min = max(a, Z_d*e^2/(4*pi*eps0*k_B*T_i)).
 *
 * Complexity: O(1).
 */
double dust_ion_collision_freq(
    double a, double n_i, double T_i, double m_i, double Z_d);

/**
 * @brief Coulomb logarithm for dust-ion collisions.
 *
 * ln Lambda_di = ln( lambda_D / max(a, b_90) )
 *
 * where b_90 = Z_d * e^2 / (4*pi*eps0 * 3 * k_B * T_i) is the
 * 90-degree scattering impact parameter.
 */
double dust_ion_coulomb_logarithm(
    double a, double Z_d, double T_i, double lambda_D);

/* ================================================================
 * L3 — Diffusion and Mobility
 * ================================================================ */

/**
 * @brief Dust diffusion coefficient (Einstein relation).
 *
 * D_d = k_B * T_d / (m_d * nu_dn)
 *
 * Valid in the weakly-coupled, collision-dominated regime.
 * In strongly coupled plasmas, D_d is modified by caging effects.
 *
 * Ref: Einstein (1905), Ann. Phys. 17, 549
 */
double dust_diffusion_coefficient(double T_d, double m_d, double nu_dn);

/**
 * @brief Dust electrical mobility.
 *
 * mu_d = e * Z_d / (m_d * nu_tot)
 *
 * where nu_tot = nu_dn + nu_di is the total momentum transfer frequency.
 * Mobility relates drift velocity to electric field: v_drift = mu_d * E.
 */
double dust_mobility(double Z_d, double m_d, double nu_tot);

/**
 * @brief Dust drift velocity in electric field E.
 *
 * v_drift = mu_d * E
 *
 * Steady-state solution of: m_d * dv/dt = Q_d * E - m_d * nu * v
 */
double dust_drift_velocity(double mu_d, double E);

/**
 * @brief Mean free path for dust grains.
 *
 * lambda_mfp_d = v_thd / nu_tot
 */
double dust_mean_free_path(double v_thd, double nu_tot);

/* ================================================================
 * L4 — Dust Temperature Equation
 * ================================================================ */

/**
 * @brief Dust heating/cooling power balance.
 *
 * dT_d/dt = (P_heat - P_cool) / ( (3/2) * k_B )
 *
 * Heating mechanisms:
 *   - Ion bombardment heating
 *   - Electric field Joule heating
 *   - Absorption of radiation
 *
 * Cooling mechanisms:
 *   - Neutral gas cooling (dominant)
 *   - Thermal radiation
 *   - Evaporative cooling (for volatile materials)
 */
double dust_heating_rate_ion_bombardment(
    double a, double n_i, double T_i, double m_i,
    double phi_f, double T_d);

/**
 * @brief Dust cooling rate by neutral gas.
 *
 * P_cool = (3/2) * k_B * (T_d - T_n) * nu_dn
 *         * 2 * m_d * m_n / (m_d + m_n)^2
 *
 * The last factor is the energy transfer efficiency per collision.
 */
double dust_cooling_rate_neutral(
    double T_d, double T_n, double nu_dn, double m_d, double m_n);

/**
 * @brief Radiative cooling of a dust grain (Stefan-Boltzmann).
 *
 * P_rad = 4 * pi * a^2 * epsilon * sigma_SB * (T_d^4 - T_env^4)
 *
 * where epsilon is the grain emissivity and sigma_SB is Stefan-Boltzmann.
 * This is usually negligible except at very high temperatures.
 */
double dust_radiative_cooling_power(
    double a, double epsilon, double T_d, double T_env);

/**
 * @brief Compute steady-state dust temperature.
 *
 * Solves P_heat(T_d) = P_cool(T_d) for T_d using bisection method.
 *
 * Complexity: O(log((T_max - T_min) / tol)).
 */
double dust_steady_state_temperature(
    double a, double n_i, double T_i, double m_i,
    double phi_f, double T_n, double nu_dn, double m_d, double m_n,
    double T_min, double T_max, double tol);

/* ================================================================
 * L5 — Non-local Transport and Memory Effects
 * ================================================================ */

/**
 * @brief Generalized (fractional) diffusion coefficient.
 *
 * In strongly coupled dusty plasmas, diffusion can be sub-diffusive
 * (alpha < 1) due to caging, or super-diffusive (alpha > 1) due to
 * long-range Yukawa interactions.
 *
 * <r^2(t)> ∝ t^alpha
 *
 * @param alpha  Anomalous diffusion exponent (alpha = 1 = normal Fickian)
 * @param D_normal  Normal diffusion coefficient
 * @param t     Time
 * @return Effective MSD at time t
 */
double dust_anomalous_diffusion_msd(
    double alpha, double D_normal, double t);

/**
 * @brief Dust viscosity in the strongly coupled regime.
 *
 * For Yukawa fluids, the shear viscosity eta is given by:
 *
 * eta = eta_0 * (Gamma / Gamma_crit)^beta
 *
 * where eta_0 = n_d * m_d * a^2 * omega_pd is a characteristic scale.
 *
 * Ref: Saigo & Hamaguchi (2002), Phys. Plasmas 9, 1210
 */
double dust_shear_viscosity(
    double n_d, double m_d, double a, double omega_pd,
    double Gamma, double Gamma_crit, double beta);

/**
 * @brief Dust thermal conductivity in the fluid regime.
 *
 * kappa_thermal = (5/2) * n_d * k_B * D_d
 *
 * Using kinetic theory (monatomic gas approximation).
 */
double dust_thermal_conductivity(
    double n_d, double D_d);

#ifdef __cplusplus
}
#endif

#endif /* DUSTY_TRANSPORT_H */
