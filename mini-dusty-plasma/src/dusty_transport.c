/**
 * @file  dusty_transport.c
 * @brief Transport coefficients in dusty plasmas.
 *
 * Transport processes govern macroscopic dust behavior: collisions with
 * neutrals and ions determine dust mobility, diffusion, and thermal
 * relaxation. In strongly coupled systems, non-Newtonian viscosity
 * and anomalous diffusion can arise.
 *
 * L3-L5: Collision frequencies, diffusion, mobility, thermal balance.
 *
 * References:
 *   Bouchoule (1999), "Dusty Plasmas", Wiley
 *   Fortov et al. (2005), Phys. Rep. 421, 1
 *   Khrapak & Morfill (2009), Contrib. Plasma Phys. 49, 148
 */

#include "dusty_transport.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* ================================================================
 * L2 — Collision Frequencies
 * ================================================================ */

double dust_neutral_collision_freq(
    double a, double n_n, double T_n, double m_n, double m_d,
    double alpha_accommodation)
{
    /* Epstein drag for spherical grain in free-molecular regime (Kn >> 1).
     * nu_dn = (8/3)*sqrt(2/pi)*a^2*n_n*v_thn*(m_n/m_d)*(1 + pi*alpha/8)
     * Epstein (1924), Phys. Rev. 23, 710 */
    if (a <= 0.0 || n_n <= 0.0 || T_n <= 0.0 || m_n <= 0.0 || m_d <= 0.0)
        return 0.0;
    double v_thn = sqrt(DUSTY_KB * T_n / m_n);
    double nu_0 = (8.0/3.0) * sqrt(2.0/M_PI) * a*a * n_n * v_thn * (m_n/m_d);
    double correction = 1.0 + M_PI * alpha_accommodation / 8.0;
    return nu_0 * correction;
}

double dust_ion_collision_freq(
    double a, double n_i, double T_i, double m_i, double Z_d)
{
    /* Dust-ion Coulomb collision frequency.
     * nu_di = (4/3)*sqrt(2*pi)*n_i*a^2*v_thi
     *         * (Z_d*e^2/(4*pi*eps0*a*k_B*T_i))^2 * ln(Lambda)
     * Fortov et al. (2005), Phys. Rep. 421, Eq. (3.17) */
    if (a <= 0.0 || n_i <= 0.0 || T_i <= 0.0 || m_i <= 0.0) return 0.0;
    double v_thi = sqrt(DUSTY_KB * T_i / m_i);
    double b90_over_a = (Z_d * DUSTY_EC * DUSTY_EC)
                        / (4.0 * M_PI * DUSTY_EPS0 * a * DUSTY_KB * T_i);
    double ln_Lambda = dust_ion_coulomb_logarithm(a, Z_d, T_i,
                        sqrt(DUSTY_EPS0 * DUSTY_KB * T_i / (n_i * DUSTY_EC * DUSTY_EC)));
    return (4.0/3.0) * sqrt(2.0*M_PI) * n_i * a*a * v_thi
           * b90_over_a * b90_over_a * ln_Lambda;
}

double dust_ion_coulomb_logarithm(
    double a, double Z_d, double T_i, double lambda_D)
{
    /* ln Lambda = ln( lambda_D / max(a, b_90) )
     * b_90 = Z_d*e^2/(4*pi*eps0 * 3*k_B*T_i) */
    if (a <= 0.0 || lambda_D <= 0.0) return 1.0;
    double b_90 = Z_d * DUSTY_EC * DUSTY_EC
                  / (4.0 * M_PI * DUSTY_EPS0 * 3.0 * DUSTY_KB * T_i);
    double r_min = (a > b_90) ? a : b_90;
    if (r_min <= 0.0) return 1.0;
    double ratio = lambda_D / r_min;
    if (ratio <= 1.0) return 0.0;
    return log(ratio);
}

/* ================================================================
 * L3 — Diffusion and Mobility
 * ================================================================ */

double dust_diffusion_coefficient(double T_d, double m_d, double nu_dn)
{
    /* Einstein relation: D_d = k_B*T_d / (m_d*nu_dn)
     * Valid in weakly-coupled, collision-dominated regime. */
    if (T_d <= 0.0 || m_d <= 0.0 || nu_dn <= 0.0) return 0.0;
    return DUSTY_KB * T_d / (m_d * nu_dn);
}

double dust_mobility(double Z_d, double m_d, double nu_tot)
{
    /* Electrical mobility: mu_d = e*Z_d / (m_d*nu_tot)
     * Relates drift velocity to electric field. */
    if (m_d <= 0.0 || nu_tot <= 0.0) return 0.0;
    return DUSTY_EC * Z_d / (m_d * nu_tot);
}

double dust_drift_velocity(double mu_d, double E)
{
    /* Steady-state drift: v_drift = mu_d * E */
    return mu_d * E;
}

double dust_mean_free_path(double v_thd, double nu_tot)
{
    /* lambda_mfp = v_thd / nu_tot */
    if (nu_tot <= 0.0) return 0.0;
    return v_thd / nu_tot;
}

/* ================================================================
 * L4 — Dust Temperature Equation
 * ================================================================ */

double dust_heating_rate_ion_bombardment(
    double a, double n_i, double T_i, double m_i,
    double phi_f, double T_d)
{
    /* Ion bombardment heating: ions accelerated into grain deposit energy.
     * P_heat ≈ I_i * |phi_f| (simplified) */
    if (a <= 0.0 || n_i <= 0.0) return 0.0;
    double I_i = M_PI * a*a * DUSTY_EC * n_i
                 * sqrt(8.0 * DUSTY_KB * T_i / (M_PI * m_i))
                 * (1.0 - DUSTY_EC * phi_f / (DUSTY_KB * T_i));
    if (I_i < 0.0) I_i = 0.0;
    return I_i * fabs(phi_f);
}

double dust_cooling_rate_neutral(
    double T_d, double T_n, double nu_dn, double m_d, double m_n)
{
    /* Neutral gas cooling: energy transfer per elastic collision.
     * P_cool = (3/2)*k_B*(T_d-T_n)*nu_dn * 2*m_d*m_n/(m_d+m_n)^2 */
    if (m_d <= 0.0 || m_n <= 0.0 || nu_dn <= 0.0) return 0.0;
    double mass_factor = 2.0 * m_d * m_n / ((m_d + m_n) * (m_d + m_n));
    return 1.5 * DUSTY_KB * (T_d - T_n) * nu_dn * mass_factor;
}

double dust_radiative_cooling_power(
    double a, double epsilon, double T_d, double T_env)
{
    /* Stefan-Boltzmann: P_rad = 4*pi*a^2 * epsilon * sigma * (T_d^4 - T_env^4)
     * sigma = 2*pi^5*k_B^4/(15*h^3*c^2) ≈ 5.67e-8 W/(m^2 K^4) */
    if (a <= 0.0) return 0.0;
    double sigma_sb = 5.670374419e-8;
    double Td4 = T_d * T_d * T_d * T_d;
    double Te4 = T_env * T_env * T_env * T_env;
    return 4.0 * M_PI * a * a * epsilon * sigma_sb * (Td4 - Te4);
}

double dust_steady_state_temperature(
    double a, double n_i, double T_i, double m_i,
    double phi_f, double T_n, double nu_dn, double m_d, double m_n,
    double T_min, double T_max, double tol)
{
    /* Solve P_heat(T_d) = P_cool(T_d) by bisection.
     * The steady state balances ion heating against neutral cooling. */
    double T_lo = T_min;
    double T_hi = T_max;

    if (T_lo <= 0.0) T_lo = T_n;
    if (T_hi <= T_lo) T_hi = T_lo + 100.0;

    for (int iter = 0; iter < 50; iter++) {
        double T_mid = 0.5 * (T_lo + T_hi);
        double P_heat = dust_heating_rate_ion_bombardment(a, n_i, T_i, m_i, phi_f, T_mid);
        double P_cool = dust_cooling_rate_neutral(T_mid, T_n, nu_dn, m_d, m_n);
        double balance = P_heat - P_cool;

        if (fabs(balance) < tol * fmax(fabs(P_heat), fabs(P_cool) + 1e-30))
            return T_mid;
        if (balance > 0.0)
            T_lo = T_mid;
        else
            T_hi = T_mid;
    }
    return 0.5 * (T_lo + T_hi);
}

/* ================================================================
 * L5 — Anomalous Diffusion and Viscosity
 * ================================================================ */

double dust_anomalous_diffusion_msd(
    double alpha, double D_normal, double t)
{
    /* Generalized diffusion: <r^2(t)> ∝ t^alpha
     * alpha=1: normal Fickian, alpha<1: sub-diffusive, alpha>1: super-diffusive */
    if (t <= 0.0 || D_normal <= 0.0) return 0.0;
    return pow(t, alpha) * D_normal;
}

double dust_shear_viscosity(
    double n_d, double m_d, double a, double omega_pd,
    double Gamma, double Gamma_crit, double beta)
{
    /* Shear viscosity in strongly coupled Yukawa fluid.
     * eta = eta_0 * (Gamma/Gamma_crit)^beta
     * Saigo & Hamaguchi (2002), Phys. Plasmas 9, 1210 */
    if (n_d <= 0.0 || m_d <= 0.0 || omega_pd <= 0.0 || Gamma_crit <= 0.0)
        return 0.0;
    double eta_0 = n_d * m_d * a * a * omega_pd;
    return eta_0 * pow(Gamma / Gamma_crit, beta);
}

double dust_thermal_conductivity(double n_d, double D_d)
{
    /* Kinetic theory: kappa_thermal = (5/2) * n_d * k_B * D_d */
    if (n_d <= 0.0 || D_d <= 0.0) return 0.0;
    return 2.5 * n_d * DUSTY_KB * D_d;
}
