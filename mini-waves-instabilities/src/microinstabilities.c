/**
 * microinstabilities.c -- Drift Wave and Micro-instabilities
 *
 * L6: Universal drift wave, ITG (ion temperature gradient),
 *     TEM (trapped electron mode), ETG (electron temperature gradient)
 * L7: Tokamak anomalous transport (ITER), core turbulence
 *
 * These instabilities are the primary drivers of turbulent transport
 * in magnetically confined fusion plasmas. They are "universal" in
 * the sense that any plasma with density/temperature gradients
 * is subject to them.
 *
 * References:
 *   Weiland, "Collective Modes in Inhomogeneous Plasma" (2000)
 *   Kadomtsev, "Plasma Turbulence" (1965)
 *   Romanelli, Phys. Fluids B 1, 1018 (1989)
 *   Dorland et al., Phys. Rev. Lett. 85, 5579 (2000)
 */

#include "waves_instabilities.h"
#include "plasma_instabilities.h"
#include <math.h>

/* ===============================================================
 * L6: Universal Drift Wave
 * =============================================================== */

double drift_wave_frequency(double k_y, double T_e, double grad_ln_n,
                             double B)
{
    /* omega_star = k_y * T_e * grad(ln n) / (e * B)
     * This is the diamagnetic drift frequency.
     * Full drift wave (cold ion slab):
     *   omega = omega_star / (1 + k_perp^2 * rho_s^2) */
    if (B <= 0.0) return 0.0;
    return k_y * T_e * grad_ln_n / (E_CHARGE * B);
}

double drift_wave_growth_rate(double omega_star, double k_parallel,
                               double v_th_e)
{
    /* Kinetic growth from inverse electron Landau damping.
     * gamma = sqrt(pi) * omega^2/|k_par|v_th_e
     *         * exp(-omega^2/(k_par^2 * v_th_e^2))
     *
     * In collisional regime:
     * gamma ~ k_perp^2 * rho_s^2 * omega_star^2 / nu_ei */
    if (fabs(k_parallel) < 1e-30 || v_th_e <= 0.0) return 0.0;
    double arg = omega_star / (fabs(k_parallel) * v_th_e);
    if (fabs(arg) > 6.0) return 0.0;
    return fabs(omega_star) * sqrt(PLASMA_PI) * fabs(arg)
           * exp(-arg * arg);
}

/* ===============================================================
 * L6: ITG (Ion Temperature Gradient) Mode
 * =============================================================== */

double itg_mode_frequency(double omega_star, double eta_i,
                           double eta_i_crit, double k_perp_rho_s)
{
    /* ITG frequency:
     * omega = omega_star * (eta_i - eta_i_crit) / (1 + k_perp^2*rho_s^2)
     * where eta_i = d(ln T_i)/d(ln n), eta_i_crit ~ 1-2 */
    if (eta_i <= eta_i_crit) return 0.0;
    double kprs2 = k_perp_rho_s * k_perp_rho_s;
    return omega_star * (eta_i - eta_i_crit) / (1.0 + kprs2);
}

double itg_growth_rate(double omega_star, double eta_i,
                        double eta_i_crit, double k_perp_rho_s)
{
    /* ITG growth: gamma = omega_star * sqrt(eta_i - eta_i_crit)
     *                     * k_perp*rho_s / (1 + k_perp^2*rho_s^2) */
    if (eta_i <= eta_i_crit) return 0.0;
    double krs = k_perp_rho_s;
    return fabs(omega_star) * sqrt(eta_i - eta_i_crit)
           * krs / (1.0 + krs * krs);
}

/* ===============================================================
 * L6: TEM (Trapped Electron Mode)
 * =============================================================== */

double tem_frequency(double omega_star, double eta_e, double epsilon)
{
    /* TEM frequency driven by toroidal precession of trapped electrons:
     * omega ~ omega_star * sqrt(epsilon) * eta_e
     * where epsilon = r/R (inverse aspect ratio) */
    if (epsilon <= 0.0 || epsilon >= 1.0) return 0.0;
    return omega_star * sqrt(epsilon) * eta_e;
}

double tem_growth_rate(double omega_star, double eta_e, double epsilon)
{
    /* TEM growth rate:
     * gamma = omega_star * sqrt(epsilon) * (eta_e - eta_e_crit)
     * with eta_e_crit ~ 0.5-1 */
    double eta_e_crit = 0.8;
    if (eta_e <= eta_e_crit || epsilon <= 0.0 || epsilon >= 1.0) return 0.0;
    return fabs(omega_star) * sqrt(epsilon) * (eta_e - eta_e_crit);
}

/* ===============================================================
 * L6: ETG (Electron Temperature Gradient) Mode
 * =============================================================== */

double etg_growth_rate(double omega_star_e, double eta_e,
                        double eta_e_crit, double k_perp_rho_e)
{
    /* ETG is the electron-scale analogue of ITG.
     * Short wavelength: k_perp * rho_e ~ 1
     * gamma = omega_star_e * sqrt(eta_e - eta_e_crit)
     *         * k_perp*rho_e / (1 + k_perp^2*rho_e^2)
     *
     * ETG produces a "streamer" dominated turbulence
     * responsible for anomalous electron thermal transport. */
    if (eta_e <= eta_e_crit) return 0.0;
    double kre = k_perp_rho_e;
    return fabs(omega_star_e) * sqrt(eta_e - eta_e_crit)
           * kre / (1.0 + kre * kre);
}

/* ===============================================================
 * L6: EPED Model (Combined ITG+TEM stability)
 * =============================================================== */

/**
 * Estimate the critical temperature gradient for ITG/TEM stability.
 *
 * In the EPED model (Snyder et al., 2002), the pedestal width
 * is set by the most unstable peeling-ballooning mode.
 * The critical gradient is:
 *
 *   (dT/dr)_crit = (T/r) * (R/r) * k_theta * rho_s * q / s
 *
 * where s = (r/q)*(dq/dr) is the magnetic shear.
 *
 * Reference: Snyder et al., Phys. Plasmas 9, 2037 (2002)
 * L7 application: ITER pedestal prediction
 *
 * @return Critical temperature gradient scale length L_T = T/(dT/dr)
 */
double critical_temperature_gradient(double T, double r, double R,
                                      double k_theta_rho_s,
                                      double q, double magnetic_shear)
{
    if (r <= 0.0 || R <= 0.0 || q <= 0.0 || magnetic_shear <= 0.0)
        return 0.0;
    return T / (r * k_theta_rho_s * q / (R * magnetic_shear) + 1e-30);
}

/**
 * Estimate the turbulent thermal diffusivity from ITG/TEM mixing length.
 *
 * chi_turb = gamma_max / k_perp^2
 *
 * This mixing-length estimate is the basis for reduced transport models
 * (GLF23, TGLF, etc.) used in integrated tokamak modeling.
 *
 * Reference: Waltz et al., Phys. Plasmas 4, 2482 (1997)
 */
double mixing_length_diffusivity(double gamma_max, double k_perp)
{
    if (k_perp <= 0.0) return 0.0;
    return gamma_max / (k_perp * k_perp);
}

/**
 * Zonal flow shearing rate for turbulence suppression.
 *
 * ExB shearing rate:
 *   omega_ExB = (r/q) * d(Er/(r*B))/dr
 *
 * When omega_ExB > gamma_max (linear growth rate),
 * turbulence is sheared apart and transport is suppressed.
 * This is the mechanism behind L-H transitions and ITBs
 * (internal transport barriers).
 *
 * Reference: Biglari, Diamond, Terry, Phys. Fluids B 2, 1 (1990)
 * L7 application: ITER H-mode access prediction
 *
 * @return 1 if turbulence is suppressed (omega_ExB > gamma), 0 otherwise
 */
int zonal_flow_suppression(double omega_ExB, double gamma_max)
{
    return (omega_ExB > gamma_max) ? 1 : 0;
}

/**
 * Dimits shift: the nonlinear upshift of the critical gradient
 * above the linear threshold due to zonal flow self-organization.
 *
 * Delta_crit = gamma_damping / omega_star * (q * rho_s / L_n)
 *
 * Reference: Dimits et al., Phys. Plasmas 7, 969 (2000)
 */
double dimits_shift(double gamma_damping, double omega_star,
                     double q, double rho_s_over_Ln)
{
    if (omega_star <= 0.0) return 0.0;
    return gamma_damping / omega_star * q * rho_s_over_Ln;
}
