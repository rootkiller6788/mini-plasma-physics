/**
 * fluid_instabilities.c -- Fluid/MHD Macro-instabilities
 *
 * L6: Rayleigh-Taylor, Kelvin-Helmholtz, interchange,
 *     ballooning, tearing mode, kink, sausage instabilities
 * L7: Tokamak stability (ITER), astrophysical MHD
 *
 * References:
 *   Chandrasekhar, "Hydrodynamic and Hydromagnetic Stability" (1961)
 *   Freidberg, "Ideal MHD" (2014)
 *   Furth, Killeen, Rosenbluth, Phys. Fluids 6, 459 (1963)
 *   Connor, Hastie, Taylor, Proc. R. Soc. A 365, 1 (1979)
 */

#include "waves_instabilities.h"
#include "plasma_instabilities.h"
#include <math.h>

/* ===============================================================
 * L6: Rayleigh-Taylor Instability
 * =============================================================== */

double rayleigh_taylor_growth(double g_eff, double k,
                               double rho_heavy, double rho_light)
{
    if (rho_heavy <= 0.0 || rho_light < 0.0 || k <= 0.0 || g_eff <= 0.0)
        return 0.0;
    double A = (rho_heavy - rho_light) / (rho_heavy + rho_light);
    if (A <= 0.0) return 0.0; /* Stable if light on top of heavy */
    return sqrt(g_eff * k * A);
}

double magnetic_rt_growth(double g_eff, double k, double rho_heavy,
                           double rho_light, double B_perp)
{
    double A = (rho_heavy > rho_light && rho_heavy + rho_light > 0.0)
               ? (rho_heavy - rho_light) / (rho_heavy + rho_light) : 0.0;
    if (A <= 0.0 || k <= 0.0) return 0.0;

    double gamma2_hydro = g_eff * k * A;
    double mag_stab = (k * B_perp) * (k * B_perp)
                      / (MU0 * (rho_heavy + rho_light));
    double gamma2 = gamma2_hydro - mag_stab;
    return (gamma2 > 0.0) ? sqrt(gamma2) : 0.0;
}

/* ===============================================================
 * L6: Kelvin-Helmholtz Instability
 * =============================================================== */

double kelvin_helmholtz_growth(double k, double v0,
                                double rho1, double rho2)
{
    if (k <= 0.0 || rho1 <= 0.0 || rho2 <= 0.0) return 0.0;
    double v0_half = fabs(v0) * 0.5;
    return k * v0_half * sqrt(rho1 * rho2) / (rho1 + rho2);
}

double magnetic_kh_growth(double k, double v0, double rho1, double rho2,
                           double B_parallel)
{
    if (k <= 0.0 || rho1 <= 0.0 || rho2 <= 0.0) return 0.0;

    double rho_sum = rho1 + rho2;
    double gamma2_hydro = k*k * v0*v0 * 0.25
                          * rho1*rho2 / (rho_sum*rho_sum);
    double mag_stab = k*k * B_parallel*B_parallel / (MU0 * rho_sum);
    double gamma2 = gamma2_hydro - mag_stab;
    return (gamma2 > 0.0) ? sqrt(gamma2) : 0.0;
}

/* ===============================================================
 * L6: Interchange (Flute) Instability
 * =============================================================== */

double interchange_growth_rate(double grad_ln_n, double cs,
                                double R_curvature)
{
    /* gamma^2 = -grad_ln_n * 2*cs^2 / R_c (for unfavorable curvature)
     * grad_ln_n is negative for outward-decreasing density in tokamak */
    if (R_curvature <= 0.0) return 0.0;
    double driving = -grad_ln_n * 2.0 * cs * cs / R_curvature;
    return (driving > 0.0) ? sqrt(driving) : 0.0;
}

/* ===============================================================
 * L6: Ballooning Instability
 * =============================================================== */

double ballooning_growth_rate(double cs, double grad_p_over_p,
                               double R_curvature, double q)
{
    /* gamma^2 = 2*cs^2 * grad_p/p * 1/R_c * (unfavorable curvature)
     * q is the safety factor; higher q -> more unstable */
    if (R_curvature <= 0.0 || q <= 0.0) return 0.0;
    double drive = 2.0 * cs * cs * grad_p_over_p / R_curvature;
    /* Ballooning is stabilized at low q, unstable at high q */
    double threshold = 1.0;
    double effective_drive = drive * (q - threshold) / q;
    return (effective_drive > 0.0) ? sqrt(effective_drive) : 0.0;
}

/* ===============================================================
 * L6: Resistive Tearing Mode
 * =============================================================== */

double tearing_mode_growth_rate(double eta, double a, double v_A,
                                 double delta_prime)
{
    /* Furth-Killeen-Rosenbluth scaling:
     * gamma = eta^(3/5) * delta_prime^(4/5) * a^(-2/5)
     *         * (mu0*rho)^(-2/5) * v_A^(2/5)
     * Simplified: gamma ~ (eta/mu0)^(3/5) * (delta_prime*a)^(4/5)
     *                    * tau_A^(-2/5) * tau_R^(-3/5) */
    if (eta <= 0.0 || a <= 0.0 || v_A <= 0.0) return 0.0;

    double tau_A = a / v_A;
    double tau_R = MU0 * a * a / eta;

    if (tau_A <= 0.0 || tau_R <= 0.0) return 0.0;

    /* gamma = tau_R^(-3/5) * tau_A^(-2/5) * (delta_prime*a)^(4/5) */
    double dp_a = delta_prime * a;
    return pow(tau_R, -0.6) * pow(tau_A, -0.4) * pow(fabs(dp_a), 0.8);
}

/* ===============================================================
 * L6: Kink and Sausage Instabilities
 * =============================================================== */

int kink_unstable(double q, int m, int n)
{
    /* Kruskal-Shafranov: m=1, n=1 kink unstable for q < 1
     * General: mode (m,n) is unstable when q < m/n */
    if (n <= 0) return 0;
    double q_crit = (double)m / (double)n;
    return (q < q_crit) ? 1 : 0;
}

int sausage_unstable(double pressure_gradient, double btheta_gradient)
{
    /* m=0 sausage mode: unstable when radial pressure gradient
     * overcomes B_theta gradient stabilization.
     * In a pure z-pinch without axial B: always unstable */
    (void)pressure_gradient;
    (void)btheta_gradient;
    /* Simplified criterion: sausage is always unstable without axial B */
    return 1;
}

/* ===============================================================
 * L6: Additional MHD Instability Metrics
 * =============================================================== */

/**
 * Growth rate for gravitational interchange (g-mode) in stratified plasma.
 * Applied to solar prominences and magnetosphere.
 *
 * gamma^2 = g * (1/rho * d(rho)/dz) for unstable stratification
 *
 * @param g         Effective gravity [m/s^2]
 * @param drho_dz   Density gradient [kg/m^4] (positive upward = unstable)
 * @param rho       Density [kg/m^3]
 * @return          Growth rate [1/s] or 0 if stable
 */
double gravitational_interchange_growth(double g, double drho_dz, double rho)
{
    if (rho <= 0.0) return 0.0;
    double n_sq = g * drho_dz / rho;
    return (n_sq > 0.0) ? sqrt(n_sq) : 0.0;
}

/**
 * Parker instability growth rate (magnetic buoyancy).
 *
 * Magnetic field lines buoyantly rise due to cosmic ray pressure
 * or reduced density in magnetic flux tubes.
 *
 * gamma_Parker = c_s * k * sqrt(beta * (1 - 1/gamma_adiabatic))
 *
 * Reference: Parker, Astrophys. J. 145, 811 (1966)
 * L7 application: Galactic dynamo, ISM structure
 */
double parker_instability_growth(double c_s, double k, double beta,
                                  double gamma_adiabatic)
{
    if (c_s <= 0.0 || k <= 0.0 || gamma_adiabatic <= 1.0) return 0.0;
    double factor = beta * (1.0 - 1.0/gamma_adiabatic);
    return (factor > 0.0) ? c_s * k * sqrt(factor) : 0.0;
}

/**
 * Resistive interchange (g-mode with resistivity).
 *
 * gamma_RIC = (eta * g * k_perp^2 / v_A^2)^(1/3)
 *
 * Similar to RT but limited by resistivity allowing field line bending.
 *
 * Reference: Furth et al., Phys. Fluids 6, 459 (1963)
 */
double resistive_interchange_growth(double eta, double g,
                                     double k_perp, double v_A)
{
    if (eta <= 0.0 || v_A <= 0.0) return 0.0;
    return pow(eta * g * k_perp * k_perp / (v_A * v_A), 1.0/3.0);
}

/**
 * Double tearing mode growth rate.
 *
 * When two rational surfaces are close (delta_r < resonant layer width),
 * the tearing modes on each surface couple, greatly enhancing growth.
 *
 * gamma_DTM ~ gamma_single * (delta_r_wall / delta_r)^(5/4)
 *
 * Reference: Pritchett, Lee, Drake, Phys. Fluids 23, 1368 (1980)
 *
 * @param gamma_single  Growth rate of isolated tearing mode [1/s]
 * @param delta_r       Distance between rational surfaces [m]
 * @param layer_width   Tearing layer width [m]
 * @return              Enhanced growth rate
 */
double double_tearing_growth(double gamma_single, double delta_r,
                              double layer_width)
{
    if (delta_r <= 0.0 || layer_width <= 0.0) return gamma_single;
    if (delta_r > 5.0 * layer_width) return gamma_single;
    return gamma_single * pow(layer_width / delta_r, 1.25);
}
