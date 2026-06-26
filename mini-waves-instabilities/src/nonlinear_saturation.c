/**
 * nonlinear_saturation.c -- Nonlinear Saturation Mechanisms
 *
 * L8: Quasilinear saturation, particle trapping, pump depletion,
 *     wave steepening, Langmuir collapse, nonlinear Landau damping,
 *     weak turbulence cascade phenomenology
 *
 * References:
 *   Vedenov, Velikhov, Sagdeev, Nucl. Fusion 1, 82 (1961)
 *   O'Neil, Phys. Fluids 8, 2255 (1965)
 *   Zakharov, Sov. Phys. JETP 35, 908 (1972)
 *   Sagdeev & Galeev, "Nonlinear Plasma Theory" (1969)
 */

#include "waves_instabilities.h"
#include "nonlinear_waves.h"
#include <math.h>

double quasilinear_saturation_time(double gamma, double omega_pe,
                                    double delta_v, double v_phase,
                                    double n_b_over_n0)
{
    /* Quasilinear plateau formation timescale:
     * tau_QL = (omega_pe/gamma)^2 * (delta_v/v_phase)^3
     *          / (omega_pe * n_b/n_0)
     * When t > tau_QL, the distribution flattens near v_phase. */
    if (gamma <= 0.0 || omega_pe <= 0.0 || v_phase <= 0.0
        || n_b_over_n0 <= 0.0) return 1e30;
    double ratio = omega_pe / gamma;
    double velocity_factor = delta_v / v_phase;
    return ratio * ratio * velocity_factor * velocity_factor
           * velocity_factor / (omega_pe * n_b_over_n0);
}

double trapping_saturation_field(double gamma, double k, double m, double q)
{
    /* When a particle is trapped in the wave potential, the instability
     * saturates. The trapping condition is:
     *   omega_bounce = sqrt(e*k*E/m) ~ gamma
     * So E_sat = m*gamma^2/(k*|q|)
     *
     * Reference: O'Neil (1965); Manheimer, Phys. Fluids 14, 579 (1971) */
    if (k <= 0.0 || fabs(q) <= 0.0) return 0.0;
    return m * gamma * gamma / (k * fabs(q));
}

double pump_depletion_fraction(double gamma_PDI, double gamma_pump_damping)
{
    /* In parametric decay, the pump is depleted when the daughter
     * waves grow to significant amplitude. The fractional energy
     * transfer before nonlinear saturation:
     *   f = gamma_PDI / (gamma_PDI + gamma_pump_damping)
     * Limited to [0, 1]. */
    if (gamma_PDI <= 0.0) return 0.0;
    double f = gamma_PDI / (gamma_PDI + gamma_pump_damping);
    return (f > 1.0) ? 1.0 : f;
}

double wave_steepening_distance(double c_s, double gamma_eff,
                                 double omega, double dn_over_n0)
{
    /* Nonlinear steepening of ion acoustic waves into shocks.
     * Steepening length:
     *   L_steep = c_s * (gamma_eff + 1) / (2 * omega * dn/n0)
     * After this distance, the wave profile develops discontinuities.
     *
     * Reference: Davidson (1972), Ch. 2 */
    if (omega <= 0.0 || dn_over_n0 <= 0.0) return 1e30;
    return c_s * (gamma_eff + 1.0) / (2.0 * omega * dn_over_n0);
}

int langmuir_collapse_threshold(double W, double n0, double T_e,
                                 double k, double lambda_De)
{
    /* Zakharov collapse criterion:
     * W/(n0*T_e) > k^2 * lambda_De^2
     *
     * When the Langmuir wave energy density exceeds this threshold,
     * the ponderomotive force expels plasma, creating a density cavity
     * that further traps the wave -> runaway collapse.
     *
     * Reference: Zakharov (1972); Robinson, Rev. Mod. Phys. 69, 507 (1997) */
    if (n0 <= 0.0 || T_e <= 0.0) return 0;
    double W_norm = W / (n0 * T_e);
    double k_lambda_sq = k * k * lambda_De * lambda_De;
    return (W_norm > k_lambda_sq) ? 1 : 0;
}

double nonlinear_landau_damping(double W_spectrum, double V_coupling_sq,
                                 double delta_omega)
{
    /* Nonlinear Landau damping (induced scattering):
     * Waves scatter off plasma particles:
     *   gamma_NL ~ sum_k' |V|^2 * W_k' / delta_omega
     *
     * This transfers energy between different spectral components.
     * It is the fundamental mechanism of weak turbulence cascades.
     *
     * Reference: Sagdeev & Galeev (1969), Ch. 2-3 */
    if (delta_omega <= 0.0) return 0.0;
    return V_coupling_sq * W_spectrum / delta_omega;
}

/* ===============================================================
 * L8: Weak Turbulence Cascade Rates
 * =============================================================== */

/**
 * Kolmogorov-type cascade rate for Langmuir turbulence.
 *
 * For isotropic 3D Langmuir turbulence:
 *   tau_cascade^(-1) = omega_pe * (W/(n0*T_e))^2
 *
 * This gives the energy transfer rate through the inertial range.
 *
 * Reference: Zakharov, Lvov, Falkovich, "Kolmogorov Spectra" (1992)
 */
double langmuir_cascade_rate(double omega_pe, double W, double n0, double T_e)
{
    if (omega_pe <= 0.0 || n0 <= 0.0 || T_e <= 0.0) return 0.0;
    double W_norm = W / (n0 * T_e);
    return omega_pe * W_norm * W_norm;
}

/**
 * Ion acoustic wave cascade rate (Kadomtsev-Petviashvili).
 *
 * tau_cascade^(-1) = omega_pi * (delta_n/n0)^2 * (k*lambda_De)^2
 *
 * Ion acoustic turbulence is important in collisionless shocks
 * and anomalous resistivity.
 */
double ion_acoustic_cascade_rate(double omega_pi, double dn_over_n0,
                                  double k_lambda_De)
{
    if (omega_pi <= 0.0) return 0.0;
    return omega_pi * dn_over_n0 * dn_over_n0
           * k_lambda_De * k_lambda_De;
}

/**
 * Energy-containing scale for drift wave turbulence.
 *
 * k_theta * rho_s ~ 0.3 (empirically, from gyrokinetic simulations
 * and experimental measurements in tokamaks).
 *
 * This gives the typical turbulent eddy size.
 *
 * Reference: Garbet et al., Nucl. Fusion 44, R93 (2004)
 * L7 application: ITER transport modeling
 */
double drift_wave_energy_containing_scale(double rho_s)
{
    return 0.3 / rho_s;
}

/**
 * Zonal flow growth rate from modulational instability of drift waves.
 *
 * The Reynolds stress from drift wave turbulence drives zonal flows
 * (n=0, m=0 ExB flows) that regulate turbulence.
 *
 * gamma_ZF = k_r^2 * rho_s^2 * gamma_drift * (chi_turb/chi_neo)
 *
 * Reference: Diamond et al., Plasma Phys. Control. Fusion 47, R35 (2005)
 */
double zonal_flow_growth_rate(double k_r_rho_s, double gamma_drift,
                               double chi_turb_over_chi_neo)
{
    if (gamma_drift <= 0.0) return 0.0;
    return k_r_rho_s * k_r_rho_s * gamma_drift
           * chi_turb_over_chi_neo;
}

/**
 * Critical gradient for avalanche propagation.
 *
 * In self-organized criticality models of transport:
 *   (dT/dr)_crit = (T/r) * (1 / (q * s))
 *
 * where s = (r/q)*dq/dr is magnetic shear.
 * When gradient exceeds critical, heat pulses (avalanches)
 * propagate ballistically.
 *
 * Reference: Carreras et al., Phys. Plasmas 3, 2903 (1996)
 */
double avalanche_critical_gradient(double T, double r, double q, double s)
{
    if (r <= 0.0 || q <= 0.0 || s <= 0.0) return 0.0;
    return T / (r * q * s);
}

/**
 * Nonlinear bounce frequency for deeply trapped particles.
 *
 * When wave amplitude is large, particles bounce in the wave potential:
 *   omega_b = sqrt(k^2 * e * phi / m)
 *
 * Saturation occurs when omega_b ~ gamma (linear growth rate).
 *
 * Reference: O'Neil & Malmberg, Phys. Fluids 11, 1754 (1968)
 */
double nonlinear_bounce_frequency(double k, double phi, double m, double q)
{
    if (m <= 0.0) return 0.0;
    return sqrt(k * k * fabs(q) * fabs(phi) / m);
}
