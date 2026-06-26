/**
 * parametric_instabilities.c -- Parametric and Nonlinear Wave Coupling
 *
 * L8: Three-wave resonance, parametric decay (PDI), SRS, SBS,
 *     two-plasmon decay, modulational instability, three-wave RK4 dynamics
 * L7: ICF laser-plasma interaction (NIF, LMJ), ionospheric heating
 *
 * References:
 *   Sagdeev & Galeev, "Nonlinear Plasma Theory" (1969)
 *   Nishikawa, J. Phys. Soc. Jpn. 24, 916 (1968)
 *   Drake et al., Phys. Fluids 17, 778 (1974)
 *   Zakharov, Sov. Phys. JETP 35, 908 (1972)
 */

#include "waves_instabilities.h"
#include "nonlinear_waves.h"
#include <math.h>
#include <stdlib.h>

/* ===============================================================
 * L8: Three-Wave Resonance Conditions
 * =============================================================== */

int three_wave_resonance_check(double w1, double k1,
                                double w2, double k2,
                                double w3, double k3,
                                int sum_freq, double tolerance)
{
    double w_match, k_match;
    if (sum_freq) {
        w_match = w1 + w2;
        k_match = k1 + k2;
    } else {
        w_match = fabs(w1 - w2);
        k_match = fabs(k1 - k2);
    }

    double dw = fabs(w_match - w3);
    double dk = fabs(k_match - k3);
    double w_scale = fmax(fmax(fabs(w1), fabs(w2)), fabs(w3));
    double k_scale = fmax(fmax(fabs(k1), fabs(k2)), fabs(k3));

    if (w_scale < 1e-30) w_scale = 1.0;
    if (k_scale < 1e-30) k_scale = 1.0;

    return (dw / w_scale < tolerance && dk / k_scale < tolerance) ? 1 : 0;
}

double three_wave_coupling_coeff(double k1, double k2, double k3,
                                  double w1, double w2, double w3,
                                  double omega_pe)
{
    /* Electrostatic coupling coefficient for Langmuir waves:
     * V = (e/(2*m_e)) * k3 / (w1*w2) * |E1||E2| coupling
     * Normalized form: V ~ omega_pe * k3 / sqrt(w1*w2)
     * w3 participates in the resonance condition; included for
     * future use in frequency-dependent coupling. */
    if (w1 <= 0.0 || w2 <= 0.0 || omega_pe <= 0.0) return 0.0;
    double k_eff = fabs(k1) + fabs(k2) + fabs(k3);
    double w3_geom = sqrt(w1 * w2); /* geometric mean frequency */
    (void)w3; /* Reserved for frequency-dependent coupling */
    return omega_pe * k_eff / w3_geom * 0.1;
}

/* ===============================================================
 * L8: Parametric Instability Growth Rates
 * =============================================================== */

double parametric_decay_growth(double V_coupling, double E_pump,
                                double w1, double w2)
{
    if (w1 <= 0.0 || w2 <= 0.0) return 0.0;
    return sqrt(V_coupling * V_coupling * E_pump * E_pump
                / (4.0 * w1 * w2));
}

double srs_growth_rate(double k_L, double v_osc,
                        double omega_pe, double omega_s)
{
    /* Stimulated Raman Scattering:
     * gamma_SRS = (k_L * v_osc / 4) * sqrt(omega_pe / omega_s)
     * where v_osc = e*E0/(m_e*omega0) is quiver velocity */
    if (omega_s <= 0.0 || omega_pe <= 0.0) return 0.0;
    return (k_L * v_osc * 0.25) * sqrt(omega_pe / omega_s);
}

double sbs_growth_rate(double omega_pi, double v_osc,
                        double w0, double ws, double v_th_e)
{
    /* Stimulated Brillouin Scattering:
     * gamma_SBS ~ (omega_pi/2) * sqrt(v_osc^2 * w0 / (v_th_e * c * sqrt(w0*ws))) */
    if (w0 <= 0.0 || ws <= 0.0 || v_th_e <= 0.0) return 0.0;
    double factor = v_osc * v_osc * w0
                    / (v_th_e * C_LIGHT * sqrt(w0 * ws));
    return 0.5 * omega_pi * sqrt(fmax(factor, 0.0));
}

double two_plasmon_decay_growth(double k_L, double v_osc,
                                 double w0, double omega_pe)
{
    /* Two-plasmon decay (TPD) at n_c/4:
     * gamma_TPD = (k_L * v_osc / 4) * sqrt(w0/omega_pe - 1)
     * Occurs when w0 ~ 2*omega_pe (quarter-critical density) */
    if (omega_pe <= 0.0 || w0 <= omega_pe) return 0.0;
    return (k_L * v_osc * 0.25) * sqrt(w0 / omega_pe - 1.0);
}

double modulational_instability_growth(double k_pert, double omega_pe,
                                        double k_De, double W,
                                        double n0, double T_e)
{
    /* Modulational instability of Langmuir waves:
     * gamma_MI = omega_pe * (k_pert/k_De) * sqrt(W/(n0*T_e))
     * Leads to Langmuir collapse and strong turbulence */
    if (omega_pe <= 0.0 || k_De <= 0.0 || n0 <= 0.0 || T_e <= 0.0)
        return 0.0;
    double W_norm = W / (n0 * T_e);
    return omega_pe * (k_pert / k_De) * sqrt(fmax(W_norm, 0.0));
}

/* ===============================================================
 * L8: Three-Wave System Dynamics (RK4)
 * =============================================================== */

/**
 * Three-wave interaction equations (normalized):
 *
 * dA1/dt = gamma1*A1 + V*A2*A3*sin(phi)
 * dA2/dt = gamma2*A2 + V*A1*A3*sin(phi)
 * dA3/dt = gamma3*A3 - V*A1*A2*sin(phi)
 * dphi/dt = delta_omega - V*cos(phi)*(A2*A3/A1 + A1*A3/A2 - A1*A2/A3)
 *
 * where phi = theta3 - theta1 - theta2 (phase mismatch)
 * and delta_omega = omega3 - omega1 - omega2 (frequency mismatch)
 *
 * Reference: Weiland & Wilhelmsson (1977), Ch. 1
 */
void three_wave_advance_rk4(ThreeWaveState *s, double dt,
                             double w1, double w2, double w3)
{
    if (!s || dt <= 0.0) return;

    double delta_omega = w3 - w1 - w2;

    /* Define the derivatives function inline */
    /* We'll use simple forward Euler for brevity in this implementation,
     * as full RK4 for a 4-variable system requires substep storage.
     * Note: RK4 would be preferred for accuracy. */

    double V = s->V;
    double A1 = s->A1, A2 = s->A2, A3 = s->A3;
    double phi = s->phi;
    double g1 = s->gamma1, g2 = s->gamma2, g3 = s->gamma3;

    /* RK4 step 1 */
    double sin_phi = sin(phi), cos_phi = cos(phi);
    double k1_A1 = g1*A1 + V*A2*A3*sin_phi;
    double k1_A2 = g2*A2 + V*A1*A3*sin_phi;
    double k1_A3 = g3*A3 - V*A1*A2*sin_phi;
    double k1_phi = delta_omega
                    - V*cos_phi*(A2*A3/(A1+1e-30)
                                 + A1*A3/(A2+1e-30)
                                 - A1*A2/(A3+1e-30));

    /* RK4 step 2 */
    double A1h = A1 + 0.5*dt*k1_A1;
    double A2h = A2 + 0.5*dt*k1_A2;
    double A3h = A3 + 0.5*dt*k1_A3;
    double phih = phi + 0.5*dt*k1_phi;
    sin_phi = sin(phih); cos_phi = cos(phih);
    double k2_A1 = g1*A1h + V*A2h*A3h*sin_phi;
    double k2_A2 = g2*A2h + V*A1h*A3h*sin_phi;
    double k2_A3 = g3*A3h - V*A1h*A2h*sin_phi;
    double k2_phi = delta_omega
                    - V*cos_phi*(A2h*A3h/(A1h+1e-30)
                                 + A1h*A3h/(A2h+1e-30)
                                 - A1h*A2h/(A3h+1e-30));

    /* RK4 step 3 */
    A1h = A1 + 0.5*dt*k2_A1;
    A2h = A2 + 0.5*dt*k2_A2;
    A3h = A3 + 0.5*dt*k2_A3;
    phih = phi + 0.5*dt*k2_phi;
    sin_phi = sin(phih); cos_phi = cos(phih);
    double k3_A1 = g1*A1h + V*A2h*A3h*sin_phi;
    double k3_A2 = g2*A2h + V*A1h*A3h*sin_phi;
    double k3_A3 = g3*A3h - V*A1h*A2h*sin_phi;
    double k3_phi = delta_omega
                    - V*cos_phi*(A2h*A3h/(A1h+1e-30)
                                 + A1h*A3h/(A2h+1e-30)
                                 - A1h*A2h/(A3h+1e-30));

    /* RK4 step 4 */
    A1h = A1 + dt*k3_A1;
    A2h = A2 + dt*k3_A2;
    A3h = A3 + dt*k3_A3;
    phih = phi + dt*k3_phi;
    sin_phi = sin(phih); cos_phi = cos(phih);
    double k4_A1 = g1*A1h + V*A2h*A3h*sin_phi;
    double k4_A2 = g2*A2h + V*A1h*A3h*sin_phi;
    double k4_A3 = g3*A3h - V*A1h*A2h*sin_phi;
    double k4_phi = delta_omega
                    - V*cos_phi*(A2h*A3h/(A1h+1e-30)
                                 + A1h*A3h/(A2h+1e-30)
                                 - A1h*A2h/(A3h+1e-30));

    /* Update */
    s->A1 += dt/6.0 * (k1_A1 + 2.0*k2_A1 + 2.0*k3_A1 + k4_A1);
    s->A2 += dt/6.0 * (k1_A2 + 2.0*k2_A2 + 2.0*k3_A2 + k4_A2);
    s->A3 += dt/6.0 * (k1_A3 + 2.0*k2_A3 + 2.0*k3_A3 + k4_A3);
    s->phi += dt/6.0 * (k1_phi + 2.0*k2_phi + 2.0*k3_phi + k4_phi);

    /* Keep phi in [-pi, pi] */
    while (s->phi > PLASMA_PI) s->phi -= 2.0 * PLASMA_PI;
    while (s->phi < -PLASMA_PI) s->phi += 2.0 * PLASMA_PI;

    /* Ensure non-negative amplitudes */
    if (s->A1 < 0.0) s->A1 = 0.0;
    if (s->A2 < 0.0) s->A2 = 0.0;
    if (s->A3 < 0.0) s->A3 = 0.0;
}

void manley_rowe_invariants(const ThreeWaveState *s,
                             double *I13, double *I23)
{
    if (!s || !I13 || !I23) return;
    *I13 = s->A1 * s->A1 + s->A3 * s->A3;
    *I23 = s->A2 * s->A2 + s->A3 * s->A3;
}

double nonlinear_frequency_shift(double A_sq, double Q)
{
    /* Nonlinear Schrodinger equation frequency shift:
     * delta_omega = Q * |A|^2
     * For Langmuir: Q = omega_pe/(4*n0*T_e) */
    return Q * A_sq;
}
