/**
 * plasma_waves.c -- Plasma Waves and Instabilities
 *
 * Implements dispersion relations for all canonical plasma wave modes
 * and instability growth rates.
 *
 * References:
 *   - Stix, "Waves in Plasmas" (1992)
 *   - Swanson, "Plasma Waves" (2003)
 *   - Buneman (1958), Weibel (1959)
 *
 * Knowledge Coverage:
 *   L4: Cold plasma dispersion, CMA diagram
 *   L5: Dispersion root finding
 *   L6: Langmuir, ion-acoustic, Alfven, whistler waves
 *   L8: Three-wave coupling, parametric decay
 */
#include "plasma_waves.h"
#include "plasma_mhd.h"
#include "plasma_params.h"
#include <math.h>
#include <stdlib.h>

/* ================================================================
 * L6.1: Electrostatic Waves (Unmagnetized)
 * ================================================================ */

double langmuir_dispersion(double k, double omega_pe, double v_the) {
    double omega2 = omega_pe * omega_pe + 3.0 * v_the * v_the * k * k;
    return sqrt(omega2);
}

double ion_acoustic_dispersion(double k, double cs, double lambda_De) {
    double k2_lD2 = k * k * lambda_De * lambda_De;
    return k * cs / sqrt(1.0 + k2_lD2);
}

void ion_acoustic_damped(double k, double cs, double lambda_De,
                         double Te, double Ti, double mi,
                         double *omega_r, double *gamma) {
    double omega = ion_acoustic_dispersion(k, cs, lambda_De);
    *omega_r = omega;
    if (omega <= 0.0 || k <= 0.0) { *gamma = 0.0; return; }
    double v_phase = omega / k;
    double v_the = sqrt(2.0 * K_B * Te / M_ELECTRON);
    double Te_Ti = Te / Ti;
    (void)mi; /* used for ion mass context in full kinetic treatment */
    double gamma_e = -omega * sqrt(M_PI / 8.0)
                     * (v_phase / v_the)
                     * exp(-0.5 * v_phase * v_phase / (v_the * v_the));
    double gamma_i = -omega * sqrt(M_PI / 8.0)
                     * pow(Te_Ti, 1.5)
                     * exp(-0.5 * Te_Ti);
    *gamma = gamma_e + gamma_i * sqrt(M_ELECTRON / mi);
}

/* ================================================================
 * L6.2: Cold Plasma Dispersion (Stix Parameters)
 * ================================================================ */

void cold_plasma_dielectric(double omega,
                            const double *wp_species,
                            const double *wc_species,
                            const int *signs,
                            int n_species,
                            ColdPlasmaDielectric *diel) {
    if (omega <= 0.0) {
        diel->S = diel->D = diel->P = diel->R = diel->L = 1.0;
        return;
    }
    double S = 1.0, D = 0.0, P = 1.0;
    for (int s = 0; s < n_species; s++) {
        double wp2 = wp_species[s] * wp_species[s];
        double wc = wc_species[s];
        double sign = (double)signs[s];
        double denom = omega * omega - wc * wc;
        if (fabs(denom) < 1e-30) denom = 1e-30;
        S -= wp2 / denom;
        D += sign * (wc / omega) * wp2 / denom;
        P -= wp2 / (omega * omega);
    }
    diel->S = S;
    diel->D = D;
    diel->P = P;
    diel->R = S + D;
    diel->L = S - D;
}

void cold_plasma_refractive_index(const ColdPlasmaDielectric *diel,
                                  double cos_theta,
                                  double *n2_fast, double *n2_slow) {
    double S = diel->S, P = diel->P;
    double R = diel->R, L = diel->L;
    (void)diel->D;
    double sin2 = 1.0 - cos_theta * cos_theta;
    double cos2 = cos_theta * cos_theta;
    double A = S * sin2 + P * cos2;
    double B = R * L * sin2 + P * S * (1.0 + cos2);
    double C = P * R * L;
    if (fabs(A) < 1e-30) {
        *n2_fast = (fabs(B) > 1e-30) ? C / B : 0.0;
        *n2_slow = 0.0;
        return;
    }
    double disc = B * B - 4.0 * A * C;
    if (disc < 0.0) disc = 0.0;
    double sqrt_disc = sqrt(disc);
    *n2_fast = (B + sqrt_disc) / (2.0 * A);
    *n2_slow = (B - sqrt_disc) / (2.0 * A);
    if (*n2_slow < 0.0) *n2_slow = 0.0;
    if (*n2_fast < 0.0) *n2_fast = 0.0;
}

int cold_plasma_resonance_condition(const ColdPlasmaDielectric *diel,
                                    double cos_theta) {
    double S = diel->S, P = diel->P;
    double sin2 = 1.0 - cos_theta * cos_theta;
    double A = S * sin2 + P * cos_theta * cos_theta;
    return (fabs(A) < 1e-10) ? 1 : 0;
}

void cold_plasma_cutoffs(const ColdPlasmaDielectric *diel,
                         int *has_P_cut, int *has_R_cut, int *has_L_cut) {
    *has_P_cut = (fabs(diel->P) < 1e-10) ? 1 : 0;
    *has_R_cut = (fabs(diel->R) < 1e-10) ? 1 : 0;
    *has_L_cut = (fabs(diel->L) < 1e-10) ? 1 : 0;
}

/* ================================================================
 * L6.3: MHD Waves
 * ================================================================ */

double alfven_wave_dispersion(double k_parallel, double vA) {
    return k_parallel * vA;
}

double fast_wave_dispersion(double k, double cs, double vA, double cos_theta) {
    double vF = fast_magnetosonic_speed(cs, vA, cos_theta);
    return k * vF;
}

double slow_wave_dispersion(double k, double cs, double vA, double cos_theta) {
    double vS = slow_magnetosonic_speed(cs, vA, cos_theta);
    return k * vS;
}

double whistler_dispersion(double k, double omega_pe, double omega_ce) {
    if (omega_pe <= 0.0 || omega_ce <= 0.0) return 0.0;
    double kc_over_wpe = k * C_LIGHT / omega_pe;
    double num = omega_ce * kc_over_wpe * kc_over_wpe;
    return num / (1.0 + kc_over_wpe * kc_over_wpe);
}

double whistler_group_velocity(double k, double omega_pe, double omega_ce) {
    if (omega_pe <= 0.0 || omega_ce <= 0.0) return 0.0;
    double kc_over_wpe = k * C_LIGHT / omega_pe;
    double denom = 1.0 + kc_over_wpe * kc_over_wpe;
    double numerator = 2.0 * omega_ce * C_LIGHT * C_LIGHT * k
                       / (omega_pe * omega_pe);
    return numerator / (denom * denom);
}

double lower_hybrid_dispersion(double omega_pi, double omega_ci,
                                double omega_ce) {
    double num = omega_ci * omega_ce + omega_pi * omega_pi;
    double denom = 1.0 + omega_pi * omega_pi / (omega_ce * omega_ce);
    if (denom <= 0.0) return 0.0;
    return sqrt(num / denom);
}

/* ================================================================
 * L6.4: Instability Growth Rates
 * ================================================================ */

double two_stream_growth_rate(double n_beam, double n_plasma, double omega_pe) {
    if (n_plasma <= 0.0 || omega_pe <= 0.0) return 0.0;
    double alpha = n_beam / n_plasma;
    if (alpha <= 0.0) return 0.0;
    double cbrt_alpha = pow(alpha, 1.0/3.0);
    return (sqrt(3.0) / pow(2.0, 4.0/3.0)) * cbrt_alpha * omega_pe;
}

double weibel_growth_rate(double omega_pe, double v_th, double T_perp,
                          double T_par, double beta) {
    if (T_par <= 0.0 || omega_pe <= 0.0) return 0.0;
    double anisotropy = T_perp / T_par - 1.0;
    if (anisotropy <= 0.0) return 0.0;
    return omega_pe * (v_th / C_LIGHT) * sqrt(anisotropy * beta / 2.0);
}

int firehose_unstable(double p_parallel, double p_perp, double B) {
    double p_mag = B * B / MU_0;
    return (p_parallel > p_perp + p_mag) ? 1 : 0;
}

int mirror_unstable(double p_parallel, double p_perp, double B) {
    if (p_parallel <= 0.0) return 0;
    double p_mag = B * B / (2.0 * MU_0);
    return (p_perp / p_parallel > 1.0 + p_mag / p_parallel) ? 1 : 0;
}

double kelvin_helmholtz_growth(double k, double v0, double rho1, double rho2) {
    if (rho1 + rho2 <= 0.0) return 0.0;
    return 0.5 * k * fabs(v0) * sqrt(rho1 * rho2) / (rho1 + rho2);
}

double rayleigh_taylor_growth(double g_eff, double k,
                              double rho_heavy, double rho_light) {
    if (rho_heavy + rho_light <= 0.0) return 0.0;
    double Atwood = (rho_heavy - rho_light) / (rho_heavy + rho_light);
    if (Atwood <= 0.0) return 0.0;
    return sqrt(g_eff * k * Atwood);
}

/* ================================================================
 * L5: Dispersion Root Finder
 * ================================================================ */

int find_dispersion_roots(DispersionFunc D, const double *params,
                          double omega_min, double omega_max,
                          int n_scan, double *roots, int max_roots) {
    if (!D || !roots || max_roots <= 0 || n_scan <= 0) return 0;
    double domega = (omega_max - omega_min) / n_scan;
    int n_found = 0;
    double prev_D = D(omega_min, params);
    for (int i = 1; i <= n_scan && n_found < max_roots; i++) {
        double omega = omega_min + i * domega;
        double curr_D = D(omega, params);
        if (prev_D * curr_D <= 0.0) {
            double a = omega - domega, b = omega;
            double fa = D(a, params);
            for (int k = 0; k < 50; k++) {
                double c = (a + b) / 2.0;
                double fc = D(c, params);
                if (fabs(fc) < 1e-12 || (b - a) < 1e-14) {
                    roots[n_found++] = c;
                    break;
                }
                if (fa * fc < 0.0) { b = c; }
                else { a = c; fa = fc; }
            }
        }
        prev_D = curr_D;
    }
    return n_found;
}

/* ================================================================
 * L8: Nonlinear Three-Wave Coupling
 * ================================================================ */

int three_wave_resonance(double w1, double k1, double w2, double k2,
                         double w3, double k3, double tolerance) {
    int freq_match = (fabs(w1 + w2 - w3) < tolerance * fmax(fabs(w3), 1.0));
    int wave_match = (fabs(k1 + k2 - k3) < tolerance * fmax(fabs(k3), 1.0));
    return (freq_match && wave_match) ? 1 : 0;
}

double parametric_decay_growth(double V_coupling, double k0,
                                double omega1, double omega2) {
    if (omega1 <= 0.0 || omega2 <= 0.0) return 0.0;
    double V2 = V_coupling * V_coupling;
    double k02 = k0 * k0;
    return sqrt(V2 * k02 / (4.0 * omega1 * omega2));
}
