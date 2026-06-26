/**
 * electromagnetic_waves.c -- Electromagnetic Waves in Magnetized Plasma
 *
 * L4: Cold plasma dielectric tensor (Stix), Appleton-Hartree dispersion
 * L6: Whistler, O-mode, X-mode, Alfven, magnetosonic, ion cyclotron waves
 * L6: Lower/upper hybrid, R/L cutoffs, CMA diagram
 *
 * References:
 *   Stix, "Waves in Plasmas" (1992), Ch. 1-4
 *   Swanson, "Plasma Waves" (2003), Ch. 3-5
 *   Allis, Buchsbaum, Bers, "Waves in Anisotropic Plasmas" (1963)
 */

#include "waves_instabilities.h"
#include <math.h>
#include <string.h>
#include <stdio.h>

/* ===============================================================
 * L4: Stix Cold Plasma Dielectric Tensor
 * =============================================================== */

void stix_dielectric_tensor(double omega,
                            const double wp_species[],
                            const double wc_species[],
                            int n_species,
                            StixDielectric *diel)
{
    diel->S = 1.0; diel->D = 0.0; diel->P = 1.0;
    diel->R = 1.0; diel->L = 1.0;
    if (omega <= 0.0) return;

    for (int s = 0; s < n_species; s++) {
        double wp2 = wp_species[s] * wp_species[s];
        double wc = wc_species[s];
        double denom = omega * omega - wc * wc;

        if (fabs(denom) < 1e-30) {
            /* At exact resonance, use small offset to avoid singularity */
            denom = (denom >= 0) ? 1e-30 : -1e-30;
        }

        double eps_perp = wp2 / denom;
        double eps_cross = (wc / omega) * eps_perp;

        diel->S -= eps_perp;
        diel->D += eps_cross;
        diel->P -= wp2 / (omega * omega);
    }

    diel->R = diel->S + diel->D;
    diel->L = diel->S - diel->D;
}

void cold_plasma_n2(const StixDielectric *diel, double cos_theta,
                    double *n2_mode1, double *n2_mode2)
{
    double sin2 = 1.0 - cos_theta * cos_theta;
    if (sin2 < 0.0) sin2 = 0.0;
    double cos2 = cos_theta * cos_theta;

    double A = diel->S * sin2 + diel->P * cos2;
    double B = diel->R * diel->L * sin2
               + diel->P * diel->S * (1.0 + cos2);
    double C = diel->P * diel->R * diel->L;

    if (fabs(A) < 1e-30) {
        /* Degenerate case: resonance cone */
        *n2_mode1 = (fabs(B) > 1e-30) ? C / B : -1.0;
        *n2_mode2 = -1.0;
        return;
    }

    double disc = B * B - 4.0 * A * C;
    if (disc < 0.0) {
        *n2_mode1 = -1.0;
        *n2_mode2 = -1.0;
        return;
    }

    double sqrt_disc = sqrt(disc);
    double n2a = (B - sqrt_disc) / (2.0 * A);
    double n2b = (B + sqrt_disc) / (2.0 * A);

    /* Order so mode1 <= mode2 (mode2 is fast wave) */
    if (n2a <= n2b) {
        *n2_mode1 = (n2a >= 0) ? n2a : -1.0;
        *n2_mode2 = (n2b >= 0) ? n2b : -1.0;
    } else {
        *n2_mode1 = (n2b >= 0) ? n2b : -1.0;
        *n2_mode2 = (n2a >= 0) ? n2a : -1.0;
    }
}

int identify_wave_mode(const StixDielectric *diel,
                       double cos_theta, double omega, double omega_ce)
{
    (void)omega; /* omega used in resonance identification */
    if (omega_ce <= 0.0) {
        /* Unmagnetized: only Langmuir and EM (cutoff at omega_pe) */
        return (fabs(diel->P) < 0.1) ? WAVE_LANGMUIR : WAVE_NONE;
    }

    double sin2 = 1.0 - cos_theta * cos_theta;

    if (sin2 < 0.01) {
        /* Parallel propagation */
        if (diel->R > 0 && diel->L < 0) return WAVE_WHISTLER;
        if (diel->L > 0 && diel->R < 0) return WAVE_ION_CYCLOTRON;
        if (diel->R > 0 && diel->L > 0) return WAVE_ALFVEN;
        return WAVE_NONE;
    }

    if (sin2 > 0.99) {
        /* Perpendicular propagation */
        if (diel->P > 0) return WAVE_O_MODE;
        if (diel->R * diel->L > 0) return WAVE_X_MODE;
        return WAVE_NONE;
    }

    /* Oblique: general mode identification */
    double n2a, n2b;
    cold_plasma_n2(diel, cos_theta, &n2a, &n2b);

    if (n2b > 0 && n2a > 0) return WAVE_FAST_MAGNETOSONIC;
    if (n2b > 0) return WAVE_WHISTLER;
    return WAVE_NONE;
}

CMARegion cma_classify(double omega_pe, double omega_ce, double omega)
{
    CMARegion r; memset(&r, 0, sizeof(r));
    if (omega <= 0.0) { r.region = 0; return r; }

    r.alpha = omega_pe * omega_pe / (omega * omega);
    r.beta = omega_ce / omega;

    double a = r.alpha, b = r.beta;

    if (a < 1.0 && b < 1.0)        { r.region = 1; strcpy(r.label, "Both slow"); }
    else if (a < 1.0 && b > 1.0)   { r.region = 2; strcpy(r.label, "Whistler zone"); }
    else if (a > 1.0 && b < 1.0)   { r.region = 3; strcpy(r.label, "Evanescent zone"); }
    else if (a > 1.0 && b > 1.0 && b*b < a*(a-1))
        { r.region = 4; strcpy(r.label, "Resonance cone"); }
    else if (a > 1.0 && b > 1.0 && b*b > a) {
        if (b > 2.0*a) { r.region = 5; strcpy(r.label, "Cyclotron dominant"); }
        else           { r.region = 6; strcpy(r.label, "Mixed propagation"); }
    } else {
        r.region = 7; strcpy(r.label, "Complex topology");
    }

    return r;
}

int cold_plasma_is_resonance(double omega, double omega_pe,
                              double omega_ce, double cos_theta)
{
    double sin2 = 1.0 - cos_theta * cos_theta;
    if (sin2 < 0.0) sin2 = 0.0;
    double cos2 = cos_theta * cos_theta;

    double wp2 = omega_pe * omega_pe;
    double wc2 = omega_ce * omega_ce;
    double w2 = omega * omega;
    double denom = w2 - wc2;
    if (fabs(denom) < 1e-30) return 1; /* Cyclotron resonance */

    double S = 1.0 - wp2 / denom;
    double P = 1.0 - wp2 / w2;
    double A = S * sin2 + P * cos2;
    return (fabs(A) < 1e-8) ? 1 : 0;
}

int cold_plasma_cutoffs(double omega_pe, double omega_ce,
                        double cutoffs[3])
{
    int n = 0;
    /* P=0 cutoff: omega = omega_pe */
    cutoffs[n++] = omega_pe;

    /* R=0 cutoff: omega_R = omega_ce/2 + sqrt(omega_pe^2 + omega_ce^2/4) */
    double wc2 = omega_ce * omega_ce * 0.25;
    cutoffs[n++] = omega_ce * 0.5
                   + sqrt(omega_pe * omega_pe + wc2);

    /* L=0 cutoff: omega_L = -omega_ce/2 + sqrt(omega_pe^2 + omega_ce^2/4) */
    double wL = -omega_ce * 0.5 + sqrt(omega_pe * omega_pe + wc2);
    if (wL > 0) cutoffs[n++] = wL;

    return n;
}

/* ===============================================================
 * L6: Whistler Wave
 * =============================================================== */

double whistler_dispersion(double k, double omega_pe, double omega_ce)
{
    if (k <= 0.0 || omega_pe <= 0.0 || omega_ce <= 0.0) return 0.0;

    /* Solve n^2 = 1 + omega_pe^2/(omega*(omega_ce - omega))
     *         = (c*k/omega)^2
     * For omega << omega_ce: omega = omega_ce*(c*k/omega_pe)^2 */
    double kc = k * C_LIGHT;
    double kc2 = kc * kc;
    double wp2 = omega_pe * omega_pe;
    double wc = omega_ce;

    /* Quadratic in omega: omega^2*(kc2+wp2) - omega*wc*kc2 - wp2*wc*omega + ... */
    /* Better: solve n^2 = 1 + wp2/(omega*(wc - omega)) directly via Newton */

    /* Initial guess from low-frequency limit */
    double w = wc * kc2 / (wp2 + kc2);
    if (w <= 0.0) w = 0.1 * wc;

    for (int it = 0; it < 30; it++) {
        double n2 = kc2 / (w * w);
        double RHS = 1.0 + wp2 / (w * (wc - w + 1e-30));
        double f = n2 - RHS;

        if (fabs(f) < 1e-10) break;

        /* Newton step */
        double dn2_dw = -2.0 * kc2 / (w * w * w);
        double dRHS_dw = -wp2 * (wc - 2.0 * w)
                         / (w * w * (wc - w) * (wc - w) + 1e-30);
        double df_dw = dn2_dw - dRHS_dw;

        if (fabs(df_dw) < 1e-30) break;
        double dw = -f / df_dw;
        /* Limit step size */
        if (fabs(dw) > 0.5 * w) dw = (dw > 0) ? 0.5*w : -0.5*w;
        w += dw;
        if (w <= 0.0) w = 0.001 * wc;
    }

    return w;
}

double whistler_group_velocity(double k, double omega_pe, double omega_ce)
{
    if (k <= 0.0) return 0.0;
    double w = whistler_dispersion(k, omega_pe, omega_ce);
    if (w <= 0.0) return 0.0;

    /* For omega << omega_ce: v_g = 2*w/k = 2*omega_ce*c^2*k/omega_pe^2 */
    if (w < 0.3 * omega_ce) return 2.0 * w / k;

    /* General: numerical derivative */
    double dk = k * 0.001;
    double wp = whistler_dispersion(k + dk, omega_pe, omega_ce);
    return (wp - w) / dk;
}

double lower_hybrid_frequency(double omega_ci, double omega_ce,
                               double omega_pi)
{
    if (omega_ci <= 0.0 || omega_ce <= 0.0 || omega_pi <= 0.0) return 0.0;
    /* 1/w_LH^2 = 1/(w_ci*w_ce) + 1/w_pi^2 */
    double inv_w2 = 1.0 / (omega_ci * omega_ce)
                    + 1.0 / (omega_pi * omega_pi);
    return 1.0 / sqrt(inv_w2);
}

double upper_hybrid_frequency(double omega_pe, double omega_ce)
{
    return sqrt(omega_pe * omega_pe + omega_ce * omega_ce);
}

double r_cutoff_frequency(double omega_pe, double omega_ce)
{
    double wc_half = omega_ce * 0.5;
    return wc_half + sqrt(omega_pe * omega_pe + wc_half * wc_half);
}

double l_cutoff_frequency(double omega_pe, double omega_ce)
{
    double wc_half = omega_ce * 0.5;
    double wL = -wc_half + sqrt(omega_pe * omega_pe + wc_half * wc_half);
    return (wL > 0.0) ? wL : 0.0;
}

/* ===============================================================
 * L6: O-mode and X-mode
 * =============================================================== */

double omode_dispersion(double k, double omega_pe)
{
    if (omega_pe <= 0.0) return k * C_LIGHT;
    double n2 = 1.0 - omega_pe * omega_pe
                / (k * k * C_LIGHT * C_LIGHT + omega_pe * omega_pe);
    double w = (k > 0.0) ? k * C_LIGHT / sqrt(fmax(n2, 0.01)) : omega_pe;
    return fmax(w, omega_pe); /* Cutoff at omega_pe */
}

double xmode_dispersion(double k, double omega_pe, double omega_ce,
                         int *is_propagating)
{
    if (is_propagating) *is_propagating = 0;
    if (omega_pe <= 0.0) {
        if (is_propagating) *is_propagating = 1;
        return k * C_LIGHT;
    }
    double w_UH = sqrt(omega_pe * omega_pe + omega_ce * omega_ce);
    /* X-mode dispersion: n^2 = 1 - wp2/w2 * (w2-wp2)/(w2-w_UH^2) */
    /* For given k, solve for omega via Newton */
    if (k <= 0.0) {
        if (is_propagating) *is_propagating = 1;
        return omega_pe;
    }

    /* Initial guess: upper hybrid or above R-cutoff */
    double w = fmax(w_UH, r_cutoff_frequency(omega_pe, omega_ce)) * 1.1;
    double kc2 = k * k * C_LIGHT * C_LIGHT;

    for (int it = 0; it < 40; it++) {
        double w2 = w * w;
        double n2 = kc2 / w2;
        double factor = omega_pe*omega_pe/w2
                        * (w2 - omega_pe*omega_pe)
                        / (w2 - w_UH*w_UH + 1e-30);
        double f = n2 - 1.0 + factor;

        if (fabs(f) < 1e-10) break;

        double df_dw = -2.0*kc2/(w2*w) + 2.0*omega_pe*omega_pe/w2/w;
        if (fabs(df_dw) < 1e-30) break;
        w -= f / df_dw;
        if (w <= 0.0) w = w_UH;
    }

    double n2 = kc2 / (w * w);
    if (is_propagating) *is_propagating = (n2 > 0.0) ? 1 : 0;
    return w;
}

double ion_cyclotron_dispersion(double k_parallel, double omega_pi,
                                 double omega_ci)
{
    if (k_parallel <= 0.0 || omega_ci <= 0.0) return 0.0;
    /* L-wave: n^2 = 1 + omega_pi^2/(omega_ci*(omega_ci - omega))
     * For omega < omega_ci, solve for omega */
    double kc2 = k_parallel * k_parallel * C_LIGHT * C_LIGHT;
    double w = omega_ci * 0.5; /* Initial guess */

    for (int it = 0; it < 30; it++) {
        double n2 = kc2 / (w * w);
        double rhs = 1.0 + omega_pi*omega_pi
                     / (omega_ci * (omega_ci - w + 1e-30));
        double f = n2 - rhs;
        if (fabs(f) < 1e-10) break;
        w *= (1.0 - 0.05 * f);
        if (w <= 0.0) w = omega_ci * 0.1;
        if (w >= omega_ci * 0.99) w = omega_ci * 0.95;
    }
    return w;
}

/* ===============================================================
 * L4: MHD Waves
 * =============================================================== */

double alfven_wave_dispersion(double k_parallel, double v_A)
{
    if (v_A <= 0.0) return 0.0;
    return fabs(k_parallel) * v_A;
}

double fast_magnetosonic_speed(double c_s, double v_A, double cos_theta)
{
    double cs2 = c_s * c_s, va2 = v_A * v_A;
    double cos2 = cos_theta * cos_theta;
    double disc = (cs2 + va2) * (cs2 + va2)
                  - 4.0 * cs2 * va2 * cos2;
    if (disc < 0.0) disc = 0.0;
    return sqrt(0.5 * (cs2 + va2 + sqrt(disc)));
}

double slow_magnetosonic_speed(double c_s, double v_A, double cos_theta)
{
    double cs2 = c_s * c_s, va2 = v_A * v_A;
    double cos2 = cos_theta * cos_theta;
    double disc = (cs2 + va2) * (cs2 + va2)
                  - 4.0 * cs2 * va2 * cos2;
    if (disc < 0.0) disc = 0.0;
    return sqrt(0.5 * (cs2 + va2 - sqrt(disc)));
}

double fast_magnetosonic_dispersion(double k, double c_s, double v_A,
                                     double cos_theta)
{
    return k * fast_magnetosonic_speed(c_s, v_A, cos_theta);
}

double slow_magnetosonic_dispersion(double k, double c_s, double v_A,
                                     double cos_theta)
{
    return k * slow_magnetosonic_speed(c_s, v_A, cos_theta);
}

void mhd_wave_speeds(double c_s, double v_A, double cos_theta,
                     double *v_alfven, double *v_fast, double *v_slow)
{
    if (v_alfven) *v_alfven = v_A;
    if (v_fast) *v_fast = fast_magnetosonic_speed(c_s, v_A, cos_theta);
    if (v_slow) *v_slow = slow_magnetosonic_speed(c_s, v_A, cos_theta);
}

double kinetic_alfven_dispersion(double k_parallel, double k_perp,
                                  double v_A, double rho_i,
                                  double Te_over_Ti)
{
    if (v_A <= 0.0 || k_parallel <= 0.0) return 0.0;
    double rho_eff2 = rho_i * rho_i
                      * (0.75 + Te_over_Ti);
    double kp2 = k_perp * k_perp;
    return fabs(k_parallel) * v_A * sqrt(1.0 + kp2 * rho_eff2);
}
