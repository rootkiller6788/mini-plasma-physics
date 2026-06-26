/**
 * example_plasma_waves.c -- Plasma Wave Modes Demo
 *
 * Demonstrates:
 *   1. Dispersion relations for all major wave modes
 *   2. CMA diagram for a range of plasma parameters
 *   3. Whistler wave propagation and group velocity
 *   4. Cutoff and resonance frequencies
 *
 * Reference: Stix "Waves in Plasmas" (1992)
 *            Kivelson & Russell S5
 */
#include "../include/space_plasma.h"
#include "../include/plasma_parameters.h"
#include "../include/plasma_waves.h"
#include <stdio.h>
#include <math.h>

/* Helper: dispersion function wrapper for Newton solver
 * Using Langmuir dispersion as target */
static double D_langmuir(double omega, double k, double theta) {
    (void)theta;
    double omega_pe = 1.78e5;
    double v_the = 1.0e6;
    double omega_bg = langmuir_frequency(k, omega_pe, v_the);
    return omega - omega_bg;
}

int main(void) {
    printf("=== Plasma Wave Modes in Space ===\n\n");

    /* Space plasma parameters (solar wind at 1 AU) */
    double n_e = 5.0e6;
    double B0 = 5.0e-9;
    double T_e = 1.0e5;
    double T_i = 1.0e4;
    double omega_pe = sp_plasma_freq_e(n_e);
    double f_pe = sp_plasma_freq_e_hz(n_e);
    double omega_pi = sp_plasma_freq_i(n_e, 1.0, SP_MP);
    double omega_ce = sp_gyrofreq(B0, SP_ME, SP_EC);
    double omega_ci = sp_gyrofreq(B0, SP_MP, SP_EC);
    double v_A = sp_alfven_speed(B0, n_e * SP_MP);
    double c_s = sp_sound_speed(T_e, SP_MP, 1.0, 5.0/3.0);
    double v_the = SP_VTHERM(T_e, SP_ME);
    double lambda_D = sp_debye_length_e(n_e, T_e);

    printf("Plasma Parameters (1 AU solar wind):\n");
    printf("  n_e = %.1f cm^-3, B = %.1f nT, T_e = %.1e K\n",
           n_e/1e6, B0*1e9, T_e);
    printf("  f_pe = %.1f kHz, f_ce = %.1f Hz\n",
           f_pe/1e3, omega_ce/(2.0*M_PI));
    printf("  f_pi = %.1f Hz, f_ci = %.2f Hz\n",
           omega_pi/(2.0*M_PI), omega_ci/(2.0*M_PI));
    printf("  v_A = %.1f km/s, c_s = %.1f km/s\n",
           v_A/1e3, c_s/1e3);
    printf("  v_the = %.1e m/s, lambda_D = %.1f m\n\n",
           v_the, lambda_D);

    /* Wave mode frequencies at representative k */
    double k = 1.0e-4;  /* 1e-4 rad/m ~ lambda = 63 km */
    double theta = 0.0;  /* parallel propagation */

    printf("Wave Frequencies at k=%.1e rad/m (parallel):\n", k);

    double w_langmuir = langmuir_frequency(k, omega_pe, v_the);
    printf("  Langmuir:        f = %.2f kHz\n", w_langmuir/(2.0*M_PI*1e3));

    double w_alfven = alfven_wave_frequency(k, theta, v_A);
    printf("  Alfven (shear):  f = %.2f Hz\n", w_alfven/(2.0*M_PI));

    double w_fast = fast_magnetosonic_frequency(k, v_A, c_s, theta);
    printf("  Fast magnetosonic: f = %.2f Hz\n", w_fast/(2.0*M_PI));

    double w_slow = slow_magnetosonic_frequency(k, v_A, c_s, theta);
    printf("  Slow magnetosonic: f = %.2f Hz\n", w_slow/(2.0*M_PI));

    double w_ia = ion_acoustic_frequency(k, c_s, lambda_D);
    printf("  Ion acoustic:    f = %.2f Hz\n", w_ia/(2.0*M_PI));

    double w_whistler = whistler_frequency(k, omega_pe, omega_ce, 0.0);
    printf("  Whistler:        f = %.2f Hz\n", w_whistler/(2.0*M_PI));

    double w_lh = lower_hybrid_wave_frequency(omega_ci, omega_ce);
    printf("  Lower hybrid:    f = %.1f Hz\n", w_lh/(2.0*M_PI));

    double w_uh = sp_upper_hybrid_freq(omega_pe, omega_ce);
    printf("  Upper hybrid:    f = %.2f kHz\n\n", w_uh/(2.0*M_PI*1e3));

    /* Whistler dispersion scan */
    printf("Whistler Wave Dispersion (theta=0):\n");
    printf("  %-12s %-12s %-12s %-12s\n", "k [1/m]", "f [Hz]", "v_ph [km/s]", "v_g [km/s]");

    double k_vals[] = {1e-5, 3e-5, 1e-4, 3e-4, 1e-3};
    for (int i = 0; i < 5; i++) {
        double ki = k_vals[i];
        double w = whistler_frequency(ki, omega_pe, omega_ce, 0.0);
        double v_phase = w / ki;
        double v_g_par, v_g_perp;

        /* Compute group velocity numerically */
        group_velocity(alfven_wave_frequency, ki, 0.0, &v_g_par, &v_g_perp, ki*1e-5, 1e-6);
        /* For whistler specifically, use the whistler function */
        double w_plus  = whistler_frequency(ki*1.001, omega_pe, omega_ce, 0.0);
        double w_minus = whistler_frequency(ki*0.999, omega_pe, omega_ce, 0.0);
        double v_g_w = (ki > 0) ? (w_plus - w_minus) / (ki*0.002) : 0.0;

        printf("  %-12.1e %-12.1f %-12.1f %-12.1f\n",
               ki, w/(2.0*M_PI), v_phase/1e3, v_g_w/1e3);
    }
    printf("\n");

    /* CMA Diagram Points */
    printf("CMA Diagram Classification:\n");
    printf("  %-8s %-8s %-8s %-8s %-6s %-6s %-6s %-6s\n",
           "omega", "alpha", "beta", "bands", "R", "L", "X", "O");

    double omegas[] = {1e3, 1e4, 1e5, 1e6, 1e7};
    for (int i = 0; i < 5; i++) {
        double omega = omegas[i];
        double alpha = (omega_pe*omega_pe) / (omega*omega);
        double beta = omega_ce / omega;
        int bands = cma_propagation_bands(alpha, beta);

        printf("  %-8.0e %-8.2f %-8.3f %-8d",
               omega, alpha, beta, bands);
        printf(" %-6s %-6s %-6s %-6s\n",
               (bands & 1) ? "YES" : "no",
               (bands & 2) ? "YES" : "no",
               (bands & 4) ? "YES" : "no",
               (bands & 8) ? "YES" : "no");
    }
    printf("\n");

    /* Cutoffs and Resonances */
    printf("Cutoff Frequencies:\n");
    double cutoffs[3], resonances[3];
    wave_cutoff_frequencies(omega_pe, omega_ce, cutoffs);
    wave_resonance_frequencies(omega_pe, omega_ce, omega_pi, omega_ci, M_PI/2.0, resonances);

    printf("  P-cutoff (O-mode):  f = %.2f kHz\n", cutoffs[0]/(2.0*M_PI*1e3));
    printf("  R-cutoff:           f = %.2f kHz\n", cutoffs[1]/(2.0*M_PI*1e3));
    printf("  L-cutoff:           f = %.2f kHz\n", cutoffs[2]/(2.0*M_PI*1e3));
    printf("  Electron cyclotron resonance: f = %.1f Hz\n", resonances[0]/(2.0*M_PI));
    printf("  Lower hybrid resonance:      f = %.1f Hz\n", resonances[1]/(2.0*M_PI));
    printf("  Upper hybrid resonance:      f = %.2f kHz\n\n", resonances[2]/(2.0*M_PI*1e3));

    /* Stix tensor at whistler frequency */
    printf("Stix Dielectric Tensor at f=100 Hz (whistler band):\n");
    stix_tensor_t eps;
    stix_dielectric_tensor(2.0*M_PI*100.0, omega_pe, omega_ce, omega_pi, omega_ci, &eps);
    printf("  S = %.6e, D = %.6e, P = %.6e\n", eps.S, eps.D, eps.P);
    printf("  R = %.6e, L = %.6e\n\n", eps.R, eps.L);

    /* Newton solver test */
    printf("Dispersion Solver (Newton-Raphson):\n");
    double omega_solved = solve_dispersion_omega(D_langmuir, 1.0, 0.0,
                                                  2.0e5, 1e-3, 50);
    if (omega_solved > 0.0) {
        printf("  Solved omega = %.2f krad/s (vs %.2f krad/s from formula)\n",
               omega_solved/1e3, w_langmuir/1e3);
    } else {
        printf("  Newton solver did not converge\n");
    }

    /* Warm plasma dielectric at Langmuir resonance */
    double eps_warm = warm_plasma_dielectric_function(w_langmuir, k,
                                                       omega_pe, v_the);
    printf("  Warm plasma epsilon(omega_pe, k) = %.6e (should be ~0)\n\n",
           eps_warm);

    printf("=== Done ===\n");
    return 0;
}
