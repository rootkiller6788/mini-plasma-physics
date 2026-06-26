/**
 * example_langmuir_dispersion.c -- Langmuir Wave Dispersion and Landau Damping
 *
 * Demonstrates Bohm-Gross dispersion relation and weak Landau damping
 * for electron plasma waves at fusion-relevant parameters.
 *
 * L7 Application: Tokamak reflectometry diagnostics (ITER)
 *
 * Build: gcc -I../include -L.. -lwaves -lm -o example_langmuir_dispersion.exe
 *        example_langmuir_dispersion.c
 */
#include "../include/waves_instabilities.h"
#include <stdio.h>
#include <math.h>

int main(void)
{
    printf("=== Langmuir Wave Dispersion and Landau Damping ===\n\n");

    /* ITER-relevant parameters */
    double n_e = 1e20;          /* Density [m^-3] */
    double T_e = 10.0 * 1.602e-19; /* Temperature [J] (~10 keV) */
    double omega_pe = electron_plasma_frequency(n_e);
    double v_th_e = sqrt(2.0 * K_BOLTZMANN * T_e / M_ELECTRON);

    printf("Plasma parameters:\n");
    printf("  n_e = %.2e m^-3\n", n_e);
    printf("  T_e = %.1f keV\n", T_e / 1.602e-19 / 1000.0);
    printf("  omega_pe = %.3e rad/s (f_pe = %.3f GHz)\n",
           omega_pe, omega_pe / (2.0 * PLASMA_PI) / 1e9);
    printf("  v_th_e = %.3e m/s\n", v_th_e);

    printf("\nDispersion curve omega(k):\n");
    printf("  %12s %16s %16s %16s\n",
           "k [1/m]", "omega [rad/s]", "omega/omega_pe", "gamma/omega_pe");

    for (int i = 0; i <= 10; i++) {
        double k = pow(10.0, i - 3.0); /* k from 0.001 to 1e7 */
        ComplexOmega w = langmuir_wave_damped(k, omega_pe, v_th_e);

        printf("  %12.3e %16.6e %16.6f %16.3e\n",
               k, w.omega_r, w.omega_r / omega_pe, w.gamma / omega_pe);
    }

    printf("\nPhysical interpretation:\n");
    printf("  - For k*lambda_De << 1: omega ~ omega_pe (cold plasma limit)\n");
    printf("  - Thermal correction: omega^2 = omega_pe^2 + 3*v_th^2*k^2\n");
    printf("  - Landau damping: negligible for k*lambda_De << 1\n");
    printf("  - Strong damping: onset when k*lambda_De > 0.3\n");

    printf("\n=== Done ===\n");
    return 0;
}
