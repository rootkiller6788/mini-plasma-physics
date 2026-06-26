/**
 * example_mhd_stability.c -- MHD Instability Analysis for Tokamak
 *
 * Demonstrates kink, ballooning, tearing, and RT instabilities
 * at tokamak-relevant parameters.
 *
 * L7 Application: ITER stability analysis, ELM control
 */
#include "../include/waves_instabilities.h"
#include "../include/plasma_instabilities.h"
#include <stdio.h>
#include <math.h>

int main(void)
{
    printf("=== MHD Stability Analysis ===\n\n");

    /* ITER-like parameters */
    double B = 5.3;               /* Toroidal field [T] */
    double n = 1e20;              /* Density [m^-3] */
    double T = 10.0 * 1.602e-19;  /* 10 keV [J] */
    double a = 2.0;               /* Minor radius [m] */
    double R = 6.2;               /* Major radius [m] */
    double q95 = 3.0;             /* Safety factor */

    printf("ITER-like parameters:\n");
    printf("  B = %.1f T, n = %.1e m^-3, T = %.0f keV\n",
           B, n, T / 1.602e-19 / 1000.0);
    printf("  a = %.1f m, R = %.1f m, q95 = %.1f\n\n", a, R, q95);

    /* Kink stability */
    printf("Kink stability (Kruskal-Shafranov):\n");
    printf("  q(edge) = %.2f -> (1,1) kink: %s\n",
           q95, kink_unstable(q95, 1, 1) ? "UNSTABLE" : "stable");
    printf("  q(edge) = %.2f -> (2,1) kink: %s\n",
           q95, kink_unstable(q95, 2, 1) ? "UNSTABLE" : "stable");

    /* Ballooning stability */
    double cs = sqrt(K_BOLTZMANN * T / M_PROTON * 2.0);
    double grad_p_over_p = 1.0 / a; /* ~1/a gradient scale */
    double gamma_ball = ballooning_growth_rate(cs, grad_p_over_p,
                                                 R, q95);
    printf("\nBallooning stability:\n");
    printf("  cs = %.3e m/s\n", cs);
    printf("  gamma_ballooning = %.3e 1/s\n", gamma_ball);
    printf("  Ballooning: %s\n",
           gamma_ball > 0 ? "UNSTABLE" : "stable");

    /* Tearing mode */
    double eta = 1e-8; /* Spitzer resistivity at 10 keV */
    double v_A = B / sqrt(MU0 * n * M_PROTON);
    double delta_prime = 10.0;
    double gamma_tear = tearing_mode_growth_rate(eta, a, v_A,
                                                   delta_prime);
    printf("\nTearing mode (FKR):\n");
    printf("  v_A = %.3e m/s\n", v_A);
    printf("  eta = %.1e ohm-m\n", eta);
    printf("  gamma_tearing = %.3e 1/s\n", gamma_tear);
    printf("  Growth time = %.3e s\n", 1.0 / fmax(gamma_tear, 1e-30));

    printf("\n=== Done ===\n");
    return 0;
}
