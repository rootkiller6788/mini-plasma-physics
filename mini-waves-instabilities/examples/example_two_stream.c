/**
 * example_two_stream.c -- Two-Stream Instability Analysis
 *
 * Demonstrates the cold two-stream instability growth rate
 * as a function of beam density and wavenumber.
 *
 * L7 Application: Fast ignition fusion, electron beam transport
 *                  in inertial confinement fusion (NIF, LMJ)
 */
#include "../include/waves_instabilities.h"
#include "../include/plasma_instabilities.h"
#include <stdio.h>
#include <math.h>

int main(void)
{
    printf("=== Two-Stream Instability Analysis ===\n\n");

    double n_plasma = 1e26;     /* Solid density plasma [m^-3] */
    double omega_pe = electron_plasma_frequency(n_plasma);
    double v_b = 2.0e8;        /* Relativistic beam [m/s] */

    printf("Plasma: n = %.2e m^-3, omega_pe = %.3e rad/s\n",
           n_plasma, omega_pe);
    printf("Beam velocity: v_b = %.2e m/s (%.2f c)\n\n",
           v_b, v_b / C_LIGHT);

    printf("Growth rate vs beam density ratio:\n");
    printf("  %12s %16s %16s\n",
           "n_b/n_0", "gamma [1/s]", "gamma/omega_pe");

    for (int i = -4; i <= 0; i++) {
        double ratio = pow(10.0, i);
        double gamma = two_stream_growth_rate(ratio * n_plasma,
                                               n_plasma, omega_pe);
        printf("  %12.2e %16.6e %16.6f\n", ratio, gamma, gamma/omega_pe);
    }

    printf("\nGrowth rate vs wavenumber for n_b/n_0 = 0.01:\n");
    double omega_b = electron_plasma_frequency(0.01 * n_plasma);
    printf("  %12s %16s\n", "k [1/m]", "gamma [1/s]");

    for (int i = 0; i <= 10; i++) {
        double k = omega_pe / v_b * (0.1 + i * 0.2);
        double gamma = two_stream_growth_at_k(k, omega_pe, omega_b, v_b);
        printf("  %12.3e %16.6e\n", k, gamma);
    }

    printf("\nPhysical interpretation:\n");
    printf("  - Maximum growth at k ~ omega_pe/v_b\n");
    printf("  - gamma_max ~ (n_b/n_0)^(1/3) * omega_pe\n");
    printf("  - Saturation via quasilinear plateau formation\n");
    printf("  - Relevant to fast electron transport in ICF\n");

    printf("\n=== Done ===\n");
    return 0;
}
