/**
 * example_whistler_propagation.c -- Whistler Wave in Earth's Magnetosphere
 *
 * Demonstrates whistler wave dispersion and group velocity
 * for ionospheric/space plasma conditions.
 *
 * L7 Application: Lightning-induced whistlers, Van Allen belt wave-particle
 *                  interactions (NASA Van Allen Probes, GPS scintillation)
 */
#include "../include/waves_instabilities.h"
#include <stdio.h>
#include <math.h>

int main(void)
{
    printf("=== Whistler Wave Propagation in Magnetosphere ===\n\n");

    /* Magnetospheric parameters (L=4, equatorial) */
    double n_e = 1e8;              /* Equatorial density [m^-3] */
    double B0 = 5e-7;              /* Magnetic field [T] (~500 nT at L=4) */
    double omega_pe = electron_plasma_frequency(n_e);
    double omega_ce = cyclotron_frequency(-E_CHARGE, B0, M_ELECTRON);

    printf("Magnetospheric plasma (L=4 equator):\n");
    printf("  n_e = %.1e m^-3\n", n_e);
    printf("  B0 = %.1e T\n", B0);
    printf("  f_pe = %.1f kHz\n", omega_pe / (2.0 * PLASMA_PI) / 1e3);
    printf("  f_ce = %.1f kHz\n", omega_ce / (2.0 * PLASMA_PI) / 1e3);

    printf("\nWhistler dispersion curve:\n");
    printf("  %12s %14s %14s %14s %14s\n",
           "f [kHz]", "k [1/km]", "n=c*k/w", "v_g/c", "t(L=4) [s]");

    double L_shell = 4.0;
    double path_length = L_shell * 6371e3 * 2.0; /* Rough field line length */

    for (int i = 1; i <= 10; i++) {
        double f = 0.5 * i * 1e3;  /* 0.5 - 5 kHz */
        double omega = 2.0 * PLASMA_PI * f;
        double k = omega / (C_LIGHT * 0.1); /* Initial guess for k */
        /* Iterate to improve k estimate */
        for (int j = 0; j < 5; j++) {
            k = omega / C_LIGHT * sqrt(1.0 + omega_pe*omega_pe
                / (omega * (omega_ce - omega + 1e-30)));
        }
        double v_g = whistler_group_velocity(k, omega_pe, omega_ce);
        double travel_time = path_length / v_g;
        double n_ref = C_LIGHT * k / omega;

        printf("  %12.1f %14.3e %14.3f %14.6f %14.3f\n",
               f/1e3, k/1e3, n_ref, v_g/C_LIGHT, travel_time);
    }

    printf("\nPhysical interpretation:\n");
    printf("  - Higher frequencies travel faster (v_g increases with f)\n");
    printf("  - This produces the characteristic descending tone\n");
    printf("  - Lightning strike -> broad spectrum -> dispersion ->\n");
    printf("    high freq arrives first at conjugate point\n");
    printf("  - Used as diagnostic for plasmasphere density (Storey, 1953)\n");
    printf("  - GPS: ionospheric scintillation from whistler turbulence\n");

    printf("\n=== Done ===\n");
    return 0;
}
