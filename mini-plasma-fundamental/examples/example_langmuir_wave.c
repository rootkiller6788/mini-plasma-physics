/**
 * example_langmuir_wave.c — Langmuir Wave Dispersion
 *
 * Demonstrates the Bohm-Gross dispersion relation for Langmuir waves
 * and the Landau damping rate.
 */
#include "plasma_params.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== Langmuir Wave Dispersion (Bohm-Gross) ===\n\n");

    double ne = 1.0e20;      /* m^-3 */
    double Te = 1.16e8;      /* K (10 keV) */
    (void)M_PROTON; /* reference mass for context */

    double omega_pe = electron_plasma_frequency(ne);
    double f_pe = plasma_frequency_hz(ne, M_ELECTRON);
    double lambda_D = debye_length(Te, ne);
    double v_the = electron_thermal_velocity(Te);

    printf("Plasma parameters:\n");
    printf("  n_e = %.2e m^-3\n", ne);
    printf("  T_e = %.2e K = %.1f eV\n", Te, Te * K_B / E_CHARGE);
    printf("  omega_pe = %.3e rad/s, f_pe = %.3e Hz\n", omega_pe, f_pe);
    printf("  lambda_D = %.3e m\n", lambda_D);
    printf("  v_th,e = %.3e m/s\n\n", v_the);

    printf("  k [m^-1]        k*lambda_D    omega [rad/s]    "
           "f [Hz]          v_phi [m/s]     gamma/omega\n");
    printf("  --------------  ------------  --------------  "
           "-------------  -------------  ------------\n");

    for (int i = -2; i <= 2; i++) {
        double k_lD = pow(10.0, i);
        double k = k_lD / lambda_D;

        /* Bohm-Gross dispersion */
        double omega = sqrt(omega_pe * omega_pe
                            + 3.0 * (K_B * Te / M_ELECTRON) * k * k);
        double f = omega / (2.0 * M_PI);
        double v_phi = omega / k;

        /* Landau damping rate (simplified model) */
        double k_lD_cube = k_lD * k_lD * k_lD;
        double gamma = 0.0;
        if (k_lD > 0.05) {
            gamma = -sqrt(M_PI / 8.0) * omega_pe / k_lD_cube
                    * exp(-0.5 / (k_lD * k_lD) - 1.5);
        }
        double gamma_over_omega = fabs(gamma / omega);

        printf("  %-14.3e  %-12.3f  %-14.3e  %-13.3e  %-13.3e  %-12.3e\n",
               k, k_lD, omega, f, v_phi, gamma_over_omega);
    }

    printf("\nInterpretation:\n");
    printf("  k*lambda_D << 1: omega ~ omega_pe (no damping)\n");
    printf("  k*lambda_D ~ 1: Strong Landau damping\n");
    printf("  k*lambda_D >> 1: omega ~ k*v_th (thermal propagation)\n");

    return 0;
}
