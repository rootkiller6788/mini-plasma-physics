/*
 * example_wakefield.c -- Laser Wakefield Acceleration (LWFA) example
 *
 * Demonstrates the calculation of LWFA operating parameters for
 * a state-of-the-art GeV-class accelerator.
 *
 * Build: gcc -std=c11 -I../include example_wakefield.c
 *        ../src/plasma_params.c ../src/wakefield.c -lm -o example_wakefield
 */

#include <stdio.h>
#include <math.h>
#include "plasma_constants.h"
#include "plasma_params.h"
#include "wakefield.h"

int main(void)
{
    printf("=== Laser Wakefield Acceleration (LWFA) Calculator ===\n\n");

    /* Three representative LWFA regimes */
    double densities_cm3[] = {1.0e17, 1.0e18, 1.0e19};
    const char *labels[] = {
        "Low-density (quasi-linear)",
        "Mid-density (nonlinear)",
        "High-density (bubble)"
    };
    double a0_vals[] = {0.5, 1.5, 3.0};
    double lambda0 = 0.8e-6;  /* Ti:Sapphire */
    double tau_fwhm = 30e-15; /* 30 fs */

    for (int i = 0; i < 3; i++) {
        double ne = densities_cm3[i] * 1e6;  /* to m^-3 */
        double a0 = a0_vals[i];
        double omega = 2.0 * M_PI * PLASMA_C / lambda0;
        double wp = plasma_frequency(ne);
        double nc = critical_density(lambda0);
        double lambda_p = plasma_wavelength(wp);

        printf("--- %s ---\n", labels[i]);
        printf("  ne = %.1e cm^-3\n", densities_cm3[i]);
        printf("  a0 = %.1f\n", a0);
        printf("  nc/ne = %.1f\n", nc / ne);
        printf("  Plasma wavelength: lambda_p = %.1f um\n", lambda_p * 1e6);
        printf("  Plasma period:     tau_p   = %.1f fs\n",
               2.0 * M_PI / wp * 1e15);

        double E_wb = cold_wavebreaking_field(wp);
        printf("  Wave-breaking field: E_wb = %.1f GV/m\n", E_wb * 1e-9);

        double L_deph = dephasing_length_1D(wp, omega, nc, ne);
        printf("  Dephasing length: L_deph = %.3f mm\n", L_deph * 1e3);

        double L_pump = pump_depletion_length(omega, wp, a0, tau_fwhm);
        printf("  Pump depletion:   L_pump = %.3f mm\n", L_pump * 1e3);

        /* Operating accelerating field */
        double E_acc;
        if (a0 < 0.5) {
            E_acc = wakefield_gradient_1D_linear(a0, E_wb);
        } else if (a0 < 2.0) {
            E_acc = wakefield_gradient_1D_nonlinear(a0, E_wb);
        } else {
            E_acc = bubble_accelerating_field(wp);
        }
        printf("  Accelerating field: E_acc = %.1f GV/m\n", E_acc * 1e-9);

        double L_acc = (L_deph < L_pump) ? L_deph : L_pump;
        double W_max = maximum_energy_gain(E_acc, L_acc);
        printf("  Max energy gain:  W_max = %.1f MeV (%.3f GeV)\n\n",
               W_max * 1e-6, W_max * 1e-9);
    }

    return 0;
}
