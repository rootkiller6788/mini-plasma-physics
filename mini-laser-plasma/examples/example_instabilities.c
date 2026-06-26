/*
 * example_instabilities.c -- Parametric Instabilities Analysis
 *
 * Demonstrates growth rate and threshold calculations for SRS, SBS,
 * and TPD in ICF-relevant plasma conditions.
 *
 * Build: gcc -std=c11 -I../include example_instabilities.c
 *        ../src/plasma_params.c ../src/instabilities.c -lm
 *        -o example_instabilities
 */

#include <stdio.h>
#include <math.h>
#include "plasma_constants.h"
#include "plasma_params.h"
#include "laser_plasma.h"
#include "instabilities.h"

int main(void)
{
    printf("=== Laser-Plasma Parametric Instabilities ===\n\n");

    /* ICF-relevant plasma: ne = 0.1 nc, Te = 3 keV, L_n = 100 um */
    double lambda0 = 0.351e-6;  /* 3 omega Nd:glass */
    double ne = 0.1 * critical_density(lambda0);

    PlasmaState ps = {
        ne,
        3000.0,    /* Te = 3 keV */
        1000.0,    /* Ti = 1 keV */
        3.5,       /* Z */
        12.0,      /* A */
        0.0,       /* B */
        lambda0,
        1.0e19     /* I = 10^15 W/cm^2 */
    };
    PlasmaDerived pd;
    compute_all_derived(&ps, &pd);
    pd.a0 = normalized_vector_potential(ps.I_laser, lambda0);

    printf("Plasma conditions:\n");
    printf("  ne/nc = %.4f\n", pd.ne_over_nc);
    printf("  Te    = %.0f eV\n", ps.Te);
    printf("  a0    = %.6f\n", pd.a0);
    printf("  v_osc/c = %.6f\n", pd.a0);
    printf("  L_n   = 100 um (assumed)\n\n");

    double L_n = 100e-6;
    double v_osc = pd.a0 * PLASMA_C;

    /* Analyze each instability type */
    InstabilityType types[] = {
        INSTABILITY_SRS_BACKSCATTER,
        INSTABILITY_SRS_FORWARD,
        INSTABILITY_SBS_BACKSCATTER,
        INSTABILITY_TPD
    };

    for (int i = 0; i < 4; i++) {
        InstabilityResult result;
        int ret = analyze_instability(types[i], &ps, &pd,
                                       v_osc, L_n, &result);

        printf("[%s]\n", instability_name(types[i]));
        if (ret != 0) {
            printf("  Not active in these conditions.\n\n");
            continue;
        }

        printf("  Growth rate:   gamma0 = %.3e s^-1 ", result.gamma0);
        printf("(%.1f fs)\n", 1.0 / result.gamma0 * 1e15);

        if (result.rosenbluth_gain > 0.0) {
            printf("  Rosenbluth gain: G = %.3f", result.rosenbluth_gain);
            if (result.rosenbluth_gain > 2.0 * M_PI) {
                printf(" (SIGNIFICANT: G > 2pi)");
            }
            printf("\n");
        }

        if (types[i] == INSTABILITY_SRS_BACKSCATTER ||
            types[i] == INSTABILITY_SRS_FORWARD) {
            printf("  EPW wavenumber: k_epw = %.3e rad/m\n", result.k_epw);
            printf("  EPW frequency:  omega_epw = %.3e rad/s\n",
                   result.omega_epw);
        }

        if (types[i] == INSTABILITY_SBS_BACKSCATTER) {
            printf("  IAW wavenumber: k_iaw = %.3e rad/m\n", result.k_iaw);
            printf("  IAW frequency:  omega_iaw = %.3e rad/s\n",
                   result.omega_iaw);
        }

        printf("\n");
    }

    /* Compare with higher intensity */
    printf("=== At higher intensity (I = 10^17 W/cm^2) ===\n\n");
    ps.I_laser = 1.0e21;
    pd.a0 = normalized_vector_potential(ps.I_laser, lambda0);
    v_osc = pd.a0 * PLASMA_C;
    printf("  a0 = %.4f, v_osc/c = %.4f\n\n", pd.a0, pd.a0);

    for (int i = 0; i < 4; i++) {
        InstabilityResult result;
        int ret = analyze_instability(types[i], &ps, &pd,
                                       v_osc, L_n, &result);
        if (ret == 0) {
            printf("[%s] gamma0 = %.3e s^-1 (%.1f fs)\n",
                   instability_name(types[i]),
                   result.gamma0, 1.0 / fmax(result.gamma0, 1e-30) * 1e15);
        }
    }

    return 0;
}
