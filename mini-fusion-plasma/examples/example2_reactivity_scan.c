/**
 * example2_reactivity_scan.c - D-T/D-D/D-He3 reactivity vs temperature
 *
 * Demonstrates: Bosch-Hale cross-section, Lawson criterion
 *               temperature optimization.
 */
#include <stdio.h>
#include <math.h>
#include "../include/fusion_plasma.h"
#include "../include/fusion_equilibrium.h"

int main(void) {
    printf("=== Fusion Reactivity Temperature Scan ===\n\n");
    printf("%8s %15s %15s %15s\n", "T[keV]", "<sv>_DT", "<sv>_DD", "<sv>_DHe3");
    printf("%8s %15s %15s %15s\n", "------", "------", "------", "------");

    double T_vals[] = {0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 30.0, 50.0, 80.0, 100.0, 150.0, 200.0};
    int n = sizeof(T_vals)/sizeof(T_vals[0]);
    double sv_max = 0.0, T_max = 0.0;

    for (int i = 0; i < n; i++) {
        double T = T_vals[i];
        double sv_dt = bosch_hale_sigma_v_dt(T);
        double sv_dd = bosch_hale_sigma_v_dd(T);
        double sv_dhe3 = bosch_hale_sigma_v_dhe3(T);
        if (sv_dt > sv_max) { sv_max = sv_dt; T_max = T; }
        printf("%8.1f %15.3e %15.3e %15.3e\n", T, sv_dt, sv_dd, sv_dhe3);
    }

    printf("\n--- Key Results ---\n");
    printf("Peak D-T reactivity: %.3e m^3/s at %.0f keV\n", sv_max, T_max);
    printf("Optimal fusion temperature: %.0f keV\n", optimal_fusion_temperature_dt());

    printf("\n--- Lawson Criterion vs T ---\n");
    printf("%8s %15s %15s\n", "T[keV]", "n*tau_E[min]", "n*T*tau_E");
    printf("%8s %15s %15s\n", "------", "------", "------");
    for (int i = 0; i < n; i++) {
        double T = T_vals[i];
        if (T < 5.0) continue;
        double sv = bosch_hale_sigma_v_dt(T);
        double E_alpha_J = E_ALPHA * E_CHARGE;
        double T_J = T * 1e3 * E_CHARGE;
        double n_tau_min = 12.0 * K_BOLTZMANN * T_J / (sv * E_alpha_J);
        printf("%8.1f %15.2e %15.2f\n", T, n_tau_min, n_tau_min * T / 1e3);
    }
    printf("\n=== Scan Complete ===\n");
    return 0;
}