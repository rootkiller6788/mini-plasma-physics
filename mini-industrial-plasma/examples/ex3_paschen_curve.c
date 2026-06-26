#include "discharge.h"
#include <stdio.h>
#include <stdlib.h>

int main(void) {
    printf("=== Example 3: Paschen Curve for Argon ===\n");
    double A = 13.6, B = 235.0, gamma = 0.1;
    double Vmin = paschen_minimum_voltage(A, B, gamma);
    double pdmin = paschen_pd_at_minimum(A, B, gamma);
    printf("Gas: Argon\n");
    printf("Paschen minimum: V_min = %.1f V at pd = %.3f Pa.m\n", Vmin, pdmin);
    printf("\np.d [Pa.m]    V_breakdown [V]\n");
    printf("---------    ----------------\n");
    double pd_vals[] = {0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0};
    for (int i = 0; i < 10; i++) {
        double Vbr = paschen_breakdown_voltage(pd_vals[i], A, B, gamma);
        printf("  %6.2f        %8.1f\n", pd_vals[i], Vbr);
    }
    double Vbr_typical = paschen_breakdown_voltage(1.0, A, B, gamma);
    printf("\nAt typical RIE conditions (1 Pa.m): V_br = %.0f V\n", Vbr_typical);
    printf("\nLowest breakdown near pd_min -> easiest ignition\n");
    return 0;
}
