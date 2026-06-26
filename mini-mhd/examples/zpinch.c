/*
 * zpinch.c -- Z-Pinch (Bennett Pinch) Equilibrium Example
 * Solves the Bennett pinch equilibrium and verifies the Bennett relation.
 * Reference: Freidberg (2014) Ch.8
 */
#include "mhd_defs.h"
#include "mhd_equilibrium.h"
#include <stdio.h>
#include <math.h>

static double Jz_parabolic(double r) {
    double rmax = 0.1;
    return 1.0e6 * (1.0 - r * r / (rmax * rmax));
}

int main(void) {
    puts("=== Z-Pinch (Bennett) Equilibrium ===\n");

    int    nr     = 40;
    double r_max  = 0.1;
    double r[40], B_phi[40], p[40];

    mhd_zpinch_equilibrium(Jz_parabolic, r_max, nr, 0.0, r, B_phi, p);

    printf("Radial profiles (parabolic current profile):\n");
    printf("%10s %15s %15s\n", "r [m]", "B_phi [T]", "p [Pa]");
    for (int i = 0; i < nr; i += 5) {
        printf("%10.4f %15.4e %15.4e\n", r[i], B_phi[i], p[i]);
    }

    /* Compute total current and line density */
    double I_total = 0.0;
    double dr = r_max / (double)(nr - 1);
    for (int i = 0; i < nr; i++) {
        I_total += Jz_parabolic(r[i]) * M_PI * r[i] * dr;
    }
    double N_line  = M_PI * r_max * r_max * 1.0e20;
    double T_bennett = mhd_bennett_temperature(I_total, N_line);

    printf("\nBennett relation verification:\n");
    printf("  Total current:   I = %.3e A\n", I_total);
    printf("  Line density:    N = %.3e m^-1\n", N_line);
    printf("  Bennett T:       T = %.3e K\n", T_bennett);

    double lhs = 2.0 * N_line * MHD_KB * T_bennett;
    double rhs = MHD_MU0 * I_total * I_total / (8.0 * M_PI);
    printf("  LHS: 2*N*kB*T  = %.4e\n", lhs);
    printf("  RHS: mu0*I^2/8pi = %.4e\n", rhs);
    printf("  Relative error   = %.3e\n", fabs(lhs - rhs) / lhs);

    /* Plasma beta at center */
    double beta0 = mhd_plasma_beta(p[0], B_phi[nr-1]);
    printf("  Central beta     = %.3e\n", beta0);

    return 0;
}
