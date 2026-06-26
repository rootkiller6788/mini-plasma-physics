/*
 * tokamak.c -- Tokamak Grad-Shafranov Equilibrium (Soloviev Solution)
 * Demonstrates axisymmetric toroidal equilibrium and stability analysis.
 * Reference: Wesson "Tokamaks" (2011) Ch.3, Soloviev (1968)
 */
#include "mhd_defs.h"
#include "mhd_equilibrium.h"
#include "mhd_instability.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    puts("=== Tokamak Soloviev Equilibrium ===\n");

    double R0    = 3.0;
    double a     = 1.0;
    double kappa = 1.7;
    double psi0  = 1.0;
    double F0    = 15.0;

    printf("Parameters:\n");
    printf("  Major radius  R0 = %.1f m\n", R0);
    printf("  Minor radius   a = %.1f m\n", a);
    printf("  Elongation  kappa = %.1f\n\n", kappa);

    printf("Poloidal flux on midplane (Z=0):\n");
    printf("%8s %12s %12s %12s %12s\n",
           "R [m]", "psi", "B_R", "B_Z", "B_phi");
    for (int i = 0; i <= 8; i++) {
        double R  = R0 - a + 2.0 * a * i / 8.0;
        double Z  = 0.0;
        double psi = mhd_soloviev_psi(R, Z, R0, kappa, psi0);
        double BR, Bphi, BZ;
        mhd_soloviev_field(R, Z, R0, kappa, psi0, F0, &BR, &Bphi, &BZ);
        printf("%8.2f %12.4e %12.4e %12.4e %12.4e\n",
               R, psi, BR, BZ, Bphi);
    }

    printf("\nSafety factor q(r):\n");
    printf("%8s %12s %12s %12s\n", "r [m]", "B_pol", "B_phi", "q");
    for (int i = 1; i <= 8; i++) {
        double r = a * i / 8.0;
        double R = R0 + r;
        double Z = 0.0;
        double BR, Bphi, BZ;
        mhd_soloviev_field(R, Z, R0, kappa, psi0, F0, &BR, &Bphi, &BZ);
        double Bpol = sqrt(BR * BR + BZ * BZ);
        double q = (Bpol > 1e-10) ? mhd_safety_factor_q(r, R0, Bphi, Bpol) : 0.0;
        printf("%8.3f %12.4e %12.4e %12.4f\n", r, Bpol, Bphi, q);
    }

    /* Kruskal-Shafranov stability check */
    double BR, Bphi, BZ;
    mhd_soloviev_field(R0 + a, 0.0, R0, kappa, psi0, F0, &BR, &Bphi, &BZ);
    double Bpol = sqrt(BR * BR + BZ * BZ);
    double q_a  = mhd_safety_factor_q(a, R0, Bphi, Bpol);
    int ks = mhd_kruskal_shafranov_condition(a, R0, Bphi, Bpol);

    printf("\nKruskal-Shafranov stability:\n");
    printf("  q(a)    = %.4f\n", q_a);
    printf("  Status  = %s\n",
           ks == 1 ? "STABLE (q > 1)" :
           ks == 0 ? "MARGINAL" : "UNSTABLE (q < 1)");

    return 0;
}
