/**
 * example3_reactor_comparison.c - Compare ITER/DEMO/SPARC/JET
 *
 * Demonstrates: operating points, triple product,
 *               beta limits, bootstrap current.
 */
#include <stdio.h>
#include <math.h>
#include "../include/fusion_plasma.h"
#include "../include/fusion_confinement.h"
#include "../include/fusion_equilibrium.h"
#include "../include/fusion_mhd.h"

static void analyze(const char *name, const PlasmaParameters *p) {
    printf("\n--- %s ---\n", name);
    printf("  R=%.2f m, a=%.2f m, kappa=%.2f, B=%.1f T, Ip=%.1f MA\n",
           p->R, p->a, p->elongation, p->B, p->Ip/1e6);
    printf("  ne=%.2e, Te=%.1f keV, Vp=%.0f m^3\n",
           p->ne, p->Te/1e3, p->V_p);

    double Ti_keV = p->Ti / 1e3;
    double sv = bosch_hale_sigma_v_dt(Ti_keV);
    double nD = 0.5 * p->ni, nT = 0.5 * p->ni;
    double P_fus = fusion_power_density_dt(nD, nT, sv, E_FUSION_DT_J) * p->V_p;
    printf("  P_fus: %.0f MW\n", P_fus/1e6);

    double tauE = energy_confinement_ipb98y2(p->Ip/1e6, p->B,
        p->ne/1e20, P_fus/1e6, p->R, p->a, p->elongation);
    double tp = triple_product(p->ne/1e20, p->Ti/1e3, tauE);
    printf("  Triple product: %.2f, tau_E: %.2f s\n", tp, tauE);

    double beta = plasma_beta(p->ne, p->Te, p->B) * 100.0;
    double betaN = normalized_beta(beta, p->a, p->B, p->Ip/1e6);
    printf("  beta: %.2f%%, beta_N: %.2f\n", beta, betaN);

    double nG = greenwald_density_limit(p->Ip, p->a);
    printf("  f_GW: %.2f, W_stored: %.0f MJ\n",
           greenwald_fraction(p->ne, nG), plasma_stored_energy(p)/1e6);

    printf("  Ignition: %s\n",
           ignition_condition(p->ne/1e20, p->Ti/1e3, tauE) ? "YES" : "NO");
}

int main(void) {
    printf("=== Fusion Reactor Comparison ===\n");

    PlasmaParameters iter, demo, sparc, jet;
    iter_operating_point(&iter);
    demo_operating_point(&demo);
    sparc_operating_point(&sparc);
    jet_dt_record(&jet);

    analyze("ITER", &iter);
    analyze("DEMO", &demo);
    analyze("SPARC", &sparc);
    analyze("JET D-T Record", &jet);

    printf("\n\n=== Triple Product Summary ===\n");
    printf("%-12s %12s %12s %12s\n", "Reactor", "nTtau", "beta_N", "Ignited");
    printf("%-12s %12s %12s %12s\n", "--------", "------", "------", "-------");

    PlasmaParameters *reactors[] = {&iter, &demo, &sparc, &jet};
    const char *names[] = {"ITER", "DEMO", "SPARC", "JET"};
    for (int i = 0; i < 4; i++) {
        PlasmaParameters *p = reactors[i];
        double tauE2 = energy_confinement_ipb98y2(p->Ip/1e6, p->B,
            p->ne/1e20, 100.0, p->R, p->a, p->elongation);
        double tp2 = triple_product(p->ne/1e20, p->Ti/1e3, tauE2);
        double b2 = plasma_beta(p->ne, p->Te, p->B) * 100.0;
        double bN2 = normalized_beta(b2, p->a, p->B, p->Ip/1e6);
        printf("%-12s %12.2f %12.2f %12s\n",
               names[i], tp2, bN2,
               ignition_condition(p->ne/1e20, p->Ti/1e3, tauE2) ? "YES" : "NO");
    }
    printf("\n=== Complete ===\n");
    return 0;
}