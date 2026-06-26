/**
 * example1_iter_power_balance.c - ITER power balance calculation
 *
 * Demonstrates: Lawson criterion, fusion power density,
 *               bremsstrahlung radiation, energy confinement,
 *               ignition condition, fusion gain Q.
 */
#include <stdio.h>
#include <math.h>
#include "../include/fusion_plasma.h"
#include "../include/fusion_confinement.h"
#include "../include/fusion_equilibrium.h"
#include "../include/fusion_mhd.h"

int main(void) {
    printf("=== ITER Power Balance Analysis ===\n\n");
    PlasmaParameters iter;
    iter_operating_point(&iter);
    printf("ITER Baseline: R=%.2f m, a=%.2f m, B=%.1f T, Ip=%.1f MA\n",
           iter.R, iter.a, iter.B, iter.Ip/1e6);
    printf("ne=%.2e m^-3, Te=%.1f keV, Ti=%.1f keV, Vp=%.0f m^3\n\n",
           iter.ne, iter.Te/1e3, iter.Ti/1e3, iter.V_p);

    double Ti_keV = iter.Ti / 1e3;
    double sigma_v = bosch_hale_sigma_v_dt(Ti_keV);
    printf("D-T reactivity at Ti=%.1f keV: %.2e m^3/s\n", Ti_keV, sigma_v);

    double nD = 0.5 * iter.ni;
    double nT = 0.5 * iter.ni;
    double P_fus_dens = fusion_power_density_dt(nD, nT, sigma_v, E_FUSION_DT_J);
    double P_fus = P_fus_dens * iter.V_p;
    printf("Fusion power: %.0f MW\n", P_fus/1e6);

    double P_alpha = alpha_power_density(nD, nT, sigma_v) * iter.V_p;
    printf("Alpha heating: %.0f MW\n", P_alpha/1e6);

    double P_brem = bremsstrahlung_power_density(iter.ne, iter.Zeff, iter.Te) * iter.V_p;
    printf("Bremsstrahlung: %.0f MW\n", P_brem/1e6);

    double P_cyc = cyclotron_power_density(iter.ne, iter.Te, iter.B) * iter.V_p;
    printf("Cyclotron: %.1f MW\n", P_cyc/1e6);

    double tauE = energy_confinement_ipb98y2(iter.Ip/1e6, iter.B,
        iter.ne/1e20, P_fus/1e6, iter.R, iter.a, iter.elongation);
    printf("\nIPB98(y,2) tau_E: %.2f s\n", tauE);

    double W = plasma_stored_energy(&iter);
    printf("Stored energy: %.0f MJ\n", W/1e6);

    double n_20 = iter.ne / 1e20;
    double T_keV = iter.Te / 1e3;
    double tp = triple_product(n_20, T_keV, tauE);
    printf("Triple product: %.2f x10^20 m^-3 keV s (ignition: 3.0)\n", tp);
    printf("Ignition condition: %s\n", ignition_condition(n_20, T_keV, tauE) ? "MET" : "NOT MET");

    double Q = fusion_gain_Q(P_fus, P_fus/10.0);
    printf("Fusion gain Q: %.1f (ITER target: 10)\n", Q);

    double nG = greenwald_density_limit(iter.Ip, iter.a);
    printf("Greenwald fraction: %.2f\n", greenwald_fraction(iter.ne, nG));

    double beta = plasma_beta(iter.ne, iter.Te, iter.B) * 100.0;
    double betaN = normalized_beta(beta, iter.a, iter.B, iter.Ip/1e6);
    printf("Beta: %.2f%%, beta_N: %.2f\n", beta, betaN);

    double B_pol = edge_poloidal_field(iter.Ip, iter.a, iter.elongation);
    double beta_p = plasma_beta_poloidal(iter.ne, iter.Te, B_pol);
    double eps = 1.0 / aspect_ratio(iter.R, iter.a);
    printf("Bootstrap fraction: %.1f%%\n", bootstrap_current_fraction(eps, beta_p)*100.0);

    printf("\n=== Analysis Complete ===\n");
    return 0;
}