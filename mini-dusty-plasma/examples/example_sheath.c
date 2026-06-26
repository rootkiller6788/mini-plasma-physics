/**
 * @file  example_sheath.c
 * @brief Example: Dust levitation in plasma sheath.
 *
 * Computes the equilibrium levitation height of a dust grain
 * in an RF discharge sheath. Balances gravity, electric force,
 * and ion drag to find the stable equilibrium point.
 *
 * L7: Application — Dust in plasma processing (semiconductor manufacturing).
 * L7: Application — Dust diagnostics in tokamaks (Fusion).
 * L7: Application — Saturn's ring dusty plasma analog (Space plasma).
 */

#include "dusty_plasma.h"
#include "dusty_forces.h"
#include "dusty_transport.h"
#include <math.h>
#include <stdio.h>

/* Simple test macro for inline checks */
#define CHECK(cond, msg) do { \
    if (!(cond)) printf("  WARNING: %s\n", msg); \
} while(0)

int main(void)
{
    printf("\n=== Dust Levitation in Plasma Sheath ===\n\n");

    /* Grain parameters */
    double a = 5.0e-6;
    double rho = DUSTY_RHO_MF;    /* melamine-formaldehyde */
    double m_d = dust_grain_mass(a, rho);

    /* Plasma parameters */
    double n_e = 1.0e15;
    double T_e = 3.0 * DUSTY_EV_IN_K;
    double n_d = 1.0e10;
    double Z_d0 = 2000.0;         /* charge at sheath edge */
    double Q_d = Z_d0 * DUSTY_EC;

    printf("Grain: a=%.0f um, m=%.2e kg, rho=%.0f kg/m^3\n",
           a * 1e6, m_d, rho);
    printf("Plasma: n_e=%.1e m^-3, T_e=%.1f eV, Z_d0=%.0f\n\n",
           n_e, T_e / DUSTY_EV_IN_K, Z_d0);

    /* Sheath parameters */
    double lambda_D = dust_debye_electron(n_e, T_e);
    double E_0_sheath = 10000.0;  /* 10 kV/m at electrode */
    double g = 9.81;

    printf("Sheath: lambda_D = %.0f um, E_0 = %.0f V/m\n\n",
           lambda_D * 1e6, E_0_sheath);

    /* Compute levitation height for different grain sizes */
    printf("Grain size effect on levitation:\n");
    printf("a [um]    Z_d0     m_d [kg]       z_eq [mm]    F_E [N]      F_g [N]\n");
    printf("------    -----    ----------     ----------    ---------    ---------\n");

    double a_arr[] = {1.0e-6, 2.0e-6, 5.0e-6, 10.0e-6, 15.0e-6};

    for (int i = 0; i < 5; i++) {
        double ai = a_arr[i];
        double mi = dust_grain_mass(ai, rho);
        double Z_di = 1400.0 * (ai / 1.0e-6) * (T_e / DUSTY_EV_IN_K);
        double Q_di = Z_di * DUSTY_EC;

        double z_eq = dust_levitation_height(
            mi, Q_di, g, E_0_sheath, lambda_D, 0.0, 0.0, 0.02, 1e-6);

        double F_E = Q_di * E_0_sheath * exp(-z_eq / lambda_D);
        double F_g = mi * g;

        printf("%6.1f    %6.0f    %10.2e  %11.4f    %10.2e  %10.2e\n",
               ai * 1e6, Z_di, mi, z_eq * 1000.0, F_E, F_g);

        if (z_eq > 0.0) {
            CHECK(fabs(F_E - F_g) / F_g < 0.01, "Force balance at equilibrium");
        }
    }

    printf("\n--- Ion Drag Effect on Levitation ---\n\n");

    /* Ion drag pushes downward, reducing levitation height */
    double n_i = 1.02e15;
    double m_i = 6.64e-26;
    double u_i = 300.0;          /* ion drift ~Mach 0.1 */
    double phi_s = -3.0;
    double T_i = 300.0;

    printf("With ion drag:\n");
    double F_id = dust_ion_drag_total(
        a, n_i, m_i, u_i, phi_s, Z_d0, T_i, lambda_D);

    double z_eq_no_drag = dust_levitation_height(
        m_d, Q_d, g, E_0_sheath, lambda_D, 0.0, 0.0, 0.02, 1e-6);
    double z_eq_drag = dust_levitation_height(
        m_d, Q_d, g, E_0_sheath, lambda_D, F_id, 0.0, 0.02, 1e-6);

    printf("  Without ion drag: z_eq = %.3f mm\n", z_eq_no_drag * 1000.0);
    printf("  With ion drag:    z_eq = %.3f mm\n", z_eq_drag * 1000.0);
    printf("  F_ion_drag = %.2e N, F_gravity = %.2e N\n", F_id, m_d * g);

    printf("\n--- Applications ---\n\n");
    printf("Semiconductor manufacturing: Dust contamination control\n");
    printf("  → Levitation height determines dust position in process chamber\n");
    printf("  → Critical for <10 nm node fabrication (ISO cleanroom standards)\n\n");

    printf("Fusion (tokamak dust): ITER dust inventory management\n");
    printf("  → Thermophoretic force dominates near hot plasma-facing components\n");
    printf("  → Dust radius a=5e-6 m, T_n=1000 K, grad_T=1e5 K/m\n");
    double F_th = dust_thermophoretic_force(
        a, 1000.0, m_i, 0.017, 1.0e5, 0.8);
    printf("  → F_th = %.2e N (comparable to gravity)\n\n", F_th);

    printf("Space dusty plasma: Saturn's rings analog\n");
    printf("  → P = Z_d*n_d/n_e >> 1: dust-dominated plasma\n");
    printf("  → Havnes parameter from Cassini: P ~ 10^2-10^3\n");
    printf("  → Radiation pressure competes with gravity for sub-micron grains\n");

    printf("\n=== Done ===\n\n");
    return 0;
}
