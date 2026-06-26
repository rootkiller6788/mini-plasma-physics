/**
 * @file  example_charging.c
 * @brief Example: Dust grain charging equilibrium.
 *
 * Solves the floating potential for a dust grain using Newton-Raphson
 * and demonstrates how charge depends on plasma parameters.
 *
 * L6: Canonical problem — Dust charging equilibrium.
 */

#include "dusty_plasma.h"
#include "dusty_charging.h"
#include <math.h>
#include <stdio.h>

int main(void)
{
    printf("\n=== Dust Grain Charging Equilibrium ===\n\n");

    double a = 5.0e-6;           /* 5 micron radius */
    DustGrain grain = dust_grain_init(a, DUSTY_RHO_SILICA, DUST_MATERIAL_SILICA);
    printf("Grain: radius = %.1f um, mass = %.2e kg\n\n", a * 1e6, grain.mass);

    double n_e = 1.0e15;
    double m_i = 6.64e-26;       /* Argon */
    double T_i = 300.0;
    double n_i = 1.1e15;

    printf("T_e [eV]    phi_f [V]    Z_d\n");
    printf("--------    ---------    --------\n");

    double T_eV_arr[] = {0.5, 1.0, 2.0, 3.0, 5.0, 8.0};

    for (int i = 0; i < 6; i++) {
        double T_eV = T_eV_arr[i];
        double T_e = T_eV * DUSTY_EV_IN_K;
        double phi_init = -2.5 * DUSTY_KB * T_e / DUSTY_EC;

        double phi_f = dust_floating_potential_solve(
            a, n_e, n_i, T_e, T_i, m_i, CHARGE_MODEL_OML,
            0.0, 0.0, 0.0, 0.0, 300.0, 5.0,
            phi_init, 1e-6, 30);

        double Z_d = dust_equilibrium_charge_number(a, phi_f);

        DustChargingCurrent cc = dust_compute_charging_currents(
            a, n_e, n_i, T_e, T_i, m_i, phi_f, CHARGE_MODEL_OML,
            0.0, 0.0, 0.0, 0.0, 300.0, 5.0);

        printf("%8.1f    %9.3f    %8.0f   (I_net=%.1e A)\n",
               T_eV, phi_f, Z_d, cc.I_total);
    }

    printf("\n--- Charge Dynamics ---\n\n");

    /* Integrate charging from neutral to equilibrium */
    double T_e = 3.0 * DUSTY_EV_IN_K;
    int n_steps = 100;
    double dt = 1e-5;

    double Q_arr[100], t_arr[100];
    dust_integrate_charge_dynamics(
        a, n_e, n_i, T_e, T_i, m_i, 0.0, n_steps * dt, dt,
        Q_arr, t_arr, n_steps);

    double phi_f = dust_floating_potential_solve(
        a, n_e, n_i, T_e, T_i, m_i, CHARGE_MODEL_OML,
        0.0, 0.0, 0.0, 0.0, 300.0, 5.0,
        -2.5 * DUSTY_KB * T_e / DUSTY_EC, 1e-6, 30);
    double Q_eq = dust_equilibrium_charge(a, phi_f);

    printf("Charge evolution:\n");
    printf("  t [ms]      Q(t)/Q_eq\n");
    for (int i = 0; i < 3; i++)
        printf("  %.5f    %.4f\n", t_arr[i] * 1000.0, Q_arr[i] / Q_eq);
    printf("  ...\n");
    for (int i = n_steps - 3; i < n_steps; i++)
        printf("  %.5f    %.4f\n", t_arr[i] * 1000.0, Q_arr[i] / Q_eq);

    printf("\nCharge relaxation time: %.2e s\n",
           dust_charge_relaxation_time(a, n_e, n_i, T_e, T_i, m_i, phi_f));

    printf("\n=== Done ===\n\n");
    return 0;
}
