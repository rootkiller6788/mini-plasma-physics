/**
 * @file  example_crystal.c
 * @brief Example: Dust crystal phase diagram and melting.
 *
 * Computes the Yukawa phase diagram and demonstrates the
 * Lindemann melting criterion for a 2D dust crystal.
 *
 * L6: Canonical problem — Dust Coulomb crystallization.
 */

#include "dusty_plasma.h"
#include "dusty_crystal.h"
#include <math.h>
#include <stdio.h>

int main(void)
{
    printf("\n=== Dust Crystal Phase Diagram ===\n\n");

    /* Yukawa phase diagram: Gamma_crit vs kappa */
    printf("kappa    Gamma_crit\n");
    printf("-----    ----------\n");
    for (int i = 0; i < 15; i++) {
        double kappa = 0.1 * i;
        double Gc = dust_critical_coupling_yukawa(kappa);
        printf("%5.1f    %10.1f\n", kappa, Gc);
    }

    printf("\n--- Crystal Stability Example ---\n\n");

    /* Typical lab crystal parameters */
    double a_ws = 2.0e-4;        /* inter-particle spacing 200 um */
    double Z_d = 2000.0;
    double Q_d = Z_d * DUSTY_EC;
    double lambda_D = 4.0e-4;    /* Debye length 400 um */
    double kappa = a_ws / lambda_D;
    double m_d = 1.15e-12;       /* dust mass */
    double n_d = 1.0e10;         /* dust density */

    printf("Parameters:\n");
    printf("  a_ws = %.0f um, lambda_D = %.0f um, kappa = %.2f\n",
           a_ws * 1e6, lambda_D * 1e6, kappa);
    printf("  Z_d = %.0f, m_d = %.2e kg\n\n", Z_d, m_d);

    printf("T_d [K]    Gamma    Gamma*   Omega_E [Hz]   rms/a    Phase\n");
    printf("------    -------   ------   ------------   ------   -----\n");

    double T_d_arr[] = {100.0, 300.0, 500.0, 1000.0, 3000.0, 10000.0};

    for (int i = 0; i < 6; i++) {
        double T_d = T_d_arr[i];
        double Gamma = dust_coulomb_coupling(Q_d, a_ws, T_d);
        double Gamma_star = dust_yukawa_coupling(Gamma, kappa);
        double Omega_E = dust_einstein_frequency(Q_d, m_d, a_ws, kappa);
        double rms = dust_rms_displacement_thermal(T_d, m_d, Omega_E);
        double rms_ratio = (a_ws > 0.0) ? rms / a_ws : 0.0;

        int crystal = dust_crystal_condition(Gamma_star, kappa);
        int phase = dust_phase_determine(Gamma_star, kappa);
        const char *phase_names[] = {"fluid", "liquid", "solid"};

        printf("%6.0f    %7.1f   %7.1f   %12.2f   %6.3f   %s\n",
               T_d, Gamma, Gamma_star, Omega_E / (2.0 * M_PI),
               rms_ratio, phase_names[phase]);
        (void)crystal;
    }

    printf("\n--- Lindemann Melting Check ---\n\n");
    double Gamma_test = 500.0;
    double Omega_E_test = dust_einstein_frequency(Q_d, m_d, a_ws, kappa);
    double T_d_test = 500.0;
    double rms_test = dust_rms_displacement_thermal(T_d_test, m_d, Omega_E_test);
    int melted = dust_lindemann_melting(rms_test, a_ws, 0.18);
    printf("T_d = %.0f K, rms/a = %.3f, melted = %s\n",
           T_d_test, rms_test / a_ws, melted ? "yes" : "no");

    printf("\n--- Madelung Constants ---\n\n");
    printf("kappa    1D chain    2D hex      3D bcc\n");
    printf("-----    --------    --------    --------\n");
    for (int i = 0; i < 10; i++) {
        double kap = 0.5 * i;
        double alpha_1d = dust_yukawa_madelung(0, kap);
        double alpha_2d = dust_yukawa_madelung(1, kap);
        double alpha_3d = dust_yukawa_madelung(2, kap);
        printf("%5.1f    %8.4f    %8.4f    %8.4f\n",
               kap, alpha_1d, alpha_2d, alpha_3d);
    }

    printf("\n=== Done ===\n\n");
    return 0;
}
