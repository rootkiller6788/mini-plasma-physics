/**
 * @file  example_daw.c
 * @brief Example: Dust-Acoustic Wave dispersion and damping.
 *
 * Demonstrates the complete DAW dispersion relation, including
 * Landau and collisional damping. Outputs omega(k) for a range
 * of wavenumbers, showing the transition from sound-like to
 * constant-frequency behavior.
 *
 * L6: Canonical problem — Dust-acoustic wave propagation.
 *
 * Compile: make examples
 * Run: build/example_daw
 */

#include "dusty_plasma.h"
#include "dusty_waves.h"
#include "dusty_transport.h"
#include <math.h>
#include <stdio.h>

int main(void)
{
    printf("\n=== Dust-Acoustic Wave Dispersion Example ===\n\n");

    /* Typical laboratory parameters */
    double Z_d = 2000.0;
    double T_e = 3.0 * DUSTY_EV_IN_K;    /* 3 eV */
    double T_i = 300.0;                   /* cold ions */
    double T_d = 300.0;                   /* room-temp dust */
    double m_d = 1.15e-12;               /* 5 um SiO2 grain */
    double n_e = 1.0e15;
    double n_d = 1.0e10;
    double n_i = n_e + Z_d * n_d;         /* quasineutrality */

    /* Derived parameters */
    double lambda_De = dust_debye_electron(n_e, T_e);
    double lambda_Di = dust_debye_ion(n_i, T_i);
    double lambda_D = dust_debye_total(lambda_De, lambda_Di);
    double c_da = dust_acoustic_speed(Z_d, T_e, m_d);
    double omega_pd = dust_plasma_freq(n_d, Z_d, m_d);
    double v_thd = sqrt(DUSTY_KB * T_d / m_d);

    /* Estimate neutral collision frequency */
    double nu_dn = dust_neutral_collision_freq(
        5.0e-6, 2.4e21, 300.0, 6.63e-26, m_d, 0.8);

    printf("Plasma parameters:\n");
    printf("  n_e = %.2e m^-3, n_i = %.2e m^-3\n", n_e, n_i);
    printf("  T_e = %.1f K (%.1f eV), T_i = %.1f K\n", T_e, T_e / DUSTY_EV_IN_K, T_i);
    printf("  lambda_De = %.2e m, lambda_Di = %.2e m\n", lambda_De, lambda_Di);
    printf("  lambda_D = %.2e m (%.1f um)\n", lambda_D, lambda_D * 1e6);
    printf("  c_da = %.4f m/s\n", c_da);
    printf("  omega_pd = %.2f rad/s (f_pd = %.2f Hz)\n", omega_pd, omega_pd / (2.0 * M_PI));
    printf("  v_thd = %.4e m/s\n", v_thd);
    printf("  nu_dn = %.2f Hz\n\n", nu_dn);

    /* Sweep wavenumbers from long to short wavelength */
    printf("  k [rad/m]    omega_r [Hz]    v_phi [m/s]    gamma [1/s]    Q-factor\n");
    printf("  ----------   ------------    -----------    -----------    --------\n");

    for (int i = 0; i < 20; i++) {
        double k = 10.0 * pow(10.0, 0.2 * i);  /* log sweep 10 to 63k */
        WaveMode wm = dust_acoustic_wave_full(
            k, c_da, lambda_D, omega_pd, v_thd, nu_dn);

        double f_hz = wm.omega_r / (2.0 * M_PI);
        double Q = (fabs(wm.omega_i) > 1e-30)
                   ? fabs(wm.omega_r / (2.0 * wm.omega_i)) : 1e30;

        printf("  %10.2e  %12.3f  %11.4f  %12.3e  %9.2f\n",
               k, f_hz, wm.phase_velocity, wm.omega_i, Q);
    }

    printf("\nInterpretation:\n");
    printf("  - Low k: omega ~ k*c_da (sound-like, v_phi constant)\n");
    printf("  - High k: omega → omega_pd (constant, v_phi → 0)\n");
    printf("  - Damping dominated by collisions for this example\n");
    printf("  - Q >> 1 → waves propagate well; Q < 1 → overdamped\n");

    /* Check if DAW is observable */
    double k_test = 1000.0;
    WaveMode wm_test = dust_acoustic_wave_full(
        k_test, c_da, lambda_D, omega_pd, v_thd, nu_dn);
    double Q_test = fabs(wm_test.omega_r / (2.0 * wm_test.omega_i + 1e-30));

    printf("\nDAW observability at k=%.0e rad/m: Q = %.2f → %s\n",
           k_test, Q_test,
           Q_test > 0.5 ? "observable" : "overdamped");

    printf("\n=== Done ===\n\n");
    return 0;
}
