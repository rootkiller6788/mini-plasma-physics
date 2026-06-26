/*
 * alfven_wave.c -- Alfven Wave Propagation Example
 * Demonstrates an Alfven wave traveling along a uniform magnetic field.
 * Reference: Goldston & Rutherford (1995) Ch.19, Alfven (1942)
 */
#include "mhd_defs.h"
#include "mhd_waves.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    puts("=== Alfven Wave Propagation ===\n");

    double B0     = 1.0;
    double rho0   = 1.67e-21;
    double p0     = 1.0e5;
    double gamma  = 5.0 / 3.0;

    MHDState state;
    mhd_state_init_uniform(&state, rho0, 0.0, 0.0, 0.0, B0, 0.0, 0.0, p0);

    double Bmag = mhd_vector_magnitude(state.Bx, state.By, state.Bz);
    double va   = mhd_alfven_speed(Bmag, state.rho);
    double cs   = mhd_sound_speed(state.p, state.rho, gamma);

    printf("Magnetic field:  B    = %.3f T\n", Bmag);
    printf("Plasma density:  rho  = %.3e kg/m^3\n", state.rho);
    printf("Alfven speed:    v_A  = %.3e m/s\n", va);
    printf("Sound speed:     c_s  = %.3e m/s\n", cs);
    printf("Plasma beta:     beta = %.3e\n\n", mhd_plasma_beta(state.p, Bmag));

    /* Wave speeds along B direction */
    MHDWaveSpeeds ws;
    mhd_wavespeeds_compute(&state, 1.0, 0.0, 0.0, gamma, &ws);
    printf("Wave speeds along B (x-direction):\n");
    printf("  Alfven:  %.3e m/s\n", ws.c_alfven);
    printf("  Fast:    %.3e m/s\n", ws.c_fast);
    printf("  Slow:    %.3e m/s\n\n", ws.c_slow);

    /* Wave speeds perpendicular to B */
    mhd_wavespeeds_compute(&state, 0.0, 1.0, 0.0, gamma, &ws);
    printf("Wave speeds perpendicular to B (y-direction):\n");
    printf("  Alfven:  %.3e m/s (zero -- cannot propagate)\n", ws.c_alfven);
    printf("  Fast:    %.3e m/s\n\n", ws.c_fast);

    /* Alfven dispersion relation */
    double k_par  = 2.0 * M_PI / 1.0;
    double k_perp = 0.0;
    double omega_r, omega_i;
    mhd_alfven_dispersion(k_par, k_perp, va, 0.0, &omega_r, &omega_i);
    printf("Alfven dispersion (k_par = %.1f m^-1):\n", k_par);
    printf("  omega  = %.3e rad/s\n", omega_r);
    printf("  period = %.3e s\n", 2.0 * M_PI / omega_r);

    return 0;
}
