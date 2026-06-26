/**
 * example_solar_wind.c -- Parker Solar Wind Model Demo
 *
 * Demonstrates:
 *   1. Parker isothermal wind solution from corona to 1 AU
 *   2. Parker spiral angle and IMF field
 *   3. CIR formation for fast/slow wind streams
 *   4. Solar wind Mach numbers
 *
 * Reference: Parker (1958) ApJ 128:664
 *            Kivelson & Russell S6
 */
#include "../include/space_plasma.h"
#include "../include/plasma_parameters.h"
#include "../include/solar_wind.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== Parker Solar Wind Model ===\n\n");

    /* Coronal base conditions */
    double T_corona = 1.0e6;          /* corona temperature [K] */
    double gamma = 5.0/3.0;
    double c_s = sqrt(gamma * SP_KB * T_corona / SP_MP);
    double r_c = parker_critical_radius(c_s);
    double r0 = SP_RSUN * 2.5;        /* source surface ~2.5 R_sun */
    double n0 = 1.0e13;               /* base density [m^-3] */

    printf("Coronal conditions:\n");
    printf("  T_corona = %.1e K\n", T_corona);
    printf("  Sound speed c_s = %.1f km/s\n", c_s/1e3);
    printf("  Critical radius r_c = %.2f R_sun\n", r_c/SP_RSUN);
    printf("  Source surface = %.1f R_sun\n\n", r0/SP_RSUN);

    /* Compute Parker wind profile */
    printf("Parker Transonic Wind Profile:\n");
    printf("  %-12s %-12s %-12s %-12s %-12s\n",
           "r [R_sun]", "v [km/s]", "Mach", "n [cm^-3]", "T [K]");

    double r[6];
    for (int i = 0; i < 6; i++)
        r[i] = SP_RSUN * (2.0 + i * 30.0);

    parker_profile_t profile[6];
    parker_transonic_profile(r, 6, r_c, c_s, r0, n0, profile);

    for (int i = 0; i < 6; i++) {
        printf("  %-12.1f %-12.1f %-12.2f %-12.1e %-12.0f\n",
               profile[i].r/SP_RSUN, profile[i].v/1e3,
               profile[i].Mach, profile[i].n/1e6, profile[i].T);
    }
    printf("\n");

    /* Parker Spiral at 1 AU */
    printf("Parker Spiral at 1 AU:\n");
    double omega_sun = 2.0 * M_PI / carrington_period();
    double r_1AU = SP_AU;
    double v_sw = profile[5].v;  /* velocity at ~1 AU */
    double B_parker[3];

    parker_spiral_field(r_1AU, M_PI/2.0, 5.0e-9, SP_RSUN, omega_sun, v_sw, B_parker);
    double psi = parker_spiral_angle(r_1AU, M_PI/2.0, omega_sun, v_sw);

    printf("  B_r     = %.2f nT\n", B_parker[0]*1e9);
    printf("  B_phi   = %.2f nT\n", B_parker[2]*1e9);
    printf("  |B|     = %.2f nT\n", sqrt(B_parker[0]*B_parker[0]+B_parker[2]*B_parker[2])*1e9);
    printf("  Spiral angle psi = %.1f deg\n", psi*180.0/M_PI);
    printf("  Solar wind speed = %.0f km/s\n\n", v_sw/1e3);

    /* Mach numbers at 1 AU */
    double n_sw = profile[5].n;
    double rho_sw = n_sw * SP_MP;
    double B_sw = sqrt(B_parker[0]*B_parker[0] + B_parker[2]*B_parker[2]);

    printf("Mach Numbers at 1 AU:\n");
    double v_A = sp_alfven_speed(B_sw, rho_sw);
    printf("  v_A = %.1f km/s\n", v_A/1e3);
    printf("  M_A (Alfven)  = %.2f\n", sp_alfven_mach(v_sw, B_sw, rho_sw));
    printf("  M_S (Sonic)   = %.2f\n", sp_sonic_mach(v_sw, c_s));
    printf("  M_ms (Magneto) = %.2f\n\n", sp_magnetosonic_mach(v_sw, v_A, c_s));

    /* CIR formation */
    printf("CIR Formation:\n");
    double v_fast = 750e3, v_slow = 350e3;
    double r_fast_src = SP_RSUN * 2.5;
    double r_slow_src = SP_RSUN * 2.5;
    double r_cir = cir_formation_radius(v_fast, v_slow, r_fast_src, r_slow_src);
    printf("  Fast wind: %.0f km/s, Slow wind: %.0f km/s\n", v_fast/1e3, v_slow/1e3);
    if (isinf(r_cir))
        printf("  No CIR: fast wind never catches slow wind\n");
    else
        printf("  CIR forms at r = %.2f AU\n", r_cir/SP_AU);

    double t_slow = solar_wind_travel_time(SP_AU, r_slow_src, v_slow);
    double t_fast = solar_wind_travel_time(SP_AU, r_fast_src, v_fast);
    printf("  Travel time to 1 AU: slow=%.1f hr, fast=%.1f hr\n",
           t_slow/3600.0, t_fast/3600.0);

    printf("\n=== Done ===\n");
    return 0;
}
