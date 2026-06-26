/*
 * example_critical_density.c -- Critical density and plasma parameters
 *
 * Demonstrates computation of fundamental plasma parameters for
 * typical laser-plasma interaction scenarios: 1 um Nd:glass,
 * 0.35 um 3-omega, and 10.6 um CO2 lasers.
 *
 * Build: gcc -std=c11 -I../include example_critical_density.c
 *        ../src/plasma_params.c -lm -o example_critical_density
 */

#include <stdio.h>
#include <math.h>
#include "plasma_constants.h"
#include "plasma_params.h"

int main(void)
{
    printf("=== Laser-Plasma Critical Density Calculation ===\n\n");

    double wavelengths[] = {10.6e-6, 1.053e-6, 0.8e-6, 0.351e-6, 0.248e-6};
    const char *names[] = {"CO2", "Nd:glass", "Ti:Sapphire", "3-omega Nd:glass", "KrF"};

    printf("%-20s %12s %16s %16s\n",
           "Laser", "lambda [um]", "nc [cm^-3]", "nc [m^-3]");
    printf("-------------------- ------------ ---------------- ----------------\n");

    for (int i = 0; i < 5; i++) {
        double lam = wavelengths[i];
        double nc = critical_density(lam);
        printf("%-20s %12.3f %16.3e %16.3e\n",
               names[i], lam * 1e6, nc * 1e-6, nc);
    }

    printf("\n=== Derived Plasma Parameters ===\n\n");

    /* Example: ICF-relevant carbon plasma */
    PlasmaState ps = {
        5.0e26,     /* ne = 5e20 cm^-3 */
        2000.0,     /* Te = 2 keV */
        500.0,      /* Ti = 0.5 keV */
        6.0,        /* fully ionized carbon */
        12.0,       /* carbon-12 */
        0.0,        /* unmagnetized */
        0.351e-6,   /* 3-omega Nd:glass */
        1.0e19      /* 10^15 W/cm^2 */
    };
    PlasmaDerived pd;
    compute_all_derived(&ps, &pd);

    printf("Input plasma state:\n");
    printf("  ne = %.2e m^-3 (%.2e cm^-3)\n", ps.ne, ps.ne * 1e-6);
    printf("  Te = %.0f eV, Ti = %.0f eV, Z = %.1f, A = %.0f\n",
           ps.Te, ps.Ti, ps.Z, ps.A);
    printf("  lambda = %.3f um\n\n", ps.lam_laser * 1e6);

    printf("Derived parameters:\n");
    printf("  Plasma frequency:  omega_p = %.3e rad/s (period = %.3f fs)\n",
           pd.wp, 2.0 * M_PI / pd.wp * 1e15);
    printf("  Critical density:  nc     = %.3e cm^-3\n", pd.nc * 1e-6);
    printf("  Density ratio:     ne/nc  = %.4f\n", pd.ne_over_nc);
    printf("  Debye length:      ld     = %.3f nm\n", pd.lambda_D * 1e9);
    printf("  Debye sphere:      N_D    = %.1f\n", pd.N_D);
    printf("  Skin depth:        delta  = %.3f um\n", pd.skin_depth * 1e6);
    printf("  Thermal velocity:  v_the  = %.3e m/s (%.3f c)\n",
           pd.v_the, pd.v_the / PLASMA_C);
    printf("  Ion sound speed:   cs     = %.3e m/s\n", pd.cs);
    printf("  Collision freq:    nu_ei  = %.3e s^-1\n", pd.nu_ei);
    printf("  Coulomb log:       lnL    = %.3f\n", pd.ln_Lambda);
    printf("  In-plasma lambda:  lam_pl = %.3f um\n", pd.lambda_plasma * 1e6);
    printf("  Refractive index:  N      = %.4f\n",
           plasma_refractive_index(ps.ne, pd.nc));
    printf("  Group velocity:    vg/c   = %.4f\n",
           group_velocity(ps.ne, pd.nc) / PLASMA_C);

    return 0;
}
