/**
 * plasma_parameter_space.c — Plasma Parameter Space Explorer
 *
 * Explores the density-temperature parameter space of plasma physics,
 * demonstrating the use of all L1-L2 functions for different
 * plasma regimes.
 */
#include "plasma_params.h"
#include <stdio.h>
#include <math.h>

static void print_plasma_regime(const char *label, double n, double Te) {
    double lD = debye_length(Te, n);
    double ND = debye_sphere_count(n, Te);
    double wpe = electron_plasma_frequency(n) / (2.0 * M_PI);
    double v_the = electron_thermal_velocity(Te);
    double lnL = coulomb_logarithm(n, Te);
    double g = plasma_parameter(n, Te);
    const char *regime = plasma_regime_name(n, Te);

    printf("\n=== %s ===\n", label);
    printf("  n     = %.2e m^-3\n", n);
    printf("  Te    = %.2e K (%.1f eV)\n", Te, Te * K_B / E_CHARGE);
    printf("  lambda_D = %.2e m\n", lD);
    printf("  N_D      = %.2e\n", ND);
    printf("  f_pe     = %.2e Hz\n", wpe);
    printf("  v_th,e   = %.2e m/s (%.2f c)\n", v_the, v_the / C_LIGHT);
    printf("  ln Lambda = %.2f\n", lnL);
    printf("  g = 1/(n lambda_D^3) = %.2e\n", g);
    printf("  Regime: %s\n", regime);
    printf("  Ideal plasma: %s\n", g < 0.01 ? "YES" : "NO");
    printf("  Collisionless: %s\n", lnL > 10.0 ? "YES" : "NO");
}

int main(void) {
    printf("=== Plasma Parameter Space Explorer ===\n");
    printf("Exploring the Hiroshima diagram of plasma physics.\n");

    /* Magnetic fusion (ITER-like) */
    print_plasma_regime("Magnetic Fusion (ITER)", 1.0e20, 1.16e8);

    /* Inertial confinement fusion */
    print_plasma_regime("Inertial Fusion (NIF)", 1.0e31, 5.0e7);

    /* Solar wind at 1 AU */
    print_plasma_regime("Solar Wind (1 AU)", 5.0e6, 1.0e5);

    /* Ionosphere F-layer */
    print_plasma_regime("Ionosphere (F-layer)", 1.0e12, 1.0e3);

    /* Glow discharge */
    print_plasma_regime("Glow Discharge", 1.0e16, 1.0e4);

    /* Lightning / arc */
    print_plasma_regime("Lightning Channel", 1.0e24, 3.0e4);

    /* Solar core */
    print_plasma_regime("Solar Core", 1.0e31, 1.5e7);

    /* Interstellar medium (warm ionized) */
    print_plasma_regime("Warm Ionized Medium", 1.0e5, 8.0e3);

    printf("\n=== Summary ===\n");
    printf("The plasma state spans over 30 orders of magnitude\n");
    printf("in density and 6 orders in temperature.\n");

    return 0;
}
