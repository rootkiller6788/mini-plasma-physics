/**
 * test_plasma_params.c — Tests for plasma parameter computation (L1-L2)
 *
 * Validates fundamental plasma parameter formulas against
 * known reference values from NRL Plasma Formulary.
 */
#include "plasma_params.h"
#include "plasma_constants.h"
#include <stdio.h>
#include <math.h>
#include <assert.h>

#define TOL 1e-10

static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name) do { tests_run++; printf("  %s... ", name); } while(0)
#define PASS() do { tests_passed++; printf("PASS\n"); } while(0)
#define FAIL(msg) do { printf("FAIL: %s\n", msg); } while(0)
#define ASSERT_NEAR(a, b, tol, name) do { \
    if (fabs((a)-(b)) > (tol)) { FAIL(name); return; } \
} while(0)

static void test_debye_length(void) {
    TEST("debye_length (fusion plasma)");
    /* Fusion plasma: Te = 10 keV = 1.16e8 K, ne = 1e20 m^-3
     * lambda_D ~ 7.4e-5 m */
    double Te = 1.16e8, ne = 1.0e20;
    double lD = debye_length(Te, ne);
    assert(lD > 1e-6 && lD < 1e-3);
    ASSERT_NEAR(lD, 7.43e-5, 0.5e-5, "fusion Debye length");
    PASS();
}

static void test_plasma_frequency(void) {
    TEST("electron plasma frequency");
    double ne = 1.0e20;
    double wpe = electron_plasma_frequency(ne);
    /* omega_pe ~ sqrt(ne*e^2/(eps0*me)) ~ 5.64e11 rad/s for ne=1e20 */
    assert(wpe > 1e11 && wpe < 1e12);
    double fpe = plasma_frequency_hz(ne, M_ELECTRON);
    assert(fabs(fpe - wpe/(2.0*M_PI)) < TOL);
    PASS();
}

static void test_thermal_velocity(void) {
    TEST("electron thermal velocity");
    double Te = 1.16e8; /* 10 keV */
    double v_the = electron_thermal_velocity(Te);
    assert(v_the > 1e7 && v_the < 1e8);
    /* v_th ~ sqrt(2*kB*Te/me) ~ 5.9e7 m/s for 10 keV */
    PASS();
}

static void test_ion_sound_speed(void) {
    TEST("ion sound speed (cold ion limit)");
    double Te = 1.16e8;
    double mi = M_PROTON * 2.0; /* Deuterium */
    double cs_cold = ion_sound_speed_cold(Te, mi);
    double cs_full = ion_sound_speed(Te, 0.0, mi);
    assert(fabs(cs_cold - cs_full) < TOL);
    assert(cs_cold > 1e5 && cs_cold < 1e7);
    PASS();
}

static void test_alfven_speed(void) {
    TEST("alfven speed (fusion)");
    double B = 5.3;
    double ne = 1.0e20;
    double mi = M_PROTON * 2.0;
    double rho = ne * mi;
    double vA = alfven_speed(B, rho);
    assert(vA > 1e6 && vA < 1e8);
    PASS();
}

static void test_gyrofrequency(void) {
    TEST("electron gyrofrequency");
    double B = 5.3;
    double wce = electron_gyrofrequency(B);
    assert(wce > 1e11 && wce < 1e12);
    double rLe = gyroradius(5.9e7, B, E_CHARGE, M_ELECTRON);
    assert(rLe > 1e-5 && rLe < 1e-3);
    PASS();
}

static void test_coulomb_logarithm(void) {
    TEST("coulomb logarithm (fusion)");
    double ne = 1.0e20, Te = 1.16e8;
    double lnL = coulomb_logarithm(ne, Te);
    /* Fusion: ln Lambda ~ 15-20 */
    assert(lnL > 10.0 && lnL < 25.0);
    PASS();
}

static void test_debye_shielding(void) {
    TEST("debye shielding potential");
    double lD = debye_length(1.16e8, 1.0e20);
    double phi_D = debye_shielding_potential(lD, E_CHARGE, lD);
    /* At r = lambda_D, potential is 1/e of Coulomb value */
    double phi_coul = E_CHARGE / (4.0 * M_PI * EPSILON_0 * lD);
    assert(fabs(phi_D / phi_coul - 1.0/exp(1.0)) < 0.01);
    PASS();
}

static void test_quasi_neutrality(void) {
    TEST("quasi-neutrality check");
    assert(is_quasi_neutral(1.0e20, 1.0e20, 1) == 1);
    assert(is_quasi_neutral(1.0e20, 0.5e20, 1) == 0);
    PASS();
}

static void test_plasma_regime(void) {
    TEST("plasma regime classification");
    PlasmaRegime regime;
    classify_plasma(1.0e20, 1.16e8, &regime);
    assert(regime.is_ideal == 1);
    assert(regime.is_collisionless == 1);
    const char *name = plasma_regime_name(1.0e20, 1.16e8);
    assert(name != NULL);
    PASS();
}

int main(void) {
    printf("=== Test: plasma_params (L1-L2) ===\n");
    test_debye_length();
    test_plasma_frequency();
    test_thermal_velocity();
    test_ion_sound_speed();
    test_alfven_speed();
    test_gyrofrequency();
    test_coulomb_logarithm();
    test_debye_shielding();
    test_quasi_neutrality();
    test_plasma_regime();
    printf("Result: %d/%d tests passed\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
