/**
 * test_core.c — Fusion plasma physics tests
 * Uses assert() for validation. Compile with -D_DEFAULT_SOURCE
 */
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "../include/fusion_plasma.h"
#include "../include/fusion_confinement.h"

static int tests_run = 0;
static int tests_passed = 0;
#define TEST(name) do { tests_run++; printf("  %s... ", name); } while(0)
#define PASS() do { tests_passed++; printf("PASS\n"); } while(0)
#define CHECK(cond) do { if(!(cond)){printf("FAIL at %s:%d\n",__FILE__,__LINE__);assert(0);} } while(0)

int main(void) {
    printf("=== Fusion Plasma Core Tests ===\n\n");

    /* L2: Debye length */
    TEST("Debye length (Te=10keV, ne=1e20)");
    {
        double ld = debye_length(10e3, 1e20);
        CHECK(ld > 1e-6 && ld < 1e-3);
        CHECK(fabs(ld - 7.43e-5) < 1e-5);
        PASS();
    }

    /* L2: Plasma frequency */
    TEST("Plasma frequency (ne=1e20)");
    {
        double wp = plasma_frequency(1e20);
        CHECK(wp > 1e11 && wp < 1e13);
        double f_ghz = wp / (2.0 * M_PI) / 1e9;
        CHECK(fabs(f_ghz - 89.8) < 10.0);
        PASS();
    }

    /* L2: Coulomb logarithm */
    TEST("Coulomb logarithm (Te=10keV, ne=1e20)");
    {
        double lnL = coulomb_logarithm(10e3, 1e20);
        CHECK(lnL >= 2.0 && lnL <= 30.0);
        CHECK(lnL > 10.0);
        PASS();
    }

    /* L2: Thermal velocities */
    TEST("Electron thermal velocity (Te=10keV)");
    {
        double vthe = electron_thermal_velocity(10e3);
        CHECK(vthe > 1e7 && vthe < 1e8);
        PASS();
    }

    TEST("Ion thermal velocity (D, Ti=10keV)");
    {
        double vthi = ion_thermal_velocity(10e3, M_D);
        CHECK(vthi > 5e5 && vthi < 5e6);
        PASS();
    }

    /* L2: Gyroradii */
    TEST("Ion gyroradius (D, 10keV, 5T)");
    {
        double rho_i = ion_gyroradius(10e3, M_D, 5.0, 1.0);
        CHECK(rho_i > 1e-4 && rho_i < 1e-2);
        PASS();
    }

    TEST("Electron gyroradius (10keV, 5T)");
    {
        double rho_e = electron_gyroradius(10e3, 5.0);
        CHECK(rho_e > 1e-6 && rho_e < 1e-4);
        PASS();
    }

    /* L2: Cyclotron frequencies */
    TEST("Electron cyclotron frequency (5T)");
    {
        double wce = electron_cyclotron_frequency(5.0);
        double f_ghz = wce / (2.0 * M_PI) / 1e9;
        CHECK(fabs(f_ghz - 140.0) < 5.0);
        PASS();
    }

    TEST("Ion cyclotron frequency (D, 5T)");
    {
        double wci = ion_cyclotron_frequency(5.0, M_D, 1.0);
        double f_mhz = wci / (2.0 * M_PI) / 1e6;
        CHECK(fabs(f_mhz - 38.0) < 5.0);
        PASS();
    }

    /* L2: Plasma beta */
    TEST("Plasma beta (n=1e20, T=10keV, B=5T)");
    {
        double beta = plasma_beta(1e20, 10e3, 5.0);
        CHECK(beta > 0.001 && beta < 0.1);
        PASS();
    }

    /* L2: Safety factor */
    TEST("Safety factor (ITER-like)");
    {
        double q95 = safety_factor_from_current(2.0, 6.2, 5.3, 15e6);
        CHECK(q95 > 1.0 && q95 < 10.0);
        PASS();
    }

    /* L2: Aspect ratio */
    TEST("Aspect ratio (ITER: R=6.2, a=2.0)");
    {
        double A = aspect_ratio(6.2, 2.0);
        CHECK(fabs(A - 3.1) < 0.01);
        PASS();
    }

    /* L2: Plasma volume */
    TEST("Plasma volume (ITER-like)");
    {
        double Vp = plasma_volume(6.2, 2.0, 1.7);
        CHECK(Vp > 500.0 && Vp < 1500.0);
        PASS();
    }

    /* L2: Greenwald density */
    TEST("Greenwald density limit (ITER: Ip=15MA, a=2.0m)");
    {
        double nG = greenwald_density_limit(15e6, 2.0);
        CHECK(nG > 1e19 && nG < 1e21);
        PASS();
    }

    /* L3: Collision frequencies */
    TEST("Electron-ion collision frequency (ITER-like)");
    {
        double lnL = coulomb_logarithm(10e3, 1e20);
        double nuei = electron_ion_collision_freq(1e20, 1.0, 10e3, lnL);
        CHECK(nuei > 1e3 && nuei < 1e6);
        PASS();
    }

    /* L4: Lawson criterion */
    TEST("Lawson criterion (ITER target)");
    {
        double lp = lawson_criterion(1.0, 10.0, 3.7);
        CHECK(lp > 1.0);
        CHECK(fabs(lp - 37.0) < 10.0);
        PASS();
    }

    TEST("Ignition condition (sub-ignition: nTtau=0.5)");
    {
        int ign = ignition_condition(0.5, 5.0, 0.2);
        CHECK(ign == 0);
        PASS();
    }

    TEST("Ignition condition (above threshold)");
    {
        int ign = ignition_condition(3.0, 10.0, 3.0);
        CHECK(ign == 1);
        PASS();
    }

    /* L4: Fusion power density */
    TEST("D-T fusion power density (ITER-like core)");
    {
        double sv = bosch_hale_sigma_v_dt(10.0);
        double P_dt = fusion_power_density_dt(5e19, 5e19, sv, E_FUSION_DT_J);
        CHECK(sv > 1e-28 && sv < 1e-18);
        CHECK(P_dt > 1e-2 && P_dt < 1e9);
        PASS();
    }

    /* L4: Bremsstrahlung */
    TEST("Bremsstrahlung power density (ITER core)");
    {
        double Pbrem = bremsstrahlung_power_density(1e20, 1.6, 10e3);
        CHECK(Pbrem > 1e3 && Pbrem < 1e6);
        PASS();
    }

    /* L5: Bosch-Hale cross-section */
    TEST("Bosch-Hale D-T reactivity at 10 keV");
    {
        double sv = bosch_hale_sigma_v_dt(10.0);
        CHECK(sv > 1e-28 && sv < 1e-18);
        PASS();
    }

    TEST("Bosch-Hale D-T reactivity at 60 keV > 0");
    {
        double sv60 = bosch_hale_sigma_v_dt(60.0);
        CHECK(sv60 > 1e-30);
        PASS();
    }

    /* L6: ITER operating point */
    TEST("ITER operating point");
    {
        PlasmaParameters iter;
        iter_operating_point(&iter);
        CHECK(fabs(iter.R - 6.2) < 0.01);
        CHECK(fabs(iter.a - 2.0) < 0.01);
        CHECK(fabs(iter.B - 5.3) < 0.01);
        CHECK(iter.V_p > 500.0);
        PASS();
    }

    /* L6: Energy confinement IPB98 */
    TEST("IPB98(y,2) for ITER");
    {
        double tauE = energy_confinement_ipb98y2(15.0, 5.3, 1.0, 100.0, 6.2, 2.0, 1.7);
        CHECK(tauE > 1.0 && tauE < 10.0);
        PASS();
    }

    /* L4: Kruskal-Shafranov */
    TEST("Kruskal-Shafranov limit (ITER q95>1)");
    {
        double q95 = safety_factor_from_current(2.0, 6.2, 5.3, 15e6);
        int stable = kruskal_shafranov_limit(q95);
        CHECK(stable == 1);
        PASS();
    }

    /* L4: Troyon beta limit */
    TEST("Troyon beta limit (ITER-like)");
    {
        double betaN = troyon_beta_limit(2.0, 2.0, 5.3, 15.0);
        CHECK(betaN > 0.5 && betaN < 5.0);
        PASS();
    }

    /* L4: Spitzer resistivity */
    TEST("Spitzer resistivity (10keV)");
    {
        double lnL = coulomb_logarithm(10e3, 1e20);
        double eta = spitzer_resistivity(1.0, 10e3, lnL);
        CHECK(eta > 1e-10 && eta < 1e-7);
        PASS();
    }

    /* Summary */
    printf("\n=== Results: %d/%d tests passed ===\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
