/*
 * test_laser_plasma.c -- Assert-based unit tests for mini-laser-plasma
 *
 * Tests cover all core APIs across plasma parameters, laser propagation,
 * wakefield acceleration, instabilities, absorption, ionization, and
 * particle motion.
 *
 * Build: gcc -std=c11 -Wall -Wextra -O2 -I../include test_laser_plasma.c
 *        ../src/*.c -lm -o test_laser_plasma && ./test_laser_plasma
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <float.h>
#include "plasma_constants.h"
#include "plasma_params.h"
#include "laser_plasma.h"
#include "laser_propagation.h"
#include "wakefield.h"
#include "instabilities.h"
#include "absorption.h"
#include "ionization.h"
#include "particle_motion.h"

static int tests_run = 0;
static int tests_passed = 0;
static int tests_failed = 0;

#define TEST(name) do { \
    tests_run++; \
    printf("  TEST %s ... ", name); \
} while(0)

#define PASS() do { \
    printf("PASSED\n"); \
    tests_passed++; \
} while(0)

#define FAIL(msg) do { \
    printf("FAILED: %s\n", msg); \
    tests_failed++; \
} while(0)

#define ASSERT_NEAR(a, b, tol, msg) do { \
    if (fabs((a) - (b)) > (tol)) { \
        FAIL(msg); \
    } else { \
        PASS(); \
    } \
} while(0)

#define ASSERT_TRUE(cond, msg) do { \
    if (!(cond)) { \
        FAIL(msg); \
    } else { \
        PASS(); \
    } \
} while(0)

/* ============================================================
 *  L1 Tests: Plasma Parameters
 * ============================================================ */
void test_plasma_frequency(void)
{
    TEST("plasma_frequency (ne=1e21 cm^-3)");
    double ne = 1.0e27;  /* m^-3 */
    double wp = plasma_frequency(ne);
    double expected = 1.78e15;
    ASSERT_NEAR(wp, expected, 1e13, "plasma freq wrong");

    TEST("plasma_frequency (ne=0)");
    ASSERT_TRUE(plasma_frequency(0.0) == 0.0, "zero density should give zero");
}

void test_critical_density(void)
{
    TEST("critical_density (lambda=1 um)");
    double nc = critical_density(1.0e-6);
    double expected = 1.115e27;  /* ~1.1e21 cm^-3 */
    ASSERT_NEAR(nc, expected, 1e26, "nc wrong for 1um");

    TEST("critical_density (lambda=0.35 um)");
    double nc_3w = critical_density(0.35e-6);
    /* nc scales as 1/lambda^2, so 9x larger */
    ASSERT_NEAR(nc_3w, 9.0 * expected, 2e27, "nc wrong for 0.35um");
}

void test_debye_length(void)
{
    TEST("debye_length (ne=1e21 cm^-3, Te=1 keV)");
    double ld = debye_length(1.0e27, 1000.0);
    ASSERT_NEAR(ld, 7.4e-9, 1e-9, "Debye length wrong");
}

void test_debye_sphere(void)
{
    TEST("debye_sphere_particles (typical ICF)");
    double ne = 1.0e27;
    double ld = debye_length(ne, 1000.0);
    double ND = debye_sphere_particles(ne, ld);
    ASSERT_TRUE(ND > 1.0, "N_D should be >> 1 for collective plasma");
}

void test_electron_thermal_velocity(void)
{
    TEST("electron_thermal_velocity (Te=1 keV)");
    double vth = electron_thermal_velocity(1000.0);
    ASSERT_NEAR(vth, 1.33e7, 1e6, "v_th wrong for 1keV");
}

void test_ion_sound_speed(void)
{
    TEST("ion_sound_speed (Te=1 keV, Z=6, A=12)");
    double cs = ion_sound_speed(1000.0, 6.0, 12.0);
    ASSERT_NEAR(cs, 2.2e5, 2e4, "cs wrong for carbon plasma");
}

/* ============================================================
 *  L2 Tests: Laser-Plasma Interaction
 * ============================================================ */
void test_normalized_vector_potential(void)
{
    TEST("normalized_vector_potential (I=1.37e18, lambda=1 um)");
    double a0 = normalized_vector_potential(1.37e22, 1.0e-6);
    ASSERT_NEAR(a0, 1.0, 0.1, "a0 should be ~1 at relativistic intensity");
}

void test_relativistic_intensity(void)
{
    TEST("relativistic_intensity (lambda=1 um)");
    double I_rel = relativistic_intensity(1.0e-6);
    ASSERT_NEAR(I_rel, 1.37e22, 1e21, "I_rel wrong for 1um");
}

void test_ponderomotive_potential(void)
{
    TEST("ponderomotive_potential (I=1e18, lambda=1 um)");
    double phi = ponderomotive_potential(1.0e22, 1.0e-6);
    ASSERT_TRUE(phi > 0.0, "ponderomotive potential must be positive");
}

void test_plasma_refractive_index(void)
{
    TEST("refractive_index (ne=0.1 nc)");
    double N = plasma_refractive_index(1.0, 10.0);
    ASSERT_NEAR(N, sqrt(0.9), 1e-9, "N wrong");

    TEST("refractive_index (cutoff, ne=nc)");
    ASSERT_NEAR(plasma_refractive_index(1.0, 1.0), 0.0, 1e-9, "N should be 0");
}

/* ============================================================
 *  L4 Tests: Fundamental Laws
 * ============================================================ */
void test_cold_plasma_permittivity(void)
{
    TEST("cold_permittivity (omega=2*wp)");
    double eps = cold_plasma_permittivity(2.0, 1.0);
    ASSERT_NEAR(eps, 0.75, 1e-9, "eps wrong");

    TEST("cold_permittivity (cutoff)");
    ASSERT_NEAR(cold_plasma_permittivity(1.0, 1.0), 0.0, 1e-9, "eps should be 0");
}

void test_warm_plasma_permittivity(void)
{
    TEST("warm_plasma_permittivity (finite k)");
    double eps = warm_plasma_permittivity(2.0, 0.5, 1.0, 0.1);
    ASSERT_TRUE(fabs(eps) < 1000.0, "warm eps should be finite");
}

/* ============================================================
 *  L2 Tests: Wakefield Acceleration
 * ============================================================ */
void test_cold_wavebreaking(void)
{
    TEST("cold_wavebreaking_field (ne=1e18 cm^-3)");
    double wp = plasma_frequency(1.0e24);
    double E_wb = cold_wavebreaking_field(wp);
    /* ~96 GV/m */
    ASSERT_NEAR(E_wb, 9.6e10, 1e10, "E_wb wrong");
}

void test_dephasing_length(void)
{
    TEST("dephasing_length (ne=1e18, lambda=0.8 um)");
    double ne = 1.0e24;
    double nc = critical_density(0.8e-6);
    double wp = plasma_frequency(ne);
    double omega = 2.0 * M_PI * PLASMA_C / 0.8e-6;
    double Ld = dephasing_length_1D(wp, omega, nc, ne);
    ASSERT_TRUE(Ld > 0.001 && Ld < 1.0, "L_deph out of range");
}

void test_wakefield_params(void)
{
    TEST("compute_wakefield_params (ne=1e18, a0=2)");
    double ne = 1.0e24;
    double nc = critical_density(0.8e-6);
    double wp = plasma_frequency(ne);
    double omega = 2.0 * M_PI * PLASMA_C / 0.8e-6;
    double tau = 30e-15;
    WakefieldParams wp_out;
    int ret = compute_wakefield_params(ne, nc, omega, wp, 2.0, tau, 5e-3, &wp_out);
    ASSERT_TRUE(ret == 0, "compute_wakefield_params failed");
    ASSERT_TRUE(wp_out.W_max > 0.0, "W_max should be positive");
    ASSERT_TRUE(wp_out.E_acc > 0.0, "E_acc should be positive");
}

/* ============================================================
 *  L2 Tests: Parametric Instabilities
 * ============================================================ */
void test_srs_matching(void)
{
    TEST("SRS matching conditions");
    double omega0 = 2.0 * M_PI * PLASMA_C / 1.0e-6;
    double k0 = 2.0 * M_PI / 1.0e-6;
    double wp = 0.3 * omega0;
    double omega_s, k_epw;
    srs_matching(omega0, k0, wp, &omega_s, &k_epw);
    ASSERT_TRUE(omega_s > 0.0, "omega_s should be positive");
    ASSERT_TRUE(omega_s < omega0, "omega_s must be Stokes shifted");
    ASSERT_TRUE(k_epw > k0, "k_epw > k0 for backscatter");
}

void test_sbs_matching(void)
{
    TEST("SBS matching conditions");
    double omega0 = 2.0 * M_PI * PLASMA_C / 1.0e-6;
    double k0 = 2.0 * M_PI / 1.0e-6;
    double cs = 2.0e5;
    double omega_s, omega_iaw, k_iaw;
    sbs_matching(omega0, k0, cs, &omega_s, &omega_iaw, &k_iaw);
    ASSERT_TRUE(omega_iaw > 0.0, "omega_iaw should be positive");
    ASSERT_TRUE(omega_iaw < omega0 * 0.01, "IAW freq too high");
}

void test_analyze_instability(void)
{
    TEST("analyze_instability (SRS backscatter)");
    PlasmaState ps = {1.0e26, 500.0, 100.0, 3.0, 12.0, 0.0, 1.0e-6, 1.0e21};
    PlasmaDerived pd;
    compute_all_derived(&ps, &pd);
    pd.a0 = 0.1;
    double v_osc = pd.a0 * PLASMA_C;
    InstabilityResult result;
    int ret = analyze_instability(INSTABILITY_SRS_BACKSCATTER, &ps, &pd, v_osc, 10e-6, &result);
    ASSERT_TRUE(ret == 0, "analyze_instability failed");
    ASSERT_TRUE(result.gamma0 > 0.0, "SRS growth rate should be positive");
}

/* ============================================================
 *  L2 Tests: Absorption Mechanisms
 * ============================================================ */
void test_ib_coefficient(void)
{
    TEST("inverse_bremsstrahlung coefficient");
    double k_ib = inverse_bremsstrahlung_coefficient(1.0e26, 1.0e27, 1.0e12);
    ASSERT_TRUE(k_ib > 0.0, "k_ib should be positive");
}

void test_ginzburg_function(void)
{
    TEST("Ginzburg function maximum");
    double phi_opt = ginzburg_function(0.8);
    ASSERT_NEAR(phi_opt, 1.0, 0.3, "phi(0.8) should be near max");
}

void test_absorption_types(void)
{
    TEST("absorption_type_name valid");
    const char *name = absorption_type_name(ABS_INVERSE_BREMSSTRAHLUNG);
    ASSERT_TRUE(name != NULL, "name should not be NULL");
    
    TEST("compute_absorption (IB)");
    PlasmaState ps = {1.0e26, 500.0, 100.0, 3.0, 12.0, 0.0, 1.0e-6, 1.0e20};
    PlasmaDerived pd;
    compute_all_derived(&ps, &pd);
    pd.a0 = 0.1;
    AbsorptionResult ar;
    int ret = compute_absorption(&ps, &pd, ABS_INVERSE_BREMSSTRAHLUNG, 10e-6, 0.0, &ar);
    ASSERT_TRUE(ret == 0, "compute_absorption failed");
    ASSERT_TRUE(ar.absorbed_fraction >= 0.0 && ar.absorbed_fraction <= 1.0,
                "absorption fraction out of range");
}

/* ============================================================
 *  L4 Tests: Ionization Models
 * ============================================================ */
void test_keldysh_parameter(void)
{
    TEST("Keldysh parameter (tunnel regime)");
    double gamma = keldysh_parameter(1.0e20, 1.0e-6, 13.6);
    ASSERT_TRUE(gamma < 1.0, "should be in tunnel regime");

    TEST("Keldysh parameter (MPI regime)");
    double gamma_mpi = keldysh_parameter(1.0e16, 0.25e-6, 13.6);
    ASSERT_TRUE(gamma_mpi > 0.0, "gamma should be positive");
}

void test_barrier_suppression(void)
{
    TEST("barrier_suppression_intensity (H)");
    double I_BSI = barrier_suppression_intensity(13.6, 1);
    ASSERT_NEAR(I_BSI, 1.4e18, 5e17, "BSI threshold for H wrong");
}

void test_adk_rate(void)
{
    TEST("ADK rate at BSI threshold");
    double I_BSI = barrier_suppression_intensity(13.6, 1);
    double rate = adk_tunneling_rate(I_BSI, 1.0e-6, 13.6, 1, 0, 0);
    ASSERT_TRUE(rate > 0.0, "ADK rate should be positive at BSI");
}

/* ============================================================
 *  L3/L5 Tests: Particle Motion
 * ============================================================ */
void test_particle_init(void)
{
    TEST("particle_init and kinetic energy");
    Particle3D p;
    particle_init(&p, 0, 0, 0, 0, 0, 0, PLASMA_ME, -PLASMA_E);
    ASSERT_NEAR(p.gamma, 1.0, 1e-9, "at-rest gamma should be 1");
    ASSERT_NEAR(kinetic_energy(&p), 0.0, 1e-9, "at-rest KE should be 0");
}

void test_boris_push_vacuum(void)
{
    TEST("boris_push (no E/B field)");
    Particle3D p;
    EMField field = {0, 0, 0, 0, 0, 0};
    particle_init(&p, 0, 0, 0, 1e-22, 0, 0, PLASMA_ME, -PLASMA_E);
    double gamma0 = p.gamma;
    boris_push(&p, &field, 1e-15);
    ASSERT_NEAR(p.gamma, gamma0, 1e-6, "gamma should be conserved in free drift");
}

void test_plane_wave_orbit(void)
{
    TEST("Volkov plane wave solution");
    double px, py, pz, gamma;
    plane_wave_orbit_analytic(1e-14, 1.0, 2.0 * M_PI * PLASMA_C / 1.0e-6,
                               0, PLASMA_ME, &px, &py, &pz, &gamma);
    ASSERT_TRUE(gamma >= 1.0, "gamma >= 1 always");
}

void test_particle_track(void)
{
    TEST("particle_track alloc and record");
    ParticleTrack *track = particle_track_alloc(100);
    ASSERT_TRUE(track != NULL, "track alloc failed");
    
    Particle3D p;
    particle_init(&p, 0, 0, 0, 0, 0, 0, PLASMA_ME, -PLASMA_E);
    int ret = particle_track_record(track, 0.0, &p);
    ASSERT_TRUE(ret == 0, "track record failed");
    ASSERT_TRUE(track->n_samples == 1, "sample count wrong");
    
    particle_track_free(track);
}

/* ============================================================
 *  L4 Tests: compute_all_derived
 * ============================================================ */
void test_compute_all_derived(void)
{
    TEST("compute_all_derived (typical ICF plasma)");
    PlasmaState ps = {
        1.0e27,    /* ne = 1e21 cm^-3 */
        3000.0,    /* Te = 3 keV */
        1000.0,    /* Ti = 1 keV */
        3.5,       /* Z = 3.5 (carbon plasma) */
        12.0,      /* A = 12 */
        0.0,       /* B = 0 */
        0.351e-6,  /* lambda = 0.35 um (3 omega) */
        1.0e19     /* I = 10^15 W/cm^2 */
    };
    PlasmaDerived pd;
    int ret = compute_all_derived(&ps, &pd);
    ASSERT_TRUE(ret == 0, "compute_all_derived failed");
    ASSERT_TRUE(pd.wp > 0.0, "wp should be positive");
    ASSERT_TRUE(pd.nc > 0.0, "nc should be positive");
    ASSERT_TRUE(pd.lambda_D > 0.0, "lambda_D should be positive");
    ASSERT_TRUE(pd.N_D > 1.0, "N_D should be >> 1");
    ASSERT_TRUE(pd.ne_over_nc > 0.0, "ne/nc should be positive");
    ASSERT_TRUE(pd.v_the > 0.0, "v_the should be positive");
    ASSERT_TRUE(pd.cs > 0.0, "cs should be positive");
}

void test_invalid_inputs(void)
{
    TEST("compute_all_derived (invalid input)");
    PlasmaState ps = {0, 0, 0, 0, 0, 0, 0, 0};
    PlasmaDerived pd;
    ASSERT_TRUE(compute_all_derived(&ps, &pd) != 0, "should reject zero density");

    TEST("plasma_frequency (negative density)");
    ASSERT_TRUE(plasma_frequency(-1.0) == 0.0, "should reject negative ne");

    TEST("critical_density (zero wavelength)");
    ASSERT_TRUE(critical_density(0.0) == 0.0, "should reject zero lambda");
}

/* ============================================================
 *  Test Runner
 * ============================================================ */
int main(void)
{
    printf("=== mini-laser-plasma Unit Tests ===\n\n");

    printf("[L1] Plasma Parameters:\n");
    test_plasma_frequency();
    test_critical_density();
    test_debye_length();
    test_debye_sphere();
    test_electron_thermal_velocity();
    test_ion_sound_speed();

    printf("\n[L2] Laser-Plasma Interaction:\n");
    test_normalized_vector_potential();
    test_relativistic_intensity();
    test_ponderomotive_potential();
    test_plasma_refractive_index();

    printf("\n[L4] Fundamental Laws:\n");
    test_cold_plasma_permittivity();
    test_warm_plasma_permittivity();

    printf("\n[L2] Wakefield Acceleration:\n");
    test_cold_wavebreaking();
    test_dephasing_length();
    test_wakefield_params();

    printf("\n[L2] Parametric Instabilities:\n");
    test_srs_matching();
    test_sbs_matching();
    test_analyze_instability();

    printf("\n[L2] Absorption Mechanisms:\n");
    test_ib_coefficient();
    test_ginzburg_function();
    test_absorption_types();

    printf("\n[L4] Ionization Models:\n");
    test_keldysh_parameter();
    test_barrier_suppression();
    test_adk_rate();

    printf("\n[L3/L5] Particle Motion:\n");
    test_particle_init();
    test_boris_push_vacuum();
    test_plane_wave_orbit();
    test_particle_track();

    printf("\n[L4] Derived Parameters:\n");
    test_compute_all_derived();

    printf("\n[Guard] Invalid Inputs:\n");
    test_invalid_inputs();

    printf("\n========================================\n");
    printf("Tests run:    %d\n", tests_run);
    printf("Tests passed: %d\n", tests_passed);
    printf("Tests failed: %d\n", tests_failed);
    printf("========================================\n");

    return (tests_failed == 0) ? 0 : 1;
}
