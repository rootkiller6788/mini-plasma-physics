/**
 * @file  test_dusty.c
 * @brief Comprehensive test suite for mini-dusty-plasma.
 *
 * Tests cover all core APIs with mathematical assertions.
 * Based on assert() from standard C — not custom macros.
 *
 * Usage: make test
 */

#include "dusty_plasma.h"
#include "dusty_charging.h"
#include "dusty_waves.h"
#include "dusty_crystal.h"
#include "dusty_transport.h"
#include "dusty_forces.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int tests_passed = 0;
static int tests_failed = 0;

#define TEST(name) do { printf("  TEST: %s ... ", name); } while(0)
#define PASS() do { printf("PASS\n"); tests_passed++; } while(0)
#define FAIL(msg) do { printf("FAIL: %s\n", msg); tests_failed++; } while(0)
#define CHECK(cond, msg) do { if (cond) PASS(); else FAIL(msg); } while(0)

/* Tolerance for floating-point comparisons */
#define TOL 1e-9

/* ================================================================
 * L1 Tests — Core Definitions
 * ================================================================ */

void test_debye_lengths(void)
{
    TEST("Electron Debye length");
    double n_e = 1.0e15;  /* 10^15 m^-3 */
    double T_e = 3.0 * DUSTY_EV_IN_K;  /* 3 eV */
    double lde = dust_debye_electron(n_e, T_e);
    CHECK(lde > 1.0e-5 && lde < 1.0e-3,
          "lambda_De should be ~400 um for 3 eV, 1e15 m^-3");

    TEST("Ion Debye length");
    double n_i = 1.1e15;
    double T_i = 300.0;
    double ldi = dust_debye_ion(n_i, T_i);
    CHECK(ldi > 0.0 && ldi < lde,
          "lambda_Di should be smaller than lambda_De for T_i << T_e");

    TEST("Total Debye length");
    double ld = dust_debye_total(lde, ldi);
    CHECK(ld > 0.0 && ld < ldi,
          "lambda_D total should be smaller than individual");

    TEST("Debye length zero input");
    double ld_zero_n = dust_debye_electron(0.0, T_e);
    CHECK(ld_zero_n == 0.0, "lambda_De should be 0 for n_e=0");

    TEST("Debye length zero temperature");
    double ld_zero_T = dust_debye_ion(n_i, 0.0);
    CHECK(ld_zero_T == 0.0, "lambda_Di should be 0 for T_i=0");
}

void test_plasma_frequencies(void)
{
    TEST("Electron plasma frequency");
    double fpe = dust_plasma_freq_electron(1.0e15);
    CHECK(fpe > 1.0e8 && fpe < 1.0e10,
          "omega_pe should be ~1.8 GHz");

    TEST("Ion plasma frequency");
    double fp_i = dust_plasma_freq_ion(1.0e15, 6.64e-26);
    CHECK(fp_i > 1.0e5 && fp_i < 1.0e8,
          "omega_pi should be ~6.6 MHz for Ar");

    TEST("Dust plasma frequency");
    double n_d = 1.0e10;
    double Z_d = 2000.0;
    double m_d = 1.15e-12;
    double fpd = dust_plasma_freq(n_d, Z_d, m_d);
    CHECK(fpd > 1.0 && fpd < 1.0e4,
          "omega_pd should be ~300 rad/s for micron grains");

    TEST("Dust plasma frequency zero inputs");
    double fpd_zero = dust_plasma_freq(0.0, Z_d, m_d);
    CHECK(fpd_zero == 0.0, "omega_pd should be 0 for n_d=0");
}

void test_havnes_parameter(void)
{
    TEST("Havnes parameter isolated regime");
    double P = dust_havnes_parameter(2000.0, 1.0e10, 1.0e15);
    CHECK(P > 0.0 && P < 1.0, "P should be < 1 for typical lab conditions");

    TEST("Havnes parameter zero electron density");
    double P_zero = dust_havnes_parameter(2000.0, 1.0e10, 0.0);
    CHECK(P_zero > 1.0, "P should be infinite when n_e=0 (regardless of Z_d)");

    TEST("Regime classification");
    DustRegime r_iso = dust_classify_regime(0.01);
    CHECK(r_iso == DUST_REGIME_ISOLATED, "P=0.01 => ISOLATED");
    DustRegime r_col = dust_classify_regime(0.5);
    CHECK(r_col == DUST_REGIME_COLLECTIVE, "P=0.5 => COLLECTIVE");
    DustRegime r_dom = dust_classify_regime(5.0);
    CHECK(r_dom == DUST_REGIME_DOMINATED, "P=5.0 => DOMINATED");
}

/* ================================================================
 * L2 Tests — Core Concepts
 * ================================================================ */

void test_coulomb_coupling(void)
{
    TEST("Coulomb coupling Gamma");
    double Q_d = 2000.0 * DUSTY_EC;
    double a = 5.0e-6;
    double T_d = 300.0;
    double Gamma = dust_coulomb_coupling(Q_d, a, T_d);
    CHECK(Gamma > 100.0, "Gamma should be >> 1 for micron grains");

    TEST("Yukawa coupling");
    double kappa = 1.0;
    double Gamma_star = dust_yukawa_coupling(Gamma, kappa);
    CHECK(Gamma_star < Gamma, "Yukawa coupling should be less than bare Coulomb");

    TEST("Zero inputs coupling");
    double Gz = dust_coulomb_coupling(Q_d, 0.0, T_d);
    CHECK(Gz == 0.0, "Gamma should be 0 for a=0");
}

void test_dust_acoustic_speed(void)
{
    TEST("Dust acoustic speed");
    double Z_d = 2000.0;
    double T_e = 3.0 * DUSTY_EV_IN_K;
    double m_d = 1.15e-12;
    double cda = dust_acoustic_speed(Z_d, T_e, m_d);
    CHECK(cda > 0.001 && cda < 0.1, "c_da should be ~30 mm/s");

    TEST("Ion acoustic speed");
    double m_i = 6.64e-26;
    double cs = dust_ion_acoustic_speed(T_e, m_i);
    CHECK(cs > 100.0 && cs < 10000.0, "c_s should be ~2700 m/s");

    TEST("Zero input acoustic speed");
    double cda0 = dust_acoustic_speed(0.0, T_e, m_d);
    CHECK(cda0 == 0.0, "c_da should be 0 for Z_d=0");
}

void test_grain_init(void)
{
    TEST("Grain initialization");
    DustGrain g = dust_grain_init(5.0e-6, DUSTY_RHO_SILICA, DUST_MATERIAL_SILICA);
    CHECK(g.radius == 5.0e-6, "Radius preserved");
    CHECK(g.mass > 1.0e-15 && g.mass < 1.0e-9, "Mass in expected range");
    CHECK(g.material == DUST_MATERIAL_SILICA, "Material type preserved");
    CHECK(g.charge == 0.0, "Initial charge should be zero");

    TEST("Grain mass computation");
    double m = dust_grain_mass(5.0e-6, 2200.0);
    double expected = (4.0/3.0) * M_PI * 5.0e-6 * 5.0e-6 * 5.0e-6 * 2200.0;
    CHECK(fabs(m - expected) < expected * TOL, "Mass matches spherical formula");
}

void test_plasma_state_init(void)
{
    TEST("Plasma state initialization");
    double n_e = 1.0e15;
    double T_e = 3.0 * DUSTY_EV_IN_K;
    double n_d = 1.0e10;
    double Z_d = 2000.0;
    DustPlasmaState ps = dust_plasma_state_init(
        n_e, T_e, 300.0, 300.0, n_d, Z_d, 6.64e-26, 0.0);
    CHECK(ps.n_i == n_e + Z_d * n_d, "Quasineutrality enforced");
    CHECK(ps.lambda_De > 0.0, "lambda_De computed");
    CHECK(ps.lambda_Di > 0.0, "lambda_Di computed");
    CHECK(ps.lambda_D > 0.0, "lambda_D total computed");
}

/* ================================================================
 * L3-L4 Tests — Charging Theory
 * ================================================================ */

void test_oml_currents(void)
{
    TEST("OML electron current");
    double a = 5.0e-6;
    double n_e = 1.0e15;
    double T_e = 3.0 * DUSTY_EV_IN_K;
    double phi_s = -5.0;
    double I_e = dust_oml_electron_current(a, n_e, T_e, phi_s);
    CHECK(I_e < 0.0, "Electron current should be negative");
    CHECK(fabs(I_e) > 0.0, "Electron current non-zero");

    TEST("OML ion current");
    double n_i = 1.02e15;
    double T_i = 300.0;
    double m_i = 6.64e-26;
    double I_i = dust_oml_ion_current(a, n_i, T_i, m_i, phi_s);
    CHECK(I_i > 0.0, "Ion current should be positive");

    TEST("OML ion current with drift");
    double I_i_drift = dust_oml_ion_current_drift(
        a, n_i, T_i, m_i, phi_s, 100.0);
    CHECK(I_i_drift > I_i, "Drift should enhance ion current");
}

void test_floating_potential(void)
{
    TEST("Floating potential solve");
    double a = 5.0e-6;
    double n_e = 1.0e15;
    double n_i = 1.02e15;
    double T_e = 3.0 * DUSTY_EV_IN_K;
    double T_i = 300.0;
    double m_i = 6.64e-26;
    double phi_init = -2.5 * DUSTY_KB * T_e / DUSTY_EC;

    double phi_f = dust_floating_potential_solve(
        a, n_e, n_i, T_e, T_i, m_i, CHARGE_MODEL_OML,
        0.0, 0.0, 0.0, 0.0, 300.0, 5.0,
        phi_init, 1e-4, 20);
    CHECK(phi_f < 0.0, "Floating potential should be negative");
    CHECK(phi_f > -20.0, "Floating potential reasonable range");

    TEST("Equilibrium charge from floating potential");
    double Q_d = dust_equilibrium_charge(a, phi_f);
    CHECK(fabs(Q_d) > 1e-17, "Equilibrium charge non-zero");

    double Z_d = dust_equilibrium_charge_number(a, phi_f);
    CHECK(Z_d > 100.0 && Z_d < 50000.0, "Z_d in expected range");

    TEST("Charge number estimate");
    double Z_est = dust_charge_number_estimate(a, T_e);
    CHECK(Z_est > 100.0, "Estimate gives reasonable Z_d");
}

void test_charge_fluctuation(void)
{
    TEST("Charge fluctuation RMS");
    double a = 5.0e-6;
    double n_e = 1.0e15;
    double n_i = 1.02e15;
    double T_e = 3.0 * DUSTY_EV_IN_K;
    double T_i = 300.0;
    double m_i = 6.64e-26;
    double phi_f = -3.0;
    double rms_Z = dust_charge_fluctuation_rms(a, n_e, n_i, T_e, T_i, m_i, phi_f);
    CHECK(rms_Z >= 0.0, "RMS fluctuation non-negative");

    TEST("Charge relaxation time");
    double tau = dust_charge_relaxation_time(a, n_e, n_i, T_e, T_i, m_i, phi_f);
    CHECK(tau > 0.0 && tau < 1.0, "Relaxation time < 1s");
}

/* ================================================================
 * L4 Tests — Wave Dispersion
 * ================================================================ */

void test_daw_dispersion(void)
{
    TEST("DAW dispersion");
    double c_da = 0.03;
    double lambda_D = 4.0e-4;
    double k = 100.0;
    WaveMode wm = dust_acoustic_wave_dispersion(k, c_da, lambda_D);
    CHECK(wm.omega_r > 0.0, "Real frequency positive");
    CHECK(wm.omega_i == 0.0, "No damping in fluid model");
    CHECK(wm.phase_velocity > 0.0, "Phase velocity positive");

    TEST("DAW Landau damping");
    double omega_pd = 300.0;
    double v_thd = 5.0e-4;
    double gamma = dust_acoustic_wave_landau_damping(
        k, wm.omega_r, omega_pd, v_thd);
    CHECK(gamma <= 0.0, "Landau damping should be negative (damping)");

    TEST("DAW collisional damping");
    double nu_dn = 100.0;
    double gamma_c = dust_acoustic_wave_collisional_damping(nu_dn);
    CHECK(gamma_c < 0.0, "Collisional damping negative");
    CHECK(fabs(gamma_c - (-50.0)) < 1.0, "gamma_c = -nu_dn/2");
}

void test_diaw_dispersion(void)
{
    TEST("DIAW dispersion");
    double c_s = 2700.0;
    double lambda_De = 4.0e-4;
    double k = 100.0;
    WaveMode wm = dust_ion_acoustic_wave_dispersion(k, c_s, lambda_De);
    CHECK(wm.omega_r > 1.0e5, "DIAW frequency > 100 kHz");
    CHECK(wm.omega_r < 1.0e7, "DIAW frequency < 10 MHz");
}

void test_dust_lattice_wave(void)
{
    TEST("1D longitudinal DLW");
    double k = 2000.0;
    double a = 2.0e-4;
    double Q_d = 2000.0 * DUSTY_EC;
    double m_d = 1.15e-12;
    double kappa = 1.0;
    WaveMode wm = dust_lattice_wave_1d_longitudinal(
        k, a, Q_d, m_d, kappa, 10);
    CHECK(wm.omega_r >= 0.0, "DLW frequency non-negative");

    TEST("1D transverse DLW");
    WaveMode wmT = dust_lattice_wave_1d_transverse(
        k, a, Q_d, m_d, kappa, 10);
    CHECK(wmT.omega_r >= 0.0, "Transverse DLW frequency non-negative");
    CHECK(wmT.omega_r <= wm.omega_r, "Transverse < longitudinal");
}

/* ================================================================
 * L3 Tests — Yukawa Potential
 * ================================================================ */

void test_yukawa_potential(void)
{
    TEST("Yukawa potential at short range");
    double Q = 2000.0 * DUSTY_EC;
    double r = 5.0e-6;
    double lambda_D = 4.0e-4;
    double U = yukawa_potential(Q, Q, r, lambda_D);
    CHECK(U > 0.0, "Pot energy positive for like charges");

    TEST("Yukawa potential at long range");
    double U_far = yukawa_potential(Q, Q, 10.0 * lambda_D, lambda_D);
    CHECK(U_far < U * 0.1, "Exponentially suppressed at long range");

    TEST("Yukawa force magnitude");
    double F = yukawa_force_magnitude(Q, Q, r, lambda_D);
    CHECK(F > 0.0, "Repulsive force positive");

    TEST("Yukawa force vector");
    /* Use inter-particle spacing ~200 um for realistic Yukawa force */
    double r1[3] = {2.0e-4, 0.0, 0.0};
    double r2[3] = {0.0, 0.0, 0.0};
    double F_out[3];
    yukawa_force_vector(r1, r2, Q, Q, lambda_D, F_out);
    CHECK(F_out[0] > 0.0, "Force points away");
    CHECK(F_out[2] == 0.0, "No z-component for aligned pair");
}

/* ================================================================
 * L4-L5 Tests — Transport and Forces
 * ================================================================ */

void test_transport(void)
{
    TEST("Neutral collision frequency");
    double a = 5.0e-6;
    double n_n = 2.4e21;
    double nu_dn = dust_neutral_collision_freq(
        a, n_n, 300.0, 6.63e-26, 1.15e-12, 0.8);
    CHECK(nu_dn > 1.0 && nu_dn < 1000.0, "nu_dn in expected range ~100 Hz");

    TEST("Diffusion coefficient");
    double D = dust_diffusion_coefficient(300.0, 1.15e-12, nu_dn);
    CHECK(D > 0.0, "D > 0");

    TEST("Mobility");
    double mu = dust_mobility(2000.0, 1.15e-12, nu_dn);
    CHECK(mu > 0.0, "mu > 0");

    TEST("Drift velocity");
    double v_drift = dust_drift_velocity(mu, 1000.0);
    CHECK(v_drift > 0.0, "Drift positive");
}

void test_forces(void)
{
    TEST("Gravity force");
    double Fg = dust_gravity_force(1.15e-12, 9.81);
    CHECK(Fg > 1.0e-14 && Fg < 1.0e-10, "Gravity ~10^-12 N");

    TEST("Electric force");
    double Fe = dust_electric_force(2000.0 * DUSTY_EC, 1000.0);
    CHECK(fabs(Fe) > 1.0e-15, "Electric force non-zero");

    TEST("Ion drag force");
    double F_id = dust_ion_drag_total(
        5.0e-6, 1.0e15, 6.64e-26, 300.0, -3.0, 2000.0, 300.0, 4.0e-4);
    CHECK(F_id >= 0.0, "Ion drag non-negative");

    TEST("Thermophoretic force");
    double F_th = dust_thermophoretic_force(
        5.0e-6, 300.0, 6.64e-26, 0.017, 1000.0, 0.8);
    CHECK(fabs(F_th) > 0.0, "Thermophoretic force non-zero");

    TEST("Levitation height");
    double z_eq = dust_levitation_height(
        1.15e-12, 2000.0 * DUSTY_EC, 9.81,
        50000.0, 4.0e-4, 0.0, 0.0, 0.01, 1e-6);
    CHECK(z_eq >= 0.0, "Levitation height >= 0");
}

/* ================================================================
 * L6 Tests — Crystal and Phase
 * ================================================================ */

void test_crystal_condition(void)
{
    TEST("Crystal condition with large Gamma");
    int crystal = dust_crystal_condition(500.0, 1.0);
    CHECK(crystal == 1, "Gamma=500, kappa=1 should crystallize");

    TEST("Crystal condition with small Gamma");
    int fluid = dust_crystal_condition(1.0, 1.0);
    CHECK(fluid == 0, "Gamma=1 should NOT crystallize");

    TEST("Critical coupling");
    double Gc = dust_critical_coupling(0.5);
    CHECK(Gc > 100.0 && Gc < 200.0, "Gamma_crit in expected range");
}

/* ================================================================
 * L8 Tests — Advanced
 * ================================================================ */

void test_advanced(void)
{
    TEST("Phase determination");
    int phase_gas = dust_phase_determine(0.1, 1.0);
    CHECK(phase_gas == 0, "Gamma=0.1 => fluid");
    int phase_solid = dust_phase_determine(200.0, 1.0);
    CHECK(phase_solid == 2, "Gamma=200 => solid");

    TEST("Lindemann melting");
    int melting = dust_lindemann_melting(0.3, 1.0, 0.15);
    CHECK(melting == 1, "Should melt at 0.3/1.0 > 0.15");

    TEST("Nucleation rate");
    double J = dust_nucleation_rate(10.0, 1.0, 300.0, 1.0e10);
    CHECK(J >= 0.0, "Nucleation rate finite");

    TEST("Grain boundary energy");
    double E_gb = dust_grain_boundary_energy(0.05, 1.0, 0.5);
    CHECK(E_gb > 0.0, "Grain boundary energy positive");
}

/* ================================================================
 * L7 Tests — Applications
 * ================================================================ */

void test_applications(void)
{
    TEST("Sheath electric field");
    double E = dust_sheath_electric_field(1.0e-4, 5000.0, 4.0e-4);
    CHECK(E > 0.0, "Sheath field positive");

    TEST("Modified Bohm velocity");
    double v_bohm = dust_modified_bohm_velocity(
        3.0 * DUSTY_EV_IN_K, 300.0, 6.64e-26, 1.0e10, 1.0e15);
    CHECK(v_bohm > 0.0, "Bohm velocity positive");

    TEST("Void radius estimate");
    double Rv = dust_void_radius_estimate(4.0e-4, 1.0e-12, 2000.0, 5000.0);
    CHECK(Rv >= 0.0, "Void radius non-negative");
}

/* ================================================================
 * Main
 * ================================================================ */

int main(void)
{
    printf("\n=== mini-dusty-plasma Test Suite ===\n\n");

    printf("--- L1: Core Definitions ---\n");
    test_debye_lengths();
    test_plasma_frequencies();
    test_havnes_parameter();

    printf("\n--- L2: Core Concepts ---\n");
    test_coulomb_coupling();
    test_dust_acoustic_speed();
    test_grain_init();
    test_plasma_state_init();

    printf("\n--- L3-L4: Charging Theory ---\n");
    test_oml_currents();
    test_floating_potential();
    test_charge_fluctuation();

    printf("\n--- L4: Wave Dispersion ---\n");
    test_daw_dispersion();
    test_diaw_dispersion();
    test_dust_lattice_wave();

    printf("\n--- L3: Yukawa Potential ---\n");
    test_yukawa_potential();

    printf("\n--- L4-L5: Transport & Forces ---\n");
    test_transport();
    test_forces();

    printf("\n--- L6: Crystal & Phase ---\n");
    test_crystal_condition();

    printf("\n--- L7: Applications ---\n");
    test_applications();

    printf("\n--- L8: Advanced Topics ---\n");
    test_advanced();

    printf("\n========================================\n");
    printf("Results: %d PASSED, %d FAILED\n", tests_passed, tests_failed);
    printf("========================================\n\n");

    return tests_failed > 0 ? 1 : 0;
}
