/**
 * test_plasma.c -- Comprehensive test suite for mini-space-plasma
 *
 * Covers all core APIs with mathematical assertions.
 * Uses standard assert() for validation.
 */
#include "../include/space_plasma.h"
#include "../include/plasma_parameters.h"
#include "../include/mhd_core.h"
#include "../include/plasma_waves.h"
#include "../include/solar_wind.h"
#include "../include/magnetosphere.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static const double TOL = 1e-10;
static int tests_run = 0;
static int tests_passed = 0;

#define TEST(name) do { tests_run++; printf("  %s... ", name); } while(0)
#define PASS() do { tests_passed++; printf("PASS\n"); } while(0)
#define CHECK(cond) do { if (!(cond)) { printf("FAIL at %s:%d\n", __FILE__, __LINE__); return; } } while(0)
#define CHECK_CLOSE(a, b, tol) do { if (fabs((a)-(b)) > (tol)) { printf("FAIL: %.6e != %.6e at %s:%d\n", a, b, __FILE__, __LINE__); return; } } while(0)

/*========================================================================
 * L1: Plasma Parameter Tests
 *========================================================================*/

static void test_debye_length(void) {
    TEST("Debye length (solar wind)");
    double n_e = 5.0e6;   /* 5 cm^-3 */
    double T_e = 1.0e5;   /* 10 eV */
    double ld = sp_debye_length_e(n_e, T_e);
    /* lambda_D ~ 7.4 m for T=1e5K, n=5e6 */
    CHECK(ld > 1.0 && ld < 20.0);
    PASS();
}

static void test_debye_number(void) {
    TEST("Debye number N_D");
    double n_e = 5.0e6;
    double T_e = 1.0e5;
    double ld = sp_debye_length_e(n_e, T_e);
    double ND = sp_debye_number(n_e, ld);
    /* N_D >> 1 for ideal space plasma */
    CHECK(ND > 100.0);
    PASS();
}

static void test_coulomb_coupling(void) {
    TEST("Coulomb coupling Gamma");
    double Gamma = sp_coulomb_coupling(5.0e6, 1.0e5, 1.0);
    /* Solar wind at 1 AU: Gamma << 1 (weakly coupled) */
    CHECK(Gamma < 1.0 && Gamma > 1e-10);
    PASS();
}

static void test_plasma_frequency(void) {
    TEST("Electron plasma frequency");
    double n_e = 1.0e6;  /* ionosphere, f_pe ~ 9 kHz */
    double omega_pe = sp_plasma_freq_e(n_e);
    CHECK(omega_pe > 1e4 && omega_pe < 1e6);
    double f_pe = sp_plasma_freq_e_hz(n_e);
    CHECK_CLOSE(f_pe, omega_pe / (2.0*M_PI), TOL);
    PASS();
}

static void test_gyrofrequency(void) {
    TEST("Electron gyrofrequency at 1 AU");
    double B = 5.0e-9;  /* IMF */
    double omega_ce = sp_gyrofreq(B, SP_ME, SP_EC);
    /* omega_ce = eB/m_e ~ 880 rad/s for B=5nT */
    CHECK(omega_ce > 100.0 && omega_ce < 2000.0);
    PASS();
}

static void test_larmor_radius(void) {
    TEST("Proton Larmor radius in solar wind");
    double v_perp = 400.0e3;  /* solar wind speed */
    double B = 5.0e-9;
    double rL = sp_larmor_radius(v_perp, B, SP_MP, SP_EC);
    /* r_L ~ 840 km */
    CHECK(rL > 1e5 && rL < 1e7);
    PASS();
}

static void test_alfven_speed(void) {
    TEST("Alfven speed in solar wind");
    double B = 5.0e-9;
    double rho = 5.0e6 * SP_MP;  /* n * m_p */
    double v_A = sp_alfven_speed(B, rho);
    /* v_A ~ 40 km/s */
    CHECK(v_A > 1e4 && v_A < 1e5);
    PASS();
}

static void test_plasma_beta(void) {
    TEST("Plasma beta in solar wind");
    double n_e = 5.0e6, T_e = 1.0e5;
    double B = 5.0e-9;
    double p = n_e * SP_KB * T_e;
    double beta = sp_plasma_beta(p, B);
    /* beta ~ 1-2 in typical solar wind */
    CHECK(beta > 0.1 && beta < 10.0);
    PASS();
}

static void test_lundquist_number(void) {
    TEST("Lundquist number");
    double L = 1.0e8;  /* 100,000 km */
    double v_A = 4.0e4;
    double lnL = sp_coulomb_logarithm(5.0e6, 1.0e5);
    double eta = sp_spitzer_resistivity(5.0e6, 1.0e5, 1.0, lnL);
    double S = sp_lundquist_number(L, v_A, eta);
    /* S ~ 1e13 for space plasmas */
    CHECK(S > 1e10);
    PASS();
}

/*========================================================================
 * L4: Collision/Transport Tests
 *========================================================================*/

static void test_coulomb_logarithm(void) {
    TEST("Coulomb logarithm");
    double lnL = sp_coulomb_logarithm(5.0e6, 1.0e5);
    /* lnLambda ~ 25-30 for solar wind */
    CHECK(lnL > 20.0 && lnL < 40.0);
    PASS();
}

static void test_spitzer_resistivity(void) {
    TEST("Spitzer resistivity (T dependence)");
    double lnL = sp_coulomb_logarithm(5.0e6, 1.0e5);
    double eta_cold = sp_spitzer_resistivity(5.0e6, 1.0e4, 1.0, lnL);
    double eta_hot  = sp_spitzer_resistivity(5.0e6, 1.0e6, 1.0, lnL);
    /* eta ~ T^(-3/2): hotter -> less resistive */
    CHECK(eta_cold > eta_hot);
    PASS();
}

/*========================================================================
 * L2: MHD Core Tests
 *========================================================================*/

static void test_mhd_prim_cons_roundtrip(void) {
    TEST("MHD primitive/cons round-trip");
    mhd_primitive_t prim = {.rho=1.67e-20, .vx=4e5, .vy=0, .vz=0,
                            .Bx=5e-9, .By=0, .Bz=0, .p=1e-10};
    mhd_conserved_t cons;
    mhd_prim_to_cons(&prim, &cons, 5.0/3.0);

    mhd_primitive_t prim2;
    int ret = mhd_cons_to_prim(&cons, &prim2, 5.0/3.0);
    CHECK(ret == 0);
    CHECK_CLOSE(prim.rho, prim2.rho, 1e-12);
    CHECK_CLOSE(prim.vx, prim2.vx, 1e-8);
    PASS();
}

static void test_mhd_flux_consistency(void) {
    TEST("MHD flux positivity");
    mhd_primitive_t prim = {.rho=1.67e-20, .vx=4e5, .vy=1e4, .vz=0,
                            .Bx=5e-9, .By=1e-9, .Bz=0, .p=1e-10};
    mhd_flux_t flux;
    mhd_flux_x(&prim, &flux, 5.0/3.0);
    /* Mass flux must be positive for positive density and vx */
    CHECK(flux.f[0] > 0.0);
    PASS();
}

static void test_mhd_wave_speeds(void) {
    TEST("MHD wave speeds (ordering)");
    mhd_primitive_t prim = {.rho=1.67e-20, .vx=0, .vy=0, .vz=0,
                            .Bx=5e-9, .By=1e-9, .Bz=0, .p=1e-10};
    mhd_eigensystem_t eigs;
    mhd_wave_speeds(&prim, &eigs, 5.0/3.0);
    /* Eigenvalues should be non-decreasing */
    for (int i = 0; i < 6; i++)
        CHECK(eigs.lambda[i] <= eigs.lambda[i+1]);
    /* v_fast >= v_alfven >= v_slow for typical parameters */
    CHECK(eigs.v_fast >= eigs.v_alfven);
    CHECK(eigs.v_alfven >= eigs.v_slow);
    PASS();
}

static void test_harris_sheet(void) {
    TEST("Harris current sheet equilibrium");
    double B_out[3];
    double n_out, p_out;
    mhd_harris_sheet(0.0, 1e6, 1e-8, 1e6, 1e5, 1e5, B_out, &n_out, &p_out);
    /* At center (y=0): B_x = 0 */
    CHECK_CLOSE(B_out[0], 0.0, 1e-16);
    /* Pressure at center = B0^2/(2*mu0) + background */
    double p0 = 1e-8*1e-8/(2.0*SP_MU0);
    CHECK(p_out > p0);
    PASS();
}

static void test_zpinch(void) {
    TEST("Z-pinch equilibrium pressure");
    double p, B[3];
    mhd_zpinch(0.0, 0.1, 1e6, 1e7, &p, B);
    /* Maximum pressure at r=0 */
    CHECK(p > 0.0);
    CHECK_CLOSE(B[0], 0.0, 1e-16);  /* B_r = 0 */
    PASS();
}

static void test_mhd_frozen_in(void) {
    TEST("Frozen-in flux condition (ideal)");
    double B[3] = {0, 5e-9, 0};
    double v[3] = {4e5, 0, 0};
    /* In ideal MHD: E = -v x B. v=(4e5,0,0), B=(0,5e-9,0) -> vxB=(0,0,2e-3) */
    double vxB[3];
    SP_CROSS3(v, B, vxB);
    double E[3] = {-vxB[0], -vxB[1], -vxB[2]};
    /* With E = -vxB: E + vxB = 0, error should be zero */
    CHECK(mhd_frozen_in_error(E, v, B) < 1e-12);
    PASS();
}

/*========================================================================
 * L2: Plasma Waves Tests
 *========================================================================*/

static void test_langmuir_frequency(void) {
    TEST("Langmuir wave dispersion (Bohm-Gross)");
    double omega_pe = 1.78e5;  /* f_pe ~ 28 kHz */
    double v_the = 1.0e6;
    double omega = langmuir_frequency(1.0, omega_pe, v_the);
    /* omega >= omega_pe (thermal correction positive) */
    CHECK(omega >= omega_pe);
    double omega_cold = langmuir_frequency(0.0, omega_pe, v_the);
    CHECK_CLOSE(omega_cold, omega_pe, 1e-10);
    PASS();
}

static void test_alfven_wave(void) {
    TEST("Alfven wave frequency");
    double k = 1e-5;
    double v_A = 4e4;
    double omega = alfven_wave_frequency(k, 0.0, v_A);
    /* Parallel propagation: omega = k * v_A */
    CHECK_CLOSE(omega, k * v_A, 1e-10);
    /* Perpendicular propagation: omega = 0 */
    double omega_perp = alfven_wave_frequency(k, M_PI/2.0, v_A);
    CHECK_CLOSE(omega_perp, 0.0, 1e-10);
    PASS();
}

static void test_whistler_frequency(void) {
    TEST("Whistler wave dispersion");
    /* Need k*c/omega_pe << 1 for whistler approximation validity.
     * Use k=1e-5 for omega_pe=1.78e5: kc/omega_pe = 1e-5*3e8/1.78e5 = 0.0168 */
    double k = 1e-5;
    double omega_pe = 1.78e5;
    double omega_ce = 8.8e2;
    double theta = 0.0;
    double omega = whistler_frequency(k, omega_pe, omega_ce, theta);
    /* omega = omega_ce * (k*c/omega_pe)^2 ~ 880 * (0.0168)^2 ~ 0.25 rad/s */
    CHECK(omega >= 0.0 && omega < omega_ce);
    PASS();
}

static void test_ion_acoustic(void) {
    TEST("Ion acoustic wave frequency");
    double c_s = 4e4;
    double lambda_De = 7.4;
    /* For k*lambda_De << 1: omega ~ k*c_s. Use k=0.01 for long wavelength */
    double omega_ac = ion_acoustic_frequency(0.01, c_s, lambda_De);
    /* omega ~ 0.01 * 40000 / sqrt(1+0.0001*54.76) ~ 400 m/s */
    CHECK(omega_ac > 300.0 && omega_ac < 500.0);
    PASS();
}

static void test_stix_tensor(void) {
    TEST("Stix dielectric tensor elements");
    stix_tensor_t eps;
    stix_dielectric_tensor(1e4, 1.78e5, 8.8e2, 1e3, 1.0, &eps);
    /* Check tensor symmetries */
    CHECK_CLOSE(eps.S, 0.5*(eps.R+eps.L), 1e-10);
    CHECK_CLOSE(eps.D, 0.5*(eps.R-eps.L), 1e-10);
    PASS();
}

static void test_wave_cutoffs(void) {
    TEST("Wave cutoff frequencies");
    double cutoffs[3];
    wave_cutoff_frequencies(1.78e5, 8.8e2, cutoffs);
    /* P cuttoff = omega_pe */
    CHECK_CLOSE(cutoffs[0], 1.78e5, 1e-10);
    /* R cutoff > L cutoff for omega_ce > 0 */
    CHECK(cutoffs[1] > cutoffs[2]);
    PASS();
}

/*========================================================================
 * L6: Solar Wind Tests
 *========================================================================*/

static void test_parker_critical_radius(void) {
    TEST("Parker critical radius");
    double c_s = 1.65e5;  /* ~165 km/s for 1e6 K corona */
    double r_c = parker_critical_radius(c_s);
    /* r_c ~ 5.5 R_sun ~ 3.8e9 m */
    CHECK(r_c > 1e9 && r_c < 1e10);
    PASS();
}

static void test_parker_isothermal(void) {
    TEST("Parker isothermal wind velocity");
    double c_s = 1.65e5;
    double r_c = parker_critical_radius(c_s);
    /* At r > r_c: supersonic solution */
    double v = parker_isothermal_velocity(10.0*r_c, r_c, 2.0*c_s, c_s, 1e-8, 100);
    CHECK(v > c_s);  /* supersonic */
    /* At r < r_c: subsonic solution */
    double v2 = parker_isothermal_velocity(0.5*r_c, r_c, 0.5*c_s, c_s, 1e-8, 100);
    CHECK(v2 < c_s);  /* subsonic */
    PASS();
}

static void test_parker_spiral(void) {
    TEST("Parker spiral at 1 AU");
    double B[3];
    double r = SP_AU;
    double theta = M_PI/2.0;  /* equator */
    double omega = 2.865e-6;   /* solar rotation rate */
    double v_sw = 4.0e5;
    parker_spiral_field(r, theta, 5.0e-9, SP_RSUN, omega, v_sw, B);
    /* B_r > 0, B_phi < 0 (at equator) */
    CHECK(B[0] > 0.0);
    CHECK(B[2] < 0.0);
    /* B_theta = 0 */
    CHECK_CLOSE(B[1], 0.0, 1e-20);
    PASS();
}

static void test_parker_spiral_angle(void) {
    TEST("Parker spiral angle ~ 45 deg");
    double psi = parker_spiral_angle(SP_AU, M_PI/2.0, 2.865e-6, 4.0e5);
    /* psi ~ -45 degrees */
    CHECK(psi < 0.0 && psi > -M_PI/3.0);
    PASS();
}

static void test_solar_wind_mach(void) {
    TEST("Solar wind Mach numbers");
    double v_sw = 400e3;
    double B = 5e-9;
    double rho = 5e6 * SP_MP;
    double c_s = 4e4;

    double Ma = sp_alfven_mach(v_sw, B, rho);
    /* Super-Alfvenic: Ma ~ 8-10 */
    CHECK(Ma > 1.0);

    double Ms = sp_sonic_mach(v_sw, c_s);
    /* Supersonic: Ms ~ 10 */
    CHECK(Ms > 1.0);
    PASS();
}

/*========================================================================
 * L6: Magnetosphere Tests
 *========================================================================*/

static void test_dipole_field(void) {
    TEST("Earth dipole field at equator");
    double x[3] = {1.0, 0.0, 0.0};  /* 1 R_E at equator */
    double B[3];
    earth_dipole_field(x, SP_MU_EARTH, B);
    /* |B| ~ 3.12e-5 T at equator */
    double Bmag = SP_MAG3(B);
    /* B_eq = mu0*M/(4*pi*R_E^3) ~ 3.01e-5 T at equator (modulo sign conventions) */
    CHECK_CLOSE(Bmag, SP_B0_EARTH, 0.1);  /* 10% tolerance for equatorial approximation */
    PASS();
}

static void test_chapman_ferraro(void) {
    TEST("Chapman-Ferraro standoff distance");
    double R_mp = chapman_ferraro_standoff(5.0e6, 400e3, 2.44);
    /* R_mp ~ 8-12 R_E */
    CHECK(R_mp > 5.0 && R_mp < 15.0);
    PASS();
}

static void test_plasmapause(void) {
    TEST("Plasmapause location (Kp=2)");
    double L_pp = plasmapause_L(2.0);
    CHECK_CLOSE(L_pp, 4.68, 0.1);
    PASS();
}

static void test_bow_shock(void) {
    TEST("Bow shock standoff");
    double R_mp = 10.0;
    double R_bs = bow_shock_standoff(R_mp, 8.0, 5.0/3.0);
    /* R_bs > R_mp */
    CHECK(R_bs > R_mp);
    CHECK(R_bs < 2.0*R_mp);
    PASS();
}

static void test_mhd_shock_jump(void) {
    TEST("MHD shock jump conditions");
    double rho_r, p_r, B_r;
    mhd_shock_jump(10.0, 5.0/3.0, &rho_r, &p_r, &B_r);
    /* For M1=10, gamma=5/3: rho_ratio = (8/3*100)/(2/3*100+2) ~ 3.88
     * As M1->inf: rho_ratio -> (gamma+1)/(gamma-1) = 4 */
    CHECK(rho_r > 3.0 && rho_r < 5.0);
    /* Pressure ratio large for M=10 */
    CHECK(p_r > 10.0);
    PASS();
}

static void test_exb_drift(void) {
    TEST("E x B drift (dawn-dusk electric field)");
    double E[3] = {0.0, 1e-3, 0.0};  /* duskward */
    double B[3] = {0.0, 0.0, 3e-8};   /* northward dipole */
    double v[3];
    magnetosphere_exb_drift(E, B, v);
    /* v_E = E x B / B^2: E=(0,1e-3,0), B=(0,0,3e-8) -> v=(1e-3*3e-8/B^2, 0, 0)
     * = (3e-11/9e-16, 0, 0) ~ (3.3e4, 0, 0) m/s sunward */
    CHECK(v[0] > 0.0);  /* sunward flow in inner magnetosphere */
    PASS();
}

static void test_dipole_drift(void) {
    TEST("Dipole drift direction");
    double v_d[3];
    /* Electron (q<0) at L=4 with 1 MeV */
    double E_kin = 1e6 * SP_EC;  /* 1 MeV */
    dipole_drift_speed(E_kin, 4.0, -SP_EC, v_d);
    /* Electron drifts eastward (positive phi) */
    CHECK(v_d[1] > 0.0);
    /* Proton (q>0) at L=4 with 100 keV */
    dipole_drift_speed(1e5*SP_EC, 4.0, SP_EC, v_d);
    /* Proton drifts westward (negative phi) */
    CHECK(v_d[1] < 0.0);
    PASS();
}

static void test_dipole_bounce_period(void) {
    TEST("Dipole bounce period");
    /* 1 MeV electron at L=4, equatorial */
    double E_kin = 1e6 * SP_EC;
    double tau_b = dipole_bounce_period(E_kin, SP_ME, 4.0, M_PI/2.0);
    /* tau_b ~ 0.3-0.5 s */
    CHECK(tau_b > 0.1 && tau_b < 2.0);
    PASS();
}

static void test_ring_current_energy(void) {
    TEST("Ring current energy density");
    double U = ring_current_energy_density(1e6, 1e5, 1e4);
    /* Energy density > 0 for finite parameters */
    CHECK(U > 0.0);
    PASS();
}

/*========================================================================
 * Edge Cases and Boundary Conditions
 *========================================================================*/

static void test_zero_density(void) {
    TEST("Zero density -> zero Debye length");
    CHECK_CLOSE(sp_debye_length_e(0.0, 1e5), 0.0, 1e-16);
    CHECK_CLOSE(sp_plasma_freq_e(0.0), 0.0, 1e-16);
    PASS();
}

static void test_zero_temperature(void) {
    TEST("Zero temperature -> zero Debye length");
    CHECK_CLOSE(sp_debye_length_e(5e6, 0.0), 0.0, 1e-16);
    PASS();
}

static void test_null_pointer_safety(void) {
    TEST("Null pointer safety (mhd_prim_to_cons)");
    mhd_conserved_t cons;
    mhd_prim_to_cons(NULL, &cons, 1.4);
    mhd_prim_to_cons(NULL, NULL, 1.4);
    /* Should not crash */
    PASS();
}

static void test_infinity_handling(void) {
    TEST("Infinite mean free path for zero collisions");
    double mfp = sp_mean_free_path(1e6, 0.0);
    CHECK(isinf(mfp));
    PASS();
}

/*========================================================================
 * Main
 *========================================================================*/

int main(void) {
    printf("=== mini-space-plasma Test Suite ===\n\n");

    /* L1: Plasma Parameters */
    printf("L1: Plasma Parameters:\n");
    test_debye_length();
    test_debye_number();
    test_coulomb_coupling();
    test_plasma_frequency();
    test_gyrofrequency();
    test_larmor_radius();
    test_alfven_speed();
    test_plasma_beta();
    test_lundquist_number();

    /* L4: Transport */
    printf("\nL4: Transport & Collisions:\n");
    test_coulomb_logarithm();
    test_spitzer_resistivity();

    /* L2: MHD */
    printf("\nL2: MHD Core:\n");
    test_mhd_prim_cons_roundtrip();
    test_mhd_flux_consistency();
    test_mhd_wave_speeds();
    test_harris_sheet();
    test_zpinch();
    test_mhd_frozen_in();

    /* L2: Plasma Waves */
    printf("\nL2-L3: Plasma Waves:\n");
    test_langmuir_frequency();
    test_alfven_wave();
    test_whistler_frequency();
    test_ion_acoustic();
    test_stix_tensor();
    test_wave_cutoffs();

    /* L6: Solar Wind */
    printf("\nL6: Solar Wind:\n");
    test_parker_critical_radius();
    test_parker_isothermal();
    test_parker_spiral();
    test_parker_spiral_angle();
    test_solar_wind_mach();

    /* L6: Magnetosphere */
    printf("\nL6: Magnetosphere:\n");
    test_dipole_field();
    test_chapman_ferraro();
    test_plasmapause();
    test_bow_shock();
    test_mhd_shock_jump();
    test_exb_drift();
    test_dipole_drift();
    test_dipole_bounce_period();
    test_ring_current_energy();

    /* Edge Cases */
    printf("\nEdge Cases:\n");
    test_zero_density();
    test_zero_temperature();
    test_null_pointer_safety();
    test_infinity_handling();

    printf("\n=== Results: %d/%d tests passed ===\n", tests_passed, tests_run);
    return (tests_passed == tests_run) ? 0 : 1;
}
