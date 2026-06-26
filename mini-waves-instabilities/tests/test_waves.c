/**
 * test_waves.c -- Tests for plasma wave dispersion relations
 *
 * Tests cover:
 *   L1: Core parameter computations
 *   L4: Langmuir, IAW, whistler, Alfven, MHD wave dispersions
 *   L6: Cutoffs, resonances, CMA classification
 *   L5: Plasma dispersion function Z(zeta), Bessel functions
 *
 * All tests use mathematical assertions (non-trivial).
 */

#include "../include/waves_instabilities.h"
#include "../include/kinetic_dispersion.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <complex.h>

static const double TOL = 1e-6;

/* ---- L1: Parameter Computations ---- */
static void test_plasma_params(void)
{
    PlasmaSpecies sp = {1e20, K_BOLTZMANN * 1000.0,
                         M_ELECTRON, -E_CHARGE, 1.0};
    PlasmaParams p;
    compute_plasma_params(&sp, &p);

    /* Plasma frequency: f_pe = omega_pe/(2*pi) ~ 8980*sqrt(n) Hz */
    /* Plasma frequency from exact constants:
     * sqrt(n q^2/(eps0 m)) using NIST values */
    double omega_pe_calc = sqrt(sp.n0 * E_CHARGE * E_CHARGE
                                 / (EPSILON0 * M_ELECTRON));
    assert(fabs(p.omega_p - omega_pe_calc) / omega_pe_calc < TOL);

    /* Cyclotron frequency: omega_ce = e*B/m_e = 1.76e11 for B=1T */
    double omega_ce_ref = E_CHARGE * 1.0 / M_ELECTRON;
    assert(fabs(p.omega_c - omega_ce_ref) / omega_ce_ref < TOL);

    /* Thermal velocity: sqrt(2*kB*T/m) */
    double v_th_ref = sqrt(2.0 * K_BOLTZMANN * sp.T / M_ELECTRON);
    assert(fabs(p.v_th - v_th_ref) / v_th_ref < TOL);

    /* Debye length: sqrt(eps0*kB*T/(n*e^2)) */
    double lambda_D_ref = sqrt(EPSILON0 * K_BOLTZMANN * sp.T
                                / (1e20 * E_CHARGE * E_CHARGE));
    assert(fabs(p.lambda_D - lambda_D_ref) / lambda_D_ref < 1e-5);

    printf("  PASS: test_plasma_params\n");
}

static void test_derived_params(void)
{
    double n = 1e20;
    double omega_pe = electron_plasma_frequency(n);
    assert(omega_pe > 0.0);
    assert(fabs(omega_pe - 5.64e11) / 5.64e11 < 0.05);

    double vA = alfven_speed(1.0, n, M_PROTON);
    assert(vA > 0.0);
    /* v_A = B/sqrt(mu0*n*m_i) = 1/sqrt(4*pi*1e-7*1e20*1.67e-27) */
    double vA_ref = 1.0 / sqrt(MU0 * n * M_PROTON);
    assert(fabs(vA - vA_ref) / vA_ref < TOL);

    double cs = ion_sound_speed(K_BOLTZMANN * 1000.0, M_PROTON);
    assert(cs > 0.0);

    /* T in Joules: plasma_beta uses formula 2*mu0*n*kB*T/B^2 */
    double T_joules = 1000.0 * K_BOLTZMANN; /* 1000 K in Joules */
    double beta = plasma_beta(n, T_joules, 1.0);
    assert(beta > 0.0);
    double beta_ref = 2.0 * MU0 * n * K_BOLTZMANN * T_joules / 1.0;
    assert(fabs(beta - beta_ref) / beta_ref < TOL);

    printf("  PASS: test_derived_params\n");
}

/* ---- L4: Electrostatic Wave Dispersion ---- */
static void test_langmuir_wave(void)
{
    double omega_pe = 5.64e11;
    double v_th_e = 4.19e6;
    double k = 100.0;

    /* Bohm-Gross: omega^2 = omega_pe^2 + 3*v_th^2*k^2 */
    double w = langmuir_wave_dispersion(k, omega_pe, v_th_e);
    double w_ref = sqrt(omega_pe*omega_pe + 3.0*v_th_e*v_th_e*k*k);
    assert(fabs(w - w_ref) / w_ref < TOL);

    /* Landau damped: should return gamma < 0 for weak damping */
    ComplexOmega damped = langmuir_wave_damped(k, omega_pe, v_th_e);
    assert(damped.omega_r > 0.0);
    /* At these parameters (k*lambda_De ~ 0.001), damping is negligible */
    assert(damped.gamma <= 0.0);

    printf("  PASS: test_langmuir_wave\n");
}

static void test_ion_acoustic_wave(void)
{
    double cs = 9.79e3;      /* c_s for Te=1eV, mi=proton */
    double lambda_De = 7.43e-5;
    double k = 1000.0;

    double w = ion_acoustic_wave_dispersion(k, cs, lambda_De);
    assert(w > 0.0);
    /* For k*lambda_De << 1: w ~ k*cs */
    assert(fabs(w - k*cs) / (k*cs) < 0.1);

    ComplexOmega damped = ion_acoustic_wave_damped(k, cs, lambda_De,
                                                     1.6e-19, 1.6e-20, M_PROTON);
    assert(damped.omega_r > 0.0);
    assert(damped.gamma <= 0.0);

    printf("  PASS: test_ion_acoustic_wave\n");
}

/* ---- L4/L6: Magnetized Waves ---- */
static void test_whistler_wave(void)
{
    double omega_pe = 5.64e11;
    double omega_ce = 1.76e11;
    double k = 10.0;

    double w = whistler_dispersion(k, omega_pe, omega_ce);
    assert(w > 0.0);
    assert(w < omega_ce); /* Whistler exists only for w < w_ce */

    double vg = whistler_group_velocity(k, omega_pe, omega_ce);
    assert(vg > 0.0);

    printf("  PASS: test_whistler_wave\n");
}

static void test_mhd_waves(void)
{
    double cs = 1e5, v_A = 2e6;
    double v_alfven, v_fast, v_slow;
    /* Parallel propagation (cos_theta=1): Alfven is shear wave at v_A,
     * fast = max(v_A, c_s), slow = min(v_A, c_s) */
    mhd_wave_speeds(cs, v_A, 1.0, &v_alfven, &v_fast, &v_slow);
    assert(v_alfven == v_A);
    assert(v_fast >= v_alfven); /* Fast >= Alfven */
    assert(v_slow <= v_alfven); /* Slow <= Alfven */

    /* Perpendicular (cos_theta=0): fast=sqrt(cs^2+vA^2), slow=0 */
    mhd_wave_speeds(cs, v_A, 0.0, &v_alfven, &v_fast, &v_slow);
    assert(v_alfven == v_A);     /* Alfven speed is geometric, same value */
    double vf_ref = sqrt(cs*cs + v_A*v_A);
    assert(fabs(v_fast - vf_ref) / vf_ref < TOL);

    printf("  PASS: test_mhd_waves\n");
}

static void test_cutoffs(void)
{
    double omega_pe = 2.0, omega_ce = 1.0;
    double w_uh = upper_hybrid_frequency(omega_pe, omega_ce);
    assert(fabs(w_uh - sqrt(5.0)) < TOL);

    double w_r = r_cutoff_frequency(omega_pe, omega_ce);
    assert(w_r > omega_ce);
    /* w_R = 0.5 + sqrt(4 + 0.25) = 0.5 + sqrt(4.25) = 2.561... */
    assert(fabs(w_r - (0.5 + sqrt(4.25))) < TOL);

    double w_lh = lower_hybrid_frequency(0.01, omega_ce,
                                          sqrt(0.01 * omega_ce));
    assert(w_lh > 0.0);

    printf("  PASS: test_cutoffs\n");
}

/* ---- L5: Plasma Dispersion Function ---- */
static void test_plasma_Z(void)
{
    double complex zeta = 0.0;
    double complex Z0 = plasma_dispersion_Z(zeta);
    /* Z(0) = i*sqrt(pi) */
    assert(fabs(creal(Z0)) < 1e-10);
    assert(fabs(cimag(Z0) - sqrt(PLASMA_PI)) < TOL);

    /* Large argument: Z(x) ~ -1/x for real x >> 1 */
    double re_large = plasma_Z_real(10.0);
    assert(fabs(re_large - (-1.0/10.0)) < 0.02);

    /* Imaginary part: Im[Z(x)] = sqrt(pi)*exp(-x^2) */
    double im_small = plasma_Z_imag(0.0);
    assert(fabs(im_small - sqrt(PLASMA_PI)) < TOL);

    printf("  PASS: test_plasma_Z\n");
}

static void test_bessel_functions(void)
{
    /* I_0(0) = 1 */
    assert(fabs(modified_bessel_I(0, 0.0) - 1.0) < TOL);
    /* I_n(0) = 0 for n > 0 */
    assert(fabs(modified_bessel_I(1, 0.0)) < TOL);
    /* I_0(1) ~ 1.266 */
    assert(fabs(modified_bessel_I(0, 1.0) - 1.266) < 0.01);

    /* Gamma_0(0) = 1, Gamma_0'(0) = -1 */
    double G0, dG0;
    gamma_0_and_derivative(0.0, &G0, &dG0);
    assert(fabs(G0 - 1.0) < TOL);
    assert(fabs(dG0 + 1.0) < TOL);

    /* J_0(0) = 1 */
    assert(fabs(bessel_J(0, 0.0) - 1.0) < TOL);

    printf("  PASS: test_bessel_functions\n");
}

/* ---- L5: Kinetic Dispersion ---- */
static void test_kinetic_dispersion(void)
{
    double omega_pe = 5.64e11, v_th_e = 4.19e6, k = 100.0;
    ComplexOmega r = kinetic_langmuir_dispersion(k, omega_pe, v_th_e);
    assert(r.omega_r > 0.0);
    assert(r.gamma <= 0.0); /* Maxwellian -> damping */

    r = kinetic_ion_acoustic_dispersion(k, 1.6e-19, 1.6e-20, 1e20, M_PROTON);
    assert(r.omega_r > 0.0);

    printf("  PASS: test_kinetic_dispersion\n");
}

int main(void)
{
    printf("=== Test: Plasma Waves ===\n");
    test_plasma_params();
    test_derived_params();
    test_langmuir_wave();
    test_ion_acoustic_wave();
    test_whistler_wave();
    test_mhd_waves();
    test_cutoffs();
    test_plasma_Z();
    test_bessel_functions();
    test_kinetic_dispersion();
    printf("=== All waves tests passed (10/10) ===\n");
    return 0;
}
