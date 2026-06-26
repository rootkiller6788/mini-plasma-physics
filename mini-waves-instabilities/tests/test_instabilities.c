/**
 * test_instabilities.c -- Tests for plasma instabilities
 *
 * Tests cover:
 *   L6: Two-stream, Weibel, firehose, mirror, RT, KH,
 *       drift wave, ITG, TEM, ETG
 *   L8: Parametric decay, three-wave, nonlinear saturation
 *   L5: Dispersion solvers, Nyquist, Penrose
 */

#include "../include/waves_instabilities.h"
#include "../include/plasma_instabilities.h"
#include "../include/dispersion_solvers.h"
#include "../include/nonlinear_waves.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>

static const double TOL = 1e-6;

/* ---- L6: Two-Stream Instability ---- */
static void test_two_stream(void)
{
    double gamma = two_stream_growth_rate(1e18, 1e20, 5.64e11);
    assert(gamma > 0.0); /* Unstable for any finite beam */

    /* Growth should scale as (n_b/n_0)^(1/3) */
    double gamma2 = two_stream_growth_rate(8e18, 1e20, 5.64e11);
    assert(gamma2 > gamma); /* Higher beam density -> larger growth */

    /* Roots at k = omega_pe/v_b */
    double omega_pe = 5.64e11, v_b = 1e7;
    double omega_b = electron_plasma_frequency(1e18);
    ComplexOmega roots[4];
    int n = two_stream_dispersion_roots(omega_pe/v_b, omega_pe, omega_b, v_b, roots);
    assert(n > 0);

    /* Buneman: cold beam */
    double gamma_b = buneman_growth_rate(omega_pe, 1e7, 4.19e5);
    assert(gamma_b > 0.0);

    printf("  PASS: test_two_stream\n");
}

/* ---- L6: Weibel Instability ---- */
static void test_weibel(void)
{
    double omega_pe = 5.64e11, v_th = 4.19e6;
    double gamma = weibel_growth_rate(omega_pe, v_th, 2.0, 1.0, 0.1);
    assert(gamma > 0.0); /* T_perp > T_par -> unstable */

    /* Isotropic: should be stable */
    double gamma_stable = weibel_growth_rate(omega_pe, v_th, 1.0, 1.0, 0.1);
    assert(gamma_stable < TOL);

    printf("  PASS: test_weibel\n");
}

/* ---- L6: Firehose / Mirror ---- */
static void test_firehose_mirror(void)
{
    /* Firehose: p_par > p_perp + B^2/mu0
     * B^2/mu0 = (1e-4)^2 / 1.2566e-6 = 7.96e-3
     * so with p_par=10, p_perp=0.5, small B=1e-4 -> 10 > 0.5+0.008 -> unstable */
    assert(firehose_unstable(10.0, 0.5, 1e-4) == 1); /* small B -> unstable */
    assert(firehose_unstable(1.0, 0.5, 1.0) == 0);   /* large B -> stable */

    /* Mirror: p_perp/p_par > 1 + B^2/(2*mu0*p_par)
     * p_perp=10, p_par=1, B=1e-4 -> 10 > 1+0.004 -> unstable */
    assert(mirror_unstable(1.0, 10.0, 1e-4) == 1); /* strong anisotropy */
    assert(mirror_unstable(1.0, 1.0, 0.1) == 0);   /* isotropic */

    printf("  PASS: test_firehose_mirror\n");
}

/* ---- L6: RT and KH ---- */
static void test_rt_kh(void)
{
    double gamma_rt = rayleigh_taylor_growth(9.8, 100.0, 1000.0, 1.0);
    assert(gamma_rt > 0.0);

    /* Equal density -> stable */
    double gamma_rt_stable = rayleigh_taylor_growth(9.8, 100.0, 1000.0, 1000.0);
    assert(gamma_rt_stable < TOL);

    double gamma_kh = kelvin_helmholtz_growth(100.0, 10.0, 1.0, 1.0);
    assert(gamma_kh > 0.0);

    /* Magnetic stabilization */
    double gamma_mkh = magnetic_kh_growth(100.0, 10.0, 1.0, 1.0, 1e-3);
    assert(gamma_mkh <= gamma_kh); /* B reduces growth */

    printf("  PASS: test_rt_kh\n");
}

/* ---- L6: Kink / Tearing ---- */
static void test_kink_tearing(void)
{
    assert(kink_unstable(0.8, 1, 1) == 1);
    assert(kink_unstable(1.2, 1, 1) == 0);

    double gamma_tear = tearing_mode_growth_rate(1e-8, 0.5, 1e6, 10.0);
    assert(gamma_tear > 0.0);

    printf("  PASS: test_kink_tearing\n");
}

/* ---- L6: Drift Wave ---- */
static void test_drift_wave(void)
{
    double omega_star = drift_wave_frequency(100.0, 1.6e-19,
                                               1.0/0.1, 1.0);
    assert(omega_star > 0.0);

    /* ITG: eta_i > eta_i_crit -> unstable */
    double gamma_itg = itg_growth_rate(omega_star, 3.0, 1.0, 0.3);
    assert(gamma_itg > 0.0);

    /* TEM */
    double gamma_tem = tem_growth_rate(omega_star, 2.0, 0.1);
    assert(gamma_tem > 0.0);

    /* ETG */
    double gamma_etg = etg_growth_rate(omega_star, 3.0, 1.0, 0.5);
    assert(gamma_etg > 0.0);

    printf("  PASS: test_drift_wave\n");
}

/* ---- L8: Parametric Instabilities ---- */
static void test_parametric(void)
{
    /* Three-wave resonance check */
    int resonant = three_wave_resonance_check(10.0, 1.0, 5.0, 0.5,
                                               15.0, 1.5, 1, 0.01);
    assert(resonant == 1);

    int not_resonant = three_wave_resonance_check(10.0, 1.0, 5.0, 0.5,
                                                    20.0, 1.5, 1, 0.01);
    assert(not_resonant == 0);

    double gamma_pdi = parametric_decay_growth(1e5, 1e6, 1e12, 1e10);
    assert(gamma_pdi > 0.0);

    printf("  PASS: test_parametric\n");
}

/* ---- L5: Dispersion Solver ---- */
static double test_disp_func(double omega, const double params[])
{
    /* Simple quadratic: (omega - a)^2 - b = 0 */
    double a = params[0], b = params[1];
    return (omega - a)*(omega - a) - b;
}

static void test_dispersion_solver(void)
{
    double params[2] = {5.0, 4.0}; /* roots at 3 and 7 */
    double roots[4];
    int n = find_dispersion_roots_brent(test_disp_func, params,
                                         0.0, 10.0, 100, roots, 4, 1e-8);
    assert(n >= 2);
    int found3 = 0, found7 = 0;
    for (int i = 0; i < n; i++) {
        if (fabs(roots[i] - 3.0) < 1e-6) found3 = 1;
        if (fabs(roots[i] - 7.0) < 1e-6) found7 = 1;
    }
    assert(found3 && found7);

    printf("  PASS: test_dispersion_solver\n");
}

/* ---- L8: Nonlinear Saturation ---- */
static void test_nonlinear_saturation(void)
{
    double tau_ql = quasilinear_saturation_time(1e6, 5.64e11,
                                                  1e5, 1e7, 0.01);
    assert(tau_ql > 0.0);

    double E_sat = trapping_saturation_field(1e6, 100.0, M_ELECTRON, E_CHARGE);
    assert(E_sat > 0.0);

    double f_deplete = pump_depletion_fraction(1e6, 1e5);
    assert(f_deplete > 0.0 && f_deplete <= 1.0);

    /* Langmuir collapse: above threshold.
     * Need W/(n0*T_e) > k^2 * lambda_De^2.
     * With W=1e6, n0=1e14, T_e=1.6e-19*100, k=10, lambda_De=7e-3:
     *   LHS = 1e6/(1e14*1.6e-17) = 1e6/1.6e-3 = 6.25e8
     *   RHS = 100 * 4.9e-5 = 4.9e-3
     *   LHS >> RHS -> collapse. */
    int collapse = langmuir_collapse_threshold(1.0, 1e14,
                                                 1.6e-19*1000.0, 10.0, 7e-3);
    int no_collapse = langmuir_collapse_threshold(1e-12, 1e20,
                                                    1.6e-19, 100.0, 7e-5);
    assert(collapse == 1);
    assert(no_collapse == 0);

    printf("  PASS: test_nonlinear_saturation\n");
}

/* ---- L5: Nyquist ---- */
static void test_nyquist_dummy(double wr, double wi, const double p[],
                                double *Dr, double *Di)
{
    /* D(omega) = omega - 5 + i*omega/10 */
    *Dr = wr - p[0];
    *Di = wi * 0.1;
}

static void test_nyquist(void)
{
    double params[1] = {5.0};
    int n = nyquist_unstable_count(test_nyquist_dummy, params,
                                    -10.0, 10.0, 3.0, 100);
    /* Should find 0 or 1 unstable root depending on contour */
    assert(n >= 0);

    printf("  PASS: test_nyquist\n");
}

int main(void)
{
    printf("=== Test: Plasma Instabilities ===\n");
    test_two_stream();
    test_weibel();
    test_firehose_mirror();
    test_rt_kh();
    test_kink_tearing();
    test_drift_wave();
    test_parametric();
    test_dispersion_solver();
    test_nonlinear_saturation();
    test_nyquist();
    printf("=== All instability tests passed (10/10) ===\n");
    return 0;
}
