/**
 * beam_instabilities.c -- Beam-Plasma and Kinetic Instabilities
 *
 * L6: Two-stream (cold, warm, Buneman), bump-on-tail, Weibel,
 *     firehose, mirror, current-driven ion acoustic instabilities
 * L7: Application to laser-plasma, space physics, fusion
 *
 * References:
 *   Buneman, Phys. Rev. 115, 503 (1959)
 *   Weibel, Phys. Rev. Lett. 2, 83 (1959)
 *   Mikhailovskii, "Theory of Plasma Instabilities" (1974)
 *   Krall & Trivelpiece (1973), Ch. 9
 */

#include "waves_instabilities.h"
#include "plasma_instabilities.h"
#include <math.h>
#include <stdlib.h>

/* ===============================================================
 * L6: Two-Stream / Beam-Plasma Instabilities
 * =============================================================== */

double two_stream_growth_rate(double n_beam, double n_plasma, double omega_pe)
{
    if (n_plasma <= 0.0 || omega_pe <= 0.0) return 0.0;
    double ratio = n_beam / n_plasma;
    if (ratio <= 0.0) return 0.0;
    /* Maximum growth: gamma_max = (sqrt(3)/2^(4/3)) * ratio^(1/3) * omega_pe
     * sqrt(3)/2^(4/3) ~ 0.687 */
    return 0.687 * pow(ratio, 1.0/3.0) * omega_pe;
}

/**
 * Solve two-stream quartic dispersion analytically for cold beams.
 *
 * The quartic: omega^4 - 2*k*v_b*omega^3 + (k^2*v_b^2 - omega_pe^2 - omega_b^2)*omega^2
 *              + 2*k*v_b*omega_pe^2*omega + omega_pe^2*(omega_b^2 - k^2*v_b^2) = 0
 *
 * For each k, we solve for omega using the discriminant method for quartics.
 */
int two_stream_dispersion_roots(double k, double omega_pe, double omega_b,
                                 double v_b, ComplexOmega roots[4])
{
    if (!roots || omega_pe <= 0.0 || v_b <= 0.0) return 0;

    double kvb = k * v_b;
    double wp2 = omega_pe * omega_pe;
    double wb2 = omega_b * omega_b;

    /* Quartic coefficients: a4*omega^4 + a3*omega^3 + a2*omega^2 + a1*omega + a0 = 0 */
    double a4 = 1.0;
    double a3 = -2.0 * kvb;
    double a2 = kvb * kvb - wp2 - wb2;
    double a1 = 2.0 * kvb * wp2;
    double a0 = wp2 * (wb2 - kvb * kvb);

    /* Normalize */
    double b3 = a3 / a4, b2 = a2 / a4, b1 = a1 / a4, b0 = a0 / a4;

    /* Ferrari's method for quartic */
    double p = b2 - 3.0*b3*b3/8.0;
    double q = b1 - b3*b2/2.0 + b3*b3*b3/8.0;
    double r = b0 - b3*b1/4.0 + b3*b3*b2/16.0 - 3.0*b3*b3*b3*b3/256.0;

    if (fabs(q) < 1e-30 && fabs(r) < 1e-30) {
        /* Degenerate: (omega - b3/4)^4 = 0 */
        roots[0].omega_r = -b3/4.0; roots[0].gamma = 0.0;
        return 1;
    }

    /* Solve resolvent cubic: y^3 + 2p*y^2 + (p^2-4r)*y - q^2 = 0 */
    double c2 = 2.0 * p;
    double c1 = p * p - 4.0 * r;
    double c0 = -q * q;

    /* Cubic solver via trigonometric method (Casus irreducibilis check) */
    double Q = (c2*c2 - 3.0*c1) / 9.0;
    double R = (2.0*c2*c2*c2 - 9.0*c2*c1 + 27.0*c0) / 54.0;
    double R2_Q3 = R*R - Q*Q*Q;

    double y1;
    if (R2_Q3 <= 0.0) {
        /* Three real roots, take the largest */
        double theta = acos(R / sqrt(Q*Q*Q));
        double sqrtQ = sqrt(Q);
        y1 = -2.0*sqrtQ*cos(theta/3.0) - c2/3.0;
        double y2 = -2.0*sqrtQ*cos((theta+2.0*PLASMA_PI)/3.0) - c2/3.0;
        double y3 = -2.0*sqrtQ*cos((theta-2.0*PLASMA_PI)/3.0) - c2/3.0;
        if (y2 > y1) y1 = y2;
        if (y3 > y1) y1 = y3;
    } else {
        /* One real root */
        double A = -pow(fabs(R) + sqrt(R2_Q3), 1.0/3.0);
        if (R < 0) A = pow(fabs(R) + sqrt(R2_Q3), 1.0/3.0);
        y1 = A + (fabs(A) > 1e-30 ? Q/A : 0.0) - c2/3.0;
    }

    if (y1 < 1e-30) y1 = 1e-30;

    double sqrt_y1 = sqrt(y1);
    double R1 = sqrt(y1 - 4.0*r);
    double D1, D2;
    double q_over_sqrt_y1 = q / sqrt_y1;

    if (p - y1 - q_over_sqrt_y1 >= 0) D1 = sqrt(p - y1 - q_over_sqrt_y1);
    else D1 = 0.0;
    if (p - y1 + q_over_sqrt_y1 >= 0) D2 = sqrt(p - y1 + q_over_sqrt_y1);
    else D2 = 0.0;

    int n = 0;
    double shift = -b3/4.0;

    /* Four roots */
    double re1 = shift + 0.5*(sqrt_y1 + D1);
    double re2 = shift + 0.5*(sqrt_y1 - D1);
    double re3 = shift - 0.5*(sqrt_y1 - D2);
    double re4 = shift - 0.5*(sqrt_y1 + D2);

    /* Imaginary parts from R1 */
    double im1 = (R1 > 0) ? 0.5*R1 : 0.0;
    double im2 = (R1 > 0) ? -0.5*R1 : 0.0;

    if (fabs(im1) < 1e-30 && fabs(im2) < 1e-30) {
        /* All real roots */
        roots[n].omega_r = re1; roots[n++].gamma = 0.0;
        roots[n].omega_r = re2; roots[n++].gamma = 0.0;
        roots[n].omega_r = re3; roots[n++].gamma = 0.0;
        roots[n].omega_r = re4; roots[n++].gamma = 0.0;
    } else {
        /* Complex conjugate pairs */
        roots[n].omega_r = re1; roots[n++].gamma = im1;
        roots[n].omega_r = re2; roots[n++].gamma = im2;
        roots[n].omega_r = re3; roots[n++].gamma = im2;
        roots[n].omega_r = re4; roots[n++].gamma = im1;
    }

    return n;
}

double two_stream_growth_at_k(double k, double omega_pe, double omega_b,
                               double v_b)
{
    ComplexOmega roots[4];
    int n = two_stream_dispersion_roots(k, omega_pe, omega_b, v_b, roots);
    double gamma_max = 0.0;
    for (int i = 0; i < n; i++) {
        if (roots[i].gamma > gamma_max) gamma_max = roots[i].gamma;
    }
    return gamma_max;
}

double buneman_growth_rate(double omega_pe, double v_d, double v_th_e)
{
    if (omega_pe <= 0.0) return 0.0;
    /* Buneman: electron beam + stationary ions
     * gamma_max = (sqrt(3)/2) * (m_e/(2*m_i))^(1/3) * omega_pe
     * for cold beams (v_d >> v_th_e) */
    if (v_d > v_th_e) {
        double mass_ratio = M_ELECTRON / (2.0 * M_PROTON);
        return 0.866 * pow(mass_ratio, 1.0/3.0) * omega_pe;
    }
    /* Warm beam: growth is reduced */
    double reduction = exp(-v_th_e * v_th_e / (v_d * v_d));
    double mass_ratio = M_ELECTRON / (2.0 * M_PROTON);
    return 0.866 * pow(mass_ratio, 1.0/3.0) * omega_pe * reduction;
}

double ion_two_stream_growth(double omega_pi, double v_b, double v_th_i)
{
    if (omega_pi <= 0.0 || v_th_i <= 0.0) return 0.0;
    return omega_pi * (v_b / v_th_i) * 0.3;
}

double bump_on_tail_growth(double omega_pe, double n_b, double n_0,
                            double v_b, double delta_v)
{
    if (omega_pe <= 0.0 || n_0 <= 0.0 || delta_v <= 0.0) return 0.0;
    double ratio = n_b / n_0;
    if (ratio <= 0.0) return 0.0;
    double vb_dv = v_b / delta_v;
    return omega_pe * sqrt(PLASMA_PI/8.0) * ratio
           * vb_dv * vb_dv * exp(-vb_dv*vb_dv/2.0);
}

/* ===============================================================
 * L6: Weibel and Electromagnetic Instabilities
 * =============================================================== */

double weibel_growth_rate(double omega_pe, double v_th, double T_perp,
                           double T_par, double beta)
{
    if (omega_pe <= 0.0 || T_par <= 0.0 || v_th <= 0.0) return 0.0;
    double anisotropy = T_perp / T_par - 1.0;
    if (anisotropy <= 0.0) return 0.0;
    return omega_pe * (v_th / C_LIGHT) * sqrt(anisotropy * beta * 0.5);
}

double weibel_relativistic_growth(double omega_pe, double n_b, double n_0,
                                   double gamma_b)
{
    if (omega_pe <= 0.0 || n_0 <= 0.0 || gamma_b <= 0.0) return 0.0;
    return omega_pe * sqrt(n_b / (gamma_b * n_0));
}

int firehose_unstable(double p_parallel, double p_perp, double B)
{
    double b2_over_mu0 = B * B / MU0;
    return (p_parallel > p_perp + b2_over_mu0) ? 1 : 0;
}

double firehose_growth_rate_sq(double k_parallel, double p_parallel,
                                double p_perp, double B, double rho)
{
    if (rho <= 0.0) return 0.0;
    double delta = p_parallel - p_perp - B * B / MU0;
    return k_parallel * k_parallel * delta / rho;
}

int mirror_unstable(double p_parallel, double p_perp, double B)
{
    if (p_parallel <= 0.0) return 0;
    double b2_over_2mu0 = B * B / (2.0 * MU0);
    return (p_perp / p_parallel > 1.0 + b2_over_2mu0 / p_parallel) ? 1 : 0;
}

double mirror_growth_rate(double k_parallel, double v_th, double beta_perp,
                           double beta_par)
{
    if (v_th <= 0.0) return 0.0;
    double delta_beta = beta_perp - beta_par - 1.0;
    if (delta_beta <= 0.0) return 0.0;
    return fabs(k_parallel) * v_th * delta_beta / (2.0 * sqrt(PLASMA_PI));
}

double current_driven_ia_growth(double k, double c_s, double v_d,
                                 double m_e, double m_i)
{
    if (c_s <= 0.0 || m_i <= 0.0 || v_d <= c_s) return 0.0;
    return sqrt(PLASMA_PI/8.0) * k * c_s
           * (v_d/c_s - 1.0) * sqrt(m_e/m_i);
}

/* ===============================================================
 * L6: Loss Cone and Cyclotron Maser
 * =============================================================== */

double loss_cone_growth_rate(double omega_ce, double n_hot, double n_cold,
                              double anisotropy)
{
    if (omega_ce <= 0.0 || n_cold <= 0.0) return 0.0;
    return omega_ce * (n_hot / n_cold) * anisotropy;
}

double cyclotron_maser_growth(double omega_pe, double omega_ce,
                               double v_perp_over_c, double n_hot_over_n)
{
    if (omega_ce <= 0.0) return 0.0;
    double gamma_factor = 1.0 / sqrt(1.0 - v_perp_over_c * v_perp_over_c);
    return omega_pe * omega_pe * PLASMA_PI / omega_ce
           * n_hot_over_n / gamma_factor
           * (v_perp_over_c * v_perp_over_c)
           * 0.1;
}
