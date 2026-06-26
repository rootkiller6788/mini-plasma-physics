/**
 * plasma_waves.c -- Plasma Waves and Instabilities Implementation
 *
 * Every function implements one independent wave phenomenon.
 * References:
 *   Stix "Waves in Plasmas" (1992)
 *   Swanson "Plasma Waves" (2003)
 *   Kivelson & Russell S5
 *   MIT 22.611 Lectures 15-18
 *   Goldston & Rutherford S6
 *
 * Knowledge: L1-L6 plasma wave theory
 */
#include "../include/space_plasma.h"
#include "../include/plasma_waves.h"
#include <math.h>
#include <string.h>

/*========================================================================
 * L1: Stix Cold Plasma Dielectric Tensor
 *========================================================================*/

void stix_dielectric_tensor(double omega, double omega_pe, double omega_ce,
                            double omega_pi, double omega_ci, stix_tensor_t *eps) {
    /* Stix cold plasma dielectric tensor elements.
     * R = 1 - sum_s omega_ps^2 / (omega*(omega+omega_cs))
     * L = 1 - sum_s omega_ps^2 / (omega*(omega-omega_cs))
     * P = 1 - sum_s omega_ps^2 / omega^2
     * S = (R+L)/2, D = (R-L)/2
     * Reference: Stix (1992) "Waves in Plasmas" Ch.1-2 */
    if (!eps || omega <= 0.0) return;

    /* Electron contributions */
    double ome2 = omega_pe * omega_pe;
    double omi2 = omega_pi * omega_pi;
    double om2  = omega * omega;

    double Re = 1.0 - ome2 / (omega * (omega + omega_ce));
    double Le = 1.0 - ome2 / (omega * (omega - omega_ce));
    double Pe = 1.0 - ome2 / om2;

    /* Ion contributions */
    double Ri = (omega_pi > 0.0) ? -omi2 / (omega * (omega + omega_ci)) : 0.0;
    double Li = (omega_pi > 0.0) ? -omi2 / (omega * (omega - omega_ci)) : 0.0;
    double Pi = (omega_pi > 0.0) ? -omi2 / om2 : 0.0;

    eps->R = Re + Ri;
    eps->L = Le + Li;
    eps->P = Pe + Pi;
    eps->S = 0.5 * (eps->R + eps->L);
    eps->D = 0.5 * (eps->R - eps->L);
}

/*========================================================================
 * L1: Appleton-Hartree Dispersion Relation
 *========================================================================*/

double cold_plasma_dispersion_residual(double omega, double k, double theta,
                                       const stix_tensor_t *eps) {
    /* Cold plasma dispersion: A*n^4 - B*n^2 + C = 0, where n = ck/omega.
     * A = S*sin^2(theta) + P*cos^2(theta)
     * B = R*L*sin^2(theta) + P*S*(1+cos^2(theta))
     * C = P*R*L
     * Reference: Appleton (1932), Hartree (1931), Stix (1992) Eq. (1-48) */
    if (!eps || omega <= 0.0) return INFINITY;

    double n = SP_C * k / omega;  /* refractive index */
    double n2 = n * n;
    double n4 = n2 * n2;

    double st = sin(theta);
    double ct = cos(theta);
    double s2 = st * st;
    double c2 = ct * ct;
    double c21 = 1.0 + c2;

    double A = eps->S * s2 + eps->P * c2;
    double B = eps->R * eps->L * s2 + eps->P * eps->S * c21;
    double C = eps->P * eps->R * eps->L;

    double residual = A * n4 - B * n2 + C;

    /* Normalize by max coefficient for numerical stability */
    double norm = fabs(A*n4) + fabs(B*n2) + fabs(C);
    return (norm > 1e-30) ? residual / norm : residual;
}

int cold_plasma_refractive_index(double omega, double theta,
                                 const stix_tensor_t *eps,
                                 double *n_fast, double *n_slow) {
    (void)omega;  /* Used implicitly through eps tensor which depends on omega */
    /* Solve A*n^4 - B*n^2 + C = 0 for n^2.
     * Two branches: n_fast (smaller n^2) and n_slow (larger n^2).
     * Returns 0 on success, -1 if no propagating solution. */
    if (!eps || !n_fast || !n_slow) return -1;

    double st = sin(theta);
    double ct = cos(theta);
    double s2 = st * st;
    double c2 = ct * ct;
    double c21 = 1.0 + c2;

    double A = eps->S * s2 + eps->P * c2;
    double B = eps->R * eps->L * s2 + eps->P * eps->S * c21;
    double C = eps->P * eps->R * eps->L;

    /* Solve quadratic: A*u^2 - B*u + C = 0 where u = n^2 */
    double disc = B*B - 4.0*A*C;
    if (disc < 0.0 || fabs(A) < 1e-30) return -1;

    double sqrt_disc = sqrt(disc);
    double u1 = (B - sqrt_disc) / (2.0 * A);  /* n_fast^2 */
    double u2 = (B + sqrt_disc) / (2.0 * A);  /* n_slow^2 */

    /* Ensure physical ordering: n_fast has smaller n (faster phase speed) */
    if (u1 > u2) { double tmp = u1; u1 = u2; u2 = tmp; }

    *n_fast = (u1 > 0.0) ? sqrt(u1) : 0.0;
    *n_slow = (u2 > 0.0) ? sqrt(u2) : 0.0;

    return ((*n_fast > 0.0) || (*n_slow > 0.0)) ? 0 : -1;
}

/*========================================================================
 * L2: Specific Wave Mode Dispersion Relations
 *========================================================================*/

double langmuir_frequency(double k, double omega_pe, double v_the) {
    /* Bohm-Gross dispersion for Langmuir waves:
     * omega^2 = omega_pe^2 + 3*k^2*v_the^2
     * Thermal correction from electron pressure gradient.
     * In cold limit (k*v_the << omega_pe): omega = omega_pe.
     * Reference: Bohm & Gross (1949) Phys. Rev. 75:1851 */
    if (omega_pe <= 0.0) return 0.0;
    double k2 = k * k;
    double vth2 = v_the * v_the;
    return sqrt(omega_pe*omega_pe + 3.0 * k2 * vth2);
}

double ion_acoustic_frequency(double k, double c_s, double lambda_De) {
    /* Ion acoustic wave dispersion in long-wavelength limit:
     * omega = k*c_s / sqrt(1 + k^2*lambda_De^2)
     * For k*lambda_De << 1: omega ~ k*c_s (sound-like)
     * For k*lambda_De >> 1: omega ~ omega_pi (ion oscillations)
     * Reference: Tonks & Langmuir (1929) Phys. Rev. 33:195 */
    if (k <= 0.0 || c_s <= 0.0) return 0.0;
    double k2ld2 = k * k * lambda_De * lambda_De;
    return k * c_s / sqrt(1.0 + k2ld2);
}

double ion_acoustic_landau_damping(double omega, double c_s, double v_the) {
    /* Landau damping rate for ion acoustic waves (T_e >> T_i):
     * gamma/omega ~ -sqrt(pi/8) * (c_s/v_the)^3 * exp(-c_s^2/(2*v_the^2))
     * Strong damping when phase speed ~ thermal speed.
     * Reference: Landau (1946) J. Phys. USSR 10:25 */
    if (omega <= 0.0 || c_s <= 0.0 || v_the <= 0.0) return 0.0;
    double ratio = c_s / v_the;
    double ratio2 = ratio * ratio;
    double damping = -sqrt(M_PI / 8.0) * ratio2 * ratio * exp(-0.5 * ratio2);
    return damping * omega;
}

double alfven_wave_frequency(double k, double theta, double v_A) {
    /* Shear Alfven wave: omega = k_parallel * v_A = k*cos(theta)*v_A
     * Non-dispersive in ideal MHD. Transverse wave, incompressible.
     * Reference: Alfven (1942) Nature 150:405 */
    if (k <= 0.0 || v_A <= 0.0) return 0.0;
    return k * cos(theta) * v_A;
}

double fast_magnetosonic_frequency(double k, double v_A, double c_s, double theta) {
    /* Fast magnetosonic wave:
     * omega^2 = k^2 * v_f^2 where v_f^2 solves v^4 - (vA^2+cs^2)v^2 + vA^2*cs^2*cos^2(theta) = 0
     * Takes the + sign (fast mode). Isotropic in perpendicular limit.
     * Reference: Kivelson & Russell S5.2 */
    if (k <= 0.0) return 0.0;
    double vA2 = v_A * v_A;
    double cs2 = c_s * c_s;
    double sum_sq = vA2 + cs2;
    double disc = sum_sq*sum_sq - 4.0 * vA2 * cs2 * cos(theta)*cos(theta);
    if (disc < 0.0) disc = 0.0;
    double vf2 = 0.5 * (sum_sq + sqrt(disc));
    return k * sqrt(vf2);
}

double slow_magnetosonic_frequency(double k, double v_A, double c_s, double theta) {
    /* Slow magnetosonic wave: same quadratic as fast mode but with - sign.
     * omega^2 = k^2 * v_s^2. Propagates mainly along B in low-beta.
     * Reference: Kivelson & Russell S5.2 */
    if (k <= 0.0) return 0.0;
    double vA2 = v_A * v_A;
    double cs2 = c_s * c_s;
    double sum_sq = vA2 + cs2;
    double disc = sum_sq*sum_sq - 4.0 * vA2 * cs2 * cos(theta)*cos(theta);
    if (disc < 0.0) disc = 0.0;
    double vs2 = 0.5 * (sum_sq - sqrt(disc));
    if (vs2 < 0.0) vs2 = 0.0;
    return k * sqrt(vs2);
}

double whistler_frequency(double k, double omega_pe, double omega_ce, double theta) {
    /* Whistler mode (R-wave, electron cyclotron range):
     * omega = omega_ce * (k*c/omega_pe)^2 * cos(theta)
     * Valid for omega_ci << omega << omega_ce.
     * Group velocity proportional to sqrt(omega) -> dispersion in chorus/whistlers.
     * Reference: Storey (1953) Phil. Trans. Roy. Soc. A246:113
     *            Helliwell "Whistlers and Related Ionospheric Phenomena" (1965) */
    if (k <= 0.0 || omega_pe <= 0.0 || omega_ce <= 0.0) return 0.0;
    double kc_over_omega_pe = k * SP_C / omega_pe;
    return omega_ce * kc_over_omega_pe * kc_over_omega_pe * cos(theta);
}

double lower_hybrid_wave_frequency(double omega_ci, double omega_ce) {
    /* Lower hybrid frequency:
     * omega_LH = omega_ci * omega_ce / sqrt(omega_ci^2 + omega_ce^2)
     * Approximately sqrt(omega_ci * omega_ce) when omega_ci << omega_ce.
     * Important for lower hybrid current drive and wave heating.
     * Reference: Stix (1992) Ch.3 */
    if (omega_ci <= 0.0 || omega_ce <= 0.0) return 0.0;
    return omega_ci * omega_ce / sqrt(omega_ci*omega_ci + omega_ce*omega_ce);
}

double kinetic_alfven_frequency(double k_par, double k_perp, double v_A,
                                double rho_i, double rho_s) {
    /* Kinetic Alfven wave with finite Larmor radius corrections:
     * omega^2 = k_par^2 * v_A^2 * (1 + k_perp^2*rho_s^2) / (1 + k_perp^2*rho_i^2)
     * At k_perp*rho_i << 1: reduces to MHD Alfven wave.
     * At k_perp*rho_i >> 1: omega ~ k_par*v_A * (rho_s/rho_i)
     * Reference: Hasegawa (1976) J. Geophys. Res. 81:5083 */
    if (k_par <= 0.0 || v_A <= 0.0) return 0.0;
    double kp2 = k_perp * k_perp;
    double num = 1.0 + kp2 * rho_s * rho_s;
    double den = 1.0 + kp2 * rho_i * rho_i;
    return k_par * v_A * sqrt(num / den);
}

/*========================================================================
 * L3: Wave Polarization
 *========================================================================*/

int wave_polarization(double omega, double k, double theta,
                      const stix_tensor_t *eps, wave_polarization_t *pol) {
    /* Compute wave polarization from cold plasma theory.
     * The electric field eigenvector E satisfies D.E = 0 where
     * D is the dispersion tensor. For given (omega,k,theta),
     * we solve for E_x/E_y from the cold plasma wave equation.
     *
     * Polarization ratio i*E_x/E_y determines wave character:
     *   >0: R-wave (whistler/helicon)
     *   <0: L-wave (ion cyclotron)
     *
     * Reference: Stix (1992) Ch.1, Goertz & Strangeway (1995) */
    if (!eps || !pol) return -1;
    memset(pol, 0, sizeof(*pol));

    double n = SP_C * k / omega;
    double n2 = n * n;
    double st = sin(theta);
    double ct = cos(theta);
    double s2 = st * st;
    double c2 = ct * ct;

    /* Dispersion matrix elements D_ij */
    double Dxx = eps->S - n2 * c2;
    double Dxy = -1.0 * eps->D;
    double Dyy = eps->S - n2;
    (void)n2; (void)s2; (void)st; (void)ct;  /* used in full matrix, simplified here */

    /* Solve D_{ij} E_j = 0 for E_x/E_y using the first two equations */
    double det_minor = Dxx*Dyy - Dxy*Dxy;
    if (fabs(det_minor) < 1e-30) return -1;

    double Ex_Ey = -Dxy / Dxx;  /* ratio from Dxx*Ex + Dxy*Ey = 0 (neglecting Dxz*Ez) */

    /* Polarization characterization */
    pol->pol_ratio   = fabs(Ex_Ey);
    pol->helicity    = (Ex_Ey > 0.0) ? 1.0 : -1.0;  /* sign of i*Ex/Ey */
    pol->ellipticity = 1.0 / pol->pol_ratio;  /* b/a ~ Ey/Ex for near-circular */

    /* Electric field eigenvector (normalized) */
    double norm = sqrt(1.0 + Ex_Ey*Ex_Ey);
    if (norm > 1e-30) {
        pol->E_vec[0] = Ex_Ey / norm;
        pol->E_vec[1] = 1.0 / norm;
        pol->E_vec[2] = 0.0;
    }

    return 0;
}

/*========================================================================
 * L3: Resonance and Cutoff Conditions
 *========================================================================*/

void wave_cutoff_frequencies(double omega_pe, double omega_ce, double cutoffs[3]) {
    /* Cutoff frequencies where refractive index n -> 0.
     * P = 0: omega = omega_pe (plasma cutoff, O-mode)
     * R = 0: omega = (omega_ce + sqrt(omega_ce^2 + 4*omega_pe^2))/2 (R-cutoff)
     * L = 0: omega = (-omega_ce + sqrt(omega_ce^2 + 4*omega_pe^2))/2 (L-cutoff)
     * Reference: Stix (1992) S1-4 */
    if (!cutoffs) return;
    double omc2 = omega_ce * omega_ce;
    double omp2 = omega_pe * omega_pe;

    cutoffs[0] = omega_pe;  /* P = 0 */
    cutoffs[1] = 0.5 * (omega_ce + sqrt(omc2 + 4.0*omp2));  /* R = 0 */
    cutoffs[2] = 0.5 * (-omega_ce + sqrt(omc2 + 4.0*omp2)); /* L = 0 */
    if (cutoffs[2] < 0.0) cutoffs[2] = 0.0;
}

void wave_resonance_frequencies(double omega_pe, double omega_ce,
                                double omega_pi, double omega_ci,
                                double theta, double res[3]) {
    /* Resonance frequencies where n -> infinity.
     * R->infinity: omega = omega_ce (electron cyclotron resonance)
     * S = 0: omega = omega_LH (lower hybrid resonance)
     * cos^2(theta)/R + sin^2(theta)/S -> infinity: omega ~ omega_UH (upper hybrid)
     *
     * Reference: Stix (1992) S1-4, Swanson (2003) S2.3 */
    (void)theta;  /* theta used for angle-dependent resonances in full treatment */
    if (!res) return;

    res[0] = omega_ce;  /* electron cyclotron resonance */

    /* Lower hybrid resonance: simplifies to omega_LH formula */
    double ome2 = omega_pe * omega_pe;
    double omc2_e = omega_ce * omega_ce;

    /* S(omega) = 0 approximation for lower hybrid:
     * Formula from Stix: omega_LH^-2 = omega_ci^-2 + omega_pi^-2 */
    if (omega_ci > 0.0 && omega_pi > 0.0) {
        res[1] = 1.0 / sqrt(1.0/(omega_ci*omega_ci) + 1.0/(omega_pi*omega_pi));
    } else {
        res[1] = 0.0;
    }

    /* Upper hybrid resonance (for perpendicular propagation theta=pi/2):
     * omega_UH^2 = omega_pe^2 + omega_ce^2 */
    res[2] = sqrt(ome2 + omc2_e);
}

/*========================================================================
 * L5: CMA Diagram Classifier
 *========================================================================*/

int cma_propagation_bands(double alpha, double beta) {
    /* CMA (Clemmow-Mullaly-Allis) diagram propagation bands.
     * alpha = omega_pe^2/omega^2, beta = omega_ce/omega
     *
     * Returns bitmask of propagating wave modes at given parameters:
     *   bit 0: R-wave (whistler/helicon)
     *   bit 1: L-wave
     *   bit 2: X-mode (extraordinary, perpendicular)
     *   bit 3: O-mode (ordinary, perpendicular)
     *
     * Bands determined by cutoff/resonance topology in (alpha, beta) plane.
     * Reference: Clemmow & Mullaly (1955), Allis (1959)
     *            Stix (1992) Ch.2 Fig. 2-1 */
    int bands = 0;

    /* R-wave: propagates below R=0 cutoff, i.e., beta > 0, alpha/(1-beta) < 1 */
    if (beta > 0.0 && alpha / (1.0 - beta) < 1.0) bands |= 1;

    /* L-wave: propagates below L=0 cutoff, i.e., alpha/(1+beta) < 1 */
    if (alpha / (1.0 + beta) < 1.0) bands |= 2;

    /* X-mode: propagates between upper hybrid and R/L cutoffs */
    /* Simplified X-mode band check */
    if (beta > 0.0) {
        double a_x_low  = 1.0 - beta;   /* lower cutoff */
        double a_x_high = 1.0 - beta*beta;  /* upper hybrid res */
        if (alpha < a_x_high && alpha > a_x_low) bands |= 4;
    }

    /* O-mode: propagates for alpha < 1 (above plasma frequency) */
    if (alpha < 1.0) bands |= 8;

    return bands;
}

/*========================================================================
 * L5: Group Velocity
 *========================================================================*/

void group_velocity(double (*omega_fn)(double, double), double k, double theta,
                    double *v_g_par, double *v_g_perp, double dk, double dtheta) {
    /* Group velocity v_g = domega/dk via central finite differences.
     * v_g_par = domega/dk (parallel direction)
     * v_g_perp = (1/k) * domega/dtheta (perpendicular direction)
     * Reference: Stix (1992) S1-5 */
    if (!omega_fn || !v_g_par || !v_g_perp) return;
    if (dk <= 0.0) dk = k * 1e-6 + 1e-10;
    if (dtheta <= 0.0) dtheta = 1e-6;

    double omega_plus_k  = omega_fn(k + dk, theta);
    double omega_minus_k = omega_fn(k - dk, theta);
    *v_g_par = (omega_plus_k - omega_minus_k) / (2.0 * dk);

    double omega_plus_th  = omega_fn(k, theta + dtheta);
    double omega_minus_th = omega_fn(k, theta - dtheta);
    *v_g_perp = (k > 1e-30) ? (omega_plus_th - omega_minus_th) / (2.0 * k * dtheta) : 0.0;
}

/*========================================================================
 * L5: Newton-Raphson Dispersion Solver
 *========================================================================*/

double solve_dispersion_omega(double (*D)(double,double,double),
                              double k, double theta, double omega0,
                              double tol, int max_iter) {
    /* Newton-Raphson root finder for D(omega,k,theta) = 0.
     * omega_{n+1} = omega_n - D(omega_n) / D'(omega_n)
     * where D' is approximated by central finite difference.
     * Reference: Press et al. "Numerical Recipes" S9.4 */
    if (!D || max_iter <= 0) return -1.0;

    double omega = omega0;
    double domega = omega * 1e-6 + 1e-10;

    for (int iter = 0; iter < max_iter; iter++) {
        double D0 = D(omega, k, theta);
        if (fabs(D0) < tol) return omega;

        /* Central difference for derivative */
        double D_plus  = D(omega + domega, k, theta);
        double D_minus = D(omega - domega, k, theta);
        double dD_domega = (D_plus - D_minus) / (2.0 * domega);

        if (fabs(dD_domega) < 1e-30) return -1.0;

        double delta = -D0 / dD_domega;

        /* Under-relaxation for robustness */
        if (fabs(delta) > 0.5 * omega) delta = 0.5 * omega * ((delta > 0) ? 1.0 : -1.0);

        omega += delta;
        if (omega <= 0.0) omega = omega0 * 0.1;
        if (omega < 1e-30) return -1.0;
    }

    return -1.0;  /* did not converge */
}

/*========================================================================
 * L4: Warm Plasma Dielectric Function (Electrostatic)
 *========================================================================*/

double warm_plasma_dielectric_function(double omega, double k,
                                       double omega_p, double v_th) {
    /* Warm plasma dielectric function (electrostatic, 1D):
     * epsilon(omega,k) = 1 - (omega_p^2/k^2*v_th^2) * Z'(omega/(sqrt(2)*k*v_th))
     *
     * Using the plasma dispersion function Z (Fried-Conte):
     * For weakly damped waves (phase speed >> v_th):
     *   epsilon ~ 1 - omega_p^2/omega^2 * (1 + 3*k^2*v_th^2/omega^2)
     *
     * This gives the Bohm-Gross dispersion when epsilon = 0.
     *
     * Reference: Fried & Conte "The Plasma Dispersion Function" (1961)
     *            Swanson (2003) S3.3 */
    if (omega <= 0.0 || k <= 0.0) return -INFINITY;

    double k2 = k * k;
    double vth2 = v_th * v_th;
    double om2 = omega * omega;

    /* Cold plasma term + thermal correction */
    double cold_term = 1.0 - (omega_p * omega_p) / om2;

    /* Thermal correction from expansion of Z' for large argument:
     * Z'(zeta) ~ -2 - 2*zeta*Z(zeta)
     * For |zeta| >> 1: Z(zeta) ~ -1/zeta - 1/(2*zeta^3) - 3/(4*zeta^5)
     * → epsilon ~ 1 - (omega_p/k*v_th)^2 * (2*zeta^2 + 1) = 1 - omega_p^2/omega^2 - 3*k^2*v_th^2*omega_p^2/omega^4
     */
    double zeta = omega / (sqrt(2.0) * k * v_th);

    if (fabs(zeta) > 3.0) {
        /* Large argument expansion: epsilon ~ 1 - omega_p^2/omega^2 * (1 + 3/(2*zeta^2)) */
        double corr = 3.0 / (2.0 * zeta * zeta);
        return cold_term - (omega_p*omega_p/om2) * corr;
    } else {
        /* Moderate argument: use full expression with approximate Z function */
        /* Z(zeta) approx using continued fraction for real zeta */
        double Z_re = -1.0/zeta - 1.0/(2.0*pow(zeta,3)) - 3.0/(4.0*pow(zeta,5));
        double Zp_re = -2.0 * (1.0 + zeta * Z_re);  /* Z'(zeta) */
        double k2vth2 = k2 * vth2;
        return 1.0 - (omega_p*omega_p / k2vth2) * Zp_re;
    }
}
