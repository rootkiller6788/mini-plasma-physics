/**
 * dispersion_solvers.c -- Numerical Dispersion Relation Solvers
 *
 * L5: Brent root finding, complex Newton iteration, continuation method,
 *     MHD stability eigenvalue solver, shooting method for ODEs,
 *     Nyquist stability analysis, Penrose criterion, CMA diagram grid
 *
 * References:
 *   Brent, "Algorithms for Minimization without Derivatives" (1973)
 *   Press et al., "Numerical Recipes" (2007), Ch. 9-11
 *   Stix, "Waves in Plasmas" (1992), Appendices
 *   Krall & Trivelpiece (1973), Appendix E
 */

#include "waves_instabilities.h"
#include "dispersion_solvers.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* ===============================================================
 * L5: Brent's Method for Real Root Finding
 * =============================================================== */

/**
 * Brent's root-finding algorithm (Brent, 1973).
 *
 * Finds a root of f(x) in [a, b] where f(a)*f(b) < 0.
 * Combines bisection (guaranteed convergence) with
 * inverse quadratic interpolation (fast convergence).
 * Superlinear convergence for well-behaved functions.
 *
 * Complexity: O(log(1/eps)) per root
 */
static double brent_root(DispersionFunc f, const double params[],
                         double a, double b, double eps)
{
    double fa = f(a, params), fb = f(b, params);
    if (fa * fb > 0.0) return a; /* No sign change, return left bound */

    double c = b, fc = fb;
    double d = b - a, e = d;

    for (int iter = 0; iter < 100; iter++) {
        if (fb * fc > 0.0) {
            c = a; fc = fa; d = b - a; e = d;
        }
        if (fabs(fc) < fabs(fb)) {
            a = b; b = c; c = a;
            fa = fb; fb = fc; fc = fa;
        }

        double tol = 2.0 * eps * fabs(b) + 0.5 * eps;
        double m = 0.5 * (c - b);

        if (fabs(m) <= tol || fabs(fb) < eps) return b;

        if (fabs(e) >= tol && fabs(fa) > fabs(fb)) {
            double s = fb / fa;
            double p, q;
            if (fabs(a - c) < 1e-30) {
                /* Linear interpolation */
                p = 2.0 * m * s;
                q = 1.0 - s;
            } else {
                /* Inverse quadratic interpolation */
                q = fa / fc;
                double r = fb / fc;
                p = s * (2.0*m*q*(q - r) - (b - a)*(r - 1.0));
                q = (q - 1.0)*(r - 1.0)*(s - 1.0);
            }
            if (p > 0.0) q = -q;
            else p = -p;
            s = e; e = d;
            if (2.0*p < 3.0*m*q - fabs(tol*q) && p < fabs(0.5*s*q)) {
                d = p / q;
            } else {
                d = m; e = m;
            }
        } else {
            d = m; e = m;
        }

        a = b; fa = fb;
        if (fabs(d) > tol) b += d;
        else if (m > 0.0) b += tol;
        else b -= tol;
        fb = f(b, params);
    }
    return b;
}

int find_dispersion_roots_brent(DispersionFunc D, const double params[],
                                 double w_min, double w_max,
                                 int n_scan, double roots[], int max_roots,
                                 double eps)
{
    if (!D || w_min >= w_max || n_scan < 2 || !roots || max_roots <= 0)
        return 0;

    double dw = (w_max - w_min) / n_scan;
    int n_found = 0;

    double f_prev = D(w_min, params);

    for (int i = 0; i < n_scan && n_found < max_roots; i++) {
        double w_i = w_min + i * dw;
        double w_next = w_i + dw;

        /* Evaluate at endpoints of subinterval */
        double f_next = D(w_next, params);

        /* Check for sign change */
        if (f_prev * f_next <= 0.0) {
            roots[n_found] = brent_root(D, params, w_i, w_next, eps);
            /* Check it is distinct from previously found roots */
            int distinct = 1;
            for (int j = 0; j < n_found; j++) {
                if (fabs(roots[n_found] - roots[j]) < 100.0 * eps)
                    { distinct = 0; break; }
            }
            if (distinct) n_found++;
        }

        f_prev = f_next;
    }

    return n_found;
}

/* ===============================================================
 * L5: Complex Root Finding (Newton in 2D)
 * =============================================================== */

int find_complex_root_newton(ComplexDispersionFunc D, const double params[],
                              double wr0, double g0,
                              double *wr_out, double *gamma_out,
                              int max_iter, double tol)
{
    if (!D || !wr_out || !gamma_out) return 2;

    double wr = wr0, gamma = g0;
    double delta = 1e-6; /* Finite difference step */

    for (int iter = 0; iter < max_iter; iter++) {
        double Dr, Di;
        D(wr, gamma, params, &Dr, &Di);

        if (sqrt(Dr*Dr + Di*Di) < tol) {
            *wr_out = wr; *gamma_out = gamma; return 0;
        }

        /* Finite-difference Jacobian */
        double Dr_wr, Di_wr, Dr_g, Di_g;
        D(wr + delta, gamma, params, &Dr_wr, &Di_wr);
        D(wr, gamma + delta, params, &Dr_g, &Di_g);

        double dDr_dwr = (Dr_wr - Dr) / delta;
        double dDi_dwr = (Di_wr - Di) / delta;
        double dDr_dg = (Dr_g - Dr) / delta;
        double dDi_dg = (Di_g - Di) / delta;

        double det = dDr_dwr * dDi_dg - dDi_dwr * dDr_dg;
        if (fabs(det) < 1e-30) return 2;

        /* Solve 2x2 linear system */
        double dwr = (-Dr * dDi_dg + Di * dDr_dg) / det;
        double dg = (-dDr_dwr * Di + dDi_dwr * Dr) / det;

        /* Damp large steps */
        double step_scale = sqrt(dwr*dwr + dg*dg);
        if (step_scale > fmax(fabs(wr), fabs(gamma))*0.5 + 1.0) {
            double s = (fmax(fabs(wr), fabs(gamma))*0.5 + 1.0) / step_scale;
            dwr *= s; dg *= s;
        }

        wr += dwr; gamma += dg;

        if (fabs(dwr) < tol && fabs(dg) < tol) {
            *wr_out = wr; *gamma_out = gamma; return 0;
        }
    }

    *wr_out = wr; *gamma_out = gamma; return 1; /* Not converged */
}

/* ===============================================================
 * L5: Dispersion Curve Continuation
 * =============================================================== */

int trace_dispersion_curve(DispersionFunc D, ParamBuilder build_params,
                            double k_min, double k_max, int n_pts,
                            double k_out[], double w_out[],
                            double w_guess)
{
    if (!D || !build_params || n_pts < 2 || !k_out || !w_out)
        return 0;

    double dk = (k_max - k_min) / (n_pts - 1);
    double params[16]; /* Generic parameter array */
    int n_success = 0;

    for (int i = 0; i < n_pts; i++) {
        double k_i = k_min + i * dk;
        k_out[i] = k_i;

        build_params(k_i, params);

        /* Use previous root as initial guess for next step */
        double roots[4];
        int n = find_dispersion_roots_brent(D, params,
                                             0.5*w_guess, 2.0*w_guess,
                                             50, roots, 1, 1e-8);
        if (n > 0) {
            w_out[i] = roots[0];
            w_guess = roots[0];
            n_success++;
        } else {
            w_out[i] = w_guess; /* Keep previous guess */
        }
    }

    return n_success;
}

/* ===============================================================
 * L5: MHD Stability Eigenvalue Solver (Inverse Iteration)
 * =============================================================== */

int mhd_stability_eigenvalue(const double rho[], const double F_diag[],
                              const double F_off[], int n,
                              double *w_sq, double v[], int max_iter)
{
    if (!rho || !F_diag || !F_off || !w_sq || !v || n < 2)
        return 1;

    /* Initialize guess */
    for (int i = 0; i < n; i++) v[i] = 1.0;
    double lambda = 0.0;

    for (int iter = 0; iter < max_iter; iter++) {
        /* Solve tridiagonal system: F*v_new = rho*v_old */
        double *a = (double*)malloc(n * sizeof(double));
        double *b = (double*)malloc(n * sizeof(double));
        double *c = (double*)malloc(n * sizeof(double));
        if (!a || !b || !c) { free(a); free(b); free(c); return 1; }

        for (int i = 0; i < n; i++) b[i] = F_diag[i];
        for (int i = 0; i < n-1; i++) {
            a[i+1] = F_off[i];
            c[i] = F_off[i];
        }
        a[0] = 0.0; c[n-1] = 0.0;

        /* Thomas algorithm (tridiagonal solver) */
        for (int i = 1; i < n; i++) {
            double w = a[i] / b[i-1];
            b[i] -= w * c[i-1];
            v[i] -= w * v[i-1];
        }
        v[n-1] /= b[n-1];
        for (int i = n-2; i >= 0; i--)
            v[i] = (v[i] - c[i]*v[i+1]) / b[i];

        free(a); free(b); free(c);

        /* Rayleigh quotient: lambda = (v^T F v) / (v^T rho v) */
        double num = 0.0, den = 0.0;
        for (int i = 0; i < n; i++) {
            double Fv = F_diag[i]*v[i]
                        + (i>0 ? F_off[i-1]*v[i-1] : 0.0)
                        + (i<n-1 ? F_off[i]*v[i+1] : 0.0);
            num += v[i] * Fv;
            den += rho[i] * v[i] * v[i];
        }
        double lambda_new = (den > 0.0) ? num / den : lambda;

        /* Normalize eigenvector */
        double norm = 0.0;
        for (int i = 0; i < n; i++) norm += v[i]*v[i];
        norm = sqrt(norm);
        if (norm > 0.0) for (int i = 0; i < n; i++) v[i] /= norm;

        if (fabs(lambda_new - lambda) < 1e-8 * (fabs(lambda_new) + 1e-30))
            { lambda = lambda_new; break; }
        lambda = lambda_new;
    }

    *w_sq = lambda;
    return 0;
}

/* ===============================================================
 * L5: Shooting Method for ODE Eigenvalues
 * =============================================================== */

int shooting_eigenvalue(ODE_RHS rhs, const double params[],
                         int n, const double r_grid[],
                         double *w_sq, double psi[])
{
    if (!rhs || !r_grid || !w_sq || !psi || n < 3) return 1;

    /* Simple shooting: integrate from r_min to r_max,
     * adjust omega^2 until BC at r_max is satisfied.
     * BC: psi(r_min) = r_min, dpsi/dr(r_min) = 1 (regular at origin)
     * BC: psi(r_max) = 0 (fixed boundary) */
    double w_lo = 0.0, w_hi = *w_sq * 10.0 + 1.0;

    for (int iter = 0; iter < 40; iter++) {
        double w_mid = 0.5 * (w_lo + w_hi);
        double dr = r_grid[1] - r_grid[0];

        psi[0] = r_grid[0];
        double dpsi = 1.0;

        /* Leapfrog integration */
        for (int i = 0; i < n-1; i++) {
            psi[i+1] = psi[i] + dpsi * dr;
            double d2psi = rhs(r_grid[i], psi[i], w_mid, params);
            dpsi += d2psi * dr;
        }

        if (psi[n-1] * psi[0] > 0.0) {
            /* No zero crossing; shift w_lo up */
            w_lo = w_mid;
        } else {
            w_hi = w_mid;
        }

        if (fabs(w_hi - w_lo) < 1e-8) break;
    }

    *w_sq = 0.5 * (w_lo + w_hi);
    return 0;
}

/* ===============================================================
 * L5: Nyquist Stability Analysis
 * =============================================================== */

int nyquist_unstable_count(ComplexDispersionFunc D, const double params[],
                            double w_min, double w_max, double R, int n_pts)
{
    if (!D || n_pts < 4) return -1;

    double dw = (w_max - w_min) / n_pts;
    double phase_total = 0.0;
    double Dr_prev, Di_prev;
    D(w_min, 0.0, params, &Dr_prev, &Di_prev);
    double phase_prev = atan2(Di_prev, Dr_prev);

    /* Real axis */
    for (int i = 1; i <= n_pts; i++) {
        double w = w_min + i * dw;
        double Dr, Di;
        D(w, 0.0, params, &Dr, &Di);
        double phase = atan2(Di, Dr);
        double dphase = phase - phase_prev;
        /* Unwrap phase */
        if (dphase > PLASMA_PI) dphase -= 2.0*PLASMA_PI;
        if (dphase < -PLASMA_PI) dphase += 2.0*PLASMA_PI;
        phase_total += dphase;
        phase_prev = phase;
        Dr_prev = Dr; Di_prev = Di;
    }

    /* Semicircle at large |omega| = R */
    for (int i = 1; i <= n_pts; i++) {
        double theta = PLASMA_PI * (0.5 - (double)i/n_pts);
        double wr = R * cos(theta);
        double wi = R * sin(theta);
        double Dr, Di;
        D(wr, wi, params, &Dr, &Di);
        double phase = atan2(Di, Dr);
        double dphase = phase - phase_prev;
        if (dphase > PLASMA_PI) dphase -= 2.0*PLASMA_PI;
        if (dphase < -PLASMA_PI) dphase += 2.0*PLASMA_PI;
        phase_total += dphase;
        phase_prev = phase;
    }

    /* N_unstable = phase_total / (2*pi) */
    return (int)round(phase_total / (2.0 * PLASMA_PI));
}

/* ===============================================================
 * L5: Penrose Criterion
 * =============================================================== */

int penrose_criterion(const double f[], const double v[], int n_v,
                       double omega_p)
{
    if (!f || !v || n_v < 4) return 0;
    (void)omega_p;

    /* Check if there exists v0 such that:
     * P integral (df/dv) / (v-v0) dv > 0
     * Approximate with trapezoidal rule and scan v0 candidates */
    double dv = v[1] - v[0];

    for (int j = 1; j < n_v - 1; j++) {
        double v0 = v[j];
        double integral = 0.0;

        for (int i = 1; i < n_v; i++) {
            double dfdv = (f[i] - f[i-1]) / dv;
            double vm = 0.5 * (v[i] + v[i-1]);
            if (fabs(vm - v0) > 1e-6 * fmax(fabs(vm), 1.0))
                integral += dfdv / (vm - v0) * dv;
        }

        if (integral > 0.01) return 1; /* Unstable */
    }

    return 0; /* Stable */
}

/* ===============================================================
 * L5: CMA Diagram Grid Generation
 * =============================================================== */

void cma_diagram_grid(double a_min, double a_max, int n_a,
                       double b_min, double b_max, int n_b,
                       int regions[], int n_prop[])
{
    if (!regions || !n_prop || n_a <= 0 || n_b <= 0) return;

    double da = (a_max - a_min) / n_a;
    double db = (b_max - b_min) / n_b;
    double w = 1.0; /* Normalized: omega is virtual */

    for (int i = 0; i < n_a; i++) {
        for (int j = 0; j < n_b; j++) {
            int idx = i * n_b + j;
            double alpha = a_min + i * da;
            double beta = b_min + j * db;

            /* Determine CMA region */
            CMARegion cr;
            /* Back-calculate omega_pe, omega_ce from alpha, beta */
            double omega_pe_calc = sqrt(alpha) * w;
            double omega_ce_calc = beta * w;
            cr = cma_classify(omega_pe_calc, omega_ce_calc, w);
            regions[idx] = cr.region;

            /* Count propagating modes using Stix tensor */
            double wp_arr[] = { omega_pe_calc };
            double wc_arr[] = { omega_ce_calc };
            StixDielectric diel;
            stix_dielectric_tensor(w, wp_arr, wc_arr, 1, &diel);
            double n2a, n2b;
            cold_plasma_n2(&diel, 0.0, &n2a, &n2b); /* Perp: count O/X */
            int count = 0;
            if (n2a > 0) count++;
            if (n2b > 0) count++;
            cold_plasma_n2(&diel, 1.0, &n2a, &n2b); /* Par: count R/L */
            if (n2a > 0) count++;
            if (n2b > 0 && n2b != n2a) count++;
            n_prop[idx] = (count > 2) ? 2 : count;
        }
    }
}

/* ===============================================================
 * L5: Special Frequencies Identification
 * =============================================================== */

int special_frequencies(double omega_pe, double omega_ce,
                         int n_sp, const double wp_ions[],
                         const double wc_ions[],
                         double special[], char labels[][32])
{
    int n = 0;
    /* Plasma frequency cutoff */
    special[n] = omega_pe;
    strcpy(labels[n++], "P cutoff (omega_pe)");

    /* R and L cutoffs */
    double wc_half = omega_ce * 0.5;
    special[n] = wc_half + sqrt(omega_pe*omega_pe + wc_half*wc_half);
    strcpy(labels[n++], "R cutoff");

    double wL = -wc_half + sqrt(omega_pe*omega_pe + wc_half*wc_half);
    if (wL > 0) {
        special[n] = wL;
        strcpy(labels[n++], "L cutoff");
    }

    /* Upper hybrid */
    special[n] = sqrt(omega_pe*omega_pe + omega_ce*omega_ce);
    strcpy(labels[n++], "Upper hybrid");

    /* Lower hybrid (uses ion contributions) */
    double wpi_tot = 0.0, wci_tot = 0.0;
    for (int s = 0; s < n_sp; s++) {
        wpi_tot += wp_ions[s] * wp_ions[s];
        wci_tot += wc_ions[s];
    }
    if (wpi_tot > 0 && omega_ce > 0) {
        double wlh2_inv = 1.0/(wci_tot*omega_ce) + 1.0/wpi_tot;
        special[n] = 1.0/sqrt(wlh2_inv);
        strcpy(labels[n++], "Lower hybrid");
    }

    return n;
}
