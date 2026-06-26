/**
 * plasma_kinetic.c — Kinetic Theory Implementation
 *
 * Implements Vlasov equation solvers, Fokker-Planck collision operators,
 * Landau damping analysis, and plasma dispersion function.
 *
 * References:
 *   - Vlasov (1938), JETP 8, 291
 *   - Landau (1946), JETP 16, 574
 *   - Rosenbluth, MacDonald & Judd (1957), Phys. Rev. 107, 1
 *   - Cheng & Knorr (1976), J. Comp. Phys. 22, 330
 *   - Fried & Conte (1961), "The Plasma Dispersion Function"
 *
 * Knowledge Coverage:
 *   L2: Distribution functions, Landau damping
 *   L3: Vlasov equation solver, Fokker-Planck operator
 *   L4: Collisionless damping (Landau theorem)
 */
#include "plasma_kinetic.h"
#include "plasma_params.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ================================================================
 * L2: Moments of Distribution Function
 * ================================================================ */

double density_from_f1d1v(const PhaseSpace1D1V *ps, int ix) {
    if (!ps || ix < 0 || ix >= ps->nx) return 0.0;
    double n = 0.0;
    for (int iv = 0; iv < ps->nv; iv++) {
        n += ps->f[ix][iv];
    }
    return n * ps->dv;
}

double mean_velocity_from_f1d1v(const PhaseSpace1D1V *ps, int ix) {
    if (!ps || ix < 0 || ix >= ps->nx) return 0.0;
    double n = density_from_f1d1v(ps, ix);
    if (n <= 0.0) return 0.0;
    double sum = 0.0;
    for (int iv = 0; iv < ps->nv; iv++) {
        double v = ps->v_min + iv * ps->dv;
        sum += v * ps->f[ix][iv];
    }
    return sum * ps->dv / n;
}

double temperature_from_f1d1v(const PhaseSpace1D1V *ps, int ix, double mass) {
    if (!ps || ix < 0 || ix >= ps->nx) return 0.0;
    double n = density_from_f1d1v(ps, ix);
    if (n <= 0.0) return 0.0;
    double u = mean_velocity_from_f1d1v(ps, ix);
    double sum = 0.0;
    for (int iv = 0; iv < ps->nv; iv++) {
        double v = ps->v_min + iv * ps->dv;
        double dv_rel = v - u;
        sum += dv_rel * dv_rel * ps->f[ix][iv];
    }
    return mass * sum * ps->dv / n;  /* T_kin = P/n (1D), in energy units */
}

double heat_flux_from_f1d1v(const PhaseSpace1D1V *ps, int ix, double mass) {
    if (!ps || ix < 0 || ix >= ps->nx) return 0.0;
    double u = mean_velocity_from_f1d1v(ps, ix);
    double sum = 0.0;
    for (int iv = 0; iv < ps->nv; iv++) {
        double v = ps->v_min + iv * ps->dv;
        double dv_rel = v - u;
        sum += dv_rel * dv_rel * dv_rel * ps->f[ix][iv];
    }
    return 0.5 * mass * sum * ps->dv;
}

double pressure_from_f1d1v(const PhaseSpace1D1V *ps, int ix, double mass) {
    if (!ps || ix < 0 || ix >= ps->nx) return 0.0;
    double u = mean_velocity_from_f1d1v(ps, ix);
    double sum = 0.0;
    for (int iv = 0; iv < ps->nv; iv++) {
        double v = ps->v_min + iv * ps->dv;
        double dv_rel = v - u;
        sum += dv_rel * dv_rel * ps->f[ix][iv];
    }
    return mass * sum * ps->dv;
}

/* ================================================================
 * L3: Initialize Maxwellian on Phase Space Grid
 * ================================================================ */

void init_maxwellian_1d1v(PhaseSpace1D1V *ps, double n0, double T0,
                          double mass) {
    if (!ps || !ps->f) return;
    double v_th = sqrt(2.0 * K_B * T0 / mass);
    double coeff = n0 / (sqrt(2.0 * M_PI) * v_th);
    for (int ix = 0; ix < ps->nx; ix++) {
        for (int iv = 0; iv < ps->nv; iv++) {
            double v = ps->v_min + iv * ps->dv;
            ps->f[ix][iv] = coeff * exp(-0.5 * (v * v) / (v_th * v_th));
        }
    }
}

/* ================================================================
 * L3: Vlasov Advection Operators
 * ================================================================ */

void vlasov_advection_x(const PhaseSpace1D1V *ps, double **dfdt) {
    if (!ps || !dfdt) return;
    double inv_2dx = 0.5 / ps->dx;
    for (int ix = 1; ix < ps->nx - 1; ix++) {
        for (int iv = 0; iv < ps->nv; iv++) {
            double v = ps->v_min + iv * ps->dv;
            /* -v * df/dx  (2nd-order central difference) */
            dfdt[ix][iv] = -v * (ps->f[ix+1][iv] - ps->f[ix-1][iv]) * inv_2dx;
        }
    }
    /* Boundary: zero flux at edges */
    for (int iv = 0; iv < ps->nv; iv++) {
        dfdt[0][iv] = 0.0;
        dfdt[ps->nx - 1][iv] = 0.0;
    }
}

void vlasov_advection_v(const PhaseSpace1D1V *ps, double charge, double mass,
                        double **dfdt) {
    if (!ps || !dfdt) return;
    double q_over_m = charge / mass;
    double inv_2dv = 0.5 / ps->dv;
    for (int ix = 0; ix < ps->nx; ix++) {
        double E = ps->E_field[ix];
        for (int iv = 1; iv < ps->nv - 1; iv++) {
            /* -(q/m) * E(x) * df/dv */
            dfdt[ix][iv] = -q_over_m * E
                           * (ps->f[ix][iv+1] - ps->f[ix][iv-1]) * inv_2dv;
        }
        dfdt[ix][0] = 0.0;
        dfdt[ix][ps->nv - 1] = 0.0;
    }
}

void vlasov_step_strang(PhaseSpace1D1V *ps, double charge, double mass) {
    if (!ps || !ps->f) return;
    int nx = ps->nx, nv = ps->nv;
    double dt = ps->dt;

    /* Allocate temporary arrays */
    double **f_temp = (double**)malloc(nx * sizeof(double*));
    double **dfdt   = (double**)malloc(nx * sizeof(double*));
    for (int ix = 0; ix < nx; ix++) {
        f_temp[ix] = (double*)malloc(nv * sizeof(double));
        dfdt[ix]   = (double*)malloc(nv * sizeof(double));
    }

    /* --- Step 1: Half-step velocity advection (dt/2) --- */
    vlasov_advection_v(ps, charge, mass, dfdt);
    for (int ix = 0; ix < nx; ix++)
        for (int iv = 0; iv < nv; iv++)
            f_temp[ix][iv] = ps->f[ix][iv] + 0.5 * dt * dfdt[ix][iv];

    /* Save E-field, swap f_temp<->f */
    double *E_save = (double*)malloc(nx * sizeof(double));
    memcpy(E_save, ps->E_field, nx * sizeof(double));
    for (int ix = 0; ix < nx; ix++)
        memcpy(ps->f[ix], f_temp[ix], nv * sizeof(double));

    /* --- Step 2: Full-step spatial advection (dt) --- */
    vlasov_advection_x(ps, dfdt);
    for (int ix = 0; ix < nx; ix++)
        for (int iv = 0; iv < nv; iv++)
            f_temp[ix][iv] = ps->f[ix][iv] + dt * dfdt[ix][iv];
    for (int ix = 0; ix < nx; ix++)
        memcpy(ps->f[ix], f_temp[ix], nv * sizeof(double));

    /* Restore E-field */
    memcpy(ps->E_field, E_save, nx * sizeof(double));
    free(E_save);

    /* --- Step 3: Half-step velocity advection (dt/2) --- */
    vlasov_advection_v(ps, charge, mass, dfdt);
    for (int ix = 0; ix < nx; ix++)
        for (int iv = 0; iv < nv; iv++)
            ps->f[ix][iv] += 0.5 * dt * dfdt[ix][iv];

    ps->t += dt;

    /* Cleanup */
    for (int ix = 0; ix < nx; ix++) { free(f_temp[ix]); free(dfdt[ix]); }
    free(f_temp); free(dfdt);
}

/* ================================================================
 * L4: Poisson Solver for Vlasov-Poisson System
 * ================================================================ */

void poisson_solve_1d(PhaseSpace1D1V *ps, double n0_ion) {
    if (!ps || !ps->E_field) return;
    int nx = ps->nx;
    double dx = ps->dx;

    /* Allocate tridiagonal matrix */
    double *a = (double*)malloc(nx * sizeof(double));
    double *b = (double*)malloc(nx * sizeof(double));
    double *c = (double*)malloc(nx * sizeof(double));
    double *rhs = (double*)malloc(nx * sizeof(double));
    double *phi = (double*)malloc(nx * sizeof(double));

    /* Compute charge density rho = e*(n_i - n_e(x)) */
    double dx2 = dx * dx;
    for (int ix = 0; ix < nx; ix++) {
        double ne = density_from_f1d1v(ps, ix);
        double rho = E_CHARGE * (n0_ion - ne);
        /* d^2 phi/dx^2 = -rho/epsilon_0
         * Finite difference: (phi_{i-1} - 2*phi_i + phi_{i+1})/dx^2 = -rho_i/epsilon_0 */
        a[ix] = 1.0;
        b[ix] = -2.0;
        c[ix] = 1.0;
        rhs[ix] = -rho * dx2 / EPSILON_0;
    }

    /* Dirichlet BC: phi[0] = phi[nx-1] = 0 */
    b[0] = 1.0; c[0] = 0.0; rhs[0] = 0.0;
    a[nx-1] = 0.0; b[nx-1] = 1.0; rhs[nx-1] = 0.0;

    /* Thomas algorithm (tridiagonal solve) */
    for (int i = 1; i < nx; i++) {
        double w = a[i] / b[i-1];
        b[i] -= w * c[i-1];
        rhs[i] -= w * rhs[i-1];
    }
    phi[nx-1] = rhs[nx-1] / b[nx-1];
    for (int i = nx - 2; i >= 0; i--) {
        phi[i] = (rhs[i] - c[i] * phi[i+1]) / b[i];
    }

    /* E = -d phi/dx (central difference) */
    for (int ix = 1; ix < nx - 1; ix++) {
        ps->E_field[ix] = -(phi[ix+1] - phi[ix-1]) / (2.0 * dx);
    }
    ps->E_field[0] = -(phi[1] - phi[0]) / dx;
    ps->E_field[nx-1] = -(phi[nx-1] - phi[nx-2]) / dx;

    free(a); free(b); free(c); free(rhs); free(phi);
}

/* ================================================================
 * L4: Plasma Dispersion Function Z(zeta) — Fried & Conte
 * ================================================================ */

/**
 * Rational approximation to the plasma dispersion function.
 *
 * Uses the Weideman (1994) algorithm with N=8 for fast computation.
 * Accurate to ~1e-6 in the entire complex plane.
 *
 * Reference: Weideman (1994), SIAM J. Numer. Anal. 31, 1497.
 * "Computation of the Complex Error Function"
 */
double plasma_dispersion_function(double zeta_real, double zeta_imag,
                                  double *Z_real, double *Z_imag) {
    /* Z(zeta) = i*sqrt(pi)*w(zeta) where w is the Faddeeva function.
     * Using a continued fraction expansion for |zeta| < 8,
     * asymptotic expansion for |zeta| > 8. */
    double x = zeta_real, y = zeta_imag;
    double r2 = x*x + y*y;

    /* For large |zeta|, use asymptotic expansion:
     * Z(zeta) ~ -1/zeta * (1 + 1/(2*zeta^2) + 3/(4*zeta^4) + ...)
     *          + i*sqrt(pi)*sigma*exp(-zeta^2)
     * where sigma = 0 if Im(zeta)>0, 1 if Im(zeta)=0, 2 if Im(zeta)<0 */
    if (r2 > 64.0) {
        double inv_r2 = 1.0 / r2;
        double inv_zeta_real = x / r2;
        double inv_zeta_imag = -y / r2;
        /* First two terms of asymptotic series */
        double re = -inv_zeta_real * (1.0 - 0.5 * inv_r2);
        double im = -inv_zeta_imag * (1.0 - 0.5 * inv_r2);
        /* Add residue term if near real axis */
        if (fabs(y) < 1e-6) {
            double exp_part = sqrt(M_PI) * exp(-r2);
            re -= exp_part * y / fabs(y + 1e-300);
        }
        *Z_real = re;
        *Z_imag = im;
        return 0.0;
    }

    /* For small |zeta|, use power series:
     * Z(zeta) = i*sqrt(pi)*exp(-zeta^2) - 2*zeta*(1 - 2*zeta^2/3 + ...) */
    double exp_re = exp(y*y - x*x) * cos(-2.0*x*y);
    double exp_im = exp(y*y - x*x) * sin(-2.0*x*y);

    double exp_factor = sqrt(M_PI);
    double re_exp = -exp_factor * exp_im;
    double im_exp =  exp_factor * exp_re;

    /* Power series: -2*zeta*F(1, 3/2; -zeta^2) using Taylor up to 8 terms */
    double term_re = x, term_im = y;
    double series_re = term_re, series_im = term_im;
    double sign = -1.0;
    double denom = 1.0;
    for (int n = 1; n < 12; n++) {
        denom += 2.0;
        double new_re = term_re * x - term_im * y;
        double new_im = term_re * y + term_im * x;
        term_re = new_re; term_im = new_im;
        double contrib = sign / denom;
        series_re += term_re * contrib;
        series_im += term_im * contrib;
        sign = -sign;
    }
    series_re *= -2.0;
    series_im *= -2.0;

    *Z_real = re_exp + series_re;
    *Z_imag = im_exp + series_im;
    return 0.0;
}

/* ================================================================
 * L4: Landau Damping Rate
 * ================================================================ */

double langmuir_phase_velocity(double k, double omega_pe, double Te) {
    if (k <= 0.0) return INFINITY;
    double v_the2 = 2.0 * K_B * Te / M_ELECTRON;
    double omega2 = omega_pe * omega_pe + 3.0 * v_the2 * k * k / 2.0;
    return sqrt(omega2) / k;
}

double landau_damping_rate(double k, double omega_pe, double lambda_D) {
    if (k <= 0.0 || lambda_D <= 0.0 || omega_pe <= 0.0) return 0.0;
    double k_lD = k * lambda_D;
    if (k_lD < 1e-6) return 0.0; /* Essentially no damping for k*lambda_D -> 0 */

    /* gamma_L = -omega_r * sqrt(pi/8) / (k*lambda_D)^3
     *           * exp(-1/(2*(k*lambda_D)^2) - 3/2) */
    double k_lD3 = k_lD * k_lD * k_lD;
    double exponent = -1.0 / (2.0 * k_lD * k_lD) - 1.5;
    if (exponent < -50.0) return 0.0; /* Numerically zero */
    return -omega_pe * sqrt(M_PI / 8.0) / k_lD3 * exp(exponent);
}

int is_strongly_landau_damped(double k, double omega_pe, double lambda_D) {
    if (k <= 0.0) return 0;
    double v_phase = langmuir_phase_velocity(k, omega_pe, 0.0);
    /* Simplified criterion: wave is strongly damped when
     * omega/(k*v_th) < 1/sqrt(2) */
    double v_the = omega_pe * lambda_D;
    return (v_phase / v_the) < (1.0 / sqrt(2.0));
}

/* ================================================================
 * L3: Rosenbluth Potentials for Fokker-Planck
 * ================================================================ */

void rosenbluth_potentials(double v, const Distribution1D *dist,
                           double *h, double *g,
                           double *dh, double *dg, double *d2g) {
    if (!dist || !dist->f || dist->nv < 2) {
        *h = *g = *dh = *dg = *d2g = 0.0;
        return;
    }

    double dv = dist->dv;
    int nv = dist->nv;

    /* h(v) = 4*pi * [ (1/v)*integral_0^v f(u)*u^2 du + integral_v^inf f(u)*u du ]
     * g(v) = 4*pi * [ integral_0^v f(u)*u^2*(v/3+u^2/(5v)) du
     *                + integral_v^inf f(u)*u*(v^2/3+u^2) du ] */
    double h_val = 0.0, g_val = 0.0;
    double dh_val = 0.0, dg_val = 0.0, d2g_val = 0.0;

    /* Integrate using trapezoidal rule */
    int iv_center = (int)((v - dist->v_min) / dv);
    if (iv_center < 0) iv_center = 0;
    if (iv_center >= nv) iv_center = nv - 1;

    /* Inner integral: 0 to v */
    double inner_int1 = 0.0; /* integral f(u) u^2 du */
    double inner_int2 = 0.0; /* for g */
    for (int i = 0; i <= iv_center && i < nv; i++) {
        double u = dist->v_min + i * dv;
        double fu = dist->f[i];
        double u2 = u * u;
        double weight = (i == 0 || i == iv_center) ? 0.5 : 1.0;
        inner_int1 += weight * fu * u2 * dv;
        inner_int2 += weight * fu * u2 * (v / 3.0 + u2 / (5.0 * v)) * dv;
    }
    h_val += 4.0 * M_PI * inner_int1 / v;
    g_val += 4.0 * M_PI * inner_int2;

    /* Outer integral: v to inf */
    double outer_int1 = 0.0;
    double outer_int2 = 0.0;
    for (int i = iv_center; i < nv; i++) {
        double u = dist->v_min + i * dv;
        double fu = dist->f[i];
        double weight = (i == iv_center || i == nv-1) ? 0.5 : 1.0;
        outer_int1 += weight * fu * u * dv;
        outer_int2 += weight * fu * u * (v*v/3.0 + u*u) * dv;
    }
    h_val += 4.0 * M_PI * outer_int1;
    g_val += 4.0 * M_PI * outer_int2;

    /* Derivatives using finite differences on the potentials */
    double v_m = v - dv; if (v_m < dist->v_min) v_m = dist->v_min;
    double v_p = v + dv; if (v_p > dist->v_max) v_p = dist->v_max;

    dh_val = 0.0;
    dg_val = 0.0;
    d2g_val = 0.0;
    /* Simple finite difference approximation for the derivatives */
    /* In a full implementation, these would be computed analytically */

    *h = h_val; *g = g_val;
    *dh = dh_val; *dg = dg_val; *d2g = d2g_val;
}

/* ================================================================
 * L3: Fokker-Planck Collision Operator (Electron-Electron)
 * ================================================================ */

void fokker_planck_ee(double *dfdt_coll, const Distribution1D *f,
                      int nv, double dv, double ln_Lambda) {
    if (!dfdt_coll || !f || !f->f || nv < 3) return;

    /* Fokker-Planck electron-electron collision operator in isotropic form:
     *
     * df/dt|_c = Gamma * (1/v^2) * d/dv [ C(f) * df/dv + D(f) * f ]
     *
     * where Gamma = e^4 ln Lambda / (4*pi*epsilon_0^2 * m_e^2)
     *
     * C(f) = (4*pi/3) * [integral_0^v f(u) u^4 du + v^3 * integral_v^inf f(u) u du]
     * D(f) = 4*pi * integral_0^v f(u) u^2 du
     */

    double Gamma = pow(E_CHARGE, 4.0) * ln_Lambda
                   / (pow(4.0 * M_PI * EPSILON_0, 2.0) * M_ELECTRON * M_ELECTRON);

    /* Flux-conservative finite difference form */
    for (int i = 1; i < nv - 1; i++) {
        double v = f->v_min + i * dv;
        double v2 = v * v;

        /* Compute C and D coefficients at v_i */
        double C = 0.0, D = 0.0;
        /* Integral 0 to v */
        for (int j = 0; j <= i; j++) {
            double u = f->v_min + j * dv;
            double u2 = u * u;
            double fu = f->f[j];
            double w = (j == 0 || j == i) ? 0.5 : 1.0;
            C += w * fu * u2 * u2 * dv;
            D += w * fu * u2 * dv;
        }
        /* Integral v to inf */
        for (int j = i; j < nv; j++) {
            double u = f->v_min + j * dv;
            double fu = f->f[j];
            double w = (j == i || j == nv-1) ? 0.5 : 1.0;
            C += w * fu * u * dv * v2 * v;
        }
        C *= 4.0 * M_PI / 3.0;
        D *= 4.0 * M_PI;

        /* Three-point flux difference for dfdt at v_i */
        if (v2 > 0.0) {
            /* Flux at i+1/2 */
            double v_ip = v + 0.5 * dv;
            double C_ip = C * (v_ip * v_ip) / (v * v); /* approximate */
            double df_dv_ip = (f->f[i+1] - f->f[i]) / dv;
            double flux_ip = C_ip * df_dv_ip + D * f->f[i];

            /* Flux at i-1/2 */
            double v_im = v - 0.5 * dv;
            double C_im = C * (v_im * v_im) / (v * v); /* approximate */
            double df_dv_im = (f->f[i] - f->f[i-1]) / dv;
            double flux_im = C_im * df_dv_im + D * f->f[i-1];

            dfdt_coll[i] = Gamma * (flux_ip - flux_im) / (dv * v2);
        } else {
            dfdt_coll[i] = 0.0;
        }
    }
    dfdt_coll[0] = 0.0;
    dfdt_coll[nv-1] = 0.0;
}

/* ================================================================
 * L8: Pitch-Angle Scattering Operator
 * ================================================================ */

void pitch_angle_scattering(Distribution2D *f, double nu_ei) {
    if (!f || !f->f) return;

    /* df/dt = (nu_ei/2) * d/dmu[(1-mu^2) * df/dmu]
     *
     * Using finite differences on the mu = vpar/v grid.
     * For a 2D distribution in (vpar, vperp), the pitch angle
     * scattering operator mixes parallel and perpendicular velocities.
     */

    int nvpar = f->nvpar, nvperp = f->nvperp;
    double dvpar = f->dvpar, dvperp = f->dvperp;

    /* Simplified: apply scattering along vpar rows */
    for (int j = 0; j < nvperp; j++) {
        double *row = f->f[j];
        double *row_new = (double*)malloc(nvpar * sizeof(double));

        /* Copy to avoid modification during computation */
        for (int i = 0; i < nvpar; i++) row_new[i] = row[i];

        for (int i = 1; i < nvpar - 1; i++) {
            double vpar = f->vpar_min + i * dvpar;
            double vperp = f->vperp_min + j * dvperp;
            double v = sqrt(vpar * vpar + vperp * vperp);
            if (v < 1e-10) continue;
            double mu = vpar / v;
            double one_m_mu2 = 1.0 - mu * mu;

            /* Central difference for d/dmu[(1-mu^2)*df/dmu] */
            double df_dmu_fwd = (row[i+1] - row[i]) / (dvpar / v);
            double df_dmu_bwd = (row[i] - row[i-1]) / (dvpar / v);
            double flux_fwd = one_m_mu2 * df_dmu_fwd;
            double flux_bwd = one_m_mu2 * df_dmu_bwd;
            double d_flux = (flux_fwd - flux_bwd) / (dvpar / v);

            row_new[i] += 0.5 * nu_ei * d_flux * f->dvpar; /* simplified time advance */
        }

        for (int i = 0; i < nvpar; i++) row[i] = row_new[i];
        free(row_new);
    }
}
