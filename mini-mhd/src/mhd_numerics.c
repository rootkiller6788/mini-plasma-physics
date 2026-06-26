/*
 * mhd_numerics.c -- Numerical Methods for MHD Simulations
 *
 * Implements finite difference operators, slope limiters, MUSCL reconstruction,
 * HLL/HLLD/Roe Riemann solvers, Runge-Kutta time integration,
 * constrained transport, Powell 8-wave formulation, boundary conditions.
 *
 * Reference:
 *   Toro "Riemann Solvers and Numerical Methods for Fluid Dynamics" (2009) Ch.15-16
 *   Balsara "MHD Riemann Solvers" (1998-2016)
 *   Miyoshi & Kusano "A multi-state HLL approximate Riemann solver for MHD" JCP 208 (2005)
 *   Powell et al. "An Upwind Scheme for MHD" JCP 154 (1999)
 *   Evans & Hawley "Simulation of MHD Flows: A Constrained Transport Method" ApJ 332 (1988)
 *   Shu & Osher "Efficient implementation of TVD Runge-Kutta" JCP 77 (1988)
 *
 * Knowledge:
 *   L5 -- Numerical finite difference/volume operators
 *   L5 -- Slope limiters for TVD/MUSCL schemes
 *   L5 -- Riemann solvers for hyperbolic conservation laws
 *   L5 -- Time integration methods
 *   L5 -- Constrained transport and div(B) control
 */

#include "mhd_numerics.h"
#include "mhd_defs.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

/* ================================================================
 * L5 -- Finite Difference Operators
 * ================================================================ */

double mhd_fd_gradient_2nd(const double *f, int i, double dx) {
    if (!f || fabs(dx) < 1e-40) return 0.0;
    return (f[i+1] - f[i-1]) / (2.0 * dx);
}

double mhd_fd_gradient_4th(const double *f, int i, double dx) {
    if (!f || fabs(dx) < 1e-40) return 0.0;
    return (-f[i+2] + 8.0*f[i+1] - 8.0*f[i-1] + f[i-2]) / (12.0 * dx);
}

double mhd_fd_laplacian_2nd(const double *f, int i, double dx) {
    if (!f || fabs(dx) < 1e-40) return 0.0;
    return (f[i+1] - 2.0*f[i] + f[i-1]) / (dx * dx);
}

double mhd_fd_divergence_2nd(const double *Fx, const double *Fy, const double *Fz,
                              int i, int j, int k, int nx, int ny,
                              double dx, double dy, double dz) {
    if (!Fx || !Fy || !Fz) return 0.0;
    if (i < 1 || i >= nx-1 || j < 1 || j >= ny-1) return 0.0;

    
    int idx_ip = (i+1) + j*nx + k*nx*ny;
    int idx_im = (i-1) + j*nx + k*nx*ny;
    int idx_jp = i + (j+1)*nx + k*nx*ny;
    int idx_jm = i + (j-1)*nx + k*nx*ny;

    double dFx = (Fx[idx_ip] - Fx[idx_im]) / (2.0 * dx);
    double dFy = (Fy[idx_jp] - Fy[idx_jm]) / (2.0 * dy);
    double dFz = 0.0;
    if (dz > 1e-40) {
        int idx_kp = i + j*nx + (k+1)*nx*ny;
        int idx_km = i + j*nx + (k-1)*nx*ny;
        dFz = (Fz[idx_kp] - Fz[idx_km]) / (2.0 * dz);
    }
    return dFx + dFy + dFz;
}

void mhd_fd_curl_2nd(const double *Bx, const double *By, const double *Bz,
                      int i, int j, int k, int nx, int ny,
                      double dx, double dy, double dz,
                      double *curl_x, double *curl_y, double *curl_z) {
    if (!Bx || !By || !Bz || !curl_x || !curl_y || !curl_z) return;
    if (i < 1 || i >= nx-1 || j < 1 || j >= ny-1) {
        *curl_x = *curl_y = *curl_z = 0.0;
        return;
    }

    int idx_ip = (i+1) + j*nx + k*nx*ny;
    int idx_im = (i-1) + j*nx + k*nx*ny;
    int idx_jp = i + (j+1)*nx + k*nx*ny;
    int idx_jm = i + (j-1)*nx + k*nx*ny;

    double dBz_dy = (Bz[idx_jp] - Bz[idx_jm]) / (2.0 * dy);
    double dBy_dx = (By[idx_ip] - By[idx_im]) / (2.0 * dx);
    double dBx_dy = (Bx[idx_jp] - Bx[idx_jm]) / (2.0 * dy);

    double dBy_dz = 0.0, dBx_dz = 0.0, dBz_dx = 0.0;
    if (dz > 1e-40) {
        int idx_kp = i + j*nx + (k+1)*nx*ny;
        int idx_km = i + j*nx + (k-1)*nx*ny;
        dBy_dz = (By[idx_kp] - By[idx_km]) / (2.0 * dz);
        dBx_dz = (Bx[idx_kp] - Bx[idx_km]) / (2.0 * dz);
    }
    dBz_dx = (Bz[idx_ip] - Bz[idx_im]) / (2.0 * dx);

    *curl_x = dBz_dy - dBy_dz;
    *curl_y = dBx_dz - dBz_dx;
    *curl_z = dBy_dx - dBx_dy;
}

/* ================================================================
 * L5 -- Slope Limiters (Sweby 1984, LeVeque 2002 Ch.9)
 *
 * Limiters enforce TVD property, preventing spurious oscillations
 * near discontinuities. phi(r) with r = (u_i - u_{i-1})/(u_{i+1} - u_i).
 *
 * Sweby TVD region: 0 <= phi(r) <= min(2r, 2) for r >= 0, phi=0 for r<0.
 * ================================================================ */

double mhd_limiter_minmod(double r) {
    /* Most diffusive: phi = max(0, min(1, r)) */
    if (r <= 0.0) return 0.0;
    return (r < 1.0) ? r : 1.0;
}

double mhd_limiter_superbee(double r) {
    /* Most compressive: phi = max(0, min(2r,1), min(r,2)) */
    if (r <= 0.0) return 0.0;
    double m1 = (2.0*r < 1.0) ? 2.0*r : 1.0;
    double m2 = (r < 2.0) ? r : 2.0;
    return (m1 > m2) ? m1 : m2;
}

double mhd_limiter_van_leer(double r) {
    /* Harmonic mean: phi = (r+|r|)/(1+|r|) */
    if (r <= 0.0) return 0.0;
    return 2.0 * r / (1.0 + r);
}

double mhd_limiter_mc(double r) {
    /* Monotonized Central: phi = max(0, min(2r, (1+r)/2, 2)) */
    if (r <= 0.0) return 0.0;
    double a = 2.0 * r;
    double b = 0.5 * (1.0 + r);
    double m1 = (a < b) ? a : b;
    return (m1 < 2.0) ? m1 : 2.0;
}

/* ================================================================
 * L5 -- MUSCL-Hancock Reconstruction (van Leer 1979)
 *
 * 2nd-order in space and time. Predictor-corrector without
 * requiring full Riemann solver at sub-step.
 * ================================================================ */

void mhd_muscl_extrapolate(const double *u, int i, double dx, double dt,
                            double (*limiter)(double),
                            double *uL, double *uR) {
    if (!u || !limiter || !uL || !uR || fabs(dx) < 1e-40) {
        if (uL && u) *uL = u[i];
        if (uR && u) *uR = u[i+1];
        return;
    }

    double du_ip1 = u[i+1] - u[i];
    double du_i   = u[i] - u[i-1];
    double r_i    = (fabs(du_ip1) < 1e-40) ? 0.0 : du_i / du_ip1;
    double r_ip1  = (fabs(u[i+2] - u[i+1]) < 1e-40) ? 0.0
                   : du_ip1 / (u[i+2] - u[i+1]);

    double slope_i   = limiter(r_i)   * du_ip1 / dx;
    double slope_ip1 = limiter(r_ip1) * (u[i+2] - u[i+1]) / dx;

    /* Half-time predictor: u* = u - 0.5*(dt/dx)*slope */
    double u_star_i   = u[i]   - 0.5 * (dt/dx) * slope_i;
    double u_star_ip1 = u[i+1] - 0.5 * (dt/dx) * slope_ip1;

    /* Reconstruct left/right states at interface i+1/2 */
    *uL = u_star_i   + 0.5 * dx * slope_i;
    *uR = u_star_ip1 - 0.5 * dx * slope_ip1;
}

void mhd_muscl_extrapolate_mhd(const MHDState *states, int i,
                                double dx, double dt,
                                double (*limiter)(double),
                                MHDState *left, MHDState *right) {
    if (!states || !limiter || !left || !right || fabs(dx) < 1e-40) {
        if (states && left)  *left  = states[i];
        if (states && right) *right = states[i+1];
        return;
    }

    /* Flatten MHD 8-field state into scalar arrays for reconstruction */
    double u[8][4], uL[8], uR[8];
    for (int k = -1; k <= 2; k++) {
        int idx = i + k;
        u[0][k+1] = states[idx].rho;
        u[1][k+1] = states[idx].vx;
        u[2][k+1] = states[idx].vy;
        u[3][k+1] = states[idx].vz;
        u[4][k+1] = states[idx].Bx;
        u[5][k+1] = states[idx].By;
        u[6][k+1] = states[idx].Bz;
        u[7][k+1] = states[idx].p;
    }

    for (int f = 0; f < 8; f++) {
        double du_ip1 = u[f][2] - u[f][1];
        double du_i   = u[f][1] - u[f][0];
        double r_i    = (fabs(du_ip1) < 1e-40) ? 0.0 : du_i / du_ip1;
        double r_ip1  = (fabs(u[f][3] - u[f][2]) < 1e-40) ? 0.0
                       : du_ip1 / (u[f][3] - u[f][2]);

        double sl_i   = limiter(r_i)   * du_ip1 / dx;
        double sl_ip1 = limiter(r_ip1) * (u[f][3] - u[f][2]) / dx;

        uL[f] = u[f][1] - 0.5*(dt/dx)*sl_i   + 0.5*dx*sl_i;
        uR[f] = u[f][2] - 0.5*(dt/dx)*sl_ip1 - 0.5*dx*sl_ip1;
    }

    left->rho  = uL[0]; left->vx = uL[1]; left->vy = uL[2]; left->vz = uL[3];
    left->Bx   = uL[4]; left->By = uL[5]; left->Bz = uL[6]; left->p  = uL[7];
    right->rho = uR[0]; right->vx= uR[1];right->vy = uR[2];right->vz = uR[3];
    right->Bx  = uR[4]; right->By= uR[5];right->Bz = uR[6];right->p  = uR[7];

    if (left->rho  < 1e-40) left->rho  = 1e-40;
    if (right->rho < 1e-40) right->rho = 1e-40;
    if (left->p    < 1e-40) left->p    = 1e-40;
    if (right->p   < 1e-40) right->p   = 1e-40;
}

/* ================================================================
 * L5 -- HLL Riemann Solver
 *
 * F_HLL = (S_R*F_L - S_L*F_R + S_R*S_L*(U_R - U_L)) / (S_R - S_L)
 *
 * Uses fast magnetosonic speeds for S_L, S_R.
 * ================================================================ */

MHDFlux mhd_hllc_simple(const MHDState *left, const MHDState *right,
                         double gamma) {
    MHDFlux flux;
    memset(&flux, 0, sizeof(MHDFlux));
    if (!left || !right) return flux;

    MHDFlux FL, FR;
    mhd_flux_compute(left,  1.0, 0.0, 0.0, gamma, &FL);
    mhd_flux_compute(right, 1.0, 0.0, 0.0, gamma, &FR);

    double Bmag_L = mhd_vector_magnitude(left->Bx, left->By, left->Bz);
    double Bmag_R = mhd_vector_magnitude(right->Bx, right->By, right->Bz);
    double cs_L = mhd_sound_speed(left->p, left->rho, gamma);
    double cs_R = mhd_sound_speed(right->p, right->rho, gamma);
    double va_L = mhd_alfven_speed(Bmag_L, left->rho);
    double va_R = mhd_alfven_speed(Bmag_R, right->rho);
    double cf_L = sqrt(cs_L*cs_L + va_L*va_L);
    double cf_R = sqrt(cs_R*cs_R + va_R*va_R);

    double SL = left->vx - cf_L;
    double SR = right->vx + cf_R;

    if (fabs(SR - SL) < 1e-40) {
        return (left->vx + right->vx >= 0.0) ? FL : FR;
    }

    double inv = 1.0 / (SR - SL);
    double SS = SR * SL * inv;

    double E_L = mhd_total_energy_density(left->rho,
        mhd_vector_magnitude(left->vx,left->vy,left->vz), Bmag_L, left->p, gamma);
    double E_R = mhd_total_energy_density(right->rho,
        mhd_vector_magnitude(right->vx,right->vy,right->vz), Bmag_R, right->p, gamma);

    flux.D_flux  = (SR*FL.D_flux  - SL*FR.D_flux  + SS*(right->rho - left->rho)) * inv;
    flux.Mx_flux = (SR*FL.Mx_flux - SL*FR.Mx_flux + SS*(right->rho*right->vx - left->rho*left->vx)) * inv;
    flux.My_flux = (SR*FL.My_flux - SL*FR.My_flux + SS*(right->rho*right->vy - left->rho*left->vy)) * inv;
    flux.Mz_flux = (SR*FL.Mz_flux - SL*FR.Mz_flux + SS*(right->rho*right->vz - left->rho*left->vz)) * inv;
    flux.Bx_flux = 0.0;
    flux.By_flux = (SR*FL.By_flux - SL*FR.By_flux + SS*(right->By - left->By)) * inv;
    flux.Bz_flux = (SR*FL.Bz_flux - SL*FR.Bz_flux + SS*(right->Bz - left->Bz)) * inv;
    flux.E_flux  = (SR*FL.E_flux  - SL*FR.E_flux  + SS*(E_R - E_L)) * inv;

    return flux;
}

/* ================================================================
 * L5 -- HLLD Riemann Solver (Miyoshi & Kusano 2005)
 *
 * 5-state solver: S_L -> S*_L -> S_M -> S*_R -> S_R
 * Resolves entropy, Alfven, fast, slow waves.
 * ================================================================ */

MHDFlux mhd_hlld_flux(const MHDState *left, const MHDState *right,
                       double gamma) {
    MHDFlux flux;
    memset(&flux, 0, sizeof(MHDFlux));
    if (!left || !right) return flux;

    double Bmag_L = mhd_vector_magnitude(left->Bx, left->By, left->Bz);
    double Bmag_R = mhd_vector_magnitude(right->Bx, right->By, right->Bz);
    double cs_L = mhd_sound_speed(left->p, left->rho, gamma);
    double cs_R = mhd_sound_speed(right->p, right->rho, gamma);
    double va_L = mhd_alfven_speed(Bmag_L, left->rho);
    double va_R = mhd_alfven_speed(Bmag_R, right->rho);
    double cf_L = sqrt(cs_L*cs_L + va_L*va_L);
    double cf_R = sqrt(cs_R*cs_R + va_R*va_R);

    double SL = (left->vx - cf_L < right->vx - cf_R) ?
                 left->vx - cf_L : right->vx - cf_R;
    double SR = (left->vx + cf_L > right->vx + cf_R) ?
                 left->vx + cf_L : right->vx + cf_R;

    if (SR - SL < 1e-40) {
        if (left->vx + right->vx >= 0.0)
            mhd_flux_compute(left, 1.0, 0.0, 0.0, gamma, &flux);
        else
            mhd_flux_compute(right, 1.0, 0.0, 0.0, gamma, &flux);
        return flux;
    }

    MHDFlux FL, FR;
    mhd_flux_compute(left,  1.0, 0.0, 0.0, gamma, &FL);
    mhd_flux_compute(right, 1.0, 0.0, 0.0, gamma, &FR);

    if (SL >= 0.0) return FL;
    if (SR <= 0.0) return FR;

    double rho_L = left->rho, rho_R = right->rho;
    double rho_L_SR = rho_L * (SL - left->vx);
    double rho_R_SL = rho_R * (SR - right->vx);
    double rho_star = (rho_L_SR - rho_R_SL) / (SR - SL);
    double inv_rho_star = 1.0 / (rho_star + 1e-40);

    double Mx_star = (rho_L_SR*left->vx - rho_R_SL*right->vx
                     + (FR.Mx_flux - FL.Mx_flux)) / (SR - SL);
    double SM = Mx_star * inv_rho_star;

    double E_L = mhd_total_energy_density(rho_L,
        mhd_vector_magnitude(left->vx,left->vy,left->vz), Bmag_L, left->p, gamma);
    double E_R = mhd_total_energy_density(rho_R,
        mhd_vector_magnitude(right->vx,right->vy,right->vz), Bmag_R, right->p, gamma);
    double E_star = (rho_L_SR*E_L - rho_R_SL*E_R
                    + FR.E_flux - FL.E_flux) / (SR - SL);

    double sqL = sqrt(rho_L), sqR = sqrt(rho_R);
    double sq_sum = sqL + sqR;
    double By_star, Bz_star, vy_star, vz_star;
    if (sq_sum > 1e-40) {
        By_star = (sqR*left->By + sqL*right->By
                  + sqL*sqR*(right->vy - left->vy)*sqrt(MHD_MU0)) / sq_sum;
        Bz_star = (sqR*left->Bz + sqL*right->Bz
                  + sqL*sqR*(right->vz - left->vz)*sqrt(MHD_MU0)) / sq_sum;
        vy_star = (sqL*left->vy + sqR*right->vy
                  + (right->By - left->By)/sqrt(MHD_MU0)) / sq_sum;
        vz_star = (sqL*left->vz + sqR*right->vz
                  + (right->Bz - left->Bz)/sqrt(MHD_MU0)) / sq_sum;
    } else {
        By_star = 0.5*(left->By + right->By);
        Bz_star = 0.5*(left->Bz + right->Bz);
        vy_star = 0.5*(left->vy + right->vy);
        vz_star = 0.5*(left->vz + right->vz);
    }

    if (SM >= 0.0) {
        flux.D_flux  = FL.D_flux  + SL * (rho_star - rho_L);
        flux.Mx_flux = FL.Mx_flux + SL * (rho_star*SM - rho_L*left->vx);
        flux.My_flux = FL.My_flux + SL * (rho_star*vy_star - rho_L*left->vy);
        flux.Mz_flux = FL.Mz_flux + SL * (rho_star*vz_star - rho_L*left->vz);
        flux.By_flux = FL.By_flux + SL * (By_star - left->By);
        flux.Bz_flux = FL.Bz_flux + SL * (Bz_star - left->Bz);
        flux.Bx_flux = 0.0;
        flux.E_flux  = FL.E_flux  + SL * (E_star - E_L);
    } else {
        flux.D_flux  = FR.D_flux  + SR * (rho_star - rho_R);
        flux.Mx_flux = FR.Mx_flux + SR * (rho_star*SM - rho_R*right->vx);
        flux.My_flux = FR.My_flux + SR * (rho_star*vy_star - rho_R*right->vy);
        flux.Mz_flux = FR.Mz_flux + SR * (rho_star*vz_star - rho_R*right->vz);
        flux.By_flux = FR.By_flux + SR * (By_star - right->By);
        flux.Bz_flux = FR.Bz_flux + SR * (Bz_star - right->Bz);
        flux.Bx_flux = 0.0;
        flux.E_flux  = FR.E_flux  + SR * (E_star - E_R);
    }
    return flux;
}

/* ================================================================
 * L5 -- Roe-type Riemann Solver
 *
 * Uses the Roe-averaged MHD state. For MHD, B fields are averaged
 * with opposite density weighting: B_roe = (sqrt(rho_R)*B_L + sqrt(rho_L)*B_R)/(sqL+sqR)
 *
 * Flux: F = 0.5*(F_L+F_R) - 0.5*|A|*(U_R-U_L)
 * We use max(|lambda|) for the dissipation (LLF-type fix).
 * ================================================================ */

static void mhd_roe_average(const MHDState *L, const MHDState *R, double gamma,
                             MHDState *roe) {
    double sqL = sqrt(L->rho), sqR = sqrt(R->rho);
    double inv = 1.0 / (sqL + sqR + 1e-40);
    roe->rho = sqL * sqR;
    roe->vx  = (sqL*L->vx + sqR*R->vx) * inv;
    roe->vy  = (sqL*L->vy + sqR*R->vy) * inv;
    roe->vz  = (sqL*L->vz + sqR*R->vz) * inv;
    roe->Bx  = (sqR*L->Bx + sqL*R->Bx) * inv;
    roe->By  = (sqR*L->By + sqL*R->By) * inv;
    roe->Bz  = (sqR*L->Bz + sqL*R->Bz) * inv;
    double B2_L = L->Bx*L->Bx+L->By*L->By+L->Bz*L->Bz;
    double B2_R = R->Bx*R->Bx+R->By*R->By+R->Bz*R->Bz;
    double H_L = gamma*L->p/((gamma-1.0)*L->rho)
               + 0.5*(L->vx*L->vx+L->vy*L->vy+L->vz*L->vz) + B2_L*MHD_MU0_INV/L->rho;
    double H_R = gamma*R->p/((gamma-1.0)*R->rho)
               + 0.5*(R->vx*R->vx+R->vy*R->vy+R->vz*R->vz) + B2_R*MHD_MU0_INV/R->rho;
    double H_roe = (sqL*H_L + sqR*H_R) * inv;
    double v2 = roe->vx*roe->vx+roe->vy*roe->vy+roe->vz*roe->vz;
    double B2 = roe->Bx*roe->Bx+roe->By*roe->By+roe->Bz*roe->Bz;
    roe->p = (gamma-1.0)/gamma*roe->rho*(H_roe-0.5*v2-B2*MHD_MU0_INV/roe->rho);
    if (roe->p < 1e-40) roe->p = 1e-40;
}

MHDFlux mhd_roe_flux(const MHDState *left, const MHDState *right,
                      double gamma) {
    MHDFlux flux;
    memset(&flux, 0, sizeof(MHDFlux));
    if (!left || !right) return flux;

    MHDState roe;
    mhd_roe_average(left, right, gamma, &roe);

    MHDFlux FL, FR;
    mhd_flux_compute(left,  1.0, 0.0, 0.0, gamma, &FL);
    mhd_flux_compute(right, 1.0, 0.0, 0.0, gamma, &FR);

    double Bmag = mhd_vector_magnitude(roe.Bx, roe.By, roe.Bz);
    double cs = mhd_sound_speed(roe.p, roe.rho, gamma);
    double va = mhd_alfven_speed(Bmag, roe.rho);
    double cf = sqrt(cs*cs + va*va);
    double lam_max = fabs(roe.vx) + cf;
    double half = 0.5;
    double diss = half * lam_max;

    double E_L = mhd_total_energy_density(left->rho,
        mhd_vector_magnitude(left->vx,left->vy,left->vz),
        mhd_vector_magnitude(left->Bx,left->By,left->Bz), left->p, gamma);
    double E_R = mhd_total_energy_density(right->rho,
        mhd_vector_magnitude(right->vx,right->vy,right->vz),
        mhd_vector_magnitude(right->Bx,right->By,right->Bz), right->p, gamma);

    flux.D_flux  = half*(FL.D_flux  + FR.D_flux)  - diss*(right->rho - left->rho);
    flux.Mx_flux = half*(FL.Mx_flux + FR.Mx_flux) - diss*(right->rho*right->vx - left->rho*left->vx);
    flux.My_flux = half*(FL.My_flux + FR.My_flux) - diss*(right->rho*right->vy - left->rho*left->vy);
    flux.Mz_flux = half*(FL.Mz_flux + FR.Mz_flux) - diss*(right->rho*right->vz - left->rho*left->vz);
    flux.Bx_flux = 0.0;
    flux.By_flux = half*(FL.By_flux + FR.By_flux) - diss*(right->By - left->By);
    flux.Bz_flux = half*(FL.Bz_flux + FR.Bz_flux) - diss*(right->Bz - left->Bz);
    flux.E_flux  = half*(FL.E_flux  + FR.E_flux)  - diss*(E_R - E_L);
    return flux;
}

/* ================================================================
 * L5 -- Runge-Kutta Time Integration
 *
 * For MHD: dU/dt = R(U) where R = -div[F(U)] + S(U).
 * SSP methods preserve TVD property under CFL constraint.
 * ================================================================ */

void mhd_rk2_step(MHDConserved *U, int n, double dt,
                   void (*rhs)(const MHDConserved*, MHDConserved*, int, void*),
                   void *ctx) {
    if (!U || !rhs || n <= 0) return;
    MHDConserved *Ut = (MHDConserved*)malloc((size_t)n * sizeof(MHDConserved));
    MHDConserved *R  = (MHDConserved*)malloc((size_t)n * sizeof(MHDConserved));
    if (!Ut || !R) { free(Ut); free(R); return; }

    /* Stage 1: U* = U^n + dt*R(U^n) */
    rhs(U, R, n, ctx);
    for (int i = 0; i < n; i++) {
        Ut[i].D  = U[i].D  + dt * R[i].D;
        Ut[i].Mx = U[i].Mx + dt * R[i].Mx;
        Ut[i].My = U[i].My + dt * R[i].My;
        Ut[i].Mz = U[i].Mz + dt * R[i].Mz;
        Ut[i].Bx = U[i].Bx + dt * R[i].Bx;
        Ut[i].By = U[i].By + dt * R[i].By;
        Ut[i].Bz = U[i].Bz + dt * R[i].Bz;
        Ut[i].E  = U[i].E  + dt * R[i].E;
    }

    /* Stage 2: U^{n+1} = 0.5*U^n + 0.5*(U* + dt*R(U*)) */
    rhs(Ut, R, n, ctx);
    for (int i = 0; i < n; i++) {
        U[i].D  = 0.5*(U[i].D  + Ut[i].D  + dt*R[i].D);
        U[i].Mx = 0.5*(U[i].Mx + Ut[i].Mx + dt*R[i].Mx);
        U[i].My = 0.5*(U[i].My + Ut[i].My + dt*R[i].My);
        U[i].Mz = 0.5*(U[i].Mz + Ut[i].Mz + dt*R[i].Mz);
        U[i].Bx = 0.5*(U[i].Bx + Ut[i].Bx + dt*R[i].Bx);
        U[i].By = 0.5*(U[i].By + Ut[i].By + dt*R[i].By);
        U[i].Bz = 0.5*(U[i].Bz + Ut[i].Bz + dt*R[i].Bz);
        U[i].E  = 0.5*(U[i].E  + Ut[i].E  + dt*R[i].E);
    }
    free(Ut); free(R);
}

void mhd_rk3_step(MHDConserved *U, int n, double dt,
                   void (*rhs)(const MHDConserved*, MHDConserved*, int, void*),
                   void *ctx) {
    if (!U || !rhs || n <= 0) return;
    MHDConserved *Us = (MHDConserved*)malloc((size_t)n * sizeof(MHDConserved));
    MHDConserved *R  = (MHDConserved*)malloc((size_t)n * sizeof(MHDConserved));
    if (!Us || !R) { free(Us); free(R); return; }

    /* Shu-Osher SSP RK3(3,3):
     * U^(1)   = U^n + dt*R(U^n)
     * U^(2)   = 3/4*U^n + 1/4*(U^(1) + dt*R(U^(1)))
     * U^{n+1} = 1/3*U^n + 2/3*(U^(2) + dt*R(U^(2))) */

    rhs(U, R, n, ctx);
    for (int i = 0; i < n; i++) {
        Us[i].D  = U[i].D  + dt*R[i].D;
        Us[i].Mx = U[i].Mx + dt*R[i].Mx;
        Us[i].My = U[i].My + dt*R[i].My;
        Us[i].Mz = U[i].Mz + dt*R[i].Mz;
        Us[i].Bx = U[i].Bx + dt*R[i].Bx;
        Us[i].By = U[i].By + dt*R[i].By;
        Us[i].Bz = U[i].Bz + dt*R[i].Bz;
        Us[i].E  = U[i].E  + dt*R[i].E;
    }

    rhs(Us, R, n, ctx);
    for (int i = 0; i < n; i++) {
        Us[i].D  = 0.75*U[i].D  + 0.25*(Us[i].D  + dt*R[i].D);
        Us[i].Mx = 0.75*U[i].Mx + 0.25*(Us[i].Mx + dt*R[i].Mx);
        Us[i].My = 0.75*U[i].My + 0.25*(Us[i].My + dt*R[i].My);
        Us[i].Mz = 0.75*U[i].Mz + 0.25*(Us[i].Mz + dt*R[i].Mz);
        Us[i].Bx = 0.75*U[i].Bx + 0.25*(Us[i].Bx + dt*R[i].Bx);
        Us[i].By = 0.75*U[i].By + 0.25*(Us[i].By + dt*R[i].By);
        Us[i].Bz = 0.75*U[i].Bz + 0.25*(Us[i].Bz + dt*R[i].Bz);
        Us[i].E  = 0.75*U[i].E  + 0.25*(Us[i].E  + dt*R[i].E);
    }

    rhs(Us, R, n, ctx);
    for (int i = 0; i < n; i++) {
        U[i].D  = (1.0/3.0)*U[i].D  + (2.0/3.0)*(Us[i].D  + dt*R[i].D);
        U[i].Mx = (1.0/3.0)*U[i].Mx + (2.0/3.0)*(Us[i].Mx + dt*R[i].Mx);
        U[i].My = (1.0/3.0)*U[i].My + (2.0/3.0)*(Us[i].My + dt*R[i].My);
        U[i].Mz = (1.0/3.0)*U[i].Mz + (2.0/3.0)*(Us[i].Mz + dt*R[i].Mz);
        U[i].Bx = (1.0/3.0)*U[i].Bx + (2.0/3.0)*(Us[i].Bx + dt*R[i].Bx);
        U[i].By = (1.0/3.0)*U[i].By + (2.0/3.0)*(Us[i].By + dt*R[i].By);
        U[i].Bz = (1.0/3.0)*U[i].Bz + (2.0/3.0)*(Us[i].Bz + dt*R[i].Bz);
        U[i].E  = (1.0/3.0)*U[i].E  + (2.0/3.0)*(Us[i].E  + dt*R[i].E);
    }
    free(Us); free(R);
}

void mhd_rk4_step(MHDConserved *U, int n, double dt,
                   void (*rhs)(const MHDConserved*, MHDConserved*, int, void*),
                   void *ctx) {
    if (!U || !rhs || n <= 0) return;
    MHDConserved *k1=(MHDConserved*)malloc((size_t)n*sizeof(MHDConserved));
    MHDConserved *k2=(MHDConserved*)malloc((size_t)n*sizeof(MHDConserved));
    MHDConserved *k3=(MHDConserved*)malloc((size_t)n*sizeof(MHDConserved));
    MHDConserved *k4=(MHDConserved*)malloc((size_t)n*sizeof(MHDConserved));
    MHDConserved *Ut=(MHDConserved*)malloc((size_t)n*sizeof(MHDConserved));
    if (!k1||!k2||!k3||!k4||!Ut) { free(k1);free(k2);free(k3);free(k4);free(Ut); return; }

    /* k1 = dt*R(U^n) */
    rhs(U, k1, n, ctx);
    for (int i = 0; i < n; i++) {
        double dth = 0.5*dt;
        Ut[i].D  = U[i].D  + dth*k1[i].D;
        Ut[i].Mx = U[i].Mx + dth*k1[i].Mx;
        Ut[i].My = U[i].My + dth*k1[i].My;
        Ut[i].Mz = U[i].Mz + dth*k1[i].Mz;
        Ut[i].Bx = U[i].Bx + dth*k1[i].Bx;
        Ut[i].By = U[i].By + dth*k1[i].By;
        Ut[i].Bz = U[i].Bz + dth*k1[i].Bz;
        Ut[i].E  = U[i].E  + dth*k1[i].E;
    }

    /* k2 = dt*R(U^n + k1/2) */
    rhs(Ut, k2, n, ctx);
    for (int i = 0; i < n; i++) {
        double dth = 0.5*dt;
        Ut[i].D  = U[i].D  + dth*k2[i].D;
        Ut[i].Mx = U[i].Mx + dth*k2[i].Mx;
        Ut[i].My = U[i].My + dth*k2[i].My;
        Ut[i].Mz = U[i].Mz + dth*k2[i].Mz;
        Ut[i].Bx = U[i].Bx + dth*k2[i].Bx;
        Ut[i].By = U[i].By + dth*k2[i].By;
        Ut[i].Bz = U[i].Bz + dth*k2[i].Bz;
        Ut[i].E  = U[i].E  + dth*k2[i].E;
    }

    /* k3 = dt*R(U^n + k2/2) */
    rhs(Ut, k3, n, ctx);
    for (int i = 0; i < n; i++) {
        Ut[i].D  = U[i].D  + dt*k3[i].D;
        Ut[i].Mx = U[i].Mx + dt*k3[i].Mx;
        Ut[i].My = U[i].My + dt*k3[i].My;
        Ut[i].Mz = U[i].Mz + dt*k3[i].Mz;
        Ut[i].Bx = U[i].Bx + dt*k3[i].Bx;
        Ut[i].By = U[i].By + dt*k3[i].By;
        Ut[i].Bz = U[i].Bz + dt*k3[i].Bz;
        Ut[i].E  = U[i].E  + dt*k3[i].E;
    }

    /* k4 = dt*R(U^n + k3) */
    rhs(Ut, k4, n, ctx);
    for (int i = 0; i < n; i++) {
        U[i].D  += (1.0/6.0)*(k1[i].D  + 2.0*k2[i].D  + 2.0*k3[i].D  + k4[i].D);
        U[i].Mx += (1.0/6.0)*(k1[i].Mx + 2.0*k2[i].Mx + 2.0*k3[i].Mx + k4[i].Mx);
        U[i].My += (1.0/6.0)*(k1[i].My + 2.0*k2[i].My + 2.0*k3[i].My + k4[i].My);
        U[i].Mz += (1.0/6.0)*(k1[i].Mz + 2.0*k2[i].Mz + 2.0*k3[i].Mz + k4[i].Mz);
        U[i].Bx += (1.0/6.0)*(k1[i].Bx + 2.0*k2[i].Bx + 2.0*k3[i].Bx + k4[i].Bx);
        U[i].By += (1.0/6.0)*(k1[i].By + 2.0*k2[i].By + 2.0*k3[i].By + k4[i].By);
        U[i].Bz += (1.0/6.0)*(k1[i].Bz + 2.0*k2[i].Bz + 2.0*k3[i].Bz + k4[i].Bz);
        U[i].E  += (1.0/6.0)*(k1[i].E  + 2.0*k2[i].E  + 2.0*k3[i].E  + k4[i].E);
    }
    free(k1); free(k2); free(k3); free(k4); free(Ut);
}

/* ================================================================
 * L5 -- Constrained Transport (Evans & Hawley 1988)
 *
 * Uses Stokes theorem on cell faces to guarantee div(B)=0.
 * B stored at face centers, EMF computed at edges.
 *
 * Induction: d/dt int B.dS = -oint E.dl
 * where E = -v x B (ideal) is the electromotive force.
 * ================================================================ */

void mhd_ct_emf_2d(double vx_face, double vy_face,
                    double Bx_corner, double By_corner,
                    double *Ez_edge) {
    if (!Ez_edge) return;
    /* E_z = -(v x B)_z = -(vx*By - vy*Bx) */
    *Ez_edge = -(vx_face * By_corner - vy_face * Bx_corner);
}

void mhd_ct_update_B_2d(double Ez_pp, double Ez_pm, double Ez_mp, double Ez_mm,
                         double dx, double dy,
                         double *dBx, double *dBy) {
    (void)Ez_mm;
    if (!dBx || !dBy) return;
    double inv_dx = (fabs(dx) > 1e-40) ? 1.0/dx : 0.0;
    double inv_dy = (fabs(dy) > 1e-40) ? 1.0/dy : 0.0;
    /* dBx/dt at i+1/2 = -(Ez(i+1/2,j+1/2) - Ez(i+1/2,j-1/2)) / dy */
    *dBx = -(Ez_pp - Ez_pm) * inv_dy;
    /* dBy/dt at j+1/2 = +(Ez(i+1/2,j+1/2) - Ez(i-1/2,j+1/2)) / dx */
    *dBy =  (Ez_pp - Ez_mp) * inv_dx;
}

/* ================================================================
 * L5 -- Powell 8-Wave Source Term (Powell et al. 1999)
 *
 * Makes the MHD system Galilean invariant and advects div(B) errors.
 * S = -div(B) * (0, B, v, v.B) in conserved variable ordering.
 *
 * When added as a source term, div(B)/rho is passively advected.
 * ================================================================ */

void mhd_powell_source(const MHDState *state,
                        double divB,
                        double *S_D, double *S_Mx, double *S_My, double *S_Mz,
                        double *S_Bx, double *S_By, double *S_Bz,
                        double *S_E) {
    if (!state) return;
    double src = -divB;
    double vB = state->vx*state->Bx + state->vy*state->By + state->vz*state->Bz;
    if (S_D)  *S_D  = 0.0;
    if (S_Mx) *S_Mx = src * state->Bx;
    if (S_My) *S_My = src * state->By;
    if (S_Mz) *S_Mz = src * state->Bz;
    if (S_Bx) *S_Bx = src * state->vx;
    if (S_By) *S_By = src * state->vy;
    if (S_Bz) *S_Bz = src * state->vz;
    if (S_E)  *S_E  = src * vB;
}

/* ================================================================
 * L5 -- MHD Boundary Conditions
 *
 * Ghost cell filling for finite-volume/finite-difference.
 * Supported types: periodic, outflow, reflecting, conducting wall, fixed.
 * ================================================================ */

static void mhd_apply_bc_one(MHDState *ghost, const MHDState *interior,
                              MHDBCType bc) {
    if (!ghost || !interior) return;
    switch (bc) {
    case MHD_BC_OUTFLOW:
        *ghost = *interior;
        break;
    case MHD_BC_REFLECTING:
        ghost->rho = interior->rho;
        ghost->vx  = -interior->vx;
        ghost->vy  =  interior->vy;
        ghost->vz  =  interior->vz;
        ghost->Bx  = -interior->Bx;
        ghost->By  =  interior->By;
        ghost->Bz  =  interior->Bz;
        ghost->p   =  interior->p;
        break;
    case MHD_BC_CONDUCTING_WALL:
        ghost->rho = interior->rho;
        ghost->vx  = 0.0;
        ghost->vy  = 0.0;
        ghost->vz  = 0.0;
        ghost->Bx  = interior->Bx;
        ghost->By  = 0.0;
        ghost->Bz  = 0.0;
        ghost->p   = interior->p;
        break;
    case MHD_BC_PERIODIC:
    case MHD_BC_FIXED:
    default:
        break;
    }
}

void mhd_apply_boundary_1d(MHDState *states, int n, MHDBCType left, MHDBCType right) {
    if (!states || n < 3) return;
    if (left == MHD_BC_PERIODIC) {
        states[0]   = states[n-2];
        states[n-1] = states[1];
    } else {
        mhd_apply_bc_one(&states[0],   &states[1],   left);
        mhd_apply_bc_one(&states[n-1], &states[n-2], right);
    }
}

void mhd_apply_boundary_2d(MHDState *states, int nx, int ny,
                            MHDBCType xl, MHDBCType xr,
                            MHDBCType yl, MHDBCType yr) {
    if (!states || nx < 3 || ny < 3) return;

    /* X edges */
    for (int j = 0; j < ny; j++) {
        if (xl == MHD_BC_PERIODIC) states[0     + j*nx] = states[nx-2 + j*nx];
        else mhd_apply_bc_one(&states[0    + j*nx], &states[1    + j*nx], xl);
        if (xr == MHD_BC_PERIODIC) states[nx-1 + j*nx] = states[1 + j*nx];
        else mhd_apply_bc_one(&states[nx-1 + j*nx], &states[nx-2 + j*nx], xr);
    }
    /* Y edges */
    for (int i = 0; i < nx; i++) {
        if (yl == MHD_BC_PERIODIC) states[i + 0*nx]      = states[i + (ny-2)*nx];
        else mhd_apply_bc_one(&states[i + 0*nx],      &states[i + 1*nx],      yl);
        if (yr == MHD_BC_PERIODIC) states[i + (ny-1)*nx] = states[i + 1*nx];
        else mhd_apply_bc_one(&states[i + (ny-1)*nx], &states[i + (ny-2)*nx], yr);
    }
}

/* ================================================================
 * L5 -- CFL Timestep
 *
 * dt = CFL * min(dx,dy,dz) / max_i( |v_i| + c_{f,i} )
 *
 * c_f = sqrt(cs^2 + va^2) is the maximum MHD signal speed.
 * CFL ~ 0.3-0.5 for RK2, ~0.5-0.8 for TVD-RK3.
 * ================================================================ */

double mhd_cfl_timestep(const MHDState *states, int n,
                         double dx, double dy, double dz,
                         double gamma, double cfl) {
    if (!states || n <= 0) return 0.0;
    if (fabs(dx) < 1e-40 && fabs(dy) < 1e-40 && fabs(dz) < 1e-40)
        return INFINITY;

    double dl = dx;
    if (dy > 1e-40 && dy < dl) dl = dy;
    if (dz > 1e-40 && dz < dl) dl = dz;
    if (dl < 1e-40) dl = 1e-40;

    double max_spd = 0.0;
    for (int i = 0; i < n; i++) {
        double Bm = mhd_vector_magnitude(states[i].Bx, states[i].By, states[i].Bz);
        double cs = mhd_sound_speed(states[i].p, states[i].rho, gamma);
        double va = mhd_alfven_speed(Bm, states[i].rho);
        double cf = sqrt(cs*cs + va*va);
        double vm = mhd_vector_magnitude(states[i].vx, states[i].vy, states[i].vz);
        double sp = vm + cf;
        if (sp > max_spd) max_spd = sp;
    }
    if (max_spd < 1e-40) return INFINITY;
    return cfl * dl / max_spd;
}
