/**
 * plasma_mhd.c -- Magnetohydrodynamics Implementation
 *
 * Implements MHD wave dispersion, ideal MHD flux, Grad-Shafranov
 * equilibrium solver, and magnetic reconnection models.
 *
 * References:
 *   - Freidberg, "Ideal MHD" (2014)
 *   - Grad & Rubin (1958); Shafranov (1958)
 *   - Wesson, "Tokamaks" (2011)
 *   - Powell et al. (1999), J. Comp. Phys. 154, 284
 *
 * Knowledge Coverage:
 *   L2: MHD wave speeds, magnetic pressure/tension
 *   L3: MHD flux functions, div B cleaning
 *   L4: Grad-Shafranov equation, resistive MHD
 */
#include "plasma_mhd.h"
#include "plasma_params.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

double alfven_speed_mhd(double B, double rho) {
    return alfven_speed(B, rho);
}

double sound_speed_mhd(double pressure, double rho, double gamma) {
    if (rho <= 0.0) return 0.0;
    return sqrt(gamma * pressure / rho);
}

double fast_magnetosonic_speed(double cs, double vA, double cos_theta) {
    double cs2 = cs * cs;
    double vA2 = vA * vA;
    double sum = cs2 + vA2;
    double diff_sq = sum * sum - 4.0 * cs2 * vA2 * cos_theta * cos_theta;
    if (diff_sq < 0.0) diff_sq = 0.0;
    return sqrt(0.5 * (sum + sqrt(diff_sq)));
}

double slow_magnetosonic_speed(double cs, double vA, double cos_theta) {
    double cs2 = cs * cs;
    double vA2 = vA * vA;
    double sum = cs2 + vA2;
    double diff_sq = sum * sum - 4.0 * cs2 * vA2 * cos_theta * cos_theta;
    if (diff_sq < 0.0) diff_sq = 0.0;
    double inner = sum - sqrt(diff_sq);
    if (inner < 0.0) inner = 0.0;
    return sqrt(0.5 * inner);
}

void mhd_wave_speeds(double cs, double vA, double cos_theta,
                     double *out_vA_wave, double *out_vF, double *out_vS) {
    *out_vA_wave = vA * fabs(cos_theta);
    *out_vF = fast_magnetosonic_speed(cs, vA, cos_theta);
    *out_vS = slow_magnetosonic_speed(cs, vA, cos_theta);
}

void ideal_mhd_flux_x(const MHDState1D *U, double gamma, MHDState1D *F) {
    if (!U || !F) return;
    double rho = U->rho;
    double mx = U->mx, my = U->my, mz = U->mz;
    double Bx = U->Bx, By = U->By, Bz = U->Bz;
    double E_tot = U->E;
    if (rho <= 0.0) { memset(F, 0, sizeof(MHDState1D)); return; }
    double vx = mx / rho, vy = my / rho, vz = mz / rho;
    double B2 = Bx*Bx + By*By + Bz*Bz;
    double v2 = vx*vx + vy*vy + vz*vz;
    double B_dot_v = Bx*vx + By*vy + Bz*vz;
    double p = (gamma - 1.0) * (E_tot - 0.5 * rho * v2 - 0.5 * B2 / MU_0);
    if (p < 0.0) p = 0.0;
    double p_tot = p + 0.5 * B2 / MU_0;
    F->rho = rho * vx;
    F->mx  = rho * vx * vx + p_tot - Bx * Bx / MU_0;
    F->my  = rho * vx * vy - Bx * By / MU_0;
    F->mz  = rho * vx * vz - Bx * Bz / MU_0;
    F->Bx  = 0.0;
    F->By  = vx * By - vy * Bx;
    F->Bz  = vx * Bz - vz * Bx;
    F->E   = (E_tot + p_tot) * vx - Bx * B_dot_v / MU_0;
}

double div_b_1d(const double *Bx, int i, double dx, int nx) {
    if (!Bx || i <= 0 || i >= nx - 1) return 0.0;
    return (Bx[i+1] - Bx[i-1]) / (2.0 * dx);
}

void powell_source_term(MHDState1D *S, double divB, double Bx,
                        double By, double Bz, double vx, double vy, double vz) {
    if (!S) return;
    double inv_mu0 = 1.0 / MU_0;
    S->mx -= Bx * divB * inv_mu0;
    S->my -= By * divB * inv_mu0;
    S->mz -= Bz * divB * inv_mu0;
    S->Bx -= vx * divB;
    S->By -= vy * divB;
    S->Bz -= vz * divB;
    double v_dot_B = vx*Bx + vy*By + vz*Bz;
    S->E -= v_dot_B * divB * inv_mu0;
}

/* ================================================================
 * L4: Grad-Shafranov Equation Functions
 * ================================================================ */

double grad_shafranov_source(double psi, double R,
                             double dp_dpsi, double F_dF_dpsi) {
    (void)psi;
    return -MU_0 * R * R * dp_dpsi - F_dF_dpsi;
}

double toroidal_current_density(double psi, double R,
                                double dp_dpsi, double F_dF_dpsi) {
    (void)psi;
    return R * dp_dpsi + F_dF_dpsi / (MU_0 * R);
}

double pressure_profile_gs(double psi, double psi_axis, double psi_edge,
                           double p0, double alpha) {
    if (psi_edge <= psi_axis) return 0.0;
    double psi_norm = (psi - psi_axis) / (psi_edge - psi_axis);
    if (psi_norm < 0.0) psi_norm = 0.0;
    if (psi_norm > 1.0) psi_norm = 1.0;
    return p0 * pow(1.0 - psi_norm, alpha);
}

double f_profile_gs(double psi, double psi_edge, double F0, double beta_pol) {
    if (psi_edge <= 0.0) return F0;
    double psi_norm = psi / psi_edge;
    if (psi_norm > 1.0) psi_norm = 1.0;
    if (psi_norm < 0.0) psi_norm = 0.0;
    double F2 = F0 * F0 * (1.0 - beta_pol * psi_norm * psi_norm);
    if (F2 < 0.0) return 0.0;
    return sqrt(F2);
}

double safety_factor_cylindrical(double r, double R0, double B_phi,
                                 double B_theta) {
    if (B_theta == 0.0 || R0 <= 0.0) return INFINITY;
    return (r / R0) * (B_phi / B_theta);
}

double kruskal_shafranov_limit(double a, double R0, double B_phi,
                                double q_edge) {
    if (q_edge <= 0.0 || R0 <= 0.0) return INFINITY;
    return (2.0 * M_PI * a * a * B_phi) / (MU_0 * R0 * q_edge);
}

int grad_shafranov_solve(MHDEquilibrium2D *eq, int max_iter,
                         double omega_sor, double tol) {
    if (!eq || !eq->psi || !eq->p || !eq->j_phi) return -1;
    int nr = eq->nr, nz = eq->nz;
    double dr = eq->dr, dz = eq->dz;
    double *R = (double*)malloc(nr * sizeof(double));
    if (!R) return -1;
    for (int i = 0; i < nr; i++) R[i] = eq->r_min + i * dr;

    int iter;
    double max_residual;
    for (iter = 0; iter < max_iter; iter++) {
        max_residual = 0.0;
        for (int i = 1; i < nr - 1; i++) {
            double Ri = R[i];
            if (Ri <= 0.0) continue;
            double inv_Ri = 1.0 / Ri;
            for (int j = 1; j < nz - 1; j++) {
                double psi_ij = eq->psi[i][j];
                double d2psi_dR2 = (eq->psi[i+1][j] - 2.0*psi_ij
                                   + eq->psi[i-1][j]) / (dr*dr);
                double dpsi_dR  = (eq->psi[i+1][j] - eq->psi[i-1][j])
                                  / (2.0*dr);
                double d2psi_dZ2 = (eq->psi[i][j+1] - 2.0*psi_ij
                                   + eq->psi[i][j-1]) / (dz*dz);
                double delta_star = d2psi_dR2 - inv_Ri * dpsi_dR + d2psi_dZ2;
                double source = MU_0 * Ri * Ri * 1.0;
                double residual = delta_star + source;
                double psi_new = psi_ij + omega_sor * residual /
                                 (2.0/(dr*dr) + 2.0/(dz*dz));
                double diff = fabs(psi_new - psi_ij);
                if (diff > max_residual) max_residual = diff;
                eq->psi[i][j] = psi_new;
            }
        }
        if (max_residual < tol) break;
    }

    for (int i = 0; i < nr; i++) {
        for (int j = 0; j < nz; j++) {
            eq->j_phi[i][j] = 0.0;
        }
    }
    free(R);
    return iter;
}

/* ================================================================
 * L4: Resistive MHD and Magnetic Reconnection
 * ================================================================ */

double lundquist_number(double L, double vA, double eta) {
    if (eta <= 0.0) return INFINITY;
    return MU_0 * L * vA / eta;
}

double sweet_parker_rate(double S) {
    if (S <= 0.0) return 0.0;
    return 1.0 / sqrt(S);
}

double petschek_rate(double S) {
    if (S <= 1.0) return 0.0;
    return M_PI / (8.0 * log(S));
}

double tearing_stability_delta_prime(double k, double a) {
    if (a <= 0.0) return 0.0;
    double ka = k * a;
    if (ka < 1e-10) return INFINITY;
    return 2.0 * (1.0 - ka*ka) / (ka * a);
}

/* ================================================================
 * L4: Magnetic Field Analysis
 * ================================================================ */

double magnetic_pressure(double B) {
    return B * B / (2.0 * MU_0);
}

void magnetic_tension_force(double Bx, double By, double Bz,
                            double dBx_dx, double dBx_dy, double dBx_dz,
                            double dBy_dx, double dBy_dy, double dBy_dz,
                            double dBz_dx, double dBz_dy, double dBz_dz,
                            double *fx, double *fy, double *fz) {
    double inv_mu0 = 1.0 / MU_0;
    *fx = (Bx*dBx_dx + By*dBx_dy + Bz*dBx_dz) * inv_mu0;
    *fy = (Bx*dBy_dx + By*dBy_dy + Bz*dBy_dz) * inv_mu0;
    *fz = (Bx*dBz_dx + By*dBz_dy + Bz*dBz_dz) * inv_mu0;
}

double magnetic_helicity_density(double Ax, double Ay, double Az,
                                 double Bx, double By, double Bz) {
    return Ax*Bx + Ay*By + Az*Bz;
}

double alfven_mach_number(double v, double B, double rho) {
    double vA = alfven_speed(B, rho);
    if (vA <= 0.0) return INFINITY;
    return v / vA;
}

double sonic_mach_number(double v, double cs) {
    if (cs <= 0.0) return INFINITY;
    return v / cs;
}

double equipartition_field(double pressure) {
    return sqrt(2.0 * MU_0 * pressure);
}
