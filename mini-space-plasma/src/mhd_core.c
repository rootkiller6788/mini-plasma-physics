/**
 * mhd_core.c -- Magnetohydrodynamics Core Implementations
 *
 * Every function implements one independent MHD concept.
 * References:
 *   Priest "Solar Magnetohydrodynamics" (1982)
 *   Kivelson & Russell S4, Goldston & Rutherford S9
 *   MIT 22.615 MHD Theory
 *   Brio & Wu (1988) J. Comput. Phys. 75:400
 *   Powell et al. (1999) J. Comput. Phys. 154:284
 *
 * Knowledge: L1-L6 MHD theory and computation
 */
#include "../include/space_plasma.h"
#include "../include/mhd_core.h"
#include <math.h>
#include <string.h>

void mhd_prim_to_cons(const mhd_primitive_t *prim, mhd_conserved_t *cons, double gamma) {
    if (!prim || !cons) return;
    double rho = prim->rho;
    if (rho < 0.0) rho = 0.0;
    cons->D = rho;
    cons->Mx = rho * prim->vx;
    cons->My = rho * prim->vy;
    cons->Mz = rho * prim->vz;
    cons->Bx = prim->Bx;
    cons->By = prim->By;
    cons->Bz = prim->Bz;
    double v2 = prim->vx*prim->vx + prim->vy*prim->vy + prim->vz*prim->vz;
    double B2 = prim->Bx*prim->Bx + prim->By*prim->By + prim->Bz*prim->Bz;
    double gm1 = gamma - 1.0;
    if (gm1 <= 1e-15) gm1 = 2.0/3.0;
    cons->E = prim->p / gm1 + 0.5 * rho * v2 + 0.5 * B2 / SP_MU0;
}

int mhd_cons_to_prim(const mhd_conserved_t *cons, mhd_primitive_t *prim, double gamma) {
    if (!cons || !prim) return -1;
    double rho = cons->D;
    if (rho <= 0.0) return -1;
    prim->rho = rho;
    prim->vx = cons->Mx / rho;
    prim->vy = cons->My / rho;
    prim->vz = cons->Mz / rho;
    prim->Bx = cons->Bx;
    prim->By = cons->By;
    prim->Bz = cons->Bz;
    double v2 = prim->vx*prim->vx + prim->vy*prim->vy + prim->vz*prim->vz;
    double B2 = prim->Bx*prim->Bx + prim->By*prim->By + prim->Bz*prim->Bz;
    double gm1 = gamma - 1.0;
    if (gm1 <= 1e-15) gm1 = 2.0/3.0;
    double p = gm1 * (cons->E - 0.5 * rho * v2 - 0.5 * B2 / SP_MU0);
    if (p < 0.0) p = 0.0;
    prim->p = p;
    prim->E = cons->E;
    return (p >= 0.0) ? 0 : -1;
}

double mhd_energy_density(const mhd_primitive_t *prim, double gamma) {
    if (!prim) return 0.0;
    double gm1 = gamma - 1.0;
    if (gm1 <= 1e-15) gm1 = 2.0/3.0;
    double v2 = prim->vx*prim->vx + prim->vy*prim->vy + prim->vz*prim->vz;
    double B2 = prim->Bx*prim->Bx + prim->By*prim->By + prim->Bz*prim->Bz;
    return prim->p / gm1 + 0.5 * prim->rho * v2 + 0.5 * B2 / SP_MU0;
}

void mhd_flux_x(const mhd_primitive_t *prim, mhd_flux_t *flux, double gamma) {
    if (!prim || !flux) return;
    double rho = prim->rho;
    double vx = prim->vx, vy = prim->vy, vz = prim->vz;
    double Bx = prim->Bx, By = prim->By, Bz = prim->Bz;
    double p  = prim->p;
    double B2 = Bx*Bx + By*By + Bz*Bz;
    double Bdotv = Bx*vx + By*vy + Bz*vz;
    double ptot = p + 0.5 * B2 / SP_MU0;
    double E = mhd_energy_density(prim, gamma);
    flux->f[0] = rho * vx;
    flux->f[1] = rho * vx * vx + ptot - Bx*Bx/SP_MU0;
    flux->f[2] = rho * vx * vy - Bx*By/SP_MU0;
    flux->f[3] = rho * vx * vz - Bx*Bz/SP_MU0;
    flux->f[4] = 0.0;
    flux->f[5] = By * vx - Bx * vy;
    flux->f[6] = Bz * vx - Bx * vz;
    flux->f[7] = (E + ptot) * vx - Bx * Bdotv / SP_MU0;
}

double mhd_divB_error(const mhd_grid_t *grid) {
    if (!grid || !grid->fields || grid->nx < 2 || grid->ny < 2) return 0.0;
    size_t nx = grid->nx, ny = grid->ny;
    double dx = grid->dx, dy = grid->dy;
    if (dx <= 0.0 || dy <= 0.0) return 0.0;
    double sum_err = 0.0;
    size_t count = 0;
    for (size_t j = 1; j < ny - 1; j++) {
        for (size_t i = 1; i < nx - 1; i++) {
            size_t idx_r = j * nx + (i+1);
            size_t idx_l = j * nx + (i-1);
            size_t idx_u = (j+1) * nx + i;
            size_t idx_d = (j-1) * nx + i;
            double dBx_dx = (grid->fields[idx_r].Bx - grid->fields[idx_l].Bx) / (2.0 * dx);
            double dBy_dy = (grid->fields[idx_u].By - grid->fields[idx_d].By) / (2.0 * dy);
            sum_err += fabs(dBx_dx + dBy_dy);
            count++;
        }
    }
    return (count > 0) ? sum_err / (double)count : 0.0;
}

void mhd_wave_speeds(const mhd_primitive_t *prim, mhd_eigensystem_t *eigs, double gamma) {
    if (!prim || !eigs) return;
    double rho = prim->rho;
    if (rho <= 0.0) { memset(eigs, 0, sizeof(*eigs)); return; }
    double vx = prim->vx;
    double Bx = prim->Bx, By = prim->By, Bz = prim->Bz;
    double p  = prim->p;
    double B2 = Bx*Bx + By*By + Bz*Bz;
    double cs2 = gamma * p / rho;
    double vA2 = B2 / (SP_MU0 * rho);
    double vAx2 = Bx*Bx / (SP_MU0 * rho);
    double v_Ax = sqrt(vAx2);
    double sum_sq = cs2 + vA2;
    double disc = sum_sq*sum_sq - 4.0 * cs2 * vAx2;
    if (disc < 0.0) disc = 0.0;
    double vf2 = 0.5 * (sum_sq + sqrt(disc));
    double vs2 = 0.5 * (sum_sq - sqrt(disc));
    if (vs2 < 0.0) vs2 = 0.0;
    double v_fast = sqrt(vf2);
    double v_slow = sqrt(vs2);
    eigs->lambda[0] = vx - v_fast;
    eigs->lambda[1] = vx - v_Ax;
    eigs->lambda[2] = vx - v_slow;
    eigs->lambda[3] = vx;
    eigs->lambda[4] = vx + v_slow;
    eigs->lambda[5] = vx + v_Ax;
    eigs->lambda[6] = vx + v_fast;
    eigs->v_entropy = vx;
    eigs->v_alfven = v_Ax;
    eigs->v_slow   = v_slow;
    eigs->v_fast   = v_fast;
}

void mhd_force_free_field(double alpha, double B0, double x, double y, double B[3]) {
    if (!B) return;
    double k = alpha * 0.5;
    double a0_sq = alpha*alpha - k*k;
    double a0 = (a0_sq > 0.0) ? sqrt(a0_sq) : alpha;
    B[0] = -B0 * (k / a0) * sin(a0 * x) * sin(k * y);
    B[1] =  B0 * cos(a0 * x) * cos(k * y);
    B[2] =  B0 * sin(a0 * x) * cos(k * y);
}

void mhd_harris_sheet(double y, double delta, double B0, double n0,
                      double n_bg, double T, double B_out[3],
                      double *n_out, double *p_out) {
    if (!B_out) return;
    double tanh_y = tanh(y / delta);
    double sech2_y = 1.0 - tanh_y * tanh_y;
    B_out[0] = B0 * tanh_y;
    B_out[1] = 0.0;
    B_out[2] = 0.0;
    double p0 = B0 * B0 / (2.0 * SP_MU0);
    if (n_out) *n_out = n0 * sech2_y + n_bg;
    if (p_out) *p_out = p0 * sech2_y + n_bg * SP_KB * T;
}

void mhd_zpinch(double r, double a, double I, double T, double *p_out, double B_out[3]) {
    (void)T;  /* Available for Bennett relation: I^2 = 8*pi*N*kB*T/mu0 */
    if (!p_out && !B_out) return;
    double rc = (r < 0.0) ? 0.0 : r;
    if (B_out) {
        B_out[0] = 0.0;
        if (rc <= a && a > 0.0)
            B_out[1] = SP_MU0 * I * rc / (2.0 * M_PI * a * a);
        else if (a > 0.0)
            B_out[1] = SP_MU0 * I / (2.0 * M_PI * rc);
        else
            B_out[1] = 0.0;
        B_out[2] = 0.0;
    }
    if (p_out) {
        if (rc <= a && a > 0.0) {
            double ratio = rc / a;
            *p_out = (SP_MU0 * I * I / (4.0 * M_PI * M_PI * a * a)) * (1.0 - ratio * ratio);
            if (*p_out < 0.0) *p_out = 0.0;
        } else *p_out = 0.0;
    }
}

void mhd_grad_shafranov(const double *psi, size_t nx, size_t ny,
                        const double *R, double *rhs) {
    if (!psi || !R || !rhs || nx < 3 || ny < 3) return;
    double dR = (R[nx-1] - R[0]) / (double)(nx - 1);
    double dZ = dR;
    for (size_t j = 1; j < ny - 1; j++) {
        for (size_t i = 1; i < nx - 1; i++) {
            size_t idx = j * nx + i;
            double Ri = R[i];
            if (Ri <= 0.0) { rhs[idx] = 0.0; continue; }
            double R_inv_half_r = 2.0 / (R[i+1] + Ri);
            double R_inv_half_l = 2.0 / (Ri + R[i-1]);
            double dpsi_dR_r = (psi[j*nx + (i+1)] - psi[idx]) / dR;
            double dpsi_dR_l = (psi[idx] - psi[j*nx + (i-1)]) / dR;
            double d2psi_dZ2 = (psi[(j+1)*nx + i] - 2.0*psi[idx] + psi[(j-1)*nx + i]) / (dZ*dZ);
            double Rinv_term = (R_inv_half_r * dpsi_dR_r - R_inv_half_l * dpsi_dR_l) / dR;
            rhs[idx] = Ri * Rinv_term + d2psi_dZ2;
        }
    }
    for (size_t j = 0; j < ny; j++) {
        rhs[j*nx] = 0.0;
        rhs[j*nx + nx-1] = 0.0;
    }
    for (size_t i = 0; i < nx; i++) {
        rhs[i] = 0.0;
        rhs[(ny-1)*nx + i] = 0.0;
    }
}

double mhd_frozen_in_error(const double E[3], const double v[3], const double B[3]) {
    if (!E || !v || !B) return 1.0;
    double vxB[3];
    SP_CROSS3(v, B, vxB);
    double mag_err = sqrt(pow(E[0]+vxB[0],2) + pow(E[1]+vxB[1],2) + pow(E[2]+vxB[2],2));
    double mag_E = SP_MAG3(E);
    double mag_vxB = SP_MAG3(vxB);
    double mag_ref = (mag_E > mag_vxB) ? mag_E : mag_vxB;
    if (mag_ref < 1e-30) mag_ref = 1e-30;
    return mag_err / mag_ref;
}

void mhd_lax_friedrichs_step(const mhd_conserved_t *u_old, mhd_conserved_t *u_new,
                             size_t N, double dx, double dt, double gamma) {
    if (!u_old || !u_new || N < 3 || dx <= 0.0 || dt <= 0.0) return;
    double dtdx = dt / (2.0 * dx);
    for (size_t i = 1; i < N - 1; i++) {
        mhd_primitive_t prim_L, prim_R;
        mhd_conserved_t cons_L = u_old[i-1];
        mhd_conserved_t cons_R = u_old[i+1];
        if (mhd_cons_to_prim(&cons_L, &prim_L, gamma) < 0) {
            memcpy(&u_new[i], &u_old[i], sizeof(mhd_conserved_t));
            continue;
        }
        if (mhd_cons_to_prim(&cons_R, &prim_R, gamma) < 0) {
            memcpy(&u_new[i], &u_old[i], sizeof(mhd_conserved_t));
            continue;
        }
        mhd_flux_t flux_L, flux_R;
        mhd_flux_x(&prim_L, &flux_L, gamma);
        mhd_flux_x(&prim_R, &flux_R, gamma);
        u_new[i].D  = 0.5*(u_old[i-1].D+u_old[i+1].D)   - dtdx*(flux_R.f[0]-flux_L.f[0]);
        u_new[i].Mx = 0.5*(u_old[i-1].Mx+u_old[i+1].Mx) - dtdx*(flux_R.f[1]-flux_L.f[1]);
        u_new[i].My = 0.5*(u_old[i-1].My+u_old[i+1].My) - dtdx*(flux_R.f[2]-flux_L.f[2]);
        u_new[i].Mz = 0.5*(u_old[i-1].Mz+u_old[i+1].Mz) - dtdx*(flux_R.f[3]-flux_L.f[3]);
        u_new[i].Bx = 0.5*(u_old[i-1].Bx+u_old[i+1].Bx);
        u_new[i].By = 0.5*(u_old[i-1].By+u_old[i+1].By) - dtdx*(flux_R.f[5]-flux_L.f[5]);
        u_new[i].Bz = 0.5*(u_old[i-1].Bz+u_old[i+1].Bz) - dtdx*(flux_R.f[6]-flux_L.f[6]);
        u_new[i].E  = 0.5*(u_old[i-1].E+u_old[i+1].E)   - dtdx*(flux_R.f[7]-flux_L.f[7]);
    }
    u_new[0] = u_old[1];
    u_new[N-1] = u_old[N-2];
}

double mhd_cfl_timestep(const mhd_grid_t *grid, double cfl, double gamma) {
    if (!grid || !grid->fields || grid->nx == 0 || grid->ny == 0) return 0.0;
    double dx = grid->dx, dy = grid->dy;
    double lambda_max = 0.0;
    size_t N = grid->nx * grid->ny;
    for (size_t i = 0; i < N; i++) {
        mhd_primitive_t *p = &grid->fields[i];
        double rho = p->rho;
        if (rho <= 0.0) continue;
        double v2 = p->vx*p->vx + p->vy*p->vy + p->vz*p->vz;
        double B2 = p->Bx*p->Bx + p->By*p->By + p->Bz*p->Bz;
        double cs2 = gamma * p->p / rho;
        double vA2 = B2 / (SP_MU0 * rho);
        double sum_sq = cs2 + vA2;
        double disc = sum_sq*sum_sq - 4.0*cs2*vA2;
        if (disc < 0.0) disc = 0.0;
        double v_fast2 = 0.5 * (sum_sq + sqrt(disc));
        double v_fast = sqrt(v_fast2);
        double lambda_cell = sqrt(v2) + v_fast;
        if (lambda_cell > lambda_max) lambda_max = lambda_cell;
    }
    if (lambda_max <= 0.0) return INFINITY;
    double dl_min = (dx < dy) ? dx : dy;
    return cfl * dl_min / lambda_max;
}

double mhd_total_pressure(const mhd_primitive_t *prim) {
    if (!prim) return 0.0;
    double B2 = prim->Bx*prim->Bx + prim->By*prim->By + prim->Bz*prim->Bz;
    return prim->p + 0.5 * B2 / SP_MU0;
}

void mhd_energy_flux(const mhd_primitive_t *prim, double flux[3], double gamma) {
    if (!prim || !flux) return;
    double v[3] = {prim->vx, prim->vy, prim->vz};
    double B[3] = {prim->Bx, prim->By, prim->Bz};
    double B2 = SP_DOT3(B, B);
    double Bdotv = SP_DOT3(B, v);
    double E_tot = mhd_energy_density(prim, gamma);
    double coef = E_tot + prim->p + 0.5 * B2 / SP_MU0;
    flux[0] = coef * v[0] - B[0] * Bdotv / SP_MU0;
    flux[1] = coef * v[1] - B[1] * Bdotv / SP_MU0;
    flux[2] = coef * v[2] - B[2] * Bdotv / SP_MU0;
}

void mhd_induction_term(const double v[3], const double B[3], double curl[3]) {
    if (!v || !B || !curl) return;
    double vxB[3];
    SP_CROSS3(v, B, vxB);
    curl[0] = vxB[0];
    curl[1] = vxB[1];
    curl[2] = vxB[2];
}

void mhd_stress_tensor(const mhd_primitive_t *prim, double T[3][3]) {
    if (!prim || !T) return;
    double rho = prim->rho;
    double pp  = prim->p;
    double B[3] = {prim->Bx, prim->By, prim->Bz};
    double v[3] = {prim->vx, prim->vy, prim->vz};
    double B2 = SP_DOT3(B, B);
    double ptot = pp + 0.5 * B2 / SP_MU0;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            double delta_ij = (i == j) ? 1.0 : 0.0;
            T[i][j] = rho * v[i] * v[j] + ptot * delta_ij - B[i] * B[j] / SP_MU0;
        }
    }
}
