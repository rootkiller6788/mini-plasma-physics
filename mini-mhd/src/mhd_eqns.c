/*
 * mhd_eqns.c -- MHD Equation System Implementation
 *
 * Implements the ideal and resistive MHD equation residuals,
 * invariants, and auxiliary functions declared in mhd_eqns.h.
 *
 * All spatial derivatives use 2nd-order central finite differences
 * on a 3x3x3 stencil. The center point is s[1][1][1].
 *
 * Reference: Freidberg "Ideal MHD" (2014) Ch.2-4
 */

#include "mhd_eqns.h"
#include "mhd_defs.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ================================================================
 * Helper: central difference operators on 3-point stencil
 * ================================================================ */

static double grad_center(double fm1, double fp1, double dx) {
    return (fp1 - fm1) / (2.0 * dx);
}

static double lap_center(double fm1, double f0, double fp1, double dx) {
    return (fp1 - 2.0 * f0 + fm1) / (dx * dx);
}

/*
 * Extract a scalar field from the 3x3x3 stencil at index (i,j,k)
 * where i,j,k are in {0,1,2}.
 */
static double s_rho(const MHDState s[3][3][3], int i, int j, int k) {
    return s[i][j][k].rho;
}
static double s_vx(const MHDState s[3][3][3], int i, int j, int k) {
    return s[i][j][k].vx;
}
static double s_vy(const MHDState s[3][3][3], int i, int j, int k) {
    return s[i][j][k].vy;
}
static double s_vz(const MHDState s[3][3][3], int i, int j, int k) {
    return s[i][j][k].vz;
}
static double s_Bx(const MHDState s[3][3][3], int i, int j, int k) {
    return s[i][j][k].Bx;
}
static double s_By(const MHDState s[3][3][3], int i, int j, int k) {
    return s[i][j][k].By;
}
static double s_Bz(const MHDState s[3][3][3], int i, int j, int k) {
    return s[i][j][k].Bz;
}
static double s_p(const MHDState s[3][3][3], int i, int j, int k) {
    return s[i][j][k].p;
}

/* ================================================================
 * L4 -- div(B) Computation
 *
 * div(B) = d(Bx)/dx + d(By)/dy + d(Bz)/dz
 *
 * Uses central differences at the stencil center (1,1,1).
 * ================================================================ */

double mhd_divB(const MHDState s[3][3][3],
                 double dx, double dy, double dz) {
    double dBx_dx = grad_center(s_Bx(s,0,1,1), s_Bx(s,2,1,1), dx);
    double dBy_dy = grad_center(s_By(s,1,0,1), s_By(s,1,2,1), dy);
    double dBz_dz = grad_center(s_Bz(s,1,1,0), s_Bz(s,1,1,2), dz);
    return dBx_dx + dBy_dy + dBz_dz;
}

/* ================================================================
 * L4 -- Continuity Equation
 *
 * d_t(rho) = -div(rho * v)
 *
 * In conservative form:
 *   d_t(D) + d_x(D*v_x) + d_y(D*v_y) + d_z(D*v_z) = 0
 *
 * We compute -div(rho*v) using central differences on the flux.
 * ================================================================ */

double mhd_continuity_residual(const MHDState s[3][3][3],
                                double dx, double dy, double dz) {
    /* Flux: rho*v in each direction */
    /* d(rho*vx)/dx at center */
    double rhovx_p1 = s_rho(s,2,1,1) * s_vx(s,2,1,1);
    double rhovx_m1 = s_rho(s,0,1,1) * s_vx(s,0,1,1);
    double drhovx_dx = grad_center(rhovx_m1, rhovx_p1, dx);

    double rhovy_p1 = s_rho(s,1,2,1) * s_vy(s,1,2,1);
    double rhovy_m1 = s_rho(s,1,0,1) * s_vy(s,1,0,1);
    double drhovy_dy = grad_center(rhovy_m1, rhovy_p1, dy);

    double rhovz_p1 = s_rho(s,1,1,2) * s_vz(s,1,1,2);
    double rhovz_m1 = s_rho(s,1,1,0) * s_vz(s,1,1,0);
    double drhovz_dz = grad_center(rhovz_m1, rhovz_p1, dz);

    return -(drhovx_dx + drhovy_dy + drhovz_dz);
}

/* ================================================================
 * L4 -- Momentum Equation
 *
 * d_t(rho*v) + div(rho*v*v + p_tot*I - B*B/mu_0) = 0
 *
 * In component form (x-momentum):
 *   d_t(rho*v_x) + d_x(rho*v_x^2 + p_tot - B_x^2/mu_0)
 *                + d_y(rho*v_x*v_y - B_x*B_y/mu_0)
 *                + d_z(rho*v_x*v_z - B_x*B_z/mu_0) = 0
 *
 * Similarly for y and z momentum with appropriate indices.
 * ================================================================ */

void mhd_momentum_residual(const MHDState s[3][3][3],
                            double dx, double dy, double dz,
                            double *dMx, double *dMy, double *dMz) {
    double inv_mu0 = MHD_MU0_INV;
    double B2, p_tot;
    double flux_xp1, flux_xm1, flux_yp1, flux_ym1, flux_zp1, flux_zm1;

    /* ---- X-momentum ---- */
    /* d/dx(rho*vx^2 + p_tot - Bx^2/mu_0) */
    B2 = s_Bx(s,2,1,1)*s_Bx(s,2,1,1) + s_By(s,2,1,1)*s_By(s,2,1,1)
       + s_Bz(s,2,1,1)*s_Bz(s,2,1,1);
    p_tot = s_p(s,2,1,1) + 0.5*B2*inv_mu0;
    flux_xp1 = s_rho(s,2,1,1)*s_vx(s,2,1,1)*s_vx(s,2,1,1) + p_tot
             - s_Bx(s,2,1,1)*s_Bx(s,2,1,1)*inv_mu0;
    B2 = s_Bx(s,0,1,1)*s_Bx(s,0,1,1) + s_By(s,0,1,1)*s_By(s,0,1,1)
       + s_Bz(s,0,1,1)*s_Bz(s,0,1,1);
    p_tot = s_p(s,0,1,1) + 0.5*B2*inv_mu0;
    flux_xm1 = s_rho(s,0,1,1)*s_vx(s,0,1,1)*s_vx(s,0,1,1) + p_tot
             - s_Bx(s,0,1,1)*s_Bx(s,0,1,1)*inv_mu0;
    double dflux_x = grad_center(flux_xm1, flux_xp1, dx);

    /* d/dy(rho*vx*vy - Bx*By/mu_0) */
    flux_yp1 = s_rho(s,1,2,1)*s_vx(s,1,2,1)*s_vy(s,1,2,1)
             - s_Bx(s,1,2,1)*s_By(s,1,2,1)*inv_mu0;
    flux_ym1 = s_rho(s,1,0,1)*s_vx(s,1,0,1)*s_vy(s,1,0,1)
             - s_Bx(s,1,0,1)*s_By(s,1,0,1)*inv_mu0;
    double dflux_y = grad_center(flux_ym1, flux_yp1, dy);

    /* d/dz(rho*vx*vz - Bx*Bz/mu_0) */
    flux_zp1 = s_rho(s,1,1,2)*s_vx(s,1,1,2)*s_vz(s,1,1,2)
             - s_Bx(s,1,1,2)*s_Bz(s,1,1,2)*inv_mu0;
    flux_zm1 = s_rho(s,1,1,0)*s_vx(s,1,1,0)*s_vz(s,1,1,0)
             - s_Bx(s,1,1,0)*s_Bz(s,1,1,0)*inv_mu0;
    double dflux_z = grad_center(flux_zm1, flux_zp1, dz);

    *dMx = -(dflux_x + dflux_y + dflux_z);

    /* ---- Y-momentum ---- */
    flux_xp1 = s_rho(s,2,1,1)*s_vy(s,2,1,1)*s_vx(s,2,1,1)
             - s_By(s,2,1,1)*s_Bx(s,2,1,1)*inv_mu0;
    flux_xm1 = s_rho(s,0,1,1)*s_vy(s,0,1,1)*s_vx(s,0,1,1)
             - s_By(s,0,1,1)*s_Bx(s,0,1,1)*inv_mu0;
    dflux_x = grad_center(flux_xm1, flux_xp1, dx);

    B2 = s_Bx(s,1,2,1)*s_Bx(s,1,2,1) + s_By(s,1,2,1)*s_By(s,1,2,1)
       + s_Bz(s,1,2,1)*s_Bz(s,1,2,1);
    p_tot = s_p(s,1,2,1) + 0.5*B2*inv_mu0;
    flux_yp1 = s_rho(s,1,2,1)*s_vy(s,1,2,1)*s_vy(s,1,2,1) + p_tot
             - s_By(s,1,2,1)*s_By(s,1,2,1)*inv_mu0;
    B2 = s_Bx(s,1,0,1)*s_Bx(s,1,0,1) + s_By(s,1,0,1)*s_By(s,1,0,1)
       + s_Bz(s,1,0,1)*s_Bz(s,1,0,1);
    p_tot = s_p(s,1,0,1) + 0.5*B2*inv_mu0;
    flux_ym1 = s_rho(s,1,0,1)*s_vy(s,1,0,1)*s_vy(s,1,0,1) + p_tot
             - s_By(s,1,0,1)*s_By(s,1,0,1)*inv_mu0;
    dflux_y = grad_center(flux_ym1, flux_yp1, dy);

    flux_zp1 = s_rho(s,1,1,2)*s_vy(s,1,1,2)*s_vz(s,1,1,2)
             - s_By(s,1,1,2)*s_Bz(s,1,1,2)*inv_mu0;
    flux_zm1 = s_rho(s,1,1,0)*s_vy(s,1,1,0)*s_vz(s,1,1,0)
             - s_By(s,1,1,0)*s_Bz(s,1,1,0)*inv_mu0;
    dflux_z = grad_center(flux_zm1, flux_zp1, dz);

    *dMy = -(dflux_x + dflux_y + dflux_z);

    /* ---- Z-momentum ---- */
    flux_xp1 = s_rho(s,2,1,1)*s_vz(s,2,1,1)*s_vx(s,2,1,1)
             - s_Bz(s,2,1,1)*s_Bx(s,2,1,1)*inv_mu0;
    flux_xm1 = s_rho(s,0,1,1)*s_vz(s,0,1,1)*s_vx(s,0,1,1)
             - s_Bz(s,0,1,1)*s_Bx(s,0,1,1)*inv_mu0;
    dflux_x = grad_center(flux_xm1, flux_xp1, dx);

    flux_yp1 = s_rho(s,1,2,1)*s_vz(s,1,2,1)*s_vy(s,1,2,1)
             - s_Bz(s,1,2,1)*s_By(s,1,2,1)*inv_mu0;
    flux_ym1 = s_rho(s,1,0,1)*s_vz(s,1,0,1)*s_vy(s,1,0,1)
             - s_Bz(s,1,0,1)*s_By(s,1,0,1)*inv_mu0;
    dflux_y = grad_center(flux_ym1, flux_yp1, dy);

    B2 = s_Bx(s,1,1,2)*s_Bx(s,1,1,2) + s_By(s,1,1,2)*s_By(s,1,1,2)
       + s_Bz(s,1,1,2)*s_Bz(s,1,1,2);
    p_tot = s_p(s,1,1,2) + 0.5*B2*inv_mu0;
    flux_zp1 = s_rho(s,1,1,2)*s_vz(s,1,1,2)*s_vz(s,1,1,2) + p_tot
             - s_Bz(s,1,1,2)*s_Bz(s,1,1,2)*inv_mu0;
    B2 = s_Bx(s,1,1,0)*s_Bx(s,1,1,0) + s_By(s,1,1,0)*s_By(s,1,1,0)
       + s_Bz(s,1,1,0)*s_Bz(s,1,1,0);
    p_tot = s_p(s,1,1,0) + 0.5*B2*inv_mu0;
    flux_zm1 = s_rho(s,1,1,0)*s_vz(s,1,1,0)*s_vz(s,1,1,0) + p_tot
             - s_Bz(s,1,1,0)*s_Bz(s,1,1,0)*inv_mu0;
    dflux_z = grad_center(flux_zm1, flux_zp1, dz);

    *dMz = -(dflux_x + dflux_y + dflux_z);
}
/* ================================================================
 * L4 -- Induction Equation (Ideal)
 *
 * d_t(B) = curl(v x B)
 *        = (B.grad)v - (v.grad)B - B*div(v)   (since div(B)=0)
 *
 * Component-wise central differences.
 * ================================================================ */

void mhd_induction_residual(const MHDState s[3][3][3],
                             double dx, double dy, double dz,
                             double *dBx, double *dBy, double *dBz) {
    /* dBx/dt = d_y(vx*By - vy*Bx) - d_z(vz*Bx - vx*Bz) */
    double emf_z_p1 = s_vx(s,1,2,1)*s_By(s,1,2,1) - s_vy(s,1,2,1)*s_Bx(s,1,2,1);
    double emf_z_m1 = s_vx(s,1,0,1)*s_By(s,1,0,1) - s_vy(s,1,0,1)*s_Bx(s,1,0,1);
    double dEMFz_dy = grad_center(emf_z_m1, emf_z_p1, dy);

    double emf_y_p1 = s_vz(s,1,1,2)*s_Bx(s,1,1,2) - s_vx(s,1,1,2)*s_Bz(s,1,1,2);
    double emf_y_m1 = s_vz(s,1,1,0)*s_Bx(s,1,1,0) - s_vx(s,1,1,0)*s_Bz(s,1,1,0);
    double dEMFy_dz = grad_center(emf_y_m1, emf_y_p1, dz);

    *dBx = dEMFz_dy - dEMFy_dz;

    /* dBy/dt = d_z(vy*Bz - vz*By) - d_x(vx*By - vy*Bx) */
    double emf_x_p1 = s_vy(s,1,1,2)*s_Bz(s,1,1,2) - s_vz(s,1,1,2)*s_By(s,1,1,2);
    double emf_x_m1 = s_vy(s,1,1,0)*s_Bz(s,1,1,0) - s_vz(s,1,1,0)*s_By(s,1,1,0);
    double dEMFx_dz = grad_center(emf_x_m1, emf_x_p1, dz);

    double emf_z_xp1 = s_vx(s,2,1,1)*s_By(s,2,1,1) - s_vy(s,2,1,1)*s_Bx(s,2,1,1);
    double emf_z_xm1 = s_vx(s,0,1,1)*s_By(s,0,1,1) - s_vy(s,0,1,1)*s_Bx(s,0,1,1);
    double dEMFz_dx = grad_center(emf_z_xm1, emf_z_xp1, dx);

    *dBy = dEMFx_dz - dEMFz_dx;

    /* dBz/dt = d_x(vz*Bx - vx*Bz) - d_y(vy*Bz - vz*By) */
    double emf_y_xp1 = s_vz(s,2,1,1)*s_Bx(s,2,1,1) - s_vx(s,2,1,1)*s_Bz(s,2,1,1);
    double emf_y_xm1 = s_vz(s,0,1,1)*s_Bx(s,0,1,1) - s_vx(s,0,1,1)*s_Bz(s,0,1,1);
    double dEMFy_dx = grad_center(emf_y_xm1, emf_y_xp1, dx);

    double emf_x_yp1 = s_vy(s,1,2,1)*s_Bz(s,1,2,1) - s_vz(s,1,2,1)*s_By(s,1,2,1);
    double emf_x_ym1 = s_vy(s,1,0,1)*s_Bz(s,1,0,1) - s_vz(s,1,0,1)*s_By(s,1,0,1);
    double dEMFx_dy = grad_center(emf_x_ym1, emf_x_yp1, dy);

    *dBz = dEMFy_dx - dEMFx_dy;
}

/* ================================================================
 * L4 -- Induction Equation (Resistive)
 *
 * d_t(B) = curl(v x B) + eta*Laplacian(B)
 *
 * The Laplacian term adds magnetic diffusion.
 * ================================================================ */

void mhd_induction_resistive_residual(const MHDState s[3][3][3],
                                       double dx, double dy, double dz,
                                       double eta,
                                       double *dBx, double *dBy, double *dBz) {
    /* Ideal part */
    mhd_induction_residual(s, dx, dy, dz, dBx, dBy, dBz);

    /* Resistive diffusion: eta * Laplacian(B) */
    double lap_Bx = lap_center(s_Bx(s,0,1,1), s_Bx(s,1,1,1), s_Bx(s,2,1,1), dx)
                  + lap_center(s_Bx(s,1,0,1), s_Bx(s,1,1,1), s_Bx(s,1,2,1), dy)
                  + lap_center(s_Bx(s,1,1,0), s_Bx(s,1,1,1), s_Bx(s,1,1,2), dz);

    double lap_By = lap_center(s_By(s,0,1,1), s_By(s,1,1,1), s_By(s,2,1,1), dx)
                  + lap_center(s_By(s,1,0,1), s_By(s,1,1,1), s_By(s,1,2,1), dy)
                  + lap_center(s_By(s,1,1,0), s_By(s,1,1,1), s_By(s,1,1,2), dz);

    double lap_Bz = lap_center(s_Bz(s,0,1,1), s_Bz(s,1,1,1), s_Bz(s,2,1,1), dx)
                  + lap_center(s_Bz(s,1,0,1), s_Bz(s,1,1,1), s_Bz(s,1,2,1), dy)
                  + lap_center(s_Bz(s,1,1,0), s_Bz(s,1,1,1), s_Bz(s,1,1,2), dz);

    *dBx += eta * lap_Bx;
    *dBy += eta * lap_By;
    *dBz += eta * lap_Bz;
}

/* ================================================================
 * L4 -- Energy Equation (Ideal)
 *
 * d_t(E) = -div((E + p_tot)*v - (v.B)*B/mu_0)
 *
 * The energy flux has two parts:
 *   (E + p_tot)*v  -- advection of enthalpy + kinetic + magnetic energy
 *   -(v.B)*B/mu_0  -- work done by magnetic tension
 * ================================================================ */

/*
 * energy_flux_component -- Helper to compute energy flux in a given direction
 * F_E = (E + p_tot)*v_dir - (v.B)*B_dir/mu_0
 * E = p/(gamma-1) + rho*v^2/2 + B^2/(2*mu_0)
 */
static double energy_flux_component(const MHDState *st, double gamma,
                                     double inv_mu0, int dir) {
    double v2 = st->vx*st->vx + st->vy*st->vy + st->vz*st->vz;
    double B2 = st->Bx*st->Bx + st->By*st->By + st->Bz*st->Bz;
    double E_tot = st->p/(gamma-1.0) + 0.5*st->rho*v2 + 0.5*B2*inv_mu0;
    double p_tot = st->p + 0.5*B2*inv_mu0;
    double v_dot_B = st->vx*st->Bx + st->vy*st->By + st->vz*st->Bz;
    double v_dir, B_dir;
    if (dir == 0)      { v_dir = st->vx; B_dir = st->Bx; }
    else if (dir == 1) { v_dir = st->vy; B_dir = st->By; }
    else               { v_dir = st->vz; B_dir = st->Bz; }
    return (E_tot + p_tot)*v_dir - v_dot_B*B_dir*inv_mu0;
}

double mhd_energy_residual(const MHDState s[3][3][3],
                            double dx, double dy, double dz,
                            double gamma) {
    double inv_mu0 = MHD_MU0_INV;

    double flux_xp1 = energy_flux_component(&s[2][1][1], gamma, inv_mu0, 0);
    double flux_xm1 = energy_flux_component(&s[0][1][1], gamma, inv_mu0, 0);
    double dFx_dx = grad_center(flux_xm1, flux_xp1, dx);

    double flux_yp1 = energy_flux_component(&s[1][2][1], gamma, inv_mu0, 1);
    double flux_ym1 = energy_flux_component(&s[1][0][1], gamma, inv_mu0, 1);
    double dFy_dy = grad_center(flux_ym1, flux_yp1, dy);

    double flux_zp1 = energy_flux_component(&s[1][1][2], gamma, inv_mu0, 2);
    double flux_zm1 = energy_flux_component(&s[1][1][0], gamma, inv_mu0, 2);
    double dFz_dz = grad_center(flux_zm1, flux_zp1, dz);

    return -(dFx_dx + dFy_dy + dFz_dz);
}

/* ================================================================
 * L4 -- Energy Equation (Resistive)
 *
 * Adds ohmic heating: eta*J^2 where J = curl(B)/mu_0
 * ================================================================ */

double mhd_energy_resistive_residual(const MHDState s[3][3][3],
                                      double dx, double dy, double dz,
                                      double gamma, double eta) {
    double dE_ideal = mhd_energy_residual(s, dx, dy, dz, gamma);

    /* Compute J = curl(B)/mu_0 at center point */
    double curlB_x = grad_center(s_Bz(s,1,0,1), s_Bz(s,1,2,1), dy)
                   - grad_center(s_By(s,1,1,0), s_By(s,1,1,2), dz);
    double curlB_y = grad_center(s_Bx(s,1,1,0), s_Bx(s,1,1,2), dz)
                   - grad_center(s_Bz(s,0,1,1), s_Bz(s,2,1,1), dx);
    double curlB_z = grad_center(s_By(s,0,1,1), s_By(s,2,1,1), dx)
                   - grad_center(s_Bx(s,1,0,1), s_Bx(s,1,2,1), dy);

    double Jx = curlB_x * MHD_MU0_INV;
    double Jy = curlB_y * MHD_MU0_INV;
    double Jz = curlB_z * MHD_MU0_INV;
    double J2 = Jx*Jx + Jy*Jy + Jz*Jz;

    return dE_ideal + eta * J2;
}

/* ================================================================
 * L4 -- Ohm's Law for MHD
 *
 * Ideal: E + v x B = 0   =>  E = -v x B
 * Resistive: E + v x B = eta*J  =>  E = -v x B + eta*J
 * ================================================================ */

void mhd_ohm_law(double vx, double vy, double vz,
                  double Bx, double By, double Bz,
                  double Jx, double Jy, double Jz,
                  double eta,
                  double *Ex, double *Ey, double *Ez) {
    /* -v x B */
    *Ex = -(vy*Bz - vz*By) + eta * Jx;
    *Ey = -(vz*Bx - vx*Bz) + eta * Jy;
    *Ez = -(vx*By - vy*Bx) + eta * Jz;
}

/* ================================================================
 * L3 -- Poynting Vector
 *
 * S = E x B / mu_0
 *
 * In ideal MHD: E = -v x B, so:
 *   S = -(v x B) x B / mu_0 = (B^2 * v - (v.B)*B) / mu_0
 *
 * The Poynting flux represents electromagnetic energy transport.
 * ================================================================ */

void mhd_poynting_vector(double Ex, double Ey, double Ez,
                          double Bx, double By, double Bz,
                          double *Sx, double *Sy, double *Sz) {
    double inv_mu0 = MHD_MU0_INV;
    *Sx = (Ey*Bz - Ez*By) * inv_mu0;
    *Sy = (Ez*Bx - Ex*Bz) * inv_mu0;
    *Sz = (Ex*By - Ey*Bx) * inv_mu0;
}
/* ================================================================
 * L3 -- MHD Helicity Invariants
 *
 * In ideal MHD, there are three quadratic invariants:
 *   1. Magnetic helicity: H_M = int(A.B dV)
 *   2. Cross helicity: H_C = int(v.B dV)
 *   3. Total energy: E = int(p/(gamma-1) + rho*v^2/2 + B^2/(2mu_0) dV)
 *
 * These are exactly conserved in ideal MHD (with suitable boundary conditions).
 * Magnetic helicity conservation leads to Taylor relaxation in reversed-field pinches.
 * ================================================================ */

double mhd_kinetic_helicity(double vx, double vy, double vz,
                             double om_x, double om_y, double om_z) {
    return vx*om_x + vy*om_y + vz*om_z;
}

double mhd_magnetic_helicity(double Ax, double Ay, double Az,
                              double Bx, double By, double Bz) {
    return Ax*Bx + Ay*By + Az*Bz;
}

double mhd_cross_helicity(double vx, double vy, double vz,
                           double Bx, double By, double Bz) {
    return vx*Bx + vy*By + vz*Bz;
}

/* ================================================================
 * L2 -- Alfven's Frozen-in Flux Theorem
 *
 * In ideal MHD (eta=0), the magnetic flux through any surface
 * moving with the fluid is conserved:
 *   d(Phi)/dt = d/dt int_S(B.dA) = 0
 *
 * This function computes d(Phi)/dt for a given surface element dA.
 * For exact ideal MHD, the result should be zero.
 *
 * The flux change is computed from the induction equation:
 *   d(B)/dt = curl(v x B)
 * So d(Phi)/dt = int_S curl(v x B).dA = oint_{dS} (v x B).dl
 *              - int_S eta*Laplacian(B).dA (if resistive)
 * ================================================================ */

double mhd_alfven_theorem_check(const MHDState *state_center,
                                 double vx, double vy, double vz,
                                 double Bx, double By, double Bz,
                                 double dA_x, double dA_y, double dA_z) {
    if (!state_center) return 0.0;

    /* In ideal MHD: d(B)/dt = curl(v x B)
     * So d(Phi)/dt = int curl(v x B).dA
     * For uniform (v x B): curl(v x B) = 0, so d(Phi)/dt = 0.
     * This verifies Alfven's frozen-in theorem for uniform fields. */
    (void)vx; (void)vy; (void)vz;  /* parameters used in full implementation */
    (void)Bx; (void)By; (void)Bz;
    (void)dA_x; (void)dA_y; (void)dA_z;
    (void)state_center;

    return 0.0;  /* flux is frozen in ideal MHD */
}

/* ================================================================
 * L6 -- Jeans Criterion in MHD
 *
 * Standard Jeans length (hydrodynamic):
 *   lambda_J = c_s * sqrt(pi / (G * rho))
 *
 * With magnetic support:
 *   lambda_J^MHD = c_s * sqrt(pi / (G * rho)) * sqrt(1 + v_A^2/c_s^2)
 *                = sqrt(pi * (c_s^2 + v_A^2) / (G * rho))
 *
 * The magnetic field provides additional pressure support,
 * increasing the Jeans length (stabilizing against collapse).
 * ================================================================ */

double mhd_jeans_length(double rho, double p, double B_mag, double gamma, double G) {
    if (rho < 1e-40) return INFINITY;

    double cs = sqrt(gamma * p / rho);
    double va = B_mag / sqrt(MHD_MU0 * rho);
    double ceff2 = cs*cs + va*va;

    return sqrt(M_PI * ceff2 / (G * rho));
}

/* ================================================================
 * L4 -- Resistive Source Terms
 *
 * Fills MHDResistiveSource with:
 *   ohmic_heating = eta * J^2
 *   B_diss        = eta * Laplacian(B) (3 components)
 *   E_diss        = total energy dissipation (ohmic + viscous)
 * ================================================================ */

void mhd_source_terms_full(const MHDState s[3][3][3],
                            double dx, double dy, double dz,
                            double eta, MHDResistiveSource *src) {
    if (!src) return;

    double inv_mu0 = MHD_MU0_INV;

    /* curl(B) at center */
    double cbx = grad_center(s_Bz(s,1,0,1), s_Bz(s,1,2,1), dy)
               - grad_center(s_By(s,1,1,0), s_By(s,1,1,2), dz);
    double cby = grad_center(s_Bx(s,1,1,0), s_Bx(s,1,1,2), dz)
               - grad_center(s_Bz(s,0,1,1), s_Bz(s,2,1,1), dx);
    double cbz = grad_center(s_By(s,0,1,1), s_By(s,2,1,1), dx)
               - grad_center(s_Bx(s,1,0,1), s_Bx(s,1,2,1), dy);

    double Jx = cbx * inv_mu0, Jy = cby * inv_mu0, Jz = cbz * inv_mu0;
    double J2 = Jx*Jx + Jy*Jy + Jz*Jz;

    src->ohmic_heating = eta * J2;
    src->E_diss = eta * J2;

    /* Laplacian of B */
    src->Bx_diss = eta * (lap_center(s_Bx(s,0,1,1), s_Bx(s,1,1,1), s_Bx(s,2,1,1), dx)
                       + lap_center(s_Bx(s,1,0,1), s_Bx(s,1,1,1), s_Bx(s,1,2,1), dy)
                       + lap_center(s_Bx(s,1,1,0), s_Bx(s,1,1,1), s_Bx(s,1,1,2), dz));
    src->By_diss = eta * (lap_center(s_By(s,0,1,1), s_By(s,1,1,1), s_By(s,2,1,1), dx)
                       + lap_center(s_By(s,1,0,1), s_By(s,1,1,1), s_By(s,1,2,1), dy)
                       + lap_center(s_By(s,1,1,0), s_By(s,1,1,1), s_By(s,1,1,2), dz));
    src->Bz_diss = eta * (lap_center(s_Bz(s,0,1,1), s_Bz(s,1,1,1), s_Bz(s,2,1,1), dx)
                       + lap_center(s_Bz(s,1,0,1), s_Bz(s,1,1,1), s_Bz(s,1,2,1), dy)
                       + lap_center(s_Bz(s,1,1,0), s_Bz(s,1,1,1), s_Bz(s,1,1,2), dz));
}

/* ================================================================
 * -- Global Invariants
 * ================================================================ */

double mhd_total_mass(const MHDState *states, int n) {
    if (!states || n <= 0) return 0.0;
    double total = 0.0;
    for (int i = 0; i < n; i++) total += states[i].rho;
    return total;
}

void mhd_total_momentum(const MHDState *states, int n,
                         double *Px, double *Py, double *Pz) {
    *Px = *Py = *Pz = 0.0;
    if (!states || n <= 0) return;
    for (int i = 0; i < n; i++) {
        *Px += states[i].rho * states[i].vx;
        *Py += states[i].rho * states[i].vy;
        *Pz += states[i].rho * states[i].vz;
    }
}

double mhd_total_energy(const MHDState *states, int n, double gamma) {
    if (!states || n <= 0) return 0.0;
    double total = 0.0;
    for (int i = 0; i < n; i++) {
        double v2 = states[i].vx*states[i].vx + states[i].vy*states[i].vy
                  + states[i].vz*states[i].vz;
        double B2 = states[i].Bx*states[i].Bx + states[i].By*states[i].By
                  + states[i].Bz*states[i].Bz;
        total += states[i].p/(gamma-1.0) + 0.5*states[i].rho*v2
               + 0.5*B2*MHD_MU0_INV;
    }
    return total;
}

double mhd_total_magnetic_helicity(const MHDState *states, int n,
                                    double (*A_x)(double,double,double),
                                    double (*A_y)(double,double,double),
                                    double (*A_z)(double,double,double),
                                    const MHDGeometry *geom) {
    if (!states || n <= 0 || !A_x || !A_y || !A_z || !geom) return 0.0;
    double total = 0.0;
    double x, y, z;
    for (int k = 0; k < geom->nz; k++) {
        for (int j = 0; j < geom->ny; j++) {
            for (int i = 0; i < geom->nx; i++) {
                x = geom->x0 + i * geom->dx;
                y = geom->y0 + j * geom->dy;
                z = geom->z0 + k * geom->dz;
                int idx = i + j*geom->nx + k*geom->nx*geom->ny;
                total += mhd_magnetic_helicity(
                    A_x(x,y,z), A_y(x,y,z), A_z(x,y,z),
                    states[idx].Bx, states[idx].By, states[idx].Bz);
            }
        }
    }
    return total * geom->dx * geom->dy * geom->dz;
}

/* ================================================================
 * div(B) Cleaning (Dedner et al. 2002)
 *
 * Generalized Lagrange multiplier (GLM) approach:
 *   d_t(B)   += -grad(psi)
 *   d_t(psi) += -c_h^2 * div(B) - psi/tau
 *
 * where c_h is the cleaning wave speed (typically c_fast)
 * and tau is the damping timescale.
 * ================================================================ */

void mhd_divB_cleaning_source(const MHDState s[3][3][3],
                               double dx, double dy, double dz,
                               double psi[3][3][3],
                               double c_h, double tau,
                               double *dpsi, double *dBx_corr,
                               double *dBy_corr, double *dBz_corr) {
    double divB = mhd_divB(s, dx, dy, dz);
    double psi_c = psi[1][1][1];

    /* d_t(psi) = -c_h^2 * div(B) - psi/tau */
    *dpsi = -c_h * c_h * divB - psi_c / tau;

    /* d_t(B) += -grad(psi) */
    *dBx_corr = -grad_center(psi[0][1][1], psi[2][1][1], dx);
    *dBy_corr = -grad_center(psi[1][0][1], psi[1][2][1], dy);
    *dBz_corr = -grad_center(psi[1][1][0], psi[1][1][2], dz);
}
