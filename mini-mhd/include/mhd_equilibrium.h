/*
 * mhd_equilibrium.h -- MHD Equilibrium Configurations
 *
 * Reference: Freidberg "Ideal MHD" (2014) Ch.8
 *            Wesson "Tokamaks" (2011) Ch.3
 *            Goedbloed & Poedts "Principles of MHD" (2004) Ch.16-17
 *
 * Knowledge: L4 -- MHD force balance JxB = grad(p)
 *             L6 -- Canonical equilibria: Z-pinch, theta-pinch, screw pinch
 *             L8 -- Grad-Shafranov equation, Taylor relaxation
 *
 * MHD equilibrium: J x B = grad p, curl(B) = mu_0 J, div(B) = 0
 * Implies B.grad(p) = 0, J.grad(p) = 0 (pressure constant along B and J).
 */

#ifndef MHD_EQUILIBRIUM_H
#define MHD_EQUILIBRIUM_H

#include "mhd_defs.h"

/* ================================================================
 * L6 -- 1D Pinch Equilibria (cylindrical)
 *
 * Force balance in cylindrical geometry (r, phi, z):
 *   dp/dr = J_phi B_z - J_z B_phi
 *   mu_0 J_phi = -dB_z/dr
 *   mu_0 J_z = (1/r) d/dr(r B_phi)
 * ================================================================ */

/* Z-pinch: B_phi only, current in z-direction.
 * Bennett relation: 2 N k_B T = mu_0 I^2/(8 pi). */
void mhd_zpinch_equilibrium(double (*J_z)(double), double r_max, int nr,
                             double p_edge,
                             double *r, double *B_phi, double *p);

/* Theta-pinch: B_z only, azimuthal current.
 * Force balance: dp/dr = J_phi B_z, mu_0 J_phi = -dB_z/dr. */
void mhd_thetapinch_equilibrium(double (*B_z)(double), double r_max, int nr,
                                 double p_edge,
                                 double *r, double *Bz, double *p);

/* Screw pinch: general 1D cylindrical with both B_phi(r) and B_z(r).
 * Input safety factor q(r) = r B_z/(R B_phi). */
void mhd_screwpinch_equilibrium(double (*q)(double), double R0,
                                 double r_max, int nr,
                                 double p0, double B0,
                                 double *r, double *B_phi, double *B_z, double *p);

/* Bennett relation for uniform-T Z-pinch.
 * 2 N k_B T = (mu_0 I^2) / (8 pi).
 * Returns T given I and N. */
double mhd_bennett_temperature(double I_total, double N_line);

/* Bennett profile characteristic radius.
 * a = sqrt(8 pi N k_B T / (mu_0 I^2)). */
double mhd_bennett_radius(double I_total, double N_line, double T);

/* ================================================================
 * L8 -- Grad-Shafranov Equation (2D axisymmetric toroidal)
 *
 * In (R, phi, Z) coordinates:
 *   Delta* psi = -mu_0 R^2 dp/dpsi - F dF/dpsi
 * where Delta* = R d/dR((1/R) dpsi/dR) + d^2psi/dZ^2
 *
 * psi(R,Z): poloidal flux function
 * F(psi) = R B_phi: toroidal field function
 * p(psi): pressure profile
 * ================================================================ */

/* Apply Delta* operator at point (iR, iZ). */
double mhd_grad_shafranov_operator(const double **psi,
                                    int iR, int iZ,
                                    double dR, double dZ, double R);

/* Soloviev analytic solution: linear dp/dpsi and F dF/dpsi.
 * psi = psi_0 * [(R^2 - R_0^2)^2/(4 R_0^2) + kappa^2 Z^2 / R_0^2]. */
double mhd_soloviev_psi(double R, double Z, double R0, double kappa, double psi0);

/* B_R, B_phi, B_Z from Soloviev solution. */
void mhd_soloviev_field(double R, double Z, double R0, double kappa,
                         double psi0, double F0,
                         double *BR, double *Bphi, double *BZ);

/* Safety factor q = r B_phi/(R B_theta) (large aspect ratio). */
double mhd_safety_factor_q(double r, double R0, double B_phi, double B_theta);

/* q from flux functions: q = (F/R_0) * dV/dpsi / (2 pi). */
double mhd_safety_factor_from_flux(double psi, double F, double dV_dpsi, double R0);

/* ================================================================
 * L8 -- Force-Free Fields and Taylor Relaxation
 *
 * Force-free: J x B = 0 => curl(B) = alpha B (Beltrami field).
 * Taylor (1974): minimum energy at fixed helicity -> curl(B) = lambda B.
 * ================================================================ */

/* Taylor relaxed state in periodic cylinder.
 * B_z = B_0 J_0(lambda r), B_phi = B_0 J_1(lambda r).
 * Returns 0 on success, -1 if root finding fails. */
int mhd_taylor_relaxation(double a, double H0_target, double B0,
                           int nr, double *r,
                           double *lambda_out,
                           double *Bz, double *Bphi);

/* Force-free parameter alpha = (J.B)/B^2 = curl(B).B/(mu_0 B^2). */
static inline double mhd_force_free_alpha(double Jx, double Jy, double Jz,
                                           double Bx, double By, double Bz) {
    double B2 = Bx*Bx + By*By + Bz*Bz;
    if (B2 < 1e-40) return 0.0;
    return (Jx*Bx + Jy*By + Jz*Bz) / B2;
}

/* Beltrami condition check: returns ||curl(B) - alpha B|| / ||B||. */
double mhd_beltrami_check(double alpha,
                           double Bx, double By, double Bz,
                           double Jx, double Jy, double Jz);

#endif /* MHD_EQUILIBRIUM_H */
