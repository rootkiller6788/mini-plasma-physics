/*
 * mhd_eqns.h -- MHD Equation System Declarations
 * Reference: Freidberg "Ideal MHD" (2014) Ch.2-4, Biskamp "Nonlinear MHD" (1993)
 * Knowledge: L4 -- Fundamental MHD equations (continuity, momentum, induction, energy)
 */

#ifndef MHD_EQNS_H
#define MHD_EQNS_H

#include "mhd_defs.h"

/**********************************************************************
 * L4 -- The MHD Equation System
 *
 * Ideal MHD (8 equations for 8 unknowns: rho, vx,vy,vz, Bx,By,Bz, p):
 *
 *  (1) Continuity: d_t(rho) + div(rho v) = 0
 *  (2) Momentum:   d_t(rho v) + div(rho v v + p_tot I - B B/mu_0) = 0
 *  (3) Induction:  d_t(B) = curl(v x B)                  (ideal)
 *  (3')            d_t(B) = curl(v x B) + eta Lap(B)     (resistive)
 *  (4) Energy:     d_t(E) + div((E+p_tot)v - (v.B)B/mu_0) = 0    (ideal)
 *  (4')            d_t(E) + div(...) = eta J^2                   (resistive)
 *
 *  where p_tot = p + B^2/(2mu_0), E = p/(gamma-1) + rho v^2/2 + B^2/(2mu_0)
 *  Constraint: div(B) = 0 at all times.
 *
 * Note: d_t denotes partial derivative with respect to time.
 *       div and curl are spatial operators.
 *       Lap = Laplacian = div(grad).
 **********************************************************************/

/* ---- Ideal MHD Equation Residuals (3x3x3 stencil) ---- */

/* Continuity: d_t(rho) = -div(rho v). Returns d(rho)/dt. */
double mhd_continuity_residual(const MHDState s[3][3][3],
                                double dx, double dy, double dz);

/* Momentum: d_t(rho v) = -div(rho v v + p_tot I - B B/mu_0).
 * Returns d(rho v_x, v_y, v_z)/dt. */
void mhd_momentum_residual(const MHDState s[3][3][3],
                            double dx, double dy, double dz,
                            double *dMx, double *dMy, double *dMz);

/* Induction (ideal): d_t(B) = curl(v x B) = (B.grad)v - (v.grad)B - B div(v).
 * The identity curl(v x B) = (B.grad)v - (v.grad)B + v div(B) - B div(v)
 * simplifies with div(B)=0. */
void mhd_induction_residual(const MHDState s[3][3][3],
                             double dx, double dy, double dz,
                             double *dBx, double *dBy, double *dBz);

/* Induction (resistive): d_t(B) = curl(v x B) + eta Lap(B). */
void mhd_induction_resistive_residual(const MHDState s[3][3][3],
                                       double dx, double dy, double dz,
                                       double eta,
                                       double *dBx, double *dBy, double *dBz);

/* Energy (ideal): d_t(E) = -div((E + p_tot)v - (v.B)B/mu_0). Returns dE/dt. */
double mhd_energy_residual(const MHDState s[3][3][3],
                            double dx, double dy, double dz,
                            double gamma);

/* Energy (resistive): d_t(E) = -div(...) + eta J^2. */
double mhd_energy_resistive_residual(const MHDState s[3][3][3],
                                      double dx, double dy, double dz,
                                      double gamma, double eta);

/* ---- Derived quantities (inline) ---- */

/* Total pressure: p_tot = p + B^2/(2 mu_0). */
static inline double mhd_total_pressure(double p, double B_mag) {
    return p + (B_mag * B_mag) / (2.0 * MHD_MU0);
}

/* ---- div(B) computation ---- */

/* div(B) = dBx/dx + dBy/dy + dBz/dz. Central finite differences. */
double mhd_divB(const MHDState s[3][3][3],
                 double dx, double dy, double dz);

/* ---- div(B) cleaning (Dedner et al. 2002) ---- */
/* Generalized Lagrange multiplier: transports div(B) errors with wave speed c_h
 * and damps them with timescale tau. */
void mhd_divB_cleaning_source(const MHDState s[3][3][3],
                               double dx, double dy, double dz,
                               double psi[3][3][3],
                               double c_h, double tau,
                               double *dpsi, double *dBx_corr,
                               double *dBy_corr, double *dBz_corr);

/* ---- Ohm's law ---- */

/* Ideal: E + v x B = 0. Resistive: E + v x B = eta J. */
void mhd_ohm_law(double vx, double vy, double vz,
                  double Bx, double By, double Bz,
                  double Jx, double Jy, double Jz,
                  double eta,
                  double *Ex, double *Ey, double *Ez);

/* ---- Poynting vector ---- */
/* S = E x B / mu_0. In ideal MHD: S = -(v x B) x B / mu_0 = (B^2 v_perp)/mu_0. */
void mhd_poynting_vector(double Ex, double Ey, double Ez,
                          double Bx, double By, double Bz,
                          double *Sx, double *Sy, double *Sz);

/* ---- MHD invariants (ideal) ---- */

/* Kinetic helicity: H_K = v . omega, omega = curl(v). Invariant in ideal HD. */
double mhd_kinetic_helicity(double vx, double vy, double vz,
                             double om_x, double om_y, double om_z);

/* Magnetic helicity: H_M = A . B, B = curl(A). Invariant in ideal MHD.
 * Woltjer (1958): minimizes magnetic energy at fixed helicity.
 * Taylor (1974): relaxation theory for reversed-field pinches. */
double mhd_magnetic_helicity(double Ax, double Ay, double Az,
                              double Bx, double By, double Bz);

/* Cross helicity: H_C = v . B. Invariant in ideal MHD (barotropic).
 * Measures alignment of v and B. Important in MHD turbulence. */
double mhd_cross_helicity(double vx, double vy, double vz,
                           double Bx, double By, double Bz);

/* ---- Alfven's frozen-in flux theorem ---- */
/* d(Phi)/dt = 0 for any surface moving with the fluid in ideal MHD.
 * Returns d(Phi)/dt for a given surface element. Should be ~0 for ideal. */
double mhd_alfven_theorem_check(const MHDState *state_center,
                                 double vx, double vy, double vz,
                                 double Bx, double By, double Bz,
                                 double dA_x, double dA_y, double dA_z);

/* ---- Gravitational Jeans criterion (MHD) ---- */
/* Magnetic Jeans length: lambda_J = c_s sqrt(pi/(G rho)) sqrt(1 + v_A^2/c_s^2).
 * Magnetic fields provide additional support against gravitational collapse. */
double mhd_jeans_length(double rho, double p, double B_mag, double gamma, double G);

/* ---- Resistive source terms ---- */
/* Compute full resistive MHD source terms from a 3x3x3 stencil. */
void mhd_source_terms_full(const MHDState s[3][3][3],
                            double dx, double dy, double dz,
                            double eta, MHDResistiveSource *src);

/* ---- Global invariants ---- */
double mhd_total_mass(const MHDState *states, int n);
void   mhd_total_momentum(const MHDState *states, int n,
                           double *Px, double *Py, double *Pz);
double mhd_total_energy(const MHDState *states, int n, double gamma);
double mhd_total_magnetic_helicity(const MHDState *states, int n,
                                    double (*A_x)(double,double,double),
                                    double (*A_y)(double,double,double),
                                    double (*A_z)(double,double,double),
                                    const MHDGeometry *geom);

#endif /* MHD_EQNS_H */
