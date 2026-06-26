/*
 * mhd_numerics.h -- Numerical Methods for MHD Simulations
 *
 * Reference: Toro "Riemann Solvers" (2009) Ch.15-16
 *            Balsara "MHD Riemann Solvers" (1998-2016)
 *            Powell et al. "An Upwind Scheme for MHD" (1999)
 *            Dedner et al. "Div(B) Cleaning" (2002)
 *
 * Knowledge: L5 -- Numerical methods for hyperbolic conservation laws
 *             L5 -- Finite volume, constrained transport, div(B) control
 */

#ifndef MHD_NUMERICS_H
#define MHD_NUMERICS_H

#include "mhd_defs.h"

/* ================================================================
 * L5 -- Finite Difference / Finite Volume Operators
 * ================================================================ */

/*
 * mhd_fd_gradient_2nd -- 2nd-order central difference gradient
 * grad(f)_i = (f_{i+1} - f_{i-1}) / (2*dx) + O(dx^2)
 *
 * Complexity: O(1) per point
 */
double mhd_fd_gradient_2nd(const double *f, int i, double dx);

/*
 * mhd_fd_gradient_4th -- 4th-order central difference gradient
 * grad(f)_i = (-f_{i+2} + 8*f_{i+1} - 8*f_{i-1} + f_{i-2}) / (12*dx) + O(dx^4)
 *
 * Complexity: O(1) per point
 */
double mhd_fd_gradient_4th(const double *f, int i, double dx);

/*
 * mhd_fd_laplacian_2nd -- 2nd-order Laplacian
 * Lap(f)_i = (f_{i+1} - 2*f_i + f_{i-1}) / dx^2 + O(dx^2)
 *
 * Complexity: O(1) per point
 */
double mhd_fd_laplacian_2nd(const double *f, int i, double dx);

/*
 * mhd_fd_divergence_2nd -- 2nd-order divergence
 * div(Fx,Fy,Fz) at cell (i,j,k)
 *
 * Complexity: O(1)
 */
double mhd_fd_divergence_2nd(const double *Fx, const double *Fy, const double *Fz,
                              int i, int j, int k, int nx, int ny,
                              double dx, double dy, double dz);

/*
 * mhd_fd_curl_2nd -- 2nd-order curl
 * curl(B) at cell (i,j,k)
 *
 * Complexity: O(1)
 */
void mhd_fd_curl_2nd(const double *Bx, const double *By, const double *Bz,
                      int i, int j, int k, int nx, int ny,
                      double dx, double dy, double dz,
                      double *curl_x, double *curl_y, double *curl_z);

/* ================================================================
 * L5 -- Upwinding and Slope Limiters
 *
 * For MUSCL-Hancock / TVD schemes. Limiter phi(r) based on ratio of
 * consecutive gradients r = (u_i - u_{i-1}) / (u_{i+1} - u_i).
 * ================================================================ */

double mhd_limiter_minmod(double r);
double mhd_limiter_superbee(double r);
double mhd_limiter_van_leer(double r);
double mhd_limiter_mc(double r);        /* monotonized central */

/*
 * mhd_muscl_extrapolate -- MUSCL-Hancock reconstruction
 * Reconstructs left/right states at cell interface i+1/2.
 *
 * Complexity: O(1) per interface
 */
void mhd_muscl_extrapolate(const double *u, int i, double dx, double dt,
                            double (*limiter)(double),
                            double *uL, double *uR);

/*
 * mhd_muscl_extrapolate_mhd -- MUSCL for full MHD state
 * Reconstructs left/right MHDState at cell interface.
 *
 * Complexity: O(1) per interface
 */
void mhd_muscl_extrapolate_mhd(const MHDState *states, int i,
                                double dx, double dt,
                                double (*limiter)(double),
                                MHDState *left, MHDState *right);

/* ================================================================
 * L5 -- HLLD Riemann Solver for MHD (Miyoshi & Kusano 2005)
 *
 * The HLLD solver resolves 5 intermediate states for the MHD system:
 *   S_L -> S*_L -> S_M -> S*_R -> S_R
 * where S_M = v*_x is the contact/Alfven wave speed.
 *
 * This is the workhorse Riemann solver for modern MHD codes.
 * ================================================================ */

/*
 * mhd_hllc_simple -- Simplified 1D MHD flux (single-direction)
 *
 * For quick prototyping. Not as accurate as HLLD but faster.
 * Uses simple averaging plus dissipation term.
 *
 * Complexity: O(1)
 */
MHDFlux mhd_hllc_simple(const MHDState *left, const MHDState *right,
                         double gamma);

/*
 * mhd_hlld_flux -- HLLD Riemann solver flux
 *
 * Returns the numerical flux at the interface between left and right states.
 * Handles all 7 MHD waves (fast x2, Alfven x2, slow x2, entropy).
 *
 * Complexity: O(1)
 */
MHDFlux mhd_hlld_flux(const MHDState *left, const MHDState *right,
                       double gamma);

/*
 * mhd_roe_flux -- Roe-type approximate Riemann solver for MHD
 *
 * Uses the Roe-averaged state to compute the flux.
 * Requires eigenvector decomposition of the MHD system.
 *
 * Complexity: O(1)
 */
MHDFlux mhd_roe_flux(const MHDState *left, const MHDState *right,
                      double gamma);

/* ================================================================
 * L5 -- Time Integration
 * ================================================================ */

/*
 * mhd_rk2_step -- 2nd-order Runge-Kutta (Heun's method)
 *
 * U* = U^n + dt * R(U^n)
 * U^{n+1} = 0.5 * U^n + 0.5 * (U* + dt * R(U*))
 *
 * Complexity: O(n) for n cells
 */
void mhd_rk2_step(MHDConserved *U, int n, double dt,
                   void (*rhs)(const MHDConserved*, MHDConserved*, int, void*),
                   void *ctx);

/*
 * mhd_rk3_step -- 3rd-order TVD Runge-Kutta (Shu & Osher 1988)
 *
 * U^(1) = U^n + dt * R(U^n)
 * U^(2) = 3/4 U^n + 1/4 (U^(1) + dt * R(U^(1)))
 * U^{n+1} = 1/3 U^n + 2/3 (U^(2) + dt * R(U^(2)))
 *
 * TVD property: monotonicity preserving for CFL <= 1.
 *
 * Complexity: O(n) per substep, 3 substeps
 */
void mhd_rk3_step(MHDConserved *U, int n, double dt,
                   void (*rhs)(const MHDConserved*, MHDConserved*, int, void*),
                   void *ctx);

/*
 * mhd_rk4_step -- 4th-order classical Runge-Kutta
 *
 * Complexity: O(n) per substep, 4 substeps
 */
void mhd_rk4_step(MHDConserved *U, int n, double dt,
                   void (*rhs)(const MHDConserved*, MHDConserved*, int, void*),
                   void *ctx);

/* ================================================================
 * L5 -- Constrained Transport (Evans & Hawley 1988)
 *
 * Preserves div(B) = 0 to machine precision by updating B
 * via face-centered electric fields (Stokes theorem on each face):
 *   d_t(Phi_face) = EMF along edges
 *
 * The EMF (electromotive force) = -(v x B) at cell edges.
 * ================================================================ */

/*
 * mhd_ct_electric_field -- Compute edge-centered EMF for constrained transport
 *
 * Uses arithmetic averaging of face-centered fluxes.
 * E_z_{i+1/2, j+1/2} = -0.25*(vx_{i+1/2,j}*By + ...)
 *
 * Complexity: O(1) per edge
 */
void mhd_ct_emf_2d(double vx_face, double vy_face,
                    double Bx_corner, double By_corner,
                    double *Ez_edge);

/*
 * mhd_ct_update_B -- Update B using constrained transport
 *
 * dBx_{i+1/2,j}/dt = -(Ez_{i+1/2,j+1/2} - Ez_{i+1/2,j-1/2}) / dy
 * dBy_{i,j+1/2}/dt =  (Ez_{i+1/2,j+1/2} - Ez_{i-1/2,j+1/2}) / dx
 *
 * Complexity: O(1) per face
 */
void mhd_ct_update_B_2d(double Ez_pp, double Ez_pm, double Ez_mp, double Ez_mm,
                         double dx, double dy,
                         double *dBx, double *dBy);

/* ================================================================
 * L5 -- Powell 8-Wave Formulation
 *
 * The standard 8-wave MHD system includes a source term proportional
 * to div(B) to make the system Galilean invariant and ensure
 * div(B) advection (not just zero div(B) at t=0).
 *
 * Source term: S = -div(B) * (0, Bx, By, Bz, v.B, vx, vy, vz)
 * added to the right-hand side.
 * ================================================================ */

/*
 * mhd_powell_source -- Powell 8-wave source term
 *
 * Computes the source term to add to the MHD equations
 * when div(B) is not exactly zero.
 *
 * Complexity: O(1)
 */
void mhd_powell_source(const MHDState *state,
                        double divB,
                        double *S_D, double *S_Mx, double *S_My, double *S_Mz,
                        double *S_Bx, double *S_By, double *S_Bz,
                        double *S_E);

/* ================================================================
 * L5 -- Boundary Conditions
 * ================================================================ */

typedef enum {
    MHD_BC_PERIODIC = 0,
    MHD_BC_OUTFLOW,
    MHD_BC_REFLECTING,
    MHD_BC_CONDUCTING_WALL,
    MHD_BC_FIXED
} MHDBCType;

void mhd_apply_boundary_1d(MHDState *states, int n, MHDBCType left, MHDBCType right);
void mhd_apply_boundary_2d(MHDState *states, int nx, int ny,
                            MHDBCType xl, MHDBCType xr,
                            MHDBCType yl, MHDBCType yr);

/*
 * mhd_cfl_timestep -- Compute CFL-limited timestep
 * dt = CFL * min(dx,dy,dz) / max(v_n + c_f)
 *
 * Complexity: O(n)
 */
double mhd_cfl_timestep(const MHDState *states, int n,
                         double dx, double dy, double dz,
                         double gamma, double cfl);

#endif /* MHD_NUMERICS_H */
