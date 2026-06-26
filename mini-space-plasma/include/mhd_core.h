/**
 * mhd_core.h — Magnetohydrodynamics Core
 *
 * References:
 *   Kivelson & Russell §4 (MHD equations)
 *   Goldston & Rutherford §9 (MHD equilibrium & stability)
 *   Priest "Solar Magnetohydrodynamics" (1982)
 *   MIT 22.615 MHD Theory
 *
 * Knowledge: L1 (definitions), L2 (MHD concepts), L3 (math structures),
 *   L4 (MHD equations), L5 (numerical MHD), L6 (MHD equilibria)
 */

#ifndef MHD_CORE_H
#define MHD_CORE_H

#include "space_plasma.h"

/*──────────────────────────────────────────────────────────────────────
 * L1: MHD State Definitions
 *──────────────────────────────────────────────────────────────────────*/

/** @brief 1D MHD primitive variable state */
typedef struct {
    double rho;   /* Mass density [kg/m³]          */
    double vx;    /* Velocity x [m/s]              */
    double vy;    /* Velocity y [m/s]              */
    double vz;    /* Velocity z [m/s]              */
    double Bx;    /* Magnetic field x [T]          */
    double By;    /* Magnetic field y [T]          */
    double Bz;    /* Magnetic field z [T]          */
    double p;     /* Pressure [Pa]                 */
    double E;     /* Total energy density [J/m³]   */
} mhd_primitive_t;

/** @brief 1D MHD conserved variable state */
typedef struct {
    double D;     /* Mass density = ρ              */
    double Mx;    /* Momentum x = ρ v_x            */
    double My;    /* Momentum y = ρ v_y            */
    double Mz;    /* Momentum z = ρ v_z            */
    double Bx;    /* B_x (no ρ factor)             */
    double By;    /* B_y                           */
    double Bz;    /* B_z                           */
    double E;     /* Total energy density          */
} mhd_conserved_t;

/** @brief MHD flux vector (for Riemann solvers) */
typedef struct {
    double f[8];  /* Flux of [D, Mx, My, Mz, By, Bz, E, ψ] */
} mhd_flux_t;

/** @brief MHD wave speeds (eigenvalues) */
typedef struct {
    double lambda[7]; /* 7 eigenvalues of 1D MHD system */
    double v_entropy; /* Entropy wave: λ = v_x          */
    double v_alfven;  /* Alfven waves: λ = v_x ± v_A_x   */
    double v_slow;    /* Slow magnetosonic               */
    double v_fast;    /* Fast magnetosonic               */
} mhd_eigensystem_t;

/** @brief 2D MHD grid (for equilibrium and simple simulations) */
typedef struct {
    size_t  nx;       /* Number of x grid points      */
    size_t  ny;       /* Number of y grid points      */
    double  dx;       /* Grid spacing x [m]           */
    double  dy;       /* Grid spacing y [m]           */
    double  Lx;       /* Domain size x [m]            */
    double  Ly;       /* Domain size y [m]            */
    mhd_primitive_t *fields; /* Flat array [nx*ny]    */
} mhd_grid_t;

/*──────────────────────────────────────────────────────────────────────
 * L1/L2: MHD conserved ↔ primitive variable conversion
 *──────────────────────────────────────────────────────────────────────*/

/**
 * @brief Convert primitive MHD variables to conserved form
 *
 * D = ρ
 * M = ρ v
 * B = B (unchanged)
 * E = p/(γ-1) + ½ρv² + B²/(2μ₀)
 *
 * @param prim  Primitive MHD state
 * @param cons  Output: conserved state
 * @param gamma Polytropic index (5/3 for adiabatic, 1 for isothermal)
 */
void mhd_prim_to_cons(const mhd_primitive_t *prim, mhd_conserved_t *cons, double gamma);

/**
 * @brief Convert conserved to primitive variables
 *
 * Requires solving for v = M/D and p from energy equation.
 *
 * @param cons  Conserved state
 * @param prim  Output: primitive state
 * @param gamma Polytropic index
 * @return      0 on success, -1 if non-physical (negative pressure)
 */
int mhd_cons_to_prim(const mhd_conserved_t *cons, mhd_primitive_t *prim, double gamma);

/**
 * @brief Compute total energy density from primitives
 */
double mhd_energy_density(const mhd_primitive_t *prim, double gamma);

/*──────────────────────────────────────────────────────────────────────
 * L2: MHD Flux Computation
 *──────────────────────────────────────────────────────────────────────*/

/**
 * @brief Compute MHD flux in x-direction
 *
 * Standard ideal MHD flux:
 *   F_D   = D v_x
 *   F_Mx  = D v_x² + p + B²/(2μ₀) - B_x²/μ₀
 *   F_My  = D v_x v_y - B_x B_y/μ₀
 *   F_Mz  = D v_x v_z - B_x B_z/μ₀
 *   F_By  = B_y v_x - B_x v_y
 *   F_Bz  = B_z v_x - B_x v_z
 *   F_E   = (E + p + B²/(2μ₀))v_x - B_x(B·v)/μ₀
 *
 * @param prim  Primitive state
 * @param flux  Output: flux vector
 * @param gamma Polytropic index
 */
void mhd_flux_x(const mhd_primitive_t *prim, mhd_flux_t *flux, double gamma);

/**
 * @brief Compute ∇·B error for a 2D grid
 * @return L1 norm of div(B) over grid
 */
double mhd_divB_error(const mhd_grid_t *grid);

/*──────────────────────────────────────────────────────────────────────
 * L2: MHD Wave Speeds (Eigenvalues)
 *──────────────────────────────────────────────────────────────────────*/

/**
 * @brief Compute all 7 MHD wave speeds in x-direction
 *
 * For waves propagating in direction making angle θ with B:
 *   v_entropy = v_x
 *   v_alfven  = v_x ± v_A cos θ
 *   v_slow²   = ½(c_s² + v_A² - sqrt((c_s²+v_A²)² - 4c_s²v_A²cos²θ))
 *   v_fast²   = ½(c_s² + v_A² + sqrt((c_s²+v_A²)² - 4c_s²v_A²cos²θ))
 *
 * @param prim  Primitive state
 * @param eigs  Output: eigenvalue structure
 * @param gamma Polytropic index
 */
void mhd_wave_speeds(const mhd_primitive_t *prim, mhd_eigensystem_t *eigs, double gamma);

/*──────────────────────────────────────────────────────────────────────
 * L3: MHD Equilibrium Solutions
 *──────────────────────────────────────────────────────────────────────*/

/**
 * @brief Compute 1D force-free field: J × B = 0
 *
 * For a linear force-free field: ∇×B = α B, α = constant.
 * Solution: B_x = -B₀/k * d/dy(sin(αx) cos(ky))
 *
 * @param alpha Force-free parameter [1/m]
 * @param B0    Reference field [T]
 * @param x, y  Coordinates [m]
 * @param B     Output: B[3] field vector [T]
 */
void mhd_force_free_field(double alpha, double B0, double x, double y, double B[3]);

/**
 * @brief Harris current sheet (1D equilibrium)
 *
 * B_x(y) = B₀ tanh(y/δ), B_y = 0
 * n(y)  = n₀ sech²(y/δ) + n_bg
 * Total pressure (thermal + magnetic) = constant.
 *
 * @param y          Coordinate perpendicular to sheet [m]
 * @param delta      Sheet half-thickness [m]
 * @param B0         Asymptotic field [T]
 * @param n0         Peak density [m⁻³]
 * @param n_bg       Background density [m⁻³]
 * @param T          Temperature [K]
 * @param B_out      Output: B field [T]
 * @param n_out      Output: density [m⁻³]
 * @param p_out      Output: pressure [Pa]
 */
void mhd_harris_sheet(double y, double delta, double B0, double n0,
                      double n_bg, double T, double B_out[3],
                      double *n_out, double *p_out);

/**
 * @brief Z-pinch equilibrium (1D cylindrical)
 *
 * Pressure balance: dp/dr = -J_z B_θ, with total current I.
 * Bennett relation: (I [A])² = 8π N k_B T / μ₀
 *
 * @param r       Radial coordinate [m]
 * @param a       Pinch radius [m]
 * @param I       Total current [A]
 * @param T       Temperature (constant) [K]
 * @param p_out   Output: pressure [Pa]
 * @param B_out   Output: B_θ(r) [T]
 */
void mhd_zpinch(double r, double a, double I, double T, double *p_out, double B_out[3]);

/**
 * @brief Grad-Shafranov operator for axisymmetric equilibrium
 *
 * Δ*ψ = R² μ₀ dp/dψ + F dF/dψ
 *
 * Simplified: assumes p = p(ψ) and F = constant.
 * Computes the elliptic operator for ψ on a 2D grid.
 *
 * @param psi   Flux function ψ on grid [nx*ny]
 * @param nx    Grid points in R
 * @param ny    Grid points in Z
 * @param R     R-coordinates [m] array [nx]
 * @param rhs   Output: Δ*ψ [T·m²]
 */
void mhd_grad_shafranov(const double *psi, size_t nx, size_t ny,
                        const double *R, double *rhs);

/*──────────────────────────────────────────────────────────────────────
 * L3: Frozen-in Flux Condition
 *──────────────────────────────────────────────────────────────────────*/

/**
 * @brief Check frozen-in flux condition
 *
 * E + v × B = η J ≈ 0 in ideal MHD.
 * Returns |E + v × B| / max(|E|, |v×B|) as normalized error.
 *
 * @param E  Electric field [V/m]
 * @param v  Velocity [m/s]
 * @param B  Magnetic field [T]
 * @return   Normalized frozen-in error (0 = ideal MHD)
 */
double mhd_frozen_in_error(const double E[3], const double v[3], const double B[3]);

/*──────────────────────────────────────────────────────────────────────
 * L5: Simple Lax-Friedrichs MHD solver (1D)
 *──────────────────────────────────────────────────────────────────────*/

/**
 * @brief Single Lax-Friedrichs step for 1D MHD
 *
 * U_i^{n+1} = (U_{i-1}ⁿ + U_{i+1}ⁿ)/2 - Δt/(2Δx)(F_{i+1}ⁿ - F_{i-1}ⁿ)
 *
 * @param u_old    Conserved state array [N cells]
 * @param u_new    Output: updated state [N cells]
 * @param N        Number of cells
 * @param dx       Cell size [m]
 * @param dt       Time step [s]
 * @param gamma    Polytropic index
 */
void mhd_lax_friedrichs_step(const mhd_conserved_t *u_old, mhd_conserved_t *u_new,
                             size_t N, double dx, double dt, double gamma);

/*──────────────────────────────────────────────────────────────────────
 * L5: MHD CFL Condition
 *──────────────────────────────────────────────────────────────────────*/

/**
 * @brief Compute maximum stable time step (CFL condition)
 *
 * Δt_max = CFL * min(Δx, Δy) / max(|v| + v_fast)
 *
 * @param grid  MHD grid
 * @param cfl   CFL number (typically 0.5-0.8)
 * @param gamma Polytropic index
 * @return      Maximum stable dt [s]
 */
double mhd_cfl_timestep(const mhd_grid_t *grid, double cfl, double gamma);

/**
 * @brief Compute MHD total pressure: p_total = p + B²/(2μ₀)
 */
double mhd_total_pressure(const mhd_primitive_t *prim);

/**
 * @brief Compute MHD energy flux (Poynting + enthalpy flux)
 * S_MHD = (E + p + B²/(2μ₀))v - B(B·v)/μ₀
 */
void mhd_energy_flux(const mhd_primitive_t *prim, double flux[3], double gamma);

/**
 * @brief Compute ∇×(v×B) induction term
 * @param v    Velocity field [3]
 * @param B    Magnetic field [3]
 * @param curl Output: ∇×(v×B) [3]
 */
void mhd_induction_term(const double v[3], const double B[3], double curl[3]);

/**
 * @brief Compute MHD stress tensor component
 * T_ij = ρ v_i v_j + (p + B²/(2μ₀))δ_ij - B_i B_j/μ₀
 */
void mhd_stress_tensor(const mhd_primitive_t *prim, double T[3][3]);

#endif /* MHD_CORE_H */
