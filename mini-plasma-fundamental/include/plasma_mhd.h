/**
 * plasma_mhd.h — Magnetohydrodynamics (MHD)
 *
 * Single-fluid and two-fluid MHD models for macroscopic plasma behavior.
 *
 * References:
 *   - Freidberg, "Ideal MHD" (2014)
 *   - Biskamp, "Nonlinear Magnetohydrodynamics" (1993)
 *   - Goedbloed & Poedts, "Principles of MHD" (2004)
 *   - Wesson, "Tokamaks" (2011)
 *
 * Knowledge Levels:
 *   L2 — Core Concepts: frozen-in flux, magnetic pressure/tension
 *   L3 — Math Structures: MHD equations, flux functions
 *   L4 — Fundamental Laws: ideal MHD, resistive MHD, Grad-Shafranov
 */

#ifndef PLASMA_MHD_H
#define PLASMA_MHD_H

#include "plasma_constants.h"

/* ================================================================
 * L2/L3: MHD State Vectors
 * ================================================================ */

/** 1D MHD state vector (conservative variables) */
typedef struct {
    double rho;      /**< Mass density [kg/m^3] */
    double mx;       /**< x-momentum density rho*vx [kg/(m^2*s)] */
    double my;       /**< y-momentum density rho*vy [kg/(m^2*s)] */
    double mz;       /**< z-momentum density rho*vz [kg/(m^2*s)] */
    double Bx;       /**< x-magnetic field [T] */
    double By;       /**< y-magnetic field [T] */
    double Bz;       /**< z-magnetic field [T] */
    double E;        /**< Total energy density [J/m^3] */
} MHDState1D;

/** 2D MHD state for Grad-Shafranov equilibrium */
typedef struct {
    int     nr;      /**< Number of radial grid points */
    int     nz;      /**< Number of axial grid points */
    double  r_min, r_max, dr;
    double  z_min, z_max, dz;
    double **psi;    /**< Poloidal flux function psi(r,z) */
    double **j_phi;  /**< Toroidal current density [A/m^2] */
    double **p;      /**< Plasma pressure [Pa] */
} MHDEquilibrium2D;

/* ================================================================
 * L2: MHD Wave Speeds
 * ================================================================ */

/**
 * Alfven speed: v_A = B / sqrt(mu_0 * rho) [m/s]
 * (Re-exported from plasma_params.h for MHD context)
 */
double alfven_speed_mhd(double B, double rho);

/**
 * Sound speed for adiabatic ideal gas:
 * c_s = sqrt(gamma * p / rho) [m/s]
 *
 * For a monatomic gas, gamma = 5/3.
 */
double sound_speed_mhd(double pressure, double rho, double gamma);

/**
 * Fast magnetosonic speed [m/s]
 *
 * c_f^2 = (1/2) * (c_s^2 + v_A^2 + sqrt((c_s^2+v_A^2)^2 - 4*c_s^2*v_A^2*cos^2(theta)))
 *
 * where theta is the angle between k and B.
 *
 * Fast mode: compressional wave, magnetic and thermal pressure in phase.
 */
double fast_magnetosonic_speed(double cs, double vA, double cos_theta);

/**
 * Slow magnetosonic speed [m/s]
 *
 * c_sl^2 = (1/2) * (c_s^2 + v_A^2 - sqrt((c_s^2+v_A^2)^2 - 4*c_s^2*v_A^2*cos^2(theta)))
 *
 * Slow mode: compressional wave, magnetic and thermal pressure out of phase.
 */
double slow_magnetosonic_speed(double cs, double vA, double cos_theta);

/**
 * Compute all three MHD wave speeds for given angle.
 *
 * @param cs Sound speed [m/s]
 * @param vA Alfven speed [m/s]
 * @param cos_theta cos(angle between k and B)
 * @param out_vA_wave Output: Alfven wave speed
 * @param out_vF Output: Fast magnetosonic speed
 * @param out_vS Output: Slow magnetosonic speed
 */
void mhd_wave_speeds(double cs, double vA, double cos_theta,
                     double *out_vA_wave, double *out_vF, double *out_vS);

/* ================================================================
 * L3: MHD Flux Functions and Divergence Cleaning
 * ================================================================ */

/**
 * Compute the ideal MHD flux vector F(U) in the x-direction.
 *
 * F(U) = [ rho*vx,
 *          rho*vx^2 + p_tot - Bx^2/mu_0,
 *          rho*vx*vy - Bx*By/mu_0,
 *          rho*vx*vz - Bx*Bz/mu_0,
 *          0,
 *          vx*By - vy*Bx,
 *          vx*Bz - vz*Bx,
 *          (E + p_tot)*vx - Bx*(v·B)/mu_0 ]
 *
 * where p_tot = p + B^2/(2*mu_0) is the total pressure.
 *
 * Complexity: O(1)
 */
void ideal_mhd_flux_x(const MHDState1D *U, double gamma, MHDState1D *F);

/**
 * Compute the divergence of B at cell i for a 1D array (for cleaning).
 *
 * divB_i = (Bx_{i+1} - Bx_{i-1}) / (2*dx)
 *
 * In ideal MHD, div B = 0 must be preserved to machine precision.
 *
 * Complexity: O(1) per cell
 */
double div_b_1d(const double *Bx, int i, double dx, int nx);

/**
 * Powell 8-wave divergence cleaning source term.
 *
 * Adds source terms proportional to div B to the MHD equations
 * to control div B errors in numerical simulations.
 *
 * Reference: Powell et al. (1999), J. Comp. Phys. 154, 284.
 */
void powell_source_term(MHDState1D *S, double divB, double Bx,
                        double By, double Bz, double vx, double vy, double vz);

/* ================================================================
 * L4: Grad-Shafranov Equation (Axisymmetric MHD Equilibrium)
 * ================================================================ */

/**
 * Evaluate the Grad-Shafranov operator on psi.
 *
 * Delta* psi = R^2 div((1/R^2) grad psi)
 *            = R * d/dR((1/R) dpsi/dR) + d^2psi/dZ^2
 *            = -mu_0 * R^2 * dp/dpsi - F * dF/dpsi
 *
 * where psi(R,Z) is the poloidal flux function, p(psi) is the
 * pressure profile, and F(psi) = R*B_phi is the toroidal field
 * function.
 *
 * Reference: Grad & Rubin (1958), Shafranov (1958)
 *
 * @param psi Poloidal flux at grid point
 * @param R Major radius coordinate [m]
 * @param dp_dpsi dp/dpsi at this psi value
 * @param F_dF_dpsi F * dF/dpsi at this psi value
 * @return The source term -mu_0*R^2*dp/dpsi - F*dF/dpsi
 */
double grad_shafranov_source(double psi, double R,
                             double dp_dpsi, double F_dF_dpsi);

/**
 * Compute the toroidal current density from Grad-Shafranov:
 *
 * j_phi = R * dp/dpsi + (1/(mu_0*R)) * F * dF/dpsi
 *
 * This is the current that generates the poloidal field.
 */
double toroidal_current_density(double psi, double R,
                                double dp_dpsi, double F_dF_dpsi);

/**
 * Pressure profile: p(psi) = p0 * (1 - psi/psi_axis)^alpha
 *
 * Common parameterization for tokamak equilibria.
 *
 * @param psi Poloidal flux (0 at axis, increasing toward edge)
 * @param psi_axis Flux at magnetic axis
 * @param psi_edge Flux at plasma boundary (separatrix)
 * @param p0 Central pressure [Pa]
 * @param alpha Peaking factor (typically 1.0–2.0)
 */
double pressure_profile_gs(double psi, double psi_axis, double psi_edge,
                           double p0, double alpha);

/**
 * F = R*B_phi profile: F^2(psi) = F0^2 * (1 - beta_pol * (psi/psi_edge)^2)
 *
 * @return F(psi) = R*B_phi (approximate constant for large aspect ratio)
 */
double f_profile_gs(double psi, double psi_edge, double F0, double beta_pol);

/**
 * Compute the safety factor q(psi) for axisymmetric equilibrium.
 *
 * q = (dPsi_tor/dPsi_pol)^(-1) = (1/2*pi) * integral (1/R) * (B_phi/B_pol) dl
 *
 * For large aspect ratio circular cross-section:
 * q ~ (r/R) * (B_phi/B_theta)
 *
 * q is fundamental to MHD stability: q < 1 is unstable to kink modes.
 *
 * Reference: Wesson, "Tokamaks" (2011), Ch. 3
 */
double safety_factor_cylindrical(double r, double R0, double B_phi,
                                 double B_theta);

/**
 * Kruskal-Shafranov limit for tokamak stability.
 *
 * q(a) > 1 for stability against external kink modes.
 * I_p < (2*pi*a^2*B_phi) / (mu_0*R0*q_edge)
 *
 * @return Maximum stable plasma current [A]
 */
double kruskal_shafranov_limit(double a, double R0, double B_phi, double q_edge);

/* ================================================================
 * L4: Resistive MHD and Magnetic Reconnection
 * ================================================================ */

/**
 * Magnetic Reynolds number (Lundquist number)
 *
 * S = mu_0 * L * v_A / eta
 *
 * S >> 1: ideal MHD (frozen-in flux)
 * S ~ 1: resistive MHD dominant
 *
 * Typical values:
 *   - Solar corona: S ~ 10^12
 *   - Tokamak:      S ~ 10^8
 *   - Laboratory:   S ~ 10^3 - 10^6
 */
double lundquist_number(double L, double vA, double eta);

/**
 * Sweet-Parker reconnection rate (2D steady-state)
 *
 * v_in / v_A = 1/sqrt(S)
 *
 * This is too slow to explain solar flare timescales — motivating
 * the Petschek model with X-point geometry.
 *
 * Reference: Sweet (1958), Parker (1957)
 */
double sweet_parker_rate(double S);

/**
 * Petschek reconnection rate (fast reconnection)
 *
 * v_in / v_A ~ pi / (8 * ln S)
 *
 * Nearly independent of S for large S, explaining fast energy release
 * in solar flares and magnetospheric substorms.
 *
 * Reference: Petschek (1964), AAS-NASA Symposium on Solar Flares
 */
double petschek_rate(double S);

/**
 * Compute the tearing mode stability parameter Delta'.
 *
 * Delta' = lim_{eps->0} [psi'(r_s+eps) - psi'(r_s-eps)] / psi(r_s)
 *
 * Delta' > 0: tearing unstable (magnetic islands grow)
 * Delta' < 0: tearing stable
 *
 * For a simple slab geometry with current sheet.
 */
double tearing_stability_delta_prime(double k, double a);

#endif /* PLASMA_MHD_H */
