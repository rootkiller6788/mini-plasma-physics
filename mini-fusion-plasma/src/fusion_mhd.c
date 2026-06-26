/**
 * fusion_mhd.c — MHD Equilibrium and Stability
 *
 * Every function implements an independent physics knowledge point.
 * Refs: Freidberg "Ideal MHD" (2014), Wesson (2011),
 *       Goedbloed & Poedts "Principles of MHD" (2004),
 *       Grad & Rubin (1958), Shafranov (1966),
 *       Lao et al. (1985), Greene & Chance (1981)
 *
 * L2: Flux surface quantities, Shafranov shift,
 *     toroidal/poloidal flux, connection length
 * L3: Grad-Shafranov equation, flux coordinates,
 *     metric tensor elements, rotational transform
 * L4: Ideal MHD equilibrium force balance,
 *     virial theorem for confined plasmas
 * L5: Numerical equilibrium reconstruction (simplified),
 *     Green's function method for external measurements
 * L8: Vertical stability, resistive wall modes,
 *     error field amplification
 * L6: ITER/DEMO/SPARC equilibrium parameters
 */

#include "include/fusion_plasma.h"
#include <math.h>
#include <assert.h>
#include <string.h>

/* ================================================================
 * L3: Flux Surface Geometry — 5 functions
 * ================================================================ */

/**
 * shafranov_shift — Shafranov (outward) shift of flux surfaces [m]
 *
 * Definition: Delta(r) = R(r) - R_0 (outward shift of flux surface)
 *          For circular cross-section:
 *          Delta'(r) ~ (r/R) * (beta_p + li/2)
 * Physics: toroidicity + plasma pressure shifts flux surfaces outward.
 *          The Shafranov shift increases with beta_p and internal inductance.
 *          Large Shafranov shift -> strong Shafranov shift -> degraded stability.
 *
 * @param r       minor radius of flux surface [m]
 * @param R0      major radius of magnetic axis [m]
 * @param beta_p  poloidal beta
 * @param li      internal inductance
 * @return        Shafranov shift Delta(r) [m]
 */
double shafranov_shift(double r, double R0, double beta_p, double li) {
    assert(R0 > 0.0);
    return (r * r / (2.0 * R0)) * (beta_p + li / 2.0);
}

/**
 * connection_length — Magnetic field line connection length [m]
 *
 * Definition: L_c = 2 * pi * R * q (distance along B between
 *          inboard and outboard midplane intersections)
 * Physics: parallel connection length sets parallel transport
 *          timescale. Important for SOL physics and divertor design.
 *
 * For ITER: L_c ~ 100-200 m (SOL), ~500 m (core)
 */
double connection_length(double R, double q) {
    return 2.0 * M_PI * R * q;
}

/**
 * rotational_transform — Rotational transform iota (stellarator)
 *
 * Definition: iota = 1/q = dtheta_poloidal/dzeta_toroidal
 * Physics: in stellarators, iota is the natural quantity
 *          (field lines rotate poloidally as they go toroidally).
 *          iota = 1/q maps between tokamak and stellarator concepts.
 */
double rotational_transform(double q) {
    assert(q > 0.0);
    return 1.0 / q;
}

/**
 * toroidal_flux — Toroidal magnetic flux [Wb]
 *
 * Definition: Phi = integral B_phi * dA_poloidal
 *          ~ B_phi * pi * r^2 (circular, uniform B)
 * Physics: total toroidal magnetic flux through poloidal
 *          cross-section. Constant on each flux surface.
 */
double toroidal_flux(double B_phi, double r) {
    return B_phi * M_PI * r * r;
}

/**
 * poloidal_flux — Poloidal magnetic flux per radian [Wb/rad]
 *
 * Definition: psi = integral B_pol * dA_toroidal / (2*pi)
 *          psi = - (1/(2*pi)) integral_0^r B_pol * 2*pi*R * dr
 * Physics: poloidal flux function. psi = constant defines
 *          magnetic surfaces. Used as radial coordinate in
 *          Grad-Shafranov equation.
 */
double poloidal_flux_from_current(double r, double R, double B_phi, double q) {
    assert(q > 0.0);
    double B_pol = (r / (R * q)) * B_phi;
    return B_pol * r * R;
}

/* ================================================================
 * L4: Ideal MHD Equilibrium — 6 functions
 * ================================================================ */

/**
 * force_balance_residual — MHD force balance residual
 *
 * Definition: F = J x B - grad(p)
 * Physics: ideal MHD equilibrium: force balance between
 *          magnetic (JxB) and pressure gradient forces.
 *          At equilibrium, F = 0 everywhere.
 * This function evaluates the residual at a point.
 *
 * @param j      current density vector [A/m^2]
 * @param B      magnetic field vector [T]
 * @param grad_p pressure gradient vector [Pa/m]
 * @return       |F|, the magnitude of force imbalance [N/m^3]
 */
double force_balance_residual(Vec3 j, Vec3 B, Vec3 grad_p) {
    Vec3 jxB = vec3_cross(j, B);
    Vec3 residual = vec3_sub(jxB, grad_p);
    return vec3_norm(residual);
}

/**
 * virial_theorem_check — Virial theorem for MHD equilibrium
 *
 * Definition: integral (3p + B^2/(2*mu0)) dV = integral (p + B^2/(2*mu0)) r·dS
 * Physics: virial theorem states that a plasma cannot be confined
 *          purely by its own magnetic field; external coils are
 *          necessary for MHD equilibrium.
 * This function computes the volumetric integrand: 3p + B^2/(2*mu0)
 *
 * @param p      thermal pressure [Pa]
 * @param B      magnetic field [T]
 * @return       virial integrand [J/m^3]
 */
double virial_integrand(double p, double B) {
    return 3.0 * p + B * B / (2.0 * MU0);
}

/**
 * internal_inductance — Plasma internal inductance li
 *
 * Definition: li = (2/B_p^2(a) * V_p) * integral (B_p^2/2*mu0) dV
 *          ~ (2*V_p)/(mu0*R*I_p^2) * integral B_p^2 dV
 * Physics: measure of current profile peakedness.
 *          li = 0.5: very peaked current (skin)
 *          li = 1.0: uniform current density
 *          li ~ 0.8-1.2: typical H-mode profiles
 *
 * Simplified formula assuming parabolic current profile:
 *   li = ln(1.65 + 0.89*(q_95 - 1))
 */
double internal_inductance(double q_95) {
    if (q_95 < 1.0) q_95 = 1.01;
    return log(1.65 + 0.89 * (q_95 - 1.0));
}

/**
 * normalized_internal_inductance — Normalized internal inductance li(3)
 *
 * Definition: li(3) = <B_p^2>_V / B_p^2(a)
 *          Used in the ITER confinement database as independent variable.
 *          Typically li(3) ~ 0.7-1.0 for standard H-mode.
 */
double normalized_internal_inductance(double li) {
    return li;
}

/**
 * poloidal_beta_from_equilibrium — beta_p from equilibrium
 *
 * Definition: beta_p = (2*mu0/(B_p^2(a)*V_p)) * integral p dV
 * Physics: ratio of volume-averaged pressure to edge poloidal
 *          magnetic pressure. Key parameter in Grad-Shafranov.
 */
double poloidal_beta_from_params(double p_avg, double B_p_a) {
    assert(B_p_a > 0.0);
    return 2.0 * MU0 * p_avg / (B_p_a * B_p_a);
}

/**
 * edge_poloidal_field — Edge poloidal magnetic field [T]
 *
 * Definition: B_p(a) = mu0 * I_p / (2 * pi * a) (circular)
 * For elongated plasma: B_p(a) = mu0 * I_p / L_pol
 *          where L_pol is the poloidal circumference.
 */
double edge_poloidal_field(double Ip, double a, double kappa) {
    double L_pol = 2.0 * M_PI * a * sqrt((1.0 + kappa*kappa) / 2.0);
    if (L_pol <= 0.0) return 0.0;
    return MU0 * Ip / L_pol;
}

/* ================================================================
 * L5: Equilibrium Reconstruction — 4 functions
 * ================================================================ */

/**
 * solve_cylindrical_equilibrium — Cylindrical equilibrium solution
 *
 * Simplified 1D cylindrical force balance:
 *   dp/dr + (B_phi/(mu0*r)) * d(r*B_phi)/dr + B_z * dB_z/dr = 0
 *
 * Returns poloidal beta for given pressure and field profiles.
 * This simplified model captures the essential force balance
 * without solving the full Grad-Shafranov equation.
 *
 * @param p0         central pressure [Pa]
 * @param B_phi0     central toroidal field [T]
 * @param B_phi_a    edge toroidal field [T]
 * @param a          plasma radius [m]
 * @return           volume-averaged beta
 */
double cylindrical_equilibrium_beta(double p0, double B_phi0,
                                     double B_phi_a, double a) {
    double B_pol_a = sqrt(fabs(B_phi_a*B_phi_a - B_phi0*B_phi0));
    if (B_pol_a < 1e-10) B_pol_a = MU0 * 1e6 / (2.0 * M_PI * a);
    double p_avg = 0.5 * p0;  /* parabolic profile average */
    return 2.0 * MU0 * p_avg / (B_pol_a * B_pol_a);
}

/**
 * grad_shafranov_operator — Grad-Shafranov elliptic operator
 *
 * Definition: Delta* psi = R^2 div((1/R^2) grad psi)
 *          = R * d/dR((1/R) * dpsi/dR) + d^2psi/dZ^2
 * Physics: the GS operator represents the toroidal current
 *          density contribution. J_phi = (1/(mu0*R)) * Delta* psi
 *
 * This function computes Delta* psi using finite differences
 * on a 2D grid (R-Z plane) at a single point.
 *
 * @param psi_grid  2D array of psi values (row-major, ny*nx)
 * @param R_grid    1D array of R coordinates (nx elements)
 * @param Z_grid    1D array of Z coordinates (ny elements)
 * @param i, j      grid indices of evaluation point
 * @param nx, ny    grid dimensions
 * @param dR, dZ    grid spacing
 * @return          Delta* psi at (i,j) [Wb/m]
 */
double delta_star_psi(const double *psi_grid, const double *R_grid,
                       const double *Z_grid, int i, int j,
                       int nx, int ny, double dR, double dZ) {
    (void)Z_grid; (void)ny;  /* parameters reserved for future use */
    if (i <= 0 || i >= nx-1 || j <= 0 || j >= ny-1) return 0.0;
    int idx = j * nx + i;
    double R = R_grid[i];
    if (R < 1e-10) return 0.0;

    /* d/dR((1/R) * dpsi/dR) */
    double dpsi_dR_plus = (psi_grid[idx+1] - psi_grid[idx]) / dR;
    double dpsi_dR_minus = (psi_grid[idx] - psi_grid[idx-1]) / dR;
    double term_R = (dpsi_dR_plus/R - dpsi_dR_minus/R) / dR;

    /* d^2psi/dZ^2 */
    double term_Z = (psi_grid[idx+nx] - 2.0*psi_grid[idx] + psi_grid[idx-nx]) / (dZ * dZ);

    return R * term_R + term_Z;
}

/**
 * magnetic_axis_find — Find magnetic axis from GS solution
 *
 * Definition: magnetic axis is where grad(psi) = 0 (maximum psi)
 * Simple method: scan for maximum psi in the grid.
 *
 * @param psi_grid  2D psi array
 * @param nx, ny    grid dimensions
 * @param R_axis    output: R of magnetic axis
 * @param Z_axis    output: Z of magnetic axis
 * @param psi_axis  output: psi at magnetic axis
 */
void magnetic_axis_find(const double *psi_grid, int nx, int ny,
                         double *R_axis, double *Z_axis, double *psi_axis,
                         const double *R_grid, const double *Z_grid) {
    double psi_max = -1e100;
    int i_max = 0, j_max = 0;
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            int idx = j * nx + i;
            if (psi_grid[idx] > psi_max) {
                psi_max = psi_grid[idx];
                i_max = i;
                j_max = j;
            }
        }
    }
    *R_axis = R_grid[i_max];
    *Z_axis = Z_grid[j_max];
    *psi_axis = psi_max;
}

/**
 * boundary_flux_find — Find plasma boundary (last closed flux surface)
 *
 * Simplification: boundary is where psi = 0 (or minimum |psi| on limiter).
 * Scans outward from axis until psi changes sign or reaches domain edge.
 */
double boundary_flux_find(const double *psi_grid, int nx, int ny,
                           double psi_axis, const double *R_grid,
                           const double *Z_grid, int axis_i, int axis_j) {
    (void)Z_grid; (void)ny; (void)psi_axis; (void)R_grid;  /* reserved for full implementation */
    double psi_min = 1e100;
    for (int i = axis_i; i < nx; i++) {
        int idx = axis_j * nx + i;
        double p = psi_grid[idx];
        if (p < psi_min) psi_min = p;
        if (p <= 0.0) break;
    }
    return (psi_min < 1e99) ? psi_min : 0.0;
}

/* ================================================================
 * L8: Advanced Stability — 5 functions
 * ================================================================ */

/**
 * vertical_stability_growth_rate — Vertical instability growth rate [1/s]
 *
 * Definition: gamma_n = (1/tau_wall) * (n_index / (1 - n_index))
 * Physics: elongated plasmas are vertically unstable.
 *          Growth rate depends on wall proximity and elongation.
 *          n_index: decay index of external field n = -(R/B_z)*dB_z/dR
 * Active feedback control required for kappa > 1.2-1.5.
 *
 * @param n_index   field decay index
 * @param tau_wall  wall penetration time [s]
 * @return          growth rate gamma [1/s]
 */
double vertical_stability_growth_rate(double n_index, double tau_wall) {
    if (tau_wall <= 0.0) return 1e10;
    if (n_index >= 1.0) n_index = 0.99;
    if (n_index <= 0.0) return 0.0;
    return (1.0 / tau_wall) * (n_index / (1.0 - n_index));
}

/**
 * resistive_wall_mode_growth_rate — RWM growth rate [1/s]
 *
 * Definition: gamma_RWM = (1/tau_wall) * (beta - beta_no_wall) /
 *                         (beta_ideal_wall - beta)
 * Physics: resistive wall modes are external kinks stabilized
 *          by a nearby conducting wall. Growth rate set by
 *          wall diffusion time.
 * For ITER: tau_wall ~ 0.1-0.5 s, gamma_RWM ~ 2-10 s^-1
 */
double resistive_wall_mode_growth_rate(double beta, double beta_no_wall,
                                        double beta_ideal_wall, double tau_wall) {
    if (tau_wall <= 0.0) return 1e10;
    if (beta >= beta_ideal_wall) return 1e6;  /* wall can't stabilize */
    if (beta <= beta_no_wall) return 0.0;      /* stable without wall */
    return (1.0 / tau_wall) * (beta - beta_no_wall) / (beta_ideal_wall - beta);
}

/**
 * error_field_amplification — Error field amplification factor
 *
 * Definition: A = 1 / (1 - n/n_crit)
 * Physics: external error fields are amplified as plasma approaches
 *          the no-wall stability limit. Resonant at rational surfaces.
 *          n_crit: critical error field for mode locking.
 *
 * @param applied_field    applied error field [T]
 * @param n_index          normalized proximity to stability boundary
 * @return                 amplification factor
 */
double error_field_amplification(double n_index) {
    if (n_index >= 1.0) return 1e6;
    if (n_index <= 0.0) return 1.0;
    return 1.0 / (1.0 - n_index);
}

/**
 * plasma_inductance — Plasma self-inductance [H]
 *
 * Definition: L_p = mu0 * R * (ln(8*R/a) - 2 + li/2)
 * Physics: self-inductance of the plasma current ring.
 *          Affects the L/R time of the discharge and
 *          the flux consumption during current ramp-up.
 * For ITER: L_p ~ 3-5 micro-H
 */
double plasma_inductance(double R, double a, double li) {
    assert(a > 0.0);
    double ln_term = log(8.0 * R / a);
    return MU0 * R * (ln_term - 2.0 + li / 2.0);
}

/**
 * flux_consumption_during_startup — Volt-second consumption [V*s = Wb]
 *
 * Definition: Delta_psi = L_p * I_p (resistive)
 *          + external flux for breakdown and ramp-up
 * Physics: total flux (volt-seconds) needed to establish
 *          the plasma current. Central solenoid must provide
 *          this flux. Critical constraint for pulsed tokamaks.
 *
 * @param L_p    plasma inductance [H]
 * @param Ip     plasma current [A]
 * @param psi_bd flux needed for breakdown [Wb]
 * @return       total flux consumption [V*s = Wb]
 */
double flux_consumption(double L_p, double Ip, double psi_bd) {
    return L_p * Ip + psi_bd;
}