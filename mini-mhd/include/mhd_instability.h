/*
 * mhd_instability.h -- MHD Instabilities and Stability Analysis
 *
 * Reference: Freidberg "Ideal MHD" (2014) Ch.9-12
 *            Bateman "MHD Instabilities" (1978)
 *            Wesson "Tokamaks" (2011) Ch.6
 *
 * Knowledge: L2 -- MHD instability concepts (kink, sausage, interchange)
 *             L4 -- Energy principle for MHD stability
 *             L6 -- Kruskal-Shafranov kink limit
 *             L8 -- Ballooning modes, resistive instabilities
 */

#ifndef MHD_INSTABILITY_H
#define MHD_INSTABILITY_H

#include "mhd_defs.h"

/* ================================================================
 * L4 -- MHD Energy Principle (Bernstein et al. 1958)
 *
 * An MHD equilibrium is linearly stable iff the potential energy
 * perturbation delta_W > 0 for all admissible displacements xi(r).
 *
 * delta_W = delta_W_F + delta_W_S + delta_W_V
 *   delta_W_F: fluid (plasma) term   -- always stabilizing
 *   delta_W_S: surface term          -- can be destabilizing
 *   delta_W_V: vacuum term           -- stabilizing
 *
 * The energy principle reduces stability analysis to a variational problem.
 * ================================================================ */

/*
 * mhd_energy_principle_fluid -- Compute delta_W_F (fluid term)
 *
 * delta_W_F = 0.5 * int dV [ |Q|^2/mu_0 + gamma*p*|div(xi)|^2
 *            + B^2*|div_perp(xi) + 2*xi_perp.kappa|^2/mu_0
 *            - 2*(xi_perp.grad p)*(kappa.xi_perp) - J_par*(xi_perp x b).Q ]
 *
 * where Q = curl(xi x B) and kappa = (b.grad)b is field line curvature.
 *
 * Input: equilibrium state, displacement xi, local geometry.
 * Returns: local delta_W_F density.
 *
 * Complexity: O(1)
 */
double mhd_energy_principle_fluid(const MHDState *eq,
                                   double xi_x, double xi_y, double xi_z,
                                   double gamma);

/*
 * mhd_suydam_criterion -- Suydam criterion for localized interchange
 *
 * For a cylindrical screw pinch, stability against localized interchange
 * requires: r*B_z^2/(8*mu_0) * (dq/dr)^2/q^2 + dp/dr > 0
 *
 * Returns: Suydam parameter. Positive = stable, negative = unstable.
 *
 * Complexity: O(1)
 */
double mhd_suydam_criterion(double r, double Bz, double q, double dq_dr,
                              double dp_dr);

/* ================================================================
 * L6 -- MHD Instability Criteria
 * ================================================================ */

/*
 * mhd_kruskal_shafranov_limit -- Kink instability condition
 *
 * For a tokamak with major radius R_0, minor radius a:
 *   q(a) = a B_phi / (R_0 B_theta) > 1
 *
 * If q(a) < 1, the plasma is kink-unstable (external kink).
 * This is the fundamental tokamak stability limit.
 *
 * Returns: 1 if stable (q > 1), 0 if marginal, -1 if unstable.
 *
 * Complexity: O(1)
 */
int mhd_kruskal_shafranov_condition(double a, double R0,
                                     double B_phi, double B_theta);

/*
 * mhd_kruskal_shafranov_margin -- Safety factor margin
 * Returns q(a) - 1. Positive = stable, negative = unstable.
 *
 * Complexity: O(1)
 */
double mhd_kruskal_shafranov_margin(double q_a);

/*
 * mhd_sausage_instability_growth -- Sausage (m=0) instability growth rate
 *
 * For a Z-pinch: sausage mode is unstable when dp/dr < 0 (Bennett profile).
 * Growth rate ~ sqrt(-(2/mu_0) * (B_phi^2/r) * (dp/dr) / (rho*|B_phi|)).
 *
 * Returns: growth rate gamma > 0 if unstable, 0 if stable.
 *
 * Complexity: O(1)
 */
double mhd_sausage_growth_rate(double r, double B_phi, double dp_dr, double rho);

/*
 * mhd_kink_growth_rate -- Kink (m=1) instability growth rate
 *
 * For a screw pinch: growth rate proportional to (1 - 1/q(a)).
 * gamma ~ (v_A / R_0) * |1 - 1/q|
 *
 * Returns: growth rate.
 *
 * Complexity: O(1)
 */
double mhd_kink_growth_rate(double va, double R0, double q);

/*
 * mhd_interchange_criterion -- Interchange (flute) instability
 *
 * In a curved magnetic field, interchange instability occurs when
 * the magnetic field curvature is "unfavorable":
 *   kappa . grad p < 0  (bad curvature)
 *
 * Equivalent: pressure gradient points away from curvature center.
 *
 * Returns: positive if unstable (bad curvature), negative if stable.
 *
 * Complexity: O(1)
 */
double mhd_interchange_criterion(double kappa_x, double kappa_y, double kappa_z,
                                  double gradP_x, double gradP_y, double gradP_z);

/* ================================================================
 * L8 -- Advanced Instabilities
 * ================================================================ */

/*
 * mhd_ballooning_mode -- Ballooning mode stability
 *
 * In tokamaks, localized ballooning modes occur on the outboard side
 * where field line curvature is unfavorable.
 * Stability determined by the Mercier criterion and ballooning equation.
 *
 * Simplified criterion: alpha = -q^2 R dp/dr * (2 mu_0/B^2) < alpha_crit.
 *
 * Returns: ballooning parameter alpha. Lower is more stable.
 *
 * Complexity: O(1)
 */
double mhd_ballooning_parameter(double q, double R, double dp_dr,
                                  double B_phi);

/*
 * mhd_mercier_criterion -- Mercier criterion for localized modes
 *
 * Generalizes Suydam to toroidal geometry.
 * D_M = D_S + D_I + D_G < 1/4 for stability
 *
 * Returns: D_M value (< 0.25 is stable).
 *
 * Complexity: O(1)
 */
double mhd_mercier_criterion(double q, double dq_dr, double dp_dr,
                              double R, double r, double B0);

/*
 * mhd_resistive_tearing_growth -- Tearing mode growth rate
 *
 * Resistive instability that reconnects field lines at rational surface.
 * Growth rate: gamma ~ tau_R^{-3/5} tau_A^{-2/5}
 * where tau_R = a^2/eta (resistive), tau_A = a/v_A (Alfven).
 *
 * Most dangerous instability in tokamaks. Causes magnetic islands.
 *
 * Complexity: O(1)
 */
double mhd_tearing_growth_rate(double a, double va, double eta,
                                double Delta_prime);

/*
 * mhd_rayleigh_taylor_growth -- MHD Rayleigh-Taylor growth rate
 *
 * When a heavy fluid is supported by a light fluid in a magnetic field,
 * the growth rate is: gamma^2 = g*k*A - (k.B)^2/(mu_0*rho)
 * where A = (rho_heavy - rho_light)/(rho_heavy + rho_light) is Atwood number.
 *
 * Magnetic field along the perturbation (k.B != 0) suppresses RT.
 *
 * Complexity: O(1)
 */
double mhd_rayleigh_taylor_growth(double g, double k, double Atwood,
                                   double k_dot_B, double rho);

/*
 * mhd_kelvin_helmholtz_growth -- MHD Kelvin-Helmholtz growth rate
 *
 * Velocity shear across interface. MHD stabilization when:
 *   (k.B)^2/(mu_0*rho) > (k.delta_v)^2/4
 *
 * Returns: growth rate.
 *
 * Complexity: O(1)
 */
double mhd_kelvin_helmholtz_growth(double k, double delta_v,
                                    double k_dot_B, double rho);

/*
 * mhd_linear_stability_scan -- Scan stability across parameter space
 *
 * Evaluates all instability criteria for a given equilibrium
 * and returns a bitmask of unstable modes:
 *   bit 0: kink, bit 1: sausage, bit 2: interchange
 *   bit 3: ballooning, bit 4: tearing, bit 5: RT, bit 6: KH
 *
 * Complexity: O(1)
 */
int mhd_linear_stability_scan(const MHDState *eq,
                               double q, double R0, double a,
                               double eta, double gamma);

#endif /* MHD_INSTABILITY_H */
