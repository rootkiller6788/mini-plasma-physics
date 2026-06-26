/**
 * fusion_confinement.c — Magnetic Confinement and Transport Physics
 *
 * Every function implements an independent physics knowledge point.
 * Refs: Wesson "Tokamaks" 4th ed (2011), Freidberg (2007),
 *       ITER Physics Basis (1999), 
 *       Connor & Taylor (1977), Troyon et al. (1984),
 *       Goldston & Rutherford (1995)
 *
 * L2: Bootstrap current, Pfirsch-Schluter factor,
 *     trapped particle fraction
 * L4: Troyon beta limit, Kruskal-Shafranov limit,
 *     neoclassical transport regimes (banana/plateau/Pfirsch-Schluter)
 * L5: Bohm/Gyro-Bohm diffusion, neoclassical ion heat diffusivity
 * L6: H-mode power threshold, pedestal physics
 * L8: Tearing mode stability, ballooning stability, Mercier criterion,
 *     magnetic shear, sawtooth oscillations, disruptions,
 *     NTM threshold, halo currents
 */

#include "include/fusion_plasma.h"
#include <math.h>
#include <assert.h>

/* ================================================================
 * L4: Beta Limits and Stability Thresholds — 4 functions
 * ================================================================ */

/**
 * troyon_beta_limit — Troyon beta limit evaluation
 *
 * Definition: beta_N = beta(%) * a * B / Ip <= beta_N_crit
 *          beta_N_crit ~ 2.8 (no wall) to 3.5 (ideal wall)
 * Physics: ideal MHD pressure-driven (ballooning/kink) limit.
 *          Exceeding Troyon limit causes major disruptions.
 * Source: Troyon et al., PPCF 26, 209 (1984)
 *
 * Returns the normalized beta value (beta_N).
 * Compare to beta_N_crit (~2.8-3.5) for stability assessment.
 */
double troyon_beta_limit(double beta_percent, double a, double B, double Ip_MA) {
    assert(Ip_MA > 0.0);
    return beta_percent * a * B / Ip_MA;
}

/**
 * normalized_beta — Normalized beta beta_N
 *
 * Definition: beta_N = beta(%) * a * B / Ip
 * Same as Troyon parameter. Used for comparing confinement
 * performance across different machines.
 * ITER: beta_N ~ 1.8 (baseline), up to 3.0 (advanced)
 */
double normalized_beta(double beta_percent, double a, double B, double Ip_MA) {
    return troyon_beta_limit(beta_percent, a, B, Ip_MA);
}

/**
 * kruskal_shafranov_limit — Kruskal-Shafranov current limit
 *
 * Definition: q_edge > 1 for external kink mode stability
 * Physics: if q_edge < 1, the plasma is unstable to ideal
 *          external kink modes (m=1, n=1 dominant).
 *          This sets the maximum plasma current:
 *          Ip_max = (2*pi*a^2*B_phi) / (mu0*R)
 * Source: Kruskal & Schwarzschild (1954), Shafranov (1957)
 *
 * Returns 1 if stable (q > 1), 0 if unstable.
 */
int kruskal_shafranov_limit(double q_edge) {
    return (q_edge > 1.0) ? 1 : 0;
}

/**
 * pfirsch_schluter_factor — Pfirsch-Schluter toroidal enhancement
 *
 * Definition: enhancement ~ (1 + q^2)
 * Physics: toroidicity enhances parallel current and transport
 *          by factor (1+q^2) relative to cylindrical approximation.
 *          Arises from the divergence of the diamagnetic current
 *          in toroidal geometry.
 */
double pfirsch_schluter_factor(double q) {
    return 1.0 + q * q;
}

/* ================================================================
 * L2: Bootstrap Current — 2 functions
 * ================================================================ */

/**
 * bootstrap_current_fraction — Bootstrap current fraction
 *
 * Definition: f_BS = I_BS / I_p ~ C_BS * sqrt(epsilon) * beta_p
 *          C_BS ~ 0.5-1.3 depending on profiles
 * Physics: neoclassical effect where radial pressure gradient
 *          drives toroidal current without external loop voltage.
 *          Key for steady-state tokamak operation.
 * For ITER: f_BS ~ 0.3-0.5 (30-50% of total current)
 */
double bootstrap_current_fraction(double epsilon, double beta_p) {
    double C_BS = 0.68;  /* typical coefficient */
    return C_BS * sqrt(epsilon) * beta_p;
}

/**
 * bootstrap_current_density — Bootstrap current density [A/m^2]
 *
 * Definition: j_BS = - (sqrt(epsilon) / B_p) * dp/dr
 *          (simplified form, omits collisionality dependence)
 * Physics: bootstrap current flows where pressure gradient exists,
 *          peaking near plasma edge in H-mode pedestal.
 */
double bootstrap_current_density(double epsilon, double B_pol, double dp_dr) {
    assert(B_pol > 0.0);
    return -sqrt(epsilon) * dp_dr / B_pol;
}

/* ================================================================
 * L5: Neoclassical Transport Regimes — 4 functions
 * ================================================================ */

/**
 * banana_regime_diffusivity — Banana regime neoclassical diffusivity [m^2/s]
 *
 * Definition: D_b = epsilon^(-3/2) * q^2 * rho_i^2 * nu_ii
 * Valid for: nu* < 1 (low collisionality)
 * Physics: trapped particles execute banana orbits, enhancing
 *          radial transport by factor epsilon^(-3/2) ~ 10-30
 *          relative to cylindrical value.
 * Source: Galeev & Sagdeev (1968); Rosenbluth, Hazeltine, Hinton (1972)
 */
double banana_regime_diffusivity(double epsilon, double q, double rho_i,
                                  double nu_ii) {
    assert(epsilon > 0.0);
    double eps_inv = 1.0 / (epsilon * sqrt(epsilon));
    return eps_inv * q * q * rho_i * rho_i * nu_ii;
}

/**
 * plateau_regime_diffusivity — Plateau regime diffusivity [m^2/s]
 *
 * Definition: D_p = (sqrt(pi)/2) * (rho_i^2 * v_thi) / (R * q)
 * Valid for: nu* ~ 1
 * Physics: resonant particles with v_parallel ~ 0 dominate.
 *          Diffusivity independent of collision frequency.
 */
double plateau_regime_diffusivity(double rho_i, double v_thi, double R, double q) {
    assert(R > 0.0);
    assert(q > 0.0);
    return 0.5 * sqrt(M_PI) * rho_i * rho_i * v_thi / (R * q);
}

/**
 * pfirsch_schluter_diffusivity — Pfirsch-Schluter regime diffusivity [m^2/s]
 *
 * Definition: D_PS = q^2 * rho_i^2 * nu_ii
 * Valid for: nu* > 1 (high collisionality)
 * Physics: toroidal enhancement of collisional (Pfirsch-Schluter)
 *          transport by factor q^2.
 */
double pfirsch_schluter_diffusivity(double q, double rho_i, double nu_ii) {
    return q * q * rho_i * rho_i * nu_ii;
}

/**
 * neoclassical_ion_heat_diffusivity — Neoclassical ion heat diffusivity [m^2/s]
 *
 * Unified formula covering all three regimes.
 * Uses Chang-Hinton formula:
 *   chi_i_NC = epsilon^(1/2) * q^2 * (rho_i^2 * nu_ii) * K(epsilon, nu*)
 * where K is a transition function between regimes.
 * Simplified as maximum of banana and Pfirsch-Schluter values.
 *
 * Source: Chang & Hinton, Phys. Fluids 25, 1493 (1982)
 */
double neoclassical_ion_heat_diffusivity(double epsilon, double q, double rho_i,
                                          double nu_ii, double v_thi, double R) {
    double nu_star = collisionality_nu_star(R, nu_ii, v_thi, epsilon);
    double D_b = banana_regime_diffusivity(epsilon, q, rho_i, nu_ii);
    double D_PS = pfirsch_schluter_diffusivity(q, rho_i, nu_ii);
    /* Smooth interpolation between regimes */
    double w_banana = 1.0 / (1.0 + nu_star * nu_star);
    double w_ps = nu_star * nu_star / (1.0 + nu_star * nu_star);
    return w_banana * D_b + w_ps * D_PS;
}

/* ================================================================
 * L5: Anomalous/Turbulent Transport — 2 functions
 * ================================================================ */

/**
 * bohm_diffusion_coefficient — Bohm diffusion coefficient [m^2/s]
 *
 * Definition: D_Bohm = (1/16) * kB*Te / (e*B)
 * Physics: semi-empirical maximum turbulent diffusion.
 *          Originally from arc discharge experiments (Bohm, 1949).
 *          Represents worst-case anomalous transport.
 *          D_Bohm ~ 1-10 m^2/s for typical fusion parameters,
 *          much larger than neoclassical (~0.01 m^2/s).
 */
double bohm_diffusion_coefficient(double Te_eV, double B) {
    assert(B > 0.0);
    double Te_J = Te_eV * E_CHARGE;
    return (1.0 / 16.0) * Te_J / (E_CHARGE * B);
}

/**
 * gyrobohm_diffusion_coefficient — Gyro-Bohm diffusion [m^2/s]
 *
 * Definition: D_gB = (rho_i / a) * D_Bohm
 *          = (1/16) * (rho_i/a) * kB*Te/(e*B)
 * Physics: drift-wave turbulence scaling.
 *          Gyro-Bohm scaling matches many experimental observations
 *          better than Bohm scaling.
 *          D_gB << D_Bohm since rho_i/a ~ 10^-3.
 */
double gyrobohm_diffusion_coefficient(double Te_eV, double B, double rho_i, double a) {
    double D_Bohm = bohm_diffusion_coefficient(Te_eV, B);
    assert(a > 0.0);
    return (rho_i / a) * D_Bohm;
}

/* ================================================================
 * L6: H-mode and Pedestal Physics — 3 functions
 * ================================================================ */

/**
 * h_mode_power_threshold — L-H transition power threshold [W]
 *
 * Definition: P_LH = 0.049 * n_20^0.72 * B^0.8 * S^0.94
 * Physics: minimum heating power required for L-H transition.
 *          Below P_LH: L-mode (low confinement).
 *          Above P_LH: H-mode (high confinement, ~2x better).
 * Source: Martin et al., J. Phys. Conf. Ser. 123, 012033 (2008)
 *          (ITER multi-machine scaling)
 *
 * @param n_20  line-averaged density [10^20 m^-3]
 * @param B     toroidal field [T]
 * @param S     plasma surface area [m^2]
 * @return      L-H threshold power [W]
 */
double h_mode_power_threshold(double n_20, double B, double S) {
    return 0.049 * pow(n_20, 0.72) * pow(B, 0.8) * pow(S, 0.94);
}

/**
 * pedestal_pressure — Pedestal pressure [Pa]
 *
 * Definition: p_ped ~ (I_p/a)^(3/2) * B^(1/2)  (simplified)
 * Physics: H-mode edge transport barrier (ETB) creates a pedestal
 *          in pressure profile. Pedestal height strongly affects
 *          global confinement (stored energy ~ pedestal).
 * EPED model: pedestal width set by kinetic ballooning mode (KBM)
 *             and peeling-ballooning stability boundary.
 *
 * @param Ip_MA  plasma current [MA]
 * @param a      minor radius [m]
 * @param B      toroidal field [T]
 * @return       pedestal pressure [Pa]
 */
double pedestal_pressure(double Ip_MA, double a, double B) {
    assert(a > 0.0);
    double j_avg = Ip_MA * 1e6 / (M_PI * a * a);
    return pow(j_avg, 1.5) * sqrt(B) * 1e-6;
}

/**
 * pedestal_width — Pedestal width [m]
 *
 * Definition: Delta_ped ~ C * sqrt(beta_p_ped) * sqrt(a * rho_i)
 *          or simpler: Delta_ped ~ 0.02-0.05 * a
 * Physics: width of H-mode edge transport barrier.
 *          Set by interplay of KBM turbulence and neoclassical
 *          (banana regime) physics.
 * For ITER: Delta_ped ~ 0.04-0.08 m
 *
 * @param beta_p_ped  poloidal beta at pedestal
 * @param a           minor radius [m]
 * @return            pedestal width [m]
 */
double pedestal_width(double beta_p_ped, double a) {
    double C = 0.04;  /* typical: 4% of minor radius */
    return C * sqrt(beta_p_ped) * a;
}

/* ================================================================
 * L8: MHD Stability — 7 functions
 * ================================================================ */

/**
 * tearing_mode_delta_prime — Tearing mode stability index delta'
 *
 * Definition: delta' = (psi'(r_s+) - psi'(r_s-)) / psi(r_s)
 * Physics: classical tearing mode stability parameter.
 *          delta' > 0: linearly unstable (grows on resistive time).
 *          delta' < 0: stable.
 *          Rational surface r_s where q = m/n.
 * Source: Furth, Killeen, Rosenbluth, Phys. Fluids 6, 459 (1963)
 *
 * @param psi_prime_plus   dpsi/dr just outside resonant surface
 * @param psi_prime_minus  dpsi/dr just inside resonant surface
 * @param psi_resonant     psi at resonant surface
 * @return                 delta' [1/m]
 */
double tearing_mode_delta_prime(double psi_prime_plus, double psi_prime_minus,
                                double psi_resonant) {
    if (fabs(psi_resonant) < 1e-30) return 0.0;
    return (psi_prime_plus - psi_prime_minus) / psi_resonant;
}

/**
 * ballooning_alpha — Ballooning stability parameter alpha
 *
 * Definition: alpha = - q^2 * R * dp/dr / (B^2 / (2*mu0))
 * Physics: normalized pressure gradient driving ballooning modes.
 *          alpha > alpha_crit: ballooning unstable (pressure-driven).
 *          alpha_crit ~ 0.4-0.8 depending on magnetic shear and geometry.
 * Source: Connor, Hastie, Taylor, Phys. Rev. Lett. 40, 396 (1978)
 *
 * For ITER: alpha ~ 0.2-0.4 (below critical)
 */
double ballooning_alpha(double q, double R, double dp_dr, double B) {
    double p_mag = B * B / (2.0 * MU0);
    if (p_mag < 1e-30) return 0.0;
    return -q * q * R * dp_dr / p_mag;
}

/**
 * mercier_criterion — Mercier interchange stability criterion
 *
 * Definition: D_M = 0.25 - (q^2 * dp/dpsi * V'' / ...)
 * Physics: localized interchange stability in toroidal geometry.
 *          D_M < 0: unstable to interchange modes.
 *          Related to magnetic well depth and shear.
 * Source: Mercier, Nucl. Fusion 1, 47 (1960)
 *
 * Simplified form: D_M = 0.25 + (shear-1)^2/4 - alpha
 * Positive -> stable to interchange.
 */
double mercier_criterion(double q, double shear, double magnetic_well,
                          double pressure_gradient) {
    double stab_term = 0.25 + (shear - 1.0) * (shear - 1.0) / 4.0;
    double drive_term = q * q * pressure_gradient * (1.0 - magnetic_well);
    return stab_term - drive_term;
}

/**
 * magnetic_shear — Magnetic shear parameter
 *
 * Definition: s = (r/q) * dq/dr
 * Physics: dimensionless measure of q-profile variation.
 *          s > 0: monotonic q profile (standard)
 *          s < 0: reversed shear (advanced tokamak)
 *          s ~ 0: low shear region (ITB formation possible)
 *
 * For typical H-mode: s ~ 1-3 in core, s_surface ~ 2-5
 */
double magnetic_shear(double r, double q, double dq_dr) {
    assert(q > 0.0);
    return (r / q) * dq_dr;
}

/**
 * magnetic_well_depth — Magnetic well depth
 *
 * Definition: U = (V'(0) - V'(psi)) / V'(0)
 * Physics: negative magnetic well stabilizes interchange modes.
 *          V'(psi) = dV/dpsi, the flux-surface differential volume.
 *          Positive U (magnetic hill): destabilizing.
 *          Negative U (magnetic well): stabilizing.
 * For tokamaks with circular cross-section: U ~ -r/R (stable).
 */
double magnetic_well_depth(double V_prime_0, double V_prime_psi) {
    if (fabs(V_prime_0) < 1e-30) return 0.0;
    return (V_prime_0 - V_prime_psi) / V_prime_0;
}

/**
 * sawtooth_period — Sawtooth oscillation period [s]
 *
 * Definition: tau_saw ~ (tau_R * tau_A)^(1/2)
 *          tau_R = mu0 * a^2 / eta (resistive diffusion time)
 *          tau_A = a / v_A (Alfven transit time)
 * Physics: periodic relaxation oscillation of core temperature
 *          caused by repeated Kadomtsev reconnection at q=1 surface.
 *          tau_saw determined by resistive diffusion to q<1,
 *          crash by fast magnetic reconnection (~100 microsec).
 * Source: Kadomtsev, Sov. J. Plasma Phys. 1, 389 (1975)
 *
 * For ITER: tau_saw ~ 1-10 s
 */
double sawtooth_period(double a, double eta, double v_A) {
    assert(v_A > 0.0);
    double tau_R = MU0 * a * a / eta;
    double tau_A = a / v_A;
    return sqrt(tau_R * tau_A);
}

/**
 * sawtooth_mixing_radius — Sawtooth mixing radius
 *
 * Definition: r_mix ~ r_q=1 (radius of q=1 surface)
 * Physics: sawtooth crash mixes plasma inside q=1 surface,
 *          flattening temperature and density.
 *          Mixing radius sets the affected plasma volume fraction.
 */
double sawtooth_mixing_radius(double r_q1) {
    return r_q1;
}

/* ================================================================
 * L8: Disruptions and NTM — 4 functions
 * ================================================================ */

/**
 * disruption_density_limit — Disruption density limit [m^-3]
 *
 * Definition: n_disrupt = margin * n_GW
 *          margin ~ 0.5-1.0 (operational safety factor)
 * Physics: exceeding Greenwald density often leads to
 *          radiative collapse and disruption.
 *          Murakami limit (pre-Greenwald) and Hugill limit
 *          are earlier versions of density limits.
 */
double disruption_density_limit(double Ip, double a, double margin) {
    double n_GW = greenwald_density_limit(Ip, a);
    return margin * n_GW;
}

/**
 * ntm_threshold_island_width — NTM seed island threshold [m]
 *
 * Definition: w_seed ~ rho_i * sqrt(beta_p)
 * Physics: neoclassical tearing modes (NTMs) require a seed island
 *          larger than the critical width to grow.
 *          Seed islands can come from sawteeth, ELMs, or error fields.
 *          Once seeded, NTMs grow due to perturbed bootstrap current
 *          loss within the island.
 * For ITER: w_seed ~ 1-2 cm
 */
double ntm_threshold_island_width(double rho_i, double beta_p) {
    return rho_i * sqrt(beta_p);
}

/**
 * halo_current_fraction — Halo current fraction during VDEs
 *
 * Definition: f_halo = I_halo / I_p0
 * Physics: during vertical displacement events (VDEs),
 *          plasma contacts wall and halo currents flow.
 *          I_halo ~ 10-30% of initial I_p.
 *          Produces large JxB forces on vessel components.
 *
 * @param I_halo  halo current [A]
 * @param I_p0    initial plasma current [A]
 * @return        halo current fraction
 */
double halo_current_fraction(double I_halo, double I_p0) {
    if (I_p0 <= 0.0) return 0.0;
    return I_halo / I_p0;
}

/**
 * power_exhaust_fraction — Power exhaust fraction
 *
 * Definition: f_rad = P_rad / P_heat
 * Physics: fraction of heating power radiated before reaching
 *          the divertor. High f_rad (>0.7-0.8) needed for
 *          detached divertor operation.
 *          P_rad = P_brem + P_cyc + P_line + P_impurity
 */
double power_exhaust_fraction(double P_rad, double P_heat) {
    if (P_heat <= 0.0) return 0.0;
    return P_rad / P_heat;
}

/**
 * divertor_heat_flux — Divertor peak heat flux [W/m^2]
 *
 * Definition: q_div = P_sol * sin(alpha) / (2*pi*R * lambda_q * f_exp)
 * Physics: power entering scrape-off layer (SOL) is conducted
 *          to divertor targets along field lines.
 *          lambda_q: power fall-off length at midplane (~1-5 mm)
 *          f_exp: flux expansion factor from midplane to target (5-10x)
 *          alpha: field line incident angle at target (~2-5 degrees)
 * For ITER: q_div ~ 5-10 MW/m^2 (detached), up to 20 MW/m^2 (attached)
 */
double divertor_heat_flux(double P_sol, double alpha_deg, double R,
                           double lambda_q, double f_exp) {
    if (lambda_q <= 0.0 || f_exp <= 0.0) return 0.0;
    double alpha_rad = alpha_deg * M_PI / 180.0;
    return P_sol * sin(alpha_rad) / (2.0 * M_PI * R * lambda_q * f_exp);
}