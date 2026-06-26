#ifndef DUSTY_FORCES_H
#define DUSTY_FORCES_H
/**
 * @file  dusty_forces.h
 * @brief Forces on dust grains in a plasma environment.
 *
 * Understanding the force balance is crucial for predicting dust
 * transport, levitation, and structure formation. The dominant forces
 * vary by environment:
 *
 * Laboratory (RF discharge): electric + ion drag + neutral drag + gravity
 * Space (Saturn's rings): gravity + radiation pressure + Lorentz force
 * Fusion devices: thermophoretic + ion drag + electric
 *
 * Refs:
 *   Barnes et al. (1992), Phys. Rev. Lett. 68, 313 — ion drag
 *   Nitter (1996), Plasma Sources Sci. Technol. 5, 93 — force review
 *   Khrapak et al. (2002), Phys. Rev. E 66, 046414 — ion drag refined
 */
#include "dusty_plasma.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================
 * L2 — Electrostatic and Gravitational Forces
 * ================================================================ */

/**
 * @brief Gravitational force on a dust grain.
 *
 * F_g = m_d * g
 *
 * @param m_d  Grain mass [kg]
 * @param g    Local gravitational acceleration [m/s^2]
 * @return Force magnitude [N]
 */
double dust_gravity_force(double m_d, double g);

/**
 * @brief Buoyancy-corrected gravitational force (for large grains).
 *
 * F_g_eff = (m_d - m_displaced) * g = m_d * (1 - rho_gas/rho_dust) * g
 *
 * Buoyancy is negligible for micron-sized grains in gas at ~1 Pa,
 * but can matter at atmospheric pressure.
 */
double dust_gravity_force_buoyancy(
    double m_d, double rho_dust, double rho_gas, double g);

/**
 * @brief Electrostatic force on a charged grain in electric field.
 *
 * F_E = Q_d * E
 *
 * @param Q_d  Dust charge [C]
 * @param E    Electric field magnitude [V/m]
 * @return Force magnitude [N]
 */
double dust_electric_force(double Q_d, double E);

/**
 * @brief Electrostatic force vector (3D).
 *
 * F_E = Q_d * E_vec
 *
 * @param E_vec  Electric field vector [V/m] (length 3)
 * @param F_out  Output force vector [N] (length 3)
 */
void dust_electric_force_vector(
    double Q_d, const double E_vec[3], double F_out[3]);

/**
 * @brief Lorentz force on a moving charged grain.
 *
 * F_L = Q_d * (v × B)
 *
 * @param v_vec  Velocity vector [m/s] (length 3)
 * @param B_vec  Magnetic field vector [T] (length 3)
 * @param F_out  Output Lorentz force vector [N] (length 3)
 */
void dust_lorentz_force_vector(
    double Q_d, const double v_vec[3], const double B_vec[3],
    double F_out[3]);

/* ================================================================
 * L2 — Ion Drag Force
 * ================================================================ */

/**
 * @brief Ion drag force — collection (direct impact) component.
 *
 * F_col = pi * b_c^2 * n_i * m_i * v_i * u_i
 *
 * where b_c = a * sqrt(1 - 2*e*phi_s/(m_i*v_i^2)) is the collection
 * impact parameter, and u_i is the ion drift velocity.
 *
 * This is from ions that physically hit the grain surface.
 *
 * Ref: Barnes et al. (1992), Phys. Rev. Lett. 68, 313
 */
double dust_ion_drag_collection(
    double a, double n_i, double m_i, double u_i, double phi_s);

/**
 * @brief Ion drag force — orbit (Coulomb scattering) component.
 *
 * F_orb = 4 * pi * b_90^2 * n_i * m_i * v_i * u_i * ln(Lambda_id)
 *
 * where b_90 = Z_d * e^2 / (4*pi*eps0 * m_i * v_i^2) is the
 * 90-degree scattering impact parameter for ion-dust collisions.
 *
 * This component comes from ions that are deflected by the
 * Coulomb field of the grain but do not hit it.
 */
double dust_ion_drag_orbit(
    double a, double n_i, double m_i, double u_i,
    double Z_d, double T_i, double lambda_D);

/**
 * @brief Total ion drag force (collection + orbit).
 *
 * F_id_total = F_collection + F_orbit
 *
 * Ion drag is typically the dominant horizontal force in RF discharges,
 * responsible for dust void formation and vortex motion.
 */
double dust_ion_drag_total(
    double a, double n_i, double m_i, double u_i,
    double phi_s, double Z_d, double T_i, double lambda_D);

/**
 * @brief Ion drag force with thermal corrections (Khrapak model).
 *
 * Refined formula including thermal ion effects:
 * F_id = sqrt(2*pi) * a^2 * n_i * m_i * v_thi * u_i
 *        * [ sqrt(pi)/4 * (b_c/a)^2 + (Z_d*e^2/(4*pi*eps0*a*k_B*T_i))^2
 *        * ln(1 + (lambda_D*(u_i^2+v_thi^2/2)/(a*k_B*T_i))^2) ]
 *
 * Ref: Khrapak et al. (2002), Phys. Rev. E 66, 046414
 */
double dust_ion_drag_khrapak(
    double a, double n_i, double m_i, double T_i,
    double u_i, double Z_d, double lambda_D);

/* ================================================================
 * L2 — Neutral Drag Force
 * ================================================================ */

/**
 * @brief Neutral drag force (Epstein regime).
 *
 * F_nd = -m_d * nu_dn * (v_d - v_n)
 *
 * where nu_dn is the dust-neutral collision frequency and
 * (v_d - v_n) is the relative velocity between dust and neutrals.
 *
 * @param m_d     Dust mass [kg]
 * @param nu_dn   Dust-neutral collision frequency [1/s]
 * @param v_rel   Relative velocity magnitude [m/s]
 * @return Force magnitude [N]
 */
double dust_neutral_drag_force(
    double m_d, double nu_dn, double v_rel);

/**
 * @brief Stokes drag force (continuum regime, Kn << 1).
 *
 * F_stokes = 6 * pi * eta * a * (v_n - v_d)
 *
 * Valid when the grain size >> gas mean free path (high pressure).
 *
 * @param a    Grain radius [m]
 * @param eta  Gas dynamic viscosity [Pa·s]
 * @param v_rel Relative velocity magnitude [m/s]
 */
double dust_stokes_drag_force(double a, double eta, double v_rel);

/**
 * @brief Interpolation between Epstein and Stokes regimes.
 *
 * Uses the Cunningham slip correction factor:
 * F = F_stokes / C_c(Kn)
 * C_c(Kn) = 1 + Kn * (A1 + A2 * exp(-A3/Kn))
 *
 * where A1=1.257, A2=0.400, A3=1.10 (Millikan's values for air).
 */
double dust_general_drag_force(
    double a, double eta, double v_rel, double Kn);

/* ================================================================
 * L4 — Thermophoretic Force
 * ================================================================ */

/**
 * @brief Thermophoretic force on a dust grain.
 *
 * F_th = -(32/15) * (a^2 / v_thn) * kappa_gas * grad_T
 *        * (1 + 5*pi/32*(1-alpha))
 *
 * where kappa_gas is the gas thermal conductivity, grad_T is the
 * temperature gradient, and alpha is the accommodation coefficient.
 *
 * Thermophoresis pushes grains from hot to cold regions.
 * This is important in fusion device dust transport.
 *
 * Ref: Talbot et al. (1980), J. Fluid Mech. 101, 737
 */
double dust_thermophoretic_force(
    double a, double T_n, double m_n,
    double kappa_gas, double grad_T, double alpha);

/* ================================================================
 * L4 — Radiation Pressure Force
 * ================================================================ */

/**
 * @brief Radiation pressure force on a dust grain.
 *
 * F_rad = pi * a^2 * Q_pr * I / c
 *
 * where Q_pr is the radiation pressure efficiency factor (~1 for
 * a >> wavelength, ~(a/lambda)^4 for a << wavelength),
 * I is the radiation intensity [W/m^2], and c is the speed of light.
 *
 * In astrophysical contexts (e.g., comet tails, protoplanetary disks),
 * radiation pressure can dominate over gravity for sub-micron grains.
 */
double dust_radiation_pressure_force(
    double a, double Q_pr, double I_radiation);

/* ================================================================
 * L5 — Force Balance Solver
 * ================================================================ */

/**
 * @brief Compute all forces on a grain and fill DustForceResult.
 *
 * This is the main diagnostic for understanding grain dynamics.
 * Combines gravity, electric, ion drag, neutral drag, thermophoretic,
 * and radiation pressure contributions.
 */
DustForceResult dust_compute_all_forces(
    double m_d, double Q_d, double g,
    double E_field, const double *B_field,
    const double *v_dust, const double *v_neutral,
    double a, double n_i, double m_i, double u_i,
    double phi_s, double Z_d, double T_i, double lambda_D,
    double nu_dn, double eta_gas, double v_rel,
    double kappa_gas, double grad_T, double alpha_th,
    double I_rad, double Q_pr);

/**
 * @brief Find the equilibrium levitation height in sheath.
 *
 * Solves F_total_z(z) = 0 for z, where the electrostatic force
 * upward balances gravity and ion drag downward.
 *
 * Uses a simple exponential sheath field model:
 * E(z) = E_0 * exp(-z / lambda_D)
 *
 * @return Equilibrium height [m] in sheath, or -1 if no solution.
 */
double dust_levitation_height(
    double m_d, double Q_d, double g,
    double E_0_sheath, double lambda_D,
    double F_ion_drag_z, double z_min, double z_max, double tol);

#ifdef __cplusplus
}
#endif

#endif /* DUSTY_FORCES_H */
