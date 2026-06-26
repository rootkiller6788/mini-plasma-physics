#ifndef DUSTY_CHARGING_H
#define DUSTY_CHARGING_H
/**
 * @file  dusty_charging.h
 * @brief Dust grain charging theory — OML and extensions.
 *
 * Orbital Motion Limited (OML) theory is the standard framework for
 * dust charging in low-pressure plasmas. The grain collects electron
 * and ion currents until the net current vanishes at the floating potential.
 *
 * References:
 *   Mott-Smith & Langmuir (1926), Phys. Rev. 28, 727
 *   Allen (1992), Physica Scripta 45, 497
 *   Goree (1994), Plasma Sources Sci. Technol. 3, 400
 *
 * L1-L4: OML charging equation, floating potential, charge fluctuation.
 */
#include "dusty_plasma.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================
 * L3 — OML Charging Current Functions
 * ================================================================ */

/**
 * @brief OML electron collection current to a spherical dust grain.
 *
 * I_e = -pi * a^2 * e * n_e * sqrt(8 * k_B * T_e / (pi * m_e))
 *       * exp(e * phi_s / (k_B * T_e))   for phi_s <= 0
 *
 * The exponential factor (Boltzmann factor) means electrons are
 * repelled by the negatively charged grain. For phi_s > 0 the
 * exponential factor saturates.
 *
 * Complexity: O(1) — closed-form expression.
 */
double dust_oml_electron_current(double a, double n_e, double T_e, double phi_s);

/**
 * @brief OML ion collection current to a spherical dust grain.
 *
 * I_i = pi * a^2 * e * n_i * sqrt(8 * k_B * T_i / (pi * m_i))
 *       * (1 - e * phi_s / (k_B * T_i))   for phi_s <= 0
 *
 * Ions are attracted by the negatively charged grain.
 * The orbital motion limit accounts for ion angular momentum.
 *
 * Complexity: O(1) — closed-form expression.
 */
double dust_oml_ion_current(double a, double n_i, double T_i, double m_i, double phi_s);

/**
 * @brief OML ion current with shifted Maxwellian (ion drift).
 *
 * When ions have a directed velocity u_i relative to the grain,
 * the current is modified by the drift. Uses the shifted-Maxwellian
 * correction factor f(u_i / v_thi).
 *
 * Ref: Hutchinson (2002), Plasma Phys. Control. Fusion 44, 1953
 */
double dust_oml_ion_current_drift(
    double a, double n_i, double T_i, double m_i, double phi_s, double u_i);

/**
 * @brief Photoelectric emission current from a dust grain.
 *
 * I_ph = pi * a^2 * J_ph * f(phi_s)
 *
 * where J_ph is the photoelectron flux density, which depends on:
 *   - Photon flux (from radiation source)
 *   - Photoelectric efficiency (material-dependent)
 *   - Work function of grain material
 *
 * For phi_s <= 0: f(phi_s) = exp(e * phi_s / (k_B * T_ph))
 * For phi_s > 0:  f(phi_s) = 1 + e * phi_s / (k_B * T_ph)
 *
 * Ref: Feuerbacher & Fitton (1972), J. Appl. Phys. 43, 1563
 */
double dust_photoelectric_current(
    double a, double phi_s, double J_ph, double T_ph);

/**
 * @brief Secondary electron emission current.
 *
 * I_se = pi * a^2 * e * n_e * sqrt(8 * k_B * T_e / (pi * m_e))
 *        * delta(E_e) * f_se(phi_s)
 *
 * where delta is the secondary electron yield (depends on primary
 * electron energy E_e ~ 2 k_B T_e) and f_se is the escape fraction.
 *
 * Ref: Meyer-Vernet (1982), Astron. Astrophys. 105, 98
 */
double dust_secondary_emission_current(
    double a, double n_e, double T_e, double phi_s,
    double se_yield_max, double E_max);

/**
 * @brief Thermionic emission current (Richardson-Dushman).
 *
 * I_th = 4 * pi * a^2 * A_R * T_d^2 * exp(-W / (k_B * T_d))
 *        * exp(min(0, e * phi_s / (k_B * T_d)))
 *
 * where A_R = 4*pi*m_e*e*k_B^2/h^3 (Richardson constant)
 * and W is the work function.
 *
 * Ref: Sodha & Guha (1971), "Physics of Colloidal Plasmas"
 */
double dust_thermionic_current(
    double a, double T_d, double W, double phi_s);

/**
 * @brief Compute all current components to a dust grain.
 *
 * Fills a DustChargingCurrent struct with the complete current budget.
 * This is the central charging diagnostic function.
 */
DustChargingCurrent dust_compute_charging_currents(
    double a, double n_e, double n_i, double T_e, double T_i,
    double m_i, double phi_s, ChargeModel model,
    double J_ph, double T_ph, double se_yield, double E_max,
    double T_d, double W);

/* ================================================================
 * L4 — Floating Potential and Equilibrium Charge
 * ================================================================ */

/**
 * @brief Compute net current to a grain at given surface potential.
 *
 * I_net(phi_s) = I_e(phi_s) + I_i(phi_s) + I_photo(phi_s) + ...
 *
 * The floating potential phi_f is the root of I_net(phi_f) = 0.
 */
double dust_net_current_to_grain(
    double a, double n_e, double n_i, double T_e, double T_i,
    double m_i, double phi_s, ChargeModel model,
    double J_ph, double T_ph, double se_yield, double E_max,
    double T_d, double W);

/**
 * @brief Solve for the floating potential phi_f using Newton-Raphson.
 *
 * Finds root of I_net(phi_f) = 0 in the range [phi_min, 0].
 * The floating potential is typically negative (electrons are more
 * mobile than ions).
 *
 * @param phi_init  Initial guess [V] (e.g., -2.5 * k_B * T_e / e)
 * @param tol       Convergence tolerance [V]
 * @param max_iter  Maximum Newton iterations
 * @return Floating potential phi_f [V]
 *
 * Complexity: O(max_iter) — Newton iterations, each O(1).
 */
double dust_floating_potential_solve(
    double a, double n_e, double n_i, double T_e, double T_i,
    double m_i, ChargeModel model,
    double J_ph, double T_ph, double se_yield, double E_max,
    double T_d, double W,
    double phi_init, double tol, int max_iter);

/**
 * @brief Compute equilibrium dust charge from floating potential.
 *
 * For a spherical grain: Q_d = 4 * pi * eps0 * a * phi_f
 * Charge number: Z_d = |Q_d| / e
 *
 * This is the spherical capacitor model valid when
 * a << lambda_D (grain radius much smaller than Debye length).
 */
double dust_equilibrium_charge(double a, double phi_f);

/**
 * @brief Compute equilibrium charge number Z_d.
 */
double dust_equilibrium_charge_number(double a, double phi_f);

/**
 * @brief Quick estimate of dust charge number using OML only.
 *
 * Uses the approximate formula:
 * Z_d ≈ 1400 * a_um * T_eV
 *
 * where a_um = grain radius in microns, T_eV = T_e in eV.
 * This is a rule-of-thumb valid for argon plasmas.
 */
double dust_charge_number_estimate(double a, double T_e);

/* ================================================================
 * L5 — Charge Fluctuation and Kinetics
 * ================================================================ */

/**
 * @brief Compute RMS charge fluctuation.
 *
 * <delta Z_d^2> = sqrt( (I_e^2 + I_i^2) * tau / e^2 )
 *
 * Charge fluctuates due to discrete nature of electron/ion collection.
 * The fluctuations are Poisson-distributed in the linear regime.
 *
 * Ref: Cui & Goree (1994), IEEE Trans. Plasma Sci. 22, 151
 */
double dust_charge_fluctuation_rms(
    double a, double n_e, double n_i, double T_e, double T_i,
    double m_i, double phi_f);

/**
 * @brief Compute charge relaxation (charging) timescale.
 *
 * tau_charge = C_grain / (dI_net/dphi_s)
 *
 * where C_grain = 4 * pi * eps0 * a is the grain capacitance.
 * This is the characteristic time for the grain charge to reach
 * equilibrium after a perturbation.
 *
 * Complexity: O(1).
 */
double dust_charge_relaxation_time(
    double a, double n_e, double n_i, double T_e, double T_i,
    double m_i, double phi_f);

/**
 * @brief Integrate the charging ODE dQ_d/dt = I_net(Q_d).
 *
 * Uses forward Euler with adaptive timestep.
 *
 * @param t_span  Integration time span [s]
 * @param dt      Timestep [s]
 * @param n_steps Number of steps (output array length)
 * @param Q_out   Output charge array [n_steps]
 * @param t_out   Output time array [n_steps]
 */
void dust_integrate_charge_dynamics(
    double a, double n_e, double n_i, double T_e, double T_i,
    double m_i, double Q0, double t_span, double dt,
    double *Q_out, double *t_out, int n_steps);

#ifdef __cplusplus
}
#endif

#endif /* DUSTY_CHARGING_H */
