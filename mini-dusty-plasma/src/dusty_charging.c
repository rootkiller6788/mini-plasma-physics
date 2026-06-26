/**
 * @file  dusty_charging.c
 * @brief Dust grain charging — OML theory and extensions.
 *
 * Implements the Orbital Motion Limited (OML) theory for dust
 * charging in collisionless plasmas. The central physical picture:
 * a spherical grain immersed in a plasma collects electrons and ions
 * until the net current vanishes at the floating potential.
 *
 * L3-L5: OML theory, floating potential, charge fluctuations.
 *
 * Governing equation:
 *   I_total(phi_s) = I_e(phi_s) + I_i(phi_s) + ... = 0
 *
 * where phi_s is the grain surface potential relative to plasma.
 *
 * References:
 *   Mott-Smith & Langmuir (1926), Phys. Rev. 28, 727
 *   Allen (1992), Physica Scripta 45, 497
 *   Goree (1994), Plasma Sources Sci. Technol. 3, 400
 *   Fortov et al. (2005), Physics Reports 421, 1-103
 */

#include "dusty_charging.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

/* ================================================================
 * L3 — OML Current Functions
 * ================================================================
 *
 * OML theory assumes:
 * 1. Collisionless ions/electrons (mean free path >> grain size)
 * 2. Maxwellian velocity distributions at infinity
 * 3. Spherical grain with absorbing surface
 * 4. Grain radius a << lambda_D (isolated grain)
 *
 * Under these assumptions, the collected current is determined
 * by conservation of energy and angular momentum of each particle
 * as it orbits the grain.
 */

double dust_oml_electron_current(double a, double n_e, double T_e, double phi_s)
{
    /* I_e(phi_s) = -pi * a^2 * e * n_e * v_th_e * I_e_norm(eta)
     *
     * where:
     *   v_th_e = sqrt(8 * k_B * T_e / (pi * m_e))
     *   eta = e * phi_s / (k_B * T_e)  (normalized potential)
     *
     * For phi_s <= 0 (repelled electrons):
     *   I_e_norm = exp(eta)
     *
     * For phi_s > 0 (attracted electrons, rare for dust):
     *   I_e_norm = 1 + eta
     *
     * The minus sign indicates electron current TO the grain
     * (electrons flow from plasma to grain surface).
     *
     * Mott-Smith & Langmuir (1926) Eq. (10) */
    if (a <= 0.0 || n_e <= 0.0 || T_e <= 0.0) return 0.0;

    double v_th_e = sqrt(8.0 * DUSTY_KB * T_e / (M_PI * DUSTY_ME));
    double area = M_PI * a * a;
    double I0 = area * DUSTY_EC * n_e * v_th_e;
    double eta = DUSTY_EC * phi_s / (DUSTY_KB * T_e);

    /* Normalized current */
    double norm;
    if (eta <= 0.0) {
        norm = exp(eta);
    } else {
        norm = 1.0 + eta;
    }

    return -I0 * norm;
}

double dust_oml_ion_current(double a, double n_i, double T_i, double m_i, double phi_s)
{
    /* I_i(phi_s) = pi * a^2 * e * n_i * v_th_i * I_i_norm(eta_i)
     *
     * where:
     *   v_th_i = sqrt(8 * k_B * T_i / (pi * m_i))
     *   eta_i = -e * phi_s / (k_B * T_i)  (positive for negative phi_s)
     *
     * For phi_s <= 0 (attracted ions):
     *   I_i_norm = 1 + eta_i
     *
     * For phi_s > 0 (repelled ions):
     *   I_i_norm = exp(-eta_i)
     *
     * The OML ion current is LARGER than the geometric current
     * because ions with finite impact parameter can be focused
     * onto the grain by the attractive potential.
     *
     * Allen (1992), Physica Scripta 45, 497 */
    if (a <= 0.0 || n_i <= 0.0 || T_i <= 0.0 || m_i <= 0.0) return 0.0;

    double v_th_i = sqrt(8.0 * DUSTY_KB * T_i / (M_PI * m_i));
    double area = M_PI * a * a;
    double I0 = area * DUSTY_EC * n_i * v_th_i;
    double eta_i = -DUSTY_EC * phi_s / (DUSTY_KB * T_i);

    double norm;
    if (eta_i >= 0.0) {
        /* Attracted ions */
        norm = 1.0 + eta_i;
    } else {
        /* Repelled ions (positive grain) */
        norm = exp(eta_i);
    }

    return I0 * norm;
}

double dust_oml_ion_current_drift(
    double a, double n_i, double T_i, double m_i, double phi_s, double u_i)
{
    /* Ion current with directed drift velocity u_i relative to grain.
     *
     * When ions have a net drift (e.g., in the presheath), the
     * current is enhanced. The shifted Maxwellian correction is:
     *
     * I_i(u_i) = I_i(0) * f(M)
     *
     * where M = u_i / c_s is the ion Mach number. For subsonic flow:
     * f(M) ≈ 1 + (2/3)*M^2 + O(M^4)
     *
     * For supersonic flow, the correction saturates.
     *
     * Hutchinson (2002), Plasma Phys. Control. Fusion 44, 1953 */
    if (a <= 0.0 || n_i <= 0.0 || T_i <= 0.0 || m_i <= 0.0) return 0.0;

    double I_base = dust_oml_ion_current(a, n_i, T_i, m_i, phi_s);

    /* Compute Mach number */
    double c_s = sqrt(2.0 * DUSTY_KB * T_i / m_i);
    if (c_s <= 0.0) return I_base;

    double M = u_i / c_s;

    /* Correction factor: empirical fit to PIC simulations.
     * Valid for 0 < M < 3. */
    double M2 = M * M;
    double correction;
    if (M < 0.5) {
        correction = 1.0 + 0.67 * M2;
    } else if (M < 2.0) {
        correction = 1.0 + 0.67 * M2 + 0.11 * M2 * M2;
    } else {
        /* Saturation at high Mach number */
        correction = 1.0 + 0.67 * 4.0 + 0.11 * 16.0;
    }

    return I_base * correction;
}

double dust_photoelectric_current(
    double a, double phi_s, double J_ph, double T_ph)
{
    /* Photoelectric emission: photons eject electrons from grain surface.
     *
     * I_ph = pi * a^2 * J_ph * f_ph(phi_s)
     *
     * where J_ph [A/m^2] is the photoelectron saturation current density
     * and T_ph [K] is the photoelectron temperature (~1-3 eV typically).
     *
     * For phi_s <= 0 (electrons escape freely):
     *   f_ph = exp(e*phi_s / (k_B*T_ph))
     *
     * For phi_s > 0 (electrons partially trapped):
     *   f_ph = 1 + e*phi_s / (k_B*T_ph)
     *
     * Photoelectric charging is dominant in:
     *   - Space plasmas (solar UV on dust)
     *   - Comet tails (solar wind interaction)
     *   - Laboratory experiments with UV sources
     *
     * Feuerbacher & Fitton (1972), J. Appl. Phys. 43, 1563 */
    if (a <= 0.0 || J_ph < 0.0) return 0.0;
    if (T_ph <= 0.0) T_ph = 11604.518; /* default 1 eV */

    double area = M_PI * a * a;
    double eta = DUSTY_EC * phi_s / (DUSTY_KB * T_ph);

    double norm;
    if (eta <= 0.0) {
        norm = exp(eta);
    } else {
        norm = 1.0 + eta;
    }

    return area * J_ph * norm;
}

double dust_secondary_emission_current(
    double a, double n_e, double T_e, double phi_s,
    double se_yield_max, double E_max)
{
    /* Secondary electron emission: incident primary electrons
     * knock out additional (secondary) electrons from the grain.
     *
     * I_se = pi*a^2 * e * n_e * v_th_e * delta(E) * f_se(phi_s)
     *
     * The secondary yield delta(E) follows the universal curve:
     *   delta(E) = delta_max * (E/E_max) * exp(1 - E/E_max)
     *
     * where delta_max is max yield and E_max is the primary energy
     * at which maximum yield occurs (typically 300-500 eV).
     *
     * For Maxwellian electrons, E_avg = 2*k_B*T_e ≈ 6 eV (for 3 eV T_e),
     * which is much less than E_max, so delta << delta_max.
     *
     * SE emission can cause grains to charge POSITIVELY in
     * high-temperature plasmas (e.g., fusion edge plasma).
     *
     * Meyer-Vernet (1982), Astron. Astrophys. 105, 98 */
    if (a <= 0.0 || n_e <= 0.0 || T_e <= 0.0) return 0.0;

    double v_th_e = sqrt(8.0 * DUSTY_KB * T_e / (M_PI * DUSTY_ME));
    double area = M_PI * a * a;
    double I0 = area * DUSTY_EC * n_e * v_th_e;

    /* Electron average energy ~ 2 k_B T_e */
    double E_avg = 2.0 * DUSTY_KB * T_e / DUSTY_EC; /* in eV */
    double x = E_avg / E_max;
    double delta = se_yield_max * x * exp(1.0 - x);

    /* Escape fraction: for phi_s < 0, all SE escape.
     * For phi_s > 0, fraction is reduced by recapture. */
    double eta = DUSTY_EC * phi_s / (DUSTY_KB * T_e);
    double escape;
    if (eta <= 0.0) {
        escape = 1.0;
    } else {
        escape = exp(-eta);
    }

    return I0 * delta * escape;
}

double dust_thermionic_current(
    double a, double T_d, double W, double phi_s)
{
    /* Thermionic emission: electrons thermally emitted from hot grain.
     *
     * Richardson-Dushman equation:
     * I_th = 4*pi*a^2 * A_R * T_d^2 * exp(-W/(k_B*T_d)) * g(phi_s)
     *
     * where A_R = 4*pi*m_e*e*k_B^2/h^3 ≈ 1.20×10^6 A/(m^2 K^2)
     *
     * For phi_s <= 0: g(phi_s) = 1 (all emitted electrons escape)
     * For phi_s > 0:  g(phi_s) = exp(-e*phi_s/(k_B*T_d))
     *
     * Important for:
     *   - Incandescent dust (e.g., in flames or near hot walls)
     *   - Fusion: dust heated by plasma to >2000 K
     *   - Meteor trails
     *
     * Sodha & Guha (1971), Physics of Colloidal Plasmas */
    if (a <= 0.0 || T_d <= 0.0 || W <= 0.0) return 0.0;

    /* Richardson constant: A_R = 4*pi*m_e*e*k_B^2 / h^3 */
    double A_R = 4.0 * M_PI * DUSTY_ME * DUSTY_EC
                 * DUSTY_KB * DUSTY_KB
                 / (DUSTY_H * DUSTY_H * DUSTY_H);

    double area = 4.0 * M_PI * a * a; /* total surface area */
    double J_richardson = A_R * T_d * T_d
                          * exp(-W * DUSTY_EC / (DUSTY_KB * T_d));

    double escape;
    if (phi_s <= 0.0) {
        escape = 1.0;
    } else {
        double eta = DUSTY_EC * phi_s / (DUSTY_KB * T_d);
        escape = exp(-eta);
    }

    return area * J_richardson * escape;
}

DustChargingCurrent dust_compute_charging_currents(
    double a, double n_e, double n_i, double T_e, double T_i,
    double m_i, double phi_s, ChargeModel model,
    double J_ph, double T_ph, double se_yield, double E_max,
    double T_d, double W)
{
    /* Assemble all current components into a diagnostic struct.
     *
     * This provides a complete current budget at a given phi_s,
     * useful for understanding which charging mechanism dominates. */
    DustChargingCurrent cc;
    cc.I_e = dust_oml_electron_current(a, n_e, T_e, phi_s);
    cc.I_i = dust_oml_ion_current(a, n_i, T_i, m_i, phi_s);
    cc.I_photo = 0.0;
    cc.I_secondary = 0.0;
    cc.I_thermionic = 0.0;

    if (model == CHARGE_MODEL_OML_PLUS || model == CHARGE_MODEL_FULL) {
        cc.I_photo = dust_photoelectric_current(a, phi_s, J_ph, T_ph);
    }
    if (model == CHARGE_MODEL_FULL) {
        cc.I_secondary = dust_secondary_emission_current(
            a, n_e, T_e, phi_s, se_yield, E_max);
        cc.I_thermionic = dust_thermionic_current(a, T_d, W, phi_s);
    }

    cc.I_total = cc.I_e + cc.I_i + cc.I_photo
               + cc.I_secondary + cc.I_thermionic;

    return cc;
}

/* ================================================================
 * L4 — Floating Potential Solver (Newton-Raphson)
 * ================================================================ */

double dust_net_current_to_grain(
    double a, double n_e, double n_i, double T_e, double T_i,
    double m_i, double phi_s, ChargeModel model,
    double J_ph, double T_ph, double se_yield, double E_max,
    double T_d, double W)
{
    DustChargingCurrent cc = dust_compute_charging_currents(
        a, n_e, n_i, T_e, T_i, m_i, phi_s, model,
        J_ph, T_ph, se_yield, E_max, T_d, W);
    return cc.I_total;
}

/**
 * @brief Compute dI_dphi numerically (central difference).
 *
 * Used internally by the Newton-Raphson solver.
 *
 * @param h  Step size for numerical derivative [V]
 */
static double dust_net_current_derivative(
    double a, double n_e, double n_i, double T_e, double T_i,
    double m_i, double phi_s, ChargeModel model,
    double J_ph, double T_ph, double se_yield, double E_max,
    double T_d, double W, double h)
{
    double I_plus = dust_net_current_to_grain(
        a, n_e, n_i, T_e, T_i, m_i, phi_s + h, model,
        J_ph, T_ph, se_yield, E_max, T_d, W);
    double I_minus = dust_net_current_to_grain(
        a, n_e, n_i, T_e, T_i, m_i, phi_s - h, model,
        J_ph, T_ph, se_yield, E_max, T_d, W);
    return (I_plus - I_minus) / (2.0 * h);
}

double dust_floating_potential_solve(
    double a, double n_e, double n_i, double T_e, double T_i,
    double m_i, ChargeModel model,
    double J_ph, double T_ph, double se_yield, double E_max,
    double T_d, double W,
    double phi_init, double tol, int max_iter)
{
    /* Newton-Raphson root finding for I_net(phi_f) = 0.
     *
     * The floating potential for dust in a typical laboratory plasma
     * is negative (phi_f ~ -2.5 * k_B*T_e/e) because electron thermal
     * velocity >> ion thermal velocity.
     *
     * Algorithm:
     *   phi_{n+1} = phi_n - I_net(phi_n) / I_net'(phi_n)
     *
     * Convergence is quadratic near the root. A good initial guess
     * (phi_init = -2.5*k_B*T_e/e) ensures convergence in 5-10 iterations.
     *
     * Safety: bounds phi to [-10*k_B*T_e/e, +k_B*T_e/e] */
    if (a <= 0.0 || n_e <= 0.0 || n_i <= 0.0 || T_e <= 0.0 || T_i <= 0.0) {
        return 0.0;
    }

    double phi = phi_init;
    double phi_min = -10.0 * DUSTY_KB * T_e / DUSTY_EC;
    double phi_max = DUSTY_KB * T_e / DUSTY_EC;
    double dh = tol * 10.0;
    if (dh < 1e-6) dh = 1e-6;

    int iter;
    for (iter = 0; iter < max_iter; iter++) {
        double I_net = dust_net_current_to_grain(
            a, n_e, n_i, T_e, T_i, m_i, phi, model,
            J_ph, T_ph, se_yield, E_max, T_d, W);
        double dI = dust_net_current_derivative(
            a, n_e, n_i, T_e, T_i, m_i, phi, model,
            J_ph, T_ph, se_yield, E_max, T_d, W, dh);

        if (fabs(dI) < 1e-30) {
            /* Zero derivative — fall back to bisection step */
            if (I_net > 0.0) phi_max = phi;
            else phi_min = phi;
            phi = 0.5 * (phi_min + phi_max);
        } else {
            double dphi = -I_net / dI;
            /* Step size limiting */
            if (fabs(dphi) > fabs(phi * 0.5)) {
                dphi = (dphi > 0.0 ? 1.0 : -1.0) * fabs(phi * 0.5);
            }
            phi += dphi;
        }

        /* Enforce bounds */
        if (phi < phi_min) phi = phi_min;
        if (phi > phi_max) phi = phi_max;

        /* Check convergence */
        if (fabs(I_net) < tol * fabs(dust_oml_electron_current(a, n_e, T_e, 0.0))) {
            return phi;
        }
    }

    return phi;
}

double dust_equilibrium_charge(double a, double phi_f)
{
    /* Q_d = C_grain * phi_f
     *
     * where C_grain = 4 * pi * eps0 * a is the capacitance of a
     * spherical grain in vacuum (valid when a << lambda_D).
     *
     * For a conductor, the charge resides on the surface and
     * the potential is constant. For dielectric grains (most
     * experimental dust), this is still a good approximation
     * because the charge distribution is effectively dipolar. */
    if (a <= 0.0) return 0.0;
    return 4.0 * M_PI * DUSTY_EPS0 * a * phi_f;
}

double dust_equilibrium_charge_number(double a, double phi_f)
{
    /* Z_d = |Q_d| / e
     *
     * Z_d is dimensionless and typically in the range 10^2-10^4
     * for micron-sized grains in laboratory plasmas. */
    double Q_d = dust_equilibrium_charge(a, phi_f);
    return fabs(Q_d) / DUSTY_EC;
}

double dust_charge_number_estimate(double a, double T_e)
{
    /* Rule-of-thumb estimate for OML charging.
     *
     * Z_d ≈ 1400 * a_um * T_eV
     *
     * where a_um is grain radius in microns and T_eV is T_e in eV.
     * This comes from balancing electron and ion OML currents
     * with m_i/m_e ≈ 7.3×10^4 (argon).
     *
     * Example: a = 5 um, T_e = 3 eV → Z_d ≈ 1400*5*3 = 21000. */
    if (a <= 0.0 || T_e <= 0.0) return 0.0;
    double a_um = a / 1.0e-6;
    double T_eV = T_e / DUSTY_EV_IN_K;
    return 1400.0 * a_um * T_eV;
}

/* ================================================================
 * L5 — Charge Fluctuation Statistics
 * ================================================================ */

double dust_charge_fluctuation_rms(
    double a, double n_e, double n_i, double T_e, double T_i,
    double m_i, double phi_f)
{
    /* RMS charge fluctuation from discrete collection statistics.
     *
     * <delta Z_d^2> = (|I_e| + |I_i|) * tau_charge / e
     *
     * The fluctuations are Poissonian because electron and ion
     * collection events are independent and random (shot noise).
     *
     * In the linear regime (small fluctuations):
     * sqrt(<delta Z_d^2>) / Z_d ~ 1/sqrt(Z_d) ~ few percent.
     *
     * Cui & Goree (1994), IEEE Trans. Plasma Sci. 22, 151 */
    if (a <= 0.0 || n_e <= 0.0 || n_i <= 0.0) return 0.0;

    double I_e = fabs(dust_oml_electron_current(a, n_e, T_e, phi_f));
    double I_i = fabs(dust_oml_ion_current(a, n_i, T_i, m_i, phi_f));

    /* Total collection rate (events per second) */
    double rate = (I_e + I_i) / DUSTY_EC;

    /* Charge relaxation time */
    double tau = dust_charge_relaxation_time(
        a, n_e, n_i, T_e, T_i, m_i, phi_f);

    /* RMS fluctuation in charge number */
    /* <delta Z^2> = rate * tau */
    double var_Z = rate * tau;
    if (var_Z < 0.0) var_Z = 0.0;

    return sqrt(var_Z);
}

double dust_charge_relaxation_time(
    double a, double n_e, double n_i, double T_e, double T_i,
    double m_i, double phi_f)
{
    /* tau_charge = C_grain / |dI/dphi|_phi_f
     *
     * This is the RC time constant of the grain-plasma system.
     * C_grain = 4*pi*eps0*a is the grain capacitance.
     * dI/dphi is the differential conductivity at the floating potential.
     *
     * For typical laboratory conditions: tau_charge ~ 10^-7 s,
     * which is much shorter than the dust plasma period (~10^-2 s),
     * meaning the charge adjusts quasi-instantaneously to perturbations.
     *
     * However, for very small grains (a < 10 nm) or very low density
     * (n_e < 10^12 m^-3), tau_charge can become significant. */
    if (a <= 0.0) return 0.0;

    double C_grain = 4.0 * M_PI * DUSTY_EPS0 * a;
    double dh = 1e-3; /* mV step for derivative */
    double dI = dust_net_current_derivative(
        a, n_e, n_i, T_e, T_i, m_i, phi_f, CHARGE_MODEL_OML,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, dh);

    if (fabs(dI) < 1e-30) {
        /* Fallback: estimate from electron current */
        double I_e = fabs(dust_oml_electron_current(a, n_e, T_e, phi_f));
        double R_eff = fabs(DUSTY_KB * T_e / (DUSTY_EC * I_e));
        return C_grain * R_eff;
    }

    return C_grain / fabs(dI);
}

void dust_integrate_charge_dynamics(
    double a, double n_e, double n_i, double T_e, double T_i,
    double m_i, double Q0, double t_span, double dt,
    double *Q_out, double *t_out, int n_steps)
{
    /* Integrate dQ_d/dt = I_net(Q_d) using forward Euler.
     *
     * The charging ODE is stiff (tau_charge ~ 10^-7 s), so a small
     * dt is needed. In practice, using variable-step RK4 is better,
     * but forward Euler with adaptive dt < tau_charge/10 suffices
     * for educational purposes.
     *
     * @param Q0   Initial charge [C]
     * @param dt   Time step [s] — should be < tau_charge/10
     * @param Q_out Output array with charge at each step
     * @param t_out Output array with time at each step
     * @param n_steps Number of output points */
    if (!Q_out || !t_out || n_steps <= 0) return;

    /* Compute phi_s from Q using spherical capacitor model */
    /* phi_s = Q / (4*pi*eps0*a) */
    double C_grain = 4.0 * M_PI * DUSTY_EPS0 * a;
    if (C_grain <= 0.0) {
        for (int i = 0; i < n_steps; i++) {
            t_out[i] = i * dt;
            Q_out[i] = Q0;
        }
        return;
    }

    double Q = Q0;
    for (int i = 0; i < n_steps; i++) {
        t_out[i] = i * dt;
        Q_out[i] = Q;

        double phi_s = Q / C_grain;
        double I_net = dust_net_current_to_grain(
            a, n_e, n_i, T_e, T_i, m_i, phi_s, CHARGE_MODEL_OML,
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

        /* Forward Euler step */
        Q += I_net * dt;
    }
}
