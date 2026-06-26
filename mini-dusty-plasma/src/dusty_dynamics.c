/**
 * @file  dusty_dynamics.c
 * @brief Dust grain dynamics — Langevin equation, trajectory integration.
 *
 * Dust grain motion in a plasma is governed by the Langevin equation:
 *   m_d * dv/dt = F_det + F_rand
 *
 * where F_det includes deterministic forces (electric, ion drag,
 * gravity, neutral drag) and F_rand represents random thermal
 * kicks from collisions with neutrals.
 *
 * L5: Langevin dynamics, RK4 integration, dust heating.
 *
 * References:
 *   Lemons & Gythiel (1997), Am. J. Phys. 65, 1079 — Langevin review
 *   Hou & Piel (2008), Phys. Plasmas 15, 073707 — dust in sheath
 *   Ivlev et al. (2012), "Complex and Dusty Plasmas", CRC Press
 */

#include "dusty_plasma.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* ================================================================
 * L5 — Langevin Equation of Motion
 * ================================================================ */

/**
 * @brief Deterministic acceleration on a dust grain.
 *
 * a = (F_electric + F_ion_drag + F_neutral_drag + F_gravity) / m_d
 *
 * @param v     Grain velocity [m/s]
 * @param E     Electric field [V/m]
 * @param Q_d   Grain charge [C]
 * @param m_d   Grain mass [kg]
 * @param nu_dn Dust-neutral collision frequency [1/s]
 * @param v_n   Neutral flow velocity [m/s]
 * @param g     Gravitational acceleration [m/s^2]
 * @param a_out Output acceleration [m/s^2] (length 3)
 */
void dust_langevin_deterministic_accel(
    const double v[3], const double E[3], double Q_d, double m_d,
    double nu_dn, const double v_n[3], double g,
    double a_out[3])
{
    /* a = (Q_d*E - m_d*nu_dn*(v - v_n) + m_d*g) / m_d */
    if (!v || !E || !v_n || !a_out || m_d <= 0.0) {
        if (a_out) a_out[0] = a_out[1] = a_out[2] = 0.0;
        return;
    }
    double inv_m = 1.0 / m_d;
    a_out[0] = inv_m * (Q_d * E[0] - m_d * nu_dn * (v[0] - v_n[0]));
    a_out[1] = inv_m * (Q_d * E[1] - m_d * nu_dn * (v[1] - v_n[1]));
    a_out[2] = inv_m * (Q_d * E[2] - m_d * nu_dn * (v[2] - v_n[2]) - m_d * g);
}

/**
 * @brief Random force amplitude from fluctuation-dissipation theorem.
 *
 * <F_rand(t) F_rand(t')> = 2 * m_d * nu_dn * k_B * T_n * delta(t - t')
 *
 * The random force amplitude (standard deviation per component):
 * sigma_F = sqrt(2 * m_d * nu_dn * k_B * T_n / dt)
 *
 * @param dt   Integration timestep [s]
 * @return Random force amplitude [N]
 */
double dust_random_force_amplitude(
    double m_d, double nu_dn, double T_n, double dt)
{
    /* Fluctuation-dissipation theorem:
     * The random force must satisfy the Einstein relation.
     * For a given timestep dt, the impulse is:
     * sigma_F * sqrt(dt) = sqrt(2*m_d*nu_dn*k_B*T_n) */
    if (m_d <= 0.0 || nu_dn <= 0.0 || T_n <= 0.0 || dt <= 0.0) return 0.0;
    double D = m_d * nu_dn * DUSTY_KB * T_n;
    return sqrt(2.0 * D / dt);
}

/**
 * @brief Box-Muller transform for Gaussian random numbers.
 *
 * Generates two independent N(0,1) random numbers from two
 * uniform [0,1] random numbers.
 *
 * This is essential for Langevin dynamics where the random
 * force must be Gaussian-distributed (central limit theorem
 * for many independent collisions).
 *
 * Box & Muller (1958), Ann. Math. Statist. 29, 610
 */
void dust_box_muller(double u1, double u2, double *n1, double *n2)
{
    /* Box-Muller: n1 = sqrt(-2*ln(u1)) * cos(2*pi*u2)
     *             n2 = sqrt(-2*ln(u1)) * sin(2*pi*u2) */
    if (!n1 || !n2) return;
    if (u1 <= 0.0) u1 = 1e-30;
    double r = sqrt(-2.0 * log(u1));
    double theta = 2.0 * M_PI * u2;
    *n1 = r * cos(theta);
    *n2 = r * sin(theta);
}

/* ================================================================
 * L5 — Trajectory Integration (RK4)
 * ================================================================ */

/**
 * @brief Single RK4 step for dust trajectory.
 *
 * Integrates the equations of motion:
 * dx/dt = v
 * dv/dt = F(x, v) / m_d
 *
 * using the classic 4th-order Runge-Kutta method.
 *
 * @param dt  Timestep [s]
 * @param E   Electric field callback: E(x,t) -> E[3]
 * @param ctx User context for E field callback
 */
void dust_rk4_step(
    double x[3], double v[3], double t, double dt,
    double Q_d, double m_d, double nu_dn, double g,
    void (*E_field)(double t, const double x[3], double E[3], void *ctx),
    void *ctx)
{
    /* Standard RK4 for 6-D phase space (position + velocity).
     *
     * k1 = dt * f(x, v)
     * k2 = dt * f(x + k1_x/2, v + k1_v/2)
     * k3 = dt * f(x + k2_x/2, v + k2_v/2)
     * k4 = dt * f(x + k3_x, v + k3_v)
     *
     * x_{n+1} = x_n + (k1_x + 2*k2_x + 2*k3_x + k4_x)/6 */
    if (!x || !v || !E_field) return;

    double E[3];
    double v_n[3] = {0.0, 0.0, 0.0};
    double a[3];
    double k1x[3], k1v[3], k2x[3], k2v[3], k3x[3], k3v[3], k4x[3], k4v[3];

    /* --- k1 --- */
    E_field(t, x, E, ctx);
    dust_langevin_deterministic_accel(v, E, Q_d, m_d, nu_dn, v_n, g, a);
    for (int j = 0; j < 3; j++) {
        k1v[j] = dt * a[j];
        k1x[j] = dt * v[j];
    }

    /* --- k2 --- */
    {
        double x2[3], v2[3];
        for (int j = 0; j < 3; j++) {
            x2[j] = x[j] + 0.5 * k1x[j];
            v2[j] = v[j] + 0.5 * k1v[j];
        }
        E_field(t + 0.5*dt, x2, E, ctx);
        dust_langevin_deterministic_accel(v2, E, Q_d, m_d, nu_dn, v_n, g, a);
        for (int j = 0; j < 3; j++) {
            k2v[j] = dt * a[j];
            k2x[j] = dt * v2[j];
        }
    }

    /* --- k3 --- */
    {
        double x3[3], v3[3];
        for (int j = 0; j < 3; j++) {
            x3[j] = x[j] + 0.5 * k2x[j];
            v3[j] = v[j] + 0.5 * k2v[j];
        }
        E_field(t + 0.5*dt, x3, E, ctx);
        dust_langevin_deterministic_accel(v3, E, Q_d, m_d, nu_dn, v_n, g, a);
        for (int j = 0; j < 3; j++) {
            k3v[j] = dt * a[j];
            k3x[j] = dt * x3[j];
        }
    }

    /* --- k4 --- */
    {
        double x4[3], v4[3];
        for (int j = 0; j < 3; j++) {
            x4[j] = x[j] + k3x[j];
            v4[j] = v[j] + k3v[j];
        }
        E_field(t + dt, x4, E, ctx);
        dust_langevin_deterministic_accel(v4, E, Q_d, m_d, nu_dn, v_n, g, a);
        for (int j = 0; j < 3; j++) {
            k4v[j] = dt * a[j];
            k4x[j] = dt * x4[j];
        }
    }

    /* --- Update --- */
    double inv6 = 1.0 / 6.0;
    for (int j = 0; j < 3; j++) {
        x[j] += inv6 * (k1x[j] + 2.0*k2x[j] + 2.0*k3x[j] + k4x[j]);
        v[j] += inv6 * (k1v[j] + 2.0*k2v[j] + 2.0*k3v[j] + k4v[j]);
    }
}

/* ================================================================
 * L5 — Dust Kinetic Energy and Effective Temperature
 * ================================================================ */

/**
 * @brief Compute kinetic energy of a dust grain.
 * E_k = (1/2) * m_d * v^2
 */
double dust_kinetic_energy(double m_d, const double v[3])
{
    if (!v || m_d <= 0.0) return 0.0;
    double v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
    return 0.5 * m_d * v2;
}

/**
 * @brief Effective temperature from kinetic energy (3D equipartition).
 * T_eff = (2/3) * E_k / k_B
 */
double dust_effective_temperature(double E_k)
{
    if (E_k <= 0.0) return 0.0;
    return (2.0 / 3.0) * E_k / DUSTY_KB;
}
