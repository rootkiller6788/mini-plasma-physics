/**
 * @file  dusty_forces.c
 * @brief Forces on dust grains — ion drag, neutral drag, thermophoresis.
 *
 * Force balance determines dust equilibrium positions and dynamics.
 * In a typical RF discharge sheath:
 *   - Gravity pulls down (~10^-14 N)
 *   - Electric force pushes up (~10^-13 N, Q_d * E_sheath)
 *   - Ion drag pushes down or sideways (~10^-14 N)
 *   - Neutral drag damps motion (friction)
 *
 * L2-L5: All force components, force balance solver, levitation.
 *
 * References:
 *   Barnes et al. (1992), Phys. Rev. Lett. 68, 313
 *   Khrapak et al. (2002), Phys. Rev. E 66, 046414
 *   Nitter (1996), Plasma Sources Sci. Technol. 5, 93
 */

#include "dusty_forces.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* ================================================================
 * L2 — Gravity and Electrostatic Forces
 * ================================================================ */

double dust_gravity_force(double m_d, double g)
{
    if (m_d <= 0.0) return 0.0;
    return m_d * g;
}

double dust_gravity_force_buoyancy(
    double m_d, double rho_dust, double rho_gas, double g)
{
    /* F_g_eff = m_d * (1 - rho_gas/rho_dust) * g */
    if (m_d <= 0.0 || rho_dust <= 0.0) return 0.0;
    double buoyancy_factor = 1.0 - rho_gas / rho_dust;
    return m_d * buoyancy_factor * g;
}

double dust_electric_force(double Q_d, double E)
{
    return Q_d * E;
}

void dust_electric_force_vector(
    double Q_d, const double E_vec[3], double F_out[3])
{
    if (!E_vec || !F_out) return;
    F_out[0] = Q_d * E_vec[0];
    F_out[1] = Q_d * E_vec[1];
    F_out[2] = Q_d * E_vec[2];
}

void dust_lorentz_force_vector(
    double Q_d, const double v_vec[3], const double B_vec[3],
    double F_out[3])
{
    /* F_L = Q_d * (v x B)
     * F_x = Q_d*(v_y*B_z - v_z*B_y) etc. */
    if (!v_vec || !B_vec || !F_out) return;
    F_out[0] = Q_d * (v_vec[1]*B_vec[2] - v_vec[2]*B_vec[1]);
    F_out[1] = Q_d * (v_vec[2]*B_vec[0] - v_vec[0]*B_vec[2]);
    F_out[2] = Q_d * (v_vec[0]*B_vec[1] - v_vec[1]*B_vec[0]);
}

/* ================================================================
 * L2 — Ion Drag Force (Barnes model + Khrapak refinement)
 * ================================================================ */

double dust_ion_drag_collection(
    double a, double n_i, double m_i, double u_i, double phi_s)
{
    /* Collection ion drag: ions that physically hit the grain.
     * F_col = pi * b_c^2 * n_i * m_i * v_i * u_i
     * b_c = a * sqrt(1 - 2*e*phi_s/(m_i*u_i^2))  for phi_s < 0
     * Barnes et al. (1992), Phys. Rev. Lett. 68, 313 */
    if (a <= 0.0 || n_i <= 0.0 || m_i <= 0.0 || u_i <= 0.0) return 0.0;

    double energy_factor = 2.0 * DUSTY_EC * fabs(phi_s) / (m_i * u_i * u_i);
    double b_c2;
    if (energy_factor < 1.0) {
        b_c2 = a * a * (1.0 - energy_factor);
    } else {
        b_c2 = a * a;
    }

    return M_PI * b_c2 * n_i * m_i * u_i * u_i;
}

double dust_ion_drag_orbit(
    double a, double n_i, double m_i, double u_i,
    double Z_d, double T_i, double lambda_D)
{
    /* Orbit (Coulomb) ion drag: ions deflected without hitting.
     * F_orb = 4*pi * b_90^2 * n_i * m_i * u_i^2 * ln(Lambda_id)
     * b_90 = Z_d*e^2 / (4*pi*eps0 * m_i * u_i^2  * 3)
     * for u_i >> v_thi.
     *
     * For thermal ions (u_i ~ v_thi), the effective velocity is
     * v_eff = sqrt(u_i^2 + (8/pi)*v_thi^2).
     *
     * Khrapak et al. (2002), Phys. Rev. E 66, 046414 */
    if (a <= 0.0 || n_i <= 0.0 || m_i <= 0.0) return 0.0;

    /* Effective velocity including thermal spread */
    double v_thi = sqrt(DUSTY_KB * T_i / m_i);
    double v_eff2 = u_i * u_i + (8.0 / M_PI) * v_thi * v_thi;
    if (v_eff2 <= 0.0) return 0.0;

    /* 90-degree scattering impact parameter */
    double b_90 = Z_d * DUSTY_EC * DUSTY_EC
                  / (4.0 * M_PI * DUSTY_EPS0 * m_i * v_eff2);
    if (b_90 <= 0.0) return 0.0;

    /* Coulomb logarithm for ion-dust scattering */
    double r_min = (a > b_90) ? a : b_90;
    double ln_Lambda;
    if (lambda_D > r_min && r_min > 0.0) {
        ln_Lambda = log(lambda_D / r_min);
    } else {
        ln_Lambda = 10.0; /* typical value */
    }
    if (ln_Lambda < 0.0) ln_Lambda = 0.0;

    return 4.0 * M_PI * b_90 * b_90 * n_i * m_i * u_i * u_i * ln_Lambda;
}

double dust_ion_drag_total(
    double a, double n_i, double m_i, double u_i,
    double phi_s, double Z_d, double T_i, double lambda_D)
{
    double F_col = dust_ion_drag_collection(a, n_i, m_i, u_i, phi_s);
    double F_orb = dust_ion_drag_orbit(a, n_i, m_i, u_i, Z_d, T_i, lambda_D);
    return F_col + F_orb;
}

double dust_ion_drag_khrapak(
    double a, double n_i, double m_i, double T_i,
    double u_i, double Z_d, double lambda_D)
{
    /* Khrapak et al. (2002) refined ion drag formula.
     * F_id = sqrt(2*pi) * a^2 * n_i * m_i * v_thi * u_i
     *        * [ (pi/4)*(b_c/a)^2 + (Z_d*e^2/(4*pi*eps0*a*k_B*T_i))^2
     *        * ln(1 + (lambda_D/(a*k_B*T_i))^2 * (v_eff^2/k_B*T_i)^2 ) ]
     *
     * This formula smoothly interpolates between subthermal (u_i << v_thi)
     * and suprathermal (u_i >> v_thi) regimes. */
    if (a <= 0.0 || n_i <= 0.0 || m_i <= 0.0 || T_i <= 0.0) return 0.0;

    double v_thi = sqrt(DUSTY_KB * T_i / m_i);
    double v_eff2 = u_i*u_i + (8.0/M_PI)*v_thi*v_thi;

    /* Collection term */
    double b_c2 = a * a;
    double term_col = (M_PI/4.0) * (b_c2 / (a*a)); /* ≈ pi/4 */

    /* Orbit term */
    double Z_factor = Z_d * DUSTY_EC * DUSTY_EC
                      / (4.0 * M_PI * DUSTY_EPS0 * a * DUSTY_KB * T_i);
    double Zf2 = Z_factor * Z_factor;
    double arg = 1.0 + pow(lambda_D * v_eff2
                           / (a * DUSTY_KB * T_i / m_i), 2.0);
    double term_orb = Zf2 * log(arg);

    return sqrt(2.0 * M_PI) * a*a * n_i * m_i * v_thi * u_i
           * (term_col + term_orb);
}

/* ================================================================
 * L2 — Neutral Drag Force
 * ================================================================ */

double dust_neutral_drag_force(
    double m_d, double nu_dn, double v_rel)
{
    /* Epstein drag: F_nd = -m_d * nu_dn * v_rel */
    if (m_d <= 0.0 || nu_dn <= 0.0) return 0.0;
    return m_d * nu_dn * v_rel;
}

double dust_stokes_drag_force(double a, double eta, double v_rel)
{
    /* Stokes drag: F = 6*pi*eta*a*v_rel
     * Valid for continuum regime (Kn << 1). */
    if (a <= 0.0 || eta <= 0.0) return 0.0;
    return 6.0 * M_PI * eta * a * v_rel;
}

double dust_general_drag_force(
    double a, double eta, double v_rel, double Kn)
{
    /* Interpolation between Epstein and Stokes regimes.
     * Uses Cunningham slip correction:
     * F = F_stokes / C_c(Kn)
     * C_c = 1 + Kn*(A1 + A2*exp(-A3/Kn))
     * Millikan's values: A1=1.257, A2=0.400, A3=1.10 */
    if (a <= 0.0 || eta <= 0.0) return 0.0;
    double F_stokes = dust_stokes_drag_force(a, eta, v_rel);
    if (Kn <= 0.0) return F_stokes;

    double A1 = 1.257, A2 = 0.400, A3 = 1.10;
    double C_c = 1.0 + Kn * (A1 + A2 * exp(-A3 / Kn));
    if (C_c <= 0.0) return F_stokes;
    return F_stokes / C_c;
}

/* ================================================================
 * L4 — Thermophoretic Force
 * ================================================================ */

double dust_thermophoretic_force(
    double a, double T_n, double m_n,
    double kappa_gas, double grad_T, double alpha)
{
    /* Thermophoretic force (hot→cold):
     * F_th = -(32/15)*(a^2/v_thn)*kappa_gas*grad_T
     *        * (1 + 5*pi/32*(1-alpha))
     *
     * Talbot et al. (1980), J. Fluid Mech. 101, 737 */
    if (a <= 0.0 || T_n <= 0.0 || m_n <= 0.0 || kappa_gas <= 0.0) return 0.0;
    double v_thn = sqrt(DUSTY_KB * T_n / m_n);
    double base = (32.0/15.0) * (a*a/v_thn) * kappa_gas * grad_T;
    double correction = 1.0 + (5.0*M_PI/32.0) * (1.0 - alpha);
    return -base * correction;
}

/* ================================================================
 * L4 — Radiation Pressure
 * ================================================================ */

double dust_radiation_pressure_force(
    double a, double Q_pr, double I_radiation)
{
    /* F_rad = pi*a^2 * Q_pr * I / c
     * Q_pr ≈ 1 for a >> wavelength, ~(a/lambda)^4 for a << wavelength */
    if (a <= 0.0 || I_radiation <= 0.0) return 0.0;
    return M_PI * a * a * Q_pr * I_radiation / DUSTY_C;
}

/* ================================================================
 * L5 — Complete Force Balance & Levitation
 * ================================================================ */

DustForceResult dust_compute_all_forces(
    double m_d, double Q_d, double g,
    double E_field, const double *B_field,
    const double *v_dust, const double *v_neutral,
    double a, double n_i, double m_i, double u_i,
    double phi_s, double Z_d, double T_i, double lambda_D,
    double nu_dn, double eta_gas, double v_rel,
    double kappa_gas, double grad_T, double alpha_th,
    double I_rad, double Q_pr)
{
    /* Assemble all forces acting on a single dust grain.
     * This is the central diagnostic for understanding grain dynamics. */
    DustForceResult fr;
    fr.F_gravity = dust_gravity_force(m_d, g);
    fr.F_electric = dust_electric_force(Q_d, E_field);

    if (B_field && v_dust) {
        double F_lorentz[3];
        dust_lorentz_force_vector(Q_d, v_dust, B_field, F_lorentz);
        fr.F_electric += sqrt(F_lorentz[0]*F_lorentz[0]
                             + F_lorentz[1]*F_lorentz[1]
                             + F_lorentz[2]*F_lorentz[2]);
    }

    fr.F_ion_drag = dust_ion_drag_total(
        a, n_i, m_i, u_i, phi_s, Z_d, T_i, lambda_D);
    fr.F_neutral_drag = dust_neutral_drag_force(m_d, nu_dn, v_rel);
    fr.F_thermophoretic = dust_thermophoretic_force(
        a, (v_neutral ? 0.0 : 300.0), m_i, kappa_gas, grad_T, alpha_th);
    fr.F_yukawa = 0.0;
    fr.F_radiation = dust_radiation_pressure_force(a, Q_pr, I_rad);

    fr.F_total = fr.F_gravity + fr.F_electric
                 + fr.F_ion_drag + fr.F_neutral_drag
                 + fr.F_thermophoretic + fr.F_yukawa
                 + fr.F_radiation;

    return fr;
}

double dust_levitation_height(
    double m_d, double Q_d, double g,
    double E_0_sheath, double lambda_D,
    double F_ion_drag_z, double z_min, double z_max, double tol)
{
    /* Find height z where F_electric(z) = F_gravity - F_ion_drag_z.
     * F_electric(z) = Q_d * E_0_sheath * exp(-z/lambda_D) (upward)
     * F_gravity = m_d*g (downward)
     * F_ion_drag_z (downward, constant approximation)
     *
     * Equilibrium: Q_d*E_0*exp(-z/lambda_D) = m_d*g + F_ion_drag_z
     * Solution: z = -lambda_D * ln((m_d*g + F_ion_drag_z)/(Q_d*E_0)) */
    if (lambda_D <= 0.0 || Q_d <= 0.0 || E_0_sheath <= 0.0) return -1.0;

    double F_total_down = m_d * g + F_ion_drag_z;
    if (F_total_down <= 0.0) return z_max; /* levitates to top */

    double ratio = F_total_down / (Q_d * E_0_sheath);
    if (ratio >= 1.0) return -1.0; /* falls to bottom */
    if (ratio <= 0.0) return z_max;

    double z_eq = -lambda_D * log(ratio);
    if (z_eq < z_min) return z_min;
    if (z_eq > z_max) return z_max;
    return z_eq;
}
