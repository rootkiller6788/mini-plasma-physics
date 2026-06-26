/**
 * @file  dusty_potential.c
 * @brief Electric potential structures in dusty plasmas.
 *
 * Dust modifies the plasma potential through its charge. Key structures:
 *   - Sheath with dust levitation
 *   - Dust void (central dust-free region in RF discharges)
 *   - Double layer formation
 *   - Dusty plasma sheath (modified Bohm criterion)
 *
 * L4-L7: Poisson-Boltzmann, sheath physics, dust voids.
 *
 * References:
 *   Nitter (1996), Plasma Sources Sci. Technol. 5, 93
 *   Morfill & Tsytovich (2000), Plasma Phys. Rep. 26, 682 — void
 *   Sheridan & Goree (1991), Phys. Fluids B 3, 2796 — sheath
 */

#include "dusty_plasma.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* ================================================================
 * L4 — Poisson-Boltzmann in Dusty Plasma
 * ================================================================ */

/**
 * @brief Linearized Poisson-Boltzmann potential around a dust grain.
 *
 * phi(r) = phi_s * (a/r) * exp(-(r-a)/lambda_D)
 *
 * This is the Debye-Hückel solution for a spherical grain.
 * Valid when |e*phi/(k_B*T)| << 1 (linear regime).
 *
 * @param r  Distance from grain center [m] (r >= a)
 * @return Electrostatic potential [V]
 */
double dust_debye_huckel_potential(
    double r, double a, double phi_s, double lambda_D)
{
    if (r <= a || lambda_D <= 0.0) return phi_s;
    return phi_s * (a / r) * exp(-(r - a) / lambda_D);
}

/**
 * @brief Electric field magnitude from Debye-Hückel potential.
 *
 * E(r) = -d phi/dr = phi_s * (a/r^2) * (1 + r/lambda_D) * exp(-(r-a)/lambda_D)
 */
double dust_debye_huckel_field(
    double r, double a, double phi_s, double lambda_D)
{
    if (r <= a || lambda_D <= 0.0) return 0.0;
    double factor = phi_s * a / (r * r);
    double yukawa = (1.0 + r / lambda_D) * exp(-(r - a) / lambda_D);
    return factor * yukawa;
}

/* ================================================================
 * L4 — Dusty Plasma Sheath
 * ================================================================ */

/**
 * @brief Sheath electric field — exponential model.
 *
 * E(z) = E_0 * exp(-z / lambda_D)
 *
 * This is the simplest sheath model, valid near the wall/electrode.
 * The field decays exponentially from the electrode surface.
 *
 * @param z  Distance from electrode [m]
 * @param E_0 Peak electric field at electrode [V/m]
 * @param lambda_D Sheath screening length [m]
 * @return Electric field magnitude [V/m]
 */
double dust_sheath_electric_field(double z, double E_0, double lambda_D)
{
    if (z < 0.0 || lambda_D <= 0.0) return 0.0;
    return E_0 * exp(-z / lambda_D);
}

/**
 * @brief Sheath potential — exponential model.
 *
 * phi(z) = phi_wall * exp(-z / lambda_D)
 *
 * where phi_wall = V_bias - V_plasma is the wall bias relative to plasma.
 */
double dust_sheath_potential(double z, double phi_wall, double lambda_D)
{
    if (z < 0.0 || lambda_D <= 0.0) return 0.0;
    return phi_wall * exp(-z / lambda_D);
}

/**
 * @brief Modified Bohm criterion with dust.
 *
 * For a dusty plasma, the Bohm velocity for ion entry into the sheath
 * is modified:
 *
 * v_Bohm = c_s * sqrt( (1 + alpha_d) / (1 + alpha_d * T_e/T_i) )
 *
 * where alpha_d = n_d/n_i is the dust-to-ion density ratio.
 *
 * For alpha_d -> 0 (no dust): v_Bohm = c_s (classic Bohm criterion).
 * For alpha_d > 0: v_Bohm < c_s (dust reduces the ion flow).
 *
 * @return Modified Bohm velocity [m/s]
 */
double dust_modified_bohm_velocity(
    double T_e, double T_i, double m_i, double n_d, double n_i)
{
    if (T_e <= 0.0 || m_i <= 0.0 || n_i <= 0.0) return 0.0;
    double c_s = sqrt(DUSTY_KB * T_e / m_i);
    if (T_i <= 0.0) return c_s;
    double alpha_d = n_d / n_i;
    double numerator = 1.0 + alpha_d;
    double denominator = 1.0 + alpha_d * T_e / T_i;
    if (denominator <= 0.0) return c_s;
    return c_s * sqrt(numerator / denominator);
}

/**
 * @brief Dust charge in sheath as function of height.
 *
 * In the sheath, the electron density drops faster than ion density
 * (electrons are repelled by negative wall), so the floating potential
 * becomes LESS negative (dust charge decreases with height).
 *
 * Simplified: Z_d(z) = Z_d0 * exp(-z / (2*lambda_D))
 *
 * where Z_d0 is the charge at the sheath edge (z=0).
 */
double dust_charge_in_sheath(double z, double Z_d0, double lambda_D)
{
    if (z < 0.0 || lambda_D <= 0.0) return Z_d0;
    return Z_d0 * exp(-z / (2.0 * lambda_D));
}

/* ================================================================
 * L7 — Dust Void Model
 * ================================================================ */

/**
 * @brief Dust void size estimation (force balance model).
 *
 * In RF discharges under microgravity, a central dust-free region
 * (void) forms due to the outward ion drag force exceeding the
 * inward electric force.
 *
 * Void radius R_void satisfies:
 * F_electric(R_void) + F_ion_drag(R_void) = 0
 *
 * Simplified estimate:
 * R_void ≈ lambda_D * ln( F_id0 / (Z_d0 * e * E_0) )
 *
 * @return Void radius [m], or 0 if no void forms.
 */
double dust_void_radius_estimate(
    double lambda_D, double F_id0, double Z_d0, double E_0)
{
    /* Void forms when ion drag pushes outward stronger than electric
     * force pulls inward. The equilibrium condition:
     * Z_d*e*E_0*exp(-R/lambda_D) = F_id0*exp(-R/(2*lambda_D)) */
    if (lambda_D <= 0.0 || Z_d0 <= 0.0 || E_0 <= 0.0) return 0.0;
    double ratio = F_id0 / (Z_d0 * DUSTY_EC * E_0);
    if (ratio <= 1.0) return 0.0; /* no void — electric force dominates */
    return lambda_D * log(ratio);
}

/**
 * @brief Dust density profile in void transition region.
 *
 * n_d(r) = n_d0 / (1 + exp((r - R_void) / delta))
 *
 * where delta is the transition width (typically ~ lambda_D).
 * This is a Fermi-function-like profile observed experimentally.
 */
double dust_density_void_profile(
    double r, double n_d0, double R_void, double delta)
{
    if (n_d0 <= 0.0 || delta <= 0.0) return 0.0;
    return n_d0 / (1.0 + exp((r - R_void) / delta));
}
