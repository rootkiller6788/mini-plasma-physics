#ifndef DUSTY_CRYSTAL_H
#define DUSTY_CRYSTAL_H
/**
 * @file  dusty_crystal.h
 * @brief Dust crystal formation, structure, and phase transitions.
 *
 * When the Coulomb coupling parameter Gamma exceeds a critical value,
 * dust grains in a plasma can form crystalline structures — so-called
 * "plasma crystals" or "Coulomb crystals". These are the classical
 * analog of Wigner crystals in electron systems.
 *
 * Refs:
 *   Ikezi (1986), Phys. Fluids 29, 1764 — first prediction
 *   Thomas et al. (1994), Phys. Rev. Lett. 73, 652 — experimental discovery
 *   Melzer et al. (1994), Phys. Lett. A 191, 301 — 3D Coulomb balls
 *   Hamaguchi et al. (1997), Phys. Rev. E 56, 4671 — phase diagram
 */
#include "dusty_plasma.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================
 * L2 — Crystal Condition & Phase Diagram
 * ================================================================ */

/**
 * @brief Compute Yukawa phase diagram boundary.
 *
 * Returns the Gamma_crit(kappa) for solid-fluid phase transition
 * based on the empirical fit:
 *
 * Gamma_crit = 170 / (1 + 0.18 * kappa^2) + 106 * (1 - 1/(1 + 0.06 * kappa^2))
 *
 * This interpolates between the Coulomb limit (kappa -> 0, Gamma_crit ≈ 170)
 * and the strong screening limit (kappa >> 1, Gamma_crit ~ 106).
 *
 * Ref: Vaulina et al. (2002), Phys. Rev. Lett. 88, 035001
 */
double dust_critical_coupling_yukawa(double kappa);

/**
 * @brief Determine the phase of a dusty plasma system.
 *
 * Returns: 0 = gas, 1 = liquid, 2 = solid (crystal).
 * Uses the empirical phase boundaries from Hamaguchi et al. (1997).
 */
int dust_phase_determine(double Gamma, double kappa);

/**
 * @brief Compute dust thermal velocity.
 * v_thd = sqrt(k_B * T_d / m_d)
 */
double dust_thermal_velocity(double T_d, double m_d);

/**
 * @brief Compute dust Einstein frequency (caged oscillation).
 *
 * Omega_E = sqrt( (Q_d^2 / (4 pi eps0 m_d)) * (kappa^2 / a^3)
 *           * (1 + kappa) * exp(-kappa) )
 *
 * This is the characteristic oscillation frequency of a dust grain
 * in the cage formed by its neighbors in a crystal.
 */
double dust_einstein_frequency(double Q_d, double m_d, double a, double kappa);

/* ================================================================
 * L3 — Yukawa Lattice Sums
 * ================================================================ */

/**
 * @brief Yukawa pair potential energy between two dust grains.
 *
 * U(r) = (Q1 * Q2 / (4 * pi * eps0 * r)) * exp(-r / lambda_D)
 *
 * @return Potential energy [J]
 */
double yukawa_potential(double Q1, double Q2, double r, double lambda_D);

/**
 * @brief Yukawa force magnitude between two dust grains.
 *
 * F(r) = -dU/dr = (Q1*Q2/(4*pi*eps0*r^2)) * (1 + r/lambda_D) * exp(-r/lambda_D)
 *
 * The force is larger than bare Coulomb at short range due to the
 * (1 + r/lambda_D) factor, but exponentially suppressed at long range.
 *
 * @return Force magnitude [N]
 */
double yukawa_force_magnitude(double Q1, double Q2, double r, double lambda_D);

/**
 * @brief Yukawa force vector between two grains (3D).
 *
 * F_12 = yukawa_force_magnitude * r_hat_12
 *
 * @param r1, r2   Position vectors [3] of grain 1 and 2
 * @param F_out    Output force on grain 1 [3]
 */
void yukawa_force_vector(
    const double r1[3], const double r2[3],
    double Q1, double Q2, double lambda_D,
    double F_out[3]);

/**
 * @brief Compute total Yukawa energy of a crystal configuration.
 *
 * U_total = sum_{i<j} U_yukawa(r_ij)
 *
 * Using minimum image convention for periodic boundaries.
 */
double dust_crystal_total_energy(
    const double **positions, size_t n_grains,
    double Q_d, double lambda_D, double Lx, double Ly);

/**
 * @brief Compute Madelung constant for a Yukawa lattice.
 *
 * For a Bravais lattice, the Madelung constant alpha determines
 * the cohesive energy: U = -alpha * Q_d^2 / (4 * pi * eps0 * a)
 *
 * @param lattice_type  0=1D chain, 1=2D hexagonal, 2=3D bcc
 * @param kappa         Screening parameter
 * @return Madelung constant (dimensionless)
 *
 * Complexity: O(N^2) where N = number of lattice sites summed.
 */
double dust_yukawa_madelung(int lattice_type, double kappa);

/* ================================================================
 * L6 — Radial Distribution Function & Structure
 * ================================================================ */

/**
 * @brief Pair correlation function g(r) in Yukawa OCP model.
 *
 * Uses the Hyper-Netted Chain (HNC) approximation for the
 * Yukawa one-component plasma. This is the standard liquid-state
 * theory for strongly coupled dusty plasmas.
 *
 * g(r) = exp(-U(r)/k_B T + h(r) - c(r))
 *
 * where h(r) = g(r) - 1 and c(r) is the direct correlation function.
 *
 * Ref: Hansen & McDonald (2013), "Theory of Simple Liquids", 4th ed.
 */
void dust_pair_correlation_hnc(
    double Gamma, double kappa, double r_max,
    int n_bins, double *r, double *g);

/**
 * @brief Static structure factor S(k) from g(r).
 *
 * S(k) = 1 + n_d * integral (g(r) - 1) * exp(i k · r) d^3 r
 *
 * For isotropic systems: S(k) = 1 + 4*pi*n_d*integral_0^inf r^2*(g(r)-1)*sin(kr)/(kr)*dr
 *
 * S(k) is directly measurable via light scattering or video microscopy.
 */
void dust_structure_factor(
    const double *r, const double *g, int n_bins,
    double n_d, double k_max, int n_k,
    double *k_out, double *S_out);

/**
 * @brief Lindemann melting criterion for dust crystals.
 *
 * Melting occurs when the RMS displacement exceeds ~15% of the
 * inter-particle spacing:
 *   sqrt(<delta r^2>) / a > c_L
 *
 * where c_L ≈ 0.15-0.18 for Yukawa crystals.
 *
 * Ref: Robbins et al. (1988), Phys. Rev. A 37, 2871
 */
int dust_lindemann_melting(double rms_displacement, double a, double c_L);

/**
 * @brief Compute RMS displacement from temperature (harmonic approximation).
 *
 * <delta r^2> = 3 * k_B * T / (m_d * Omega_E^2)
 *
 * where Omega_E is the Einstein frequency.
 */
double dust_rms_displacement_thermal(double T_d, double m_d, double Omega_E);

/* ================================================================
 * L8 — Phase Transition Dynamics
 * ================================================================ */

/**
 * @brief Nucleation rate for dust crystal formation (classical theory).
 *
 * J = J_0 * exp(-Delta G* / (k_B * T))
 *
 * where Delta G* is the critical nucleation barrier.
 *
 * Ref: Kelton & Greer (2010), "Nucleation in Condensed Matter"
 */
double dust_nucleation_rate(
    double Gamma, double kappa, double T_d, double n_d);

/**
 * @brief Dust crystal grain boundary energy (Read-Shockley model).
 *
 * E_gb = E_0 * theta * (A - ln(theta))
 *
 * where theta is the misorientation angle.
 */
double dust_grain_boundary_energy(
    double theta, double E_0, double A);

#ifdef __cplusplus
}
#endif

#endif /* DUSTY_CRYSTAL_H */
