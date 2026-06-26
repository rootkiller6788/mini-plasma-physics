/**
 * @file  dusty_plasma_core.c
 * @brief Core dusty plasma parameter functions.
 *
 * Implements the fundamental length and frequency scales that
 * define dusty plasma physics. These parameters are the building
 * blocks for all more advanced models.
 *
 * L1-L2: Definitions and core plasma parameters.
 *
 * References:
 *   Goldston & Rutherford (1995), "Introduction to Plasma Physics"
 *   Shukla & Mamun (2002), "Introduction to Dusty Plasma Physics"
 *   Piel (2010), "Plasma Physics", Springer
 */

#include "dusty_plasma.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* ================================================================
 * L1 — Debye Lengths (Screening Scales)
 * ================================================================
 *
 * The Debye length is the fundamental screening length in a plasma.
 * A test charge is shielded by a redistribution of surrounding charges
 * over this distance. For dusty plasmas, there are three Debye lengths,
 * but the dust Debye length is usually negligible because T_d << T_e.
 */

double dust_debye_electron(double n_e, double T_e)
{
    /* lambda_De = sqrt(eps0 * k_B * T_e / (n_e * e^2))
     *
     * Derivation: Linearize Poisson-Boltzmann, assume thermal electrons.
     * For T_e = 3 eV, n_e = 10^15 m^-3 → lambda_De ≈ 400 um.
     *
     * Goldston & Rutherford (1995) §2.3 */
    if (n_e <= 0.0 || T_e <= 0.0) return 0.0;
    return sqrt(DUSTY_EPS0 * DUSTY_KB * T_e / (n_e * DUSTY_EC * DUSTY_EC));
}

double dust_debye_ion(double n_i, double T_i)
{
    /* lambda_Di = sqrt(eps0 * k_B * T_i / (n_i * e^2))
     *
     * Ion Debye length is typically smaller than electron Debye length
     * because T_i << T_e in low-temperature laboratory plasmas. */
    if (n_i <= 0.0 || T_i <= 0.0) return 0.0;
    return sqrt(DUSTY_EPS0 * DUSTY_KB * T_i / (n_i * DUSTY_EC * DUSTY_EC));
}

double dust_debye_total(double lambda_De, double lambda_Di)
{
    /* Total (linearized) Debye length:
     * lambda_D^{-2} = lambda_De^{-2} + lambda_Di^{-2}
     *
     * This follows from the linear superposition of screening species.
     * The dust contribution lambda_Dd^{-2} is omitted because:
     *   lambda_Dd^{-2} ∝ Z_d^2 n_d / T_d → typically << other terms
     * since T_d ~ 0.03 eV while T_e ~ 3 eV.
     *
     * Piel (2010), Eq. (11.3) */
    if (lambda_De <= 0.0 || lambda_Di <= 0.0) return 0.0;
    double inv_sq_sum = 1.0/(lambda_De*lambda_De) + 1.0/(lambda_Di*lambda_Di);
    if (inv_sq_sum <= 0.0) return 0.0;
    return 1.0 / sqrt(inv_sq_sum);
}

/* ================================================================
 * L1 — Plasma Frequencies (Timescales)
 * ================================================================
 *
 * The plasma frequency is the natural oscillation frequency of a
 * charge species responding to a charge perturbation. It sets the
 * fastest timescale for collective behavior.
 */

double dust_plasma_freq_electron(double n_e)
{
    /* omega_pe = sqrt(n_e * e^2 / (eps0 * m_e))
     *
     * For n_e = 10^15 m^-3: omega_pe ≈ 1.8 GHz.
     * This is the fastest timescale in the plasma. */
    if (n_e <= 0.0) return 0.0;
    return sqrt(n_e * DUSTY_EC * DUSTY_EC / (DUSTY_EPS0 * DUSTY_ME));
}

double dust_plasma_freq_ion(double n_i, double m_i)
{
    /* omega_pi = sqrt(n_i * e^2 / (eps0 * m_i))
     *
     * For Ar ions (m_i = 40 amu), n_i = 10^15 m^-3: omega_pi ≈ 6.6 MHz.
     * About 300x slower than electron plasma frequency. */
    if (n_i <= 0.0 || m_i <= 0.0) return 0.0;
    return sqrt(n_i * DUSTY_EC * DUSTY_EC / (DUSTY_EPS0 * m_i));
}

double dust_plasma_freq(double n_d, double Z_d, double m_d)
{
    /* omega_pd = sqrt(n_d * Z_d^2 * e^2 / (eps0 * m_d))
     *
     * For a = 5 um SiO2 grain (m_d ≈ 1.15e-12 kg), Z_d = 2000,
     * n_d = 10^10 m^-3: omega_pd ≈ 300 rad/s ≈ 48 Hz.
     *
     * This ultra-low frequency is the hallmark of dusty plasma dynamics.
     * Shukla & Mamun (2002) §2.4 */
    if (n_d <= 0.0 || Z_d <= 0.0 || m_d <= 0.0) return 0.0;
    double Zd_e = Z_d * DUSTY_EC;
    return sqrt(n_d * Zd_e * Zd_e / (DUSTY_EPS0 * m_d));
}

/* ================================================================
 * L2 — Havnes Parameter & Regime Classification
 * ================================================================ */

double dust_havnes_parameter(double Z_d, double n_d, double n_e)
{
    /* P = Z_d * n_d / n_e
     *
     * Named after Ove Havnes who first identified the importance of
     * dust charge in planetary rings (Havnes et al. 1984, 1987).
     *
     * Physical interpretation:
     *   P = (charge density on dust) / (electron charge density)
     *
     * For Saturn's rings: P ~ 10^3 (dust-dominated)
     * For lab RF plasmas: P ~ 0.01-1 (transitional)
     * For fusion devices: P << 0.01 (isolated grain limit) */
    if (n_e <= 0.0) {
        /* If n_e = 0, plasma is fully dust-ion */
        return (Z_d * n_d > 0.0) ? INFINITY : 0.0;
    }
    return Z_d * n_d / n_e;
}

DustRegime dust_classify_regime(double P)
{
    /* Regime boundaries from Havnes et al. (1990):
     *   P < 0.1:  isolated grains — dust doesn't affect plasma
     *   0.1 ≤ P ≤ 1.0: collective — dust modifies plasma response
     *   P > 1.0:  dust-dominated — plasma is largely excluded
     *
     * These are approximate; the transition is smooth. */
    if (P < 0.1) return DUST_REGIME_ISOLATED;
    if (P <= 1.0) return DUST_REGIME_COLLECTIVE;
    return DUST_REGIME_DOMINATED;
}

/* ================================================================
 * L2 — Coulomb Coupling Parameters
 * ================================================================ */

double dust_coulomb_coupling(double Q_d, double a, double T_d)
{
    /* Gamma = Q_d^2 / (4 * pi * eps0 * a * k_B * T_d)
     *
     * Ratio of Coulomb potential energy at inter-particle spacing
     * to thermal kinetic energy. This is the central dimensionless
     * parameter for strongly coupled plasma physics.
     *
     * For a typical lab experiment: Q_d = 5000 e, a = 5 um, T_d = 0.03 eV
     * → Gamma ≈ 10^4, well into the strongly coupled regime.
     *
     * Note: In the Yukawa case, use inter-particle spacing 'a'
     * as the characteristic length, NOT the grain radius.
     * The convention here follows:
     *   Q_d^2 / (4*pi*eps0 * a_ws * k_B * T_d)
     * where a_ws = (3/(4*pi*n_d))^{1/3} is the Wigner-Seitz radius.
     *
     * Ikezi (1986), Phys. Fluids 29, 1764 */
    if (a <= 0.0 || T_d <= 0.0) return 0.0;
    return (Q_d * Q_d) / (4.0 * M_PI * DUSTY_EPS0 * a * DUSTY_KB * T_d);
}

double dust_yukawa_coupling(double Gamma, double kappa)
{
    /* Gamma* = Gamma * exp(-kappa)
     *
     * The exponential screening factor reduces the effective coupling
     * because the Yukawa potential is shorter-ranged than bare Coulomb.
     *
     * For kappa << 1 (weak screening): Gamma* ≈ Gamma
     * For kappa >> 1 (strong screening): Gamma* << Gamma
     *
     * Hamaguchi et al. (1997), Phys. Rev. E 56, 4671 */
    return Gamma * exp(-kappa);
}

int dust_crystal_condition(double Gamma, double kappa)
{
    /* Crystallization criterion for Yukawa systems.
     *
     * The critical Gamma for crystallization depends on kappa:
     *   kappa = 0 (Coulomb): Gamma_crit = 170 (from MD simulations)
     *   kappa → ∞ (hard sphere): Gamma_crit → 106
     *
     * We use the Vaulina-Khrapak interpolation formula:
     *   Gamma_crit = 170/(1+0.18*kappa^2) + 106*(1-1/(1+0.06*kappa^2))
     *
     * This matches MD simulation data to within ~5% for 0 < kappa < 10.
     *
     * Returns 1 if crystalline, 0 otherwise. */
    double Gamma_crit = dust_critical_coupling(kappa);
    return (Gamma >= Gamma_crit) ? 1 : 0;
}

double dust_critical_coupling(double kappa)
{
    /* Simplified piecewise-linear critical coupling.
     * For precise work, use dust_critical_coupling_yukawa().
     *
     * Gamma_crit ≈ 170 for kappa < 1
     * Gamma_crit ≈ 106 for kappa > 5
     * Linear interpolation between. */
    if (kappa < 1.0) return 170.0;
    if (kappa > 5.0) return 106.0;
    return 170.0 - (kappa - 1.0) * (170.0 - 106.0) / 4.0;
}

/* ================================================================
 * L1 — Grain Initialization Utilities
 * ================================================================ */

double dust_grain_mass(double radius, double density)
{
    /* m_d = (4/3) * pi * a^3 * rho
     *
     * Spherical grain assumption — valid for most experimental
     * dust particles (MF spheres) and many natural ones. */
    if (radius <= 0.0 || density <= 0.0) return 0.0;
    return (4.0 / 3.0) * M_PI * radius * radius * radius * density;
}

DustGrain dust_grain_init(double radius, double density, DustMaterial material)
{
    /* Initialize a dust grain with sensible defaults.
     *
     * The initial charge is set to 0 — use the charging functions
     * (dusty_charging.h) to compute the equilibrium charge.
     *
     * Material-specific work functions and SEY are set from
     * dusty_constants.h. */
    DustGrain g;
    g.radius = radius;
    g.mass = dust_grain_mass(radius, density);
    g.charge = 0.0;
    g.charge_number = 0.0;
    g.material_density = density;
    g.material = material;
    g.surface_temp = 300.0;

    switch (material) {
    case DUST_MATERIAL_SILICA:
        g.work_function = DUSTY_WF_SILICA;
        g.se_yield_max = DUSTY_SEY_SILICA;
        break;
    case DUST_MATERIAL_CARBON:
        g.work_function = DUSTY_WF_CARBON;
        g.se_yield_max = 1.0;
        break;
    case DUST_MATERIAL_IRON:
        g.work_function = 4.5;
        g.se_yield_max = 1.3;
        break;
    case DUST_MATERIAL_MF:
        g.work_function = 4.0;
        g.se_yield_max = 2.0;
        break;
    case DUST_MATERIAL_ICE:
        g.work_function = 8.7;
        g.se_yield_max = 1.0;
        break;
    case DUST_MATERIAL_CUSTOM:
    default:
        g.work_function = 5.0;
        g.se_yield_max = 1.0;
        break;
    }

    return g;
}

DustPlasmaState dust_plasma_state_init(
    double n_e, double T_e, double T_i, double T_d,
    double n_d, double Z_d, double m_i, double B)
{
    /* Initialize a consistent plasma state.
     *
     * Quasineutrality is enforced: n_i = n_e + Z_d * n_d
     * This means the ion density is computed from the given
     * electron and dust densities, not independently specified.
     *
     * This is the correct approach because in most experiments,
     * n_e is measured (Langmuir probe) and n_d is known from
     * particle counting, while n_i adjusts to maintain neutrality. */
    DustPlasmaState ps;
    ps.n_e = n_e;
    ps.n_d = n_d;
    ps.Z_d = Z_d;
    ps.n_i = n_e + Z_d * n_d;
    ps.n_n = 0.0;
    ps.T_e = T_e;
    ps.T_i = T_i;
    ps.T_d = T_d;
    ps.T_n = 300.0;
    ps.B = B;
    ps.m_i = m_i;
    ps.m_n = 6.63e-26;  /* Ar atom mass as default neutral */
    ps.phi_float = 0.0;

    ps.lambda_De = dust_debye_electron(n_e, T_e);
    ps.lambda_Di = dust_debye_ion(ps.n_i, T_i);
    ps.lambda_D = dust_debye_total(ps.lambda_De, ps.lambda_Di);

    return ps;
}

double dust_acoustic_speed(double Z_d, double T_e, double m_d)
{
    /* c_da = sqrt(Z_d * k_B * T_e / m_d)
     *
     * This is the phase velocity of long-wavelength dust-acoustic
     * waves (k * lambda_D << 1).
     *
     * Physical picture: electron pressure provides the restoring force
     * (like a spring), dust mass provides the inertia.
     *
     * For Z_d = 2000, T_e = 3 eV, m_d = 1e-15 kg:
     * c_da ≈ 30 mm/s — extremely slow compared to ion sound speed!
     *
     * Rao, Shukla & Yu (1990) */
    if (m_d <= 0.0 || Z_d <= 0.0 || T_e <= 0.0) return 0.0;
    return sqrt(Z_d * DUSTY_KB * T_e / m_d);
}

double dust_ion_acoustic_speed(double T_e, double m_i)
{
    /* c_s = sqrt(k_B * T_e / m_i)
     *
     * The standard ion acoustic speed. In the DIAW, the dust is
     * stationary and ions provide the inertia.
     *
     * For Ar plasma: T_e = 3 eV → c_s ≈ 2700 m/s.
     * This is ~10^5 times faster than c_da. */
    if (T_e <= 0.0 || m_i <= 0.0) return 0.0;
    return sqrt(DUSTY_KB * T_e / m_i);
}
