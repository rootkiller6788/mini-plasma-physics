#ifndef ABSORPTION_H
#define ABSORPTION_H

/*
 * absorption.h -- Laser absorption mechanisms in plasma
 *
 * Models for the principal laser-energy-to-plasma coupling channels:
 * inverse bremsstrahlung (collisional), resonance absorption,
 * vacuum heating (Brunel effect), and J x B heating.
 *
 * References:
 *   - Kruer (1988), Ch. 3-5
 *   - Ginzburg (1964) "Propagation of EM Waves in Plasmas"
 *   - Brunel (1987), Phys. Rev. Lett. 59, 52
 *   - Wilks & Kruer (1997), IEEE J. Quantum Electron. 33, 1954
 *
 * Knowledge Layers:
 *   L1: Absorption coefficient, skin depth, absorption fraction
 *   L2: IB, resonance absorption, Brunel effect, JxB heating
 *   L4: Fresnel equations for stratified plasma, wave equation with damping
 *   L5: Numerical absorption profile integration
 *
 * Courses: MIT 22.611, Princeton PHY 525, ETH 402-0891
 */

#include "plasma_constants.h"
#include "plasma_params.h"

/* ============================================================
 *  L1/L2: Absorption Mechanism Types
 * ============================================================ */

typedef enum {
    ABS_INVERSE_BREMSSTRAHLUNG = 0,
    ABS_RESONANCE              = 1,
    ABS_BRUNEL                 = 2,
    ABS_JXB_HEATING            = 3,
    ABS_COUNT                  = 4
} AbsorptionType;

typedef struct {
    double absorbed_fraction;  /* fraction of incident energy absorbed */
    double absorption_length;  /* characteristic absorption length [m] */
    double heating_rate;       /* energy deposition rate [W/m^3]      */
    double albedo;             /* 1 - absorbed_fraction               */
} AbsorptionResult;

/* ============================================================
 *  L1/L4: Inverse Bremsstrahlung Absorption
 * ============================================================ */

/**
 * inverse_bremsstrahlung_coefficient -- IB absorption coefficient
 *
 *   k_ib = (nu_ei / c) * (ne / nc) / sqrt(1 - ne/nc)
 *
 * Collisional absorption coefficient due to electron-ion
 * Coulomb collisions in the laser field.
 * [Kruer Sec 3.2; Ginzburg Sec 5]
 *
 * Complexity: O(1)
 * Theorem: From complex dielectric function (Drude model)
 *          alpha = 2 (omega/c) Im(sqrt(eps))
 */
double inverse_bremsstrahlung_coefficient(double ne, double nc, double nu_ei);

/**
 * ib_absorption_fraction -- Fraction absorbed by IB over path L
 *
 *   f_abs = 1 - exp(-int_0^L k_ib(z) dz)
 *
 * For an exponential density profile ne = nc exp(-z/L_n)
 * from the critical surface outward, this can be integrated
 * analytically.
 *
 * Complexity: O(1)
 */
double ib_absorption_fraction(double nc, double L_n, double nu_ei,
                               double Te_eV, double lambda_m);

/**
 * ib_linear_profile -- IB absorption for linear density ramp
 *
 * Returns absorbed fraction for a linear density profile
 * from 0 to ne_max over length L.
 */
double ib_linear_profile(double ne_max, double nc, double L,
                          double nu_ei, double lambda_m);

/* ============================================================
 *  L2/L4: Resonance Absorption
 * ============================================================ */

/**
 * resonance_absorption_fraction -- Fraction absorbed at critical surface
 *
 *   f_RA = phi^2(tau) / 2   where phi is the Ginzburg function
 *
 * For p-polarized light, obliquely incident on a density gradient,
 * the EM wave tunnels to the critical surface and resonantly excites
 * a plasma wave, leading to efficient absorption.
 * [Kruer Sec 4.2; Ginzburg 1964]
 *
 * Complexity: O(1)
 * Theorem: Mode conversion at the critical surface:
 *   EM wave -> electron plasma wave through resonant tunnelling
 */
double resonance_absorption_fraction(double theta_inc, double L_n,
                                      double lambda0);

/**
 * ginzburg_function -- Ginzburg function for resonance absorption
 *
 *   phi(tau) ~ 2.3 tau exp(-2 tau^3 / 3)
 *   tau = (omega0 L_n / c)^{1/3} sin(theta)
 *
 * The coupling efficiency function. Maximum absorption ~50%
 * occurs at optimal angle where tau ~ 0.8.
 */
double ginzburg_function(double tau);

/**
 * optimal_resonance_angle -- Angle for maximum resonance absorption
 *
 *   theta_opt = arcsin[0.8 / (omega0 L_n / c)^{1/3}]
 */
double optimal_resonance_angle(double omega0, double L_n);

/* ============================================================
 *  L2/L4: Brunel Effect (Vacuum Heating)
 * ============================================================ */

/**
 * brunel_absorption_fraction -- Vacuum heating absorption
 *
 *   f_B = (eta / pi) * (v_osc^3 / v_the^2 c) * ...
 *
 * At high intensities (a0 v_osc > v_the), electrons pulled into
 * vacuum by the laser field are re-injected into the plasma with
 * the return field, leading to collisionless absorption.
 * [Brunel 1987; Gibbon Sec 4.3]
 *
 * Complexity: O(1)
 * Theorem: Non-adiabatic electron heating from oscillating
 *          sheath field at plasma-vacuum boundary
 */
double brunel_absorption_fraction(double v_osc, double v_the,
                                   double lambda0, double L_n);

/**
 * jxb_heating_fraction -- J x B heating for ultra-intense pulses
 *
 *   f_JxB ~ (omega_p^2 / omega0^2) * (v_osc/c)^2 * ...
 *
 * At relativistic intensities, the v x B term in the Lorentz
 * force drives electron oscillations at 2*omega0, leading to
 * heating at the critical surface.
 * [Wilks & Kruer 1997; Gibbon Sec 4.4]
 */
double jxb_heating_fraction(double v_osc, double wp, double omega0);

/* ============================================================
 *  L5: Integrated Absorption Calculation
 * ============================================================ */

/**
 * total_absorption_fraction -- Combined absorption from all mechanisms
 *
 * Sums contributions from IB, resonance absorption, Brunel effect,
 * and JxB heating, with appropriate weighting based on plasma
 * parameters and laser intensity.
 *
 * Complexity: O(1)
 */
double total_absorption_fraction(const PlasmaState *ps,
                                  const PlasmaDerived *pd,
                                  double L_n, double theta_inc);

/**
 * compute_absorption -- Detailed absorption analysis
 *
 * Fills an AbsorptionResult with per-mechanism breakdown.
 */
int compute_absorption(const PlasmaState *ps, const PlasmaDerived *pd,
                       AbsorptionType type, double L_n, double theta_inc,
                       AbsorptionResult *result);

/**
 * absorption_type_name -- Human-readable name for AbsorptionType
 */
const char *absorption_type_name(AbsorptionType type);

#endif /* ABSORPTION_H */
