/*
 * absorption.c -- Laser absorption mechanisms in plasma
 *
 * Implements the principal laser-energy-to-plasma coupling channels:
 * inverse bremsstrahlung, resonance absorption, Brunel effect
 * (vacuum heating), and JxB heating.
 *
 * References:
 *   - Kruer (1988), Ch. 3-5
 *   - Ginzburg (1964) "Propagation of EM Waves in Plasmas"
 *   - Brunel (1987), Phys. Rev. Lett. 59, 52
 *   - Wilks & Kruer (1997), IEEE J. Quantum Electron. 33, 1954
 *
 * Knowledge Layers: L1, L2, L4, L5
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "plasma_constants.h"
#include "plasma_params.h"
#include "absorption.h"

static const double c   = PLASMA_C;

/* ============================================================
 *  L1/L4: inverse_bremsstrahlung_coefficient
 *
 *  k_ib = (nu_ei / c) * (ne/nc) / sqrt(1 - ne/nc)
 *
 * This is the spatial absorption coefficient due to electron-
 * ion Coulomb collisions in the oscillating laser field.
 * The 1/sqrt(1 - ne/nc) factor represents the swelling of
 * the electric field as the wave approaches the critical
 * surface (N -> 0 => E-field resonance in cold theory).
 *
 * For 1 um light, ne = 0.1 nc, Te = 1 keV, Z = 6:
 *   nu_ei ~ 3.5e12 s^-1
 *   k_ib ~ 1.3e3 m^-1  (absorption length ~ 0.8 mm)
 *
 * High-Z, low-Te plasmas absorb more strongly via IB.
 *
 * Source: Kruer Sec 3.2; Ginzburg Sec 5
 * Theorem: From Poynting theorem with complex dielectric:
 *   div S = -omega eps0 Im(eps) |E|^2
 *   k_ib = (omega/c) Im(sqrt(eps))
 *
 * Complexity: O(1)
 * ============================================================ */
double inverse_bremsstrahlung_coefficient(double ne, double nc,
                                           double nu_ei)
{
    if (ne <= 0.0 || nc <= 0.0 || nu_ei <= 0.0) return 0.0;
    if (ne >= nc) return DBL_MAX;  /* evanescent region */

    double denom = sqrt(1.0 - ne / nc);
    if (denom <= 0.0) return DBL_MAX;

    return (nu_ei / c) * (ne / nc) / denom;
}

/* ============================================================
 *  L2/L4: ib_absorption_fraction -- Exponential profile
 *
 *  f_abs = 1 - exp(-2 * integral(k_ib(z) dz))
 *
 * For an exponential density profile ne(z) = nc * exp(-(z_cr - z)/L_n)
 * with z_cr being the critical surface position, the integral
 * can be evaluated analytically.
 *
 * This yields the fraction of incident laser energy absorbed
 * by IB for a given density scale length and collision rate.
 *
 * Source: Kruer Sec 3.2
 *
 * Complexity: O(1)
 * ============================================================ */
double ib_absorption_fraction(double nc, double L_n, double nu_ei,
                               double Te_eV, double lambda_m)
{
    if (nc <= 0.0 || L_n <= 0.0 || nu_ei <= 0.0) return 0.0;

    /* For exponential profile, the integral is:
     *   tau = nu_ei * L_n / c  (optical depth at critical) */
    double tau = nu_ei * L_n / c;

    /* Absorption fraction: f = 1 - exp(-8 tau / 3) (Ginzburg) */
    if (tau > 50.0) return 0.999;  /* saturated */

    double f_abs = 1.0 - exp(-8.0 * tau / 3.0);
    return (f_abs > 0.0) ? f_abs : 0.0;
}

/* ============================================================
 *  L5: ib_linear_profile -- Linear density ramp
 *
 * For ne(z) linear from 0 to ne_max over length L.
 * The integral of ne/sqrt(1 - ne/nc) can be done analytically.
 *
 * This is useful for targets with a well-defined linear
 * density ramp (e.g., gas jets, exploded foils).
 *
 * Complexity: O(1)
 * ============================================================ */
double ib_linear_profile(double ne_max, double nc, double L,
                          double nu_ei, double lambda_m)
{
    if (nc <= 0.0 || L <= 0.0 || nu_ei <= 0.0 || ne_max <= 0.0)
        return 0.0;

    /* Linear profile: ne(z) = ne_max * z / L  for z in [0, L] */
    /* Assume ne_max < nc (underdense) for simplicity */
    if (ne_max >= nc) return 1.0;  /* total absorption at critical */

    /* Integral of ne(z)/sqrt(1 - ne(z)/nc) dz for linear ramp */
    double u = ne_max / nc;  /* < 1 */
    /* Integral result: 2 nc L / ne_max * (1 - (1-u)^{3/2}) / 3 */
    double integral = (2.0 * nc * L / ne_max)
                      * (1.0 - pow(1.0 - u, 1.5)) / 3.0;

    double optical_depth = (nu_ei / c) * (1.0 / nc) * integral;
    double f_abs = 1.0 - exp(-2.0 * optical_depth);

    return (f_abs > 0.0 && f_abs < 1.0) ? f_abs : 0.0;
}

/* ============================================================
 *  L2/L4: resonance_absorption_fraction
 *
 *  f_RA = phi^2(tau) / 2
 *
 * For p-polarized light obliquely incident on a plasma with
 * a density gradient, EM energy tunnels from the turning
 * point to the critical surface, where it resonantly excites
 * an electron plasma wave.
 *
 * This is a collisionless absorption mechanism that can
 * dominate over IB for steep gradients and oblique incidence.
 * Maximum absorption ~50% at optimal angle.
 *
 * Source: Kruer Sec 4.2; Ginzburg (1964) Sec 20
 * Theorem: Linear mode conversion of EM wave to EPW
 *          at the critical surface
 *
 * Complexity: O(1)
 * ============================================================ */
double resonance_absorption_fraction(double theta_inc, double L_n,
                                      double lambda0)
{
    if (L_n <= 0.0 || lambda0 <= 0.0) return 0.0;

    double omega0 = 2.0 * M_PI * c / lambda0;
    double tau = pow(omega0 * L_n / c, 1.0/3.0) * sin(theta_inc);
    double phi = ginzburg_function(tau);

    return phi * phi / 2.0;
}

/* ============================================================
 *  L2: ginzburg_function -- Coupling efficiency for resonance
 *       absorption
 *
 *  phi(tau) ~ 2.3 tau exp(-2 tau^3 / 3)
 *
 * This is the Ginzburg function that describes the efficiency
 * of linear mode conversion at the critical surface.
 *
 * tau << 1 : phi ~ 2.3 tau  (linear growth)
 * tau ~ 0.8 : phi ~ 1.0    (maximum, ~50% absorption)
 * tau >> 1 : phi -> 0      (wave reflected before reaching nc)
 *
 * Source: Ginzburg (1964), Sec 20
 *
 * Complexity: O(1)
 * ============================================================ */
double ginzburg_function(double tau)
{
    if (tau <= 0.0) return 0.0;
    if (tau > 5.0) return 0.0;  /* exponentially small */

    return 1.76 * tau * exp(-2.0 * tau * tau * tau / 3.0);
}

/* ============================================================
 *  L2: optimal_resonance_angle -- Angle for maximum absorption
 *
 *  theta_opt = arcsin(0.8 / (omega0 L_n / c)^{1/3})
 *
 * This angle maximizes the Ginzburg function at tau ~ 0.8.
 *
 * For a 1 um laser with L_n = 10 um:
 *   omega0 L_n / c ~ 63
 *   theta_opt ~ arcsin(0.8 / 4.0) ~ 11.5 degrees
 *
 * Complexity: O(1)
 * ============================================================ */
double optimal_resonance_angle(double omega0, double L_n)
{
    if (omega0 <= 0.0 || L_n <= 0.0) return 0.0;
    double tau_factor = pow(omega0 * L_n / c, 1.0/3.0);
    double sin_theta = 0.8 / tau_factor;
    if (sin_theta >= 1.0) return M_PI / 2.0;
    return asin(sin_theta);
}

/* ============================================================
 *  L2: brunel_absorption_fraction -- Vacuum heating
 *
 *  f_B = (eta / pi) * (v_osc^3 / (v_the^2 c)) * ...
 *
 * For high-intensity, short-pulse interactions, electrons
 * pulled into vacuum during one half-cycle are re-injected
 * into the overdense plasma with the return field, leading
 * to efficient collisionless heating.
 *
 * This "Brunel effect" or "vacuum heating" dominates for
 * a0 v_osc > v_the (i.e., when the quiver amplitude exceeds
 * the Debye length).
 *
 * Source: Brunel (1987); Gibbon Sec 4.3
 *
 * Complexity: O(1)
 * ============================================================ */
double brunel_absorption_fraction(double v_osc, double v_the,
                                   double lambda0, double L_n)
{
    if (v_the <= 0.0 || lambda0 <= 0.0) return 0.0;

    double eta = 1.0;  /* efficiency factor, ~1 for step profile */

    /* Brunel scaling */
    double f_abs = (eta / M_PI) * (v_osc * v_osc * v_osc)
                   / (v_the * v_the * c);

    /* Apply geometric factor for finite gradient */
    if (L_n > 0.0) {
        double omega0 = 2.0 * M_PI * c / lambda0;
        double xi = v_osc / (omega0 * L_n);
        f_abs *= (xi < 1.0) ? xi : 1.0;
    }

    return (f_abs < 1.0) ? f_abs : 0.5;
}

/* ============================================================
 *  L2: jxb_heating_fraction -- JxB heating
 *
 *  f_JxB ~ (omega_p^2 / omega0^2) * (v_osc/c)^2
 *
 * At relativistic intensities, the v x B term in the Lorentz
 * force drives electron oscillations at 2*omega0 along the
 * laser propagation direction.  This "JxB heating" is most
 * effective at the critical surface where the E-field is
 * enhanced.
 *
 * Source: Wilks & Kruer (1997); Gibbon Sec 4.4
 *
 * Complexity: O(1)
 * ============================================================ */
double jxb_heating_fraction(double v_osc, double wp, double omega0)
{
    if (omega0 <= 0.0) return 0.0;
    double frac = (wp * wp) / (omega0 * omega0)
                  * (v_osc / c) * (v_osc / c);
    return (frac < 1.0) ? frac : 0.5;
}

/* ============================================================
 *  L5: total_absorption_fraction -- All mechanisms combined
 *
 * Sums contributions from IB, resonance absorption, Brunel
 * effect, and JxB heating with appropriate cross-coupling.
 *
 * The total is not simply additive because the mechanisms
 * compete for the same incident energy.  We use:
 *   f_total = 1 - (1-f_IB)(1-f_RA)(1-f_B)(1-f_JxB)
 *
 * which treats each mechanism as absorbing an independent
 * fraction of the remaining flux.
 *
 * Complexity: O(1)
 * ============================================================ */
double total_absorption_fraction(const PlasmaState *ps,
                                  const PlasmaDerived *pd,
                                  double L_n, double theta_inc)
{
    if (!ps || !pd) return 0.0;

    /* Inverse bremsstrahlung */
    double f_ib = ib_absorption_fraction(pd->nc, L_n, pd->nu_ei,
                                          ps->Te, ps->lam_laser);

    /* Resonance absorption (p-polarized only) */
    double f_ra = 0.0;
    if (theta_inc > 0.0) {
        f_ra = resonance_absorption_fraction(theta_inc, L_n,
                                              ps->lam_laser);
    }

    /* Brunel effect (relevant for high intensity) */
    double v_osc = pd->a0 * c;
    double f_br = brunel_absorption_fraction(v_osc, pd->v_the,
                                              ps->lam_laser, L_n);

    /* JxB heating */
    double f_jxb = jxb_heating_fraction(v_osc, pd->wp,
                                         pd->omega_laser);

    /* Combined absorption (independent-channel model) */
    double surviving = (1.0 - f_ib) * (1.0 - f_ra)
                       * (1.0 - f_br) * (1.0 - f_jxb);
    double f_total = 1.0 - surviving;

    return (f_total >= 0.0 && f_total <= 1.0) ? f_total
           : (f_total > 1.0 ? 1.0 : 0.0);
}

/* ============================================================
 *  L5: compute_absorption -- Detailed absorption analysis
 *
 * Fills AbsorptionResult for a specific mechanism type.
 *
 * Complexity: O(1)
 * ============================================================ */
int compute_absorption(const PlasmaState *ps, const PlasmaDerived *pd,
                       AbsorptionType type, double L_n, double theta_inc,
                       AbsorptionResult *result)
{
    if (!ps || !pd || !result) return -1;

    memset(result, 0, sizeof(AbsorptionResult));

    double v_osc = pd->a0 * c;

    switch (type) {
    case ABS_INVERSE_BREMSSTRAHLUNG:
        result->absorbed_fraction = ib_absorption_fraction(
            pd->nc, L_n, pd->nu_ei, ps->Te, ps->lam_laser);
        if (pd->nu_ei > 0.0 && pd->nc > 0.0) {
            double k_ib = inverse_bremsstrahlung_coefficient(
                ps->ne, pd->nc, pd->nu_ei);
            if (k_ib > 0.0)
                result->absorption_length = 1.0 / k_ib;
        }
        break;

    case ABS_RESONANCE:
        if (theta_inc > 0.0) {
            result->absorbed_fraction = resonance_absorption_fraction(
                theta_inc, L_n, ps->lam_laser);
        }
        break;

    case ABS_BRUNEL:
        result->absorbed_fraction = brunel_absorption_fraction(
            v_osc, pd->v_the, ps->lam_laser, L_n);
        break;

    case ABS_JXB_HEATING:
        result->absorbed_fraction = jxb_heating_fraction(
            v_osc, pd->wp, pd->omega_laser);
        break;

    default:
        return -1;
    }

    result->albedo = 1.0 - result->absorbed_fraction;
    result->heating_rate = result->absorbed_fraction
                           * ps->I_laser / L_n;

    return 0;
}

/* ============================================================
 *  L1: absorption_type_name -- Human-readable names
 * ============================================================ */
const char *absorption_type_name(AbsorptionType type)
{
    static const char *names[] = {
        "Inverse Bremsstrahlung",
        "Resonance Absorption",
        "Brunel (Vacuum Heating)",
        "JxB Heating"
    };
    if (type < 0 || type >= ABS_COUNT) return "Unknown";
    return names[type];
}
