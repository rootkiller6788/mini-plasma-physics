/**
 * plasma_waves.h — Plasma Waves and Instabilities
 *
 * Dispersion relations, growth rates, and stability analysis for
 * the canonical wave modes in magnetized and unmagnetized plasmas.
 *
 * References:
 *   - Stix, "Waves in Plasmas" (1992)
 *   - Swanson, "Plasma Waves" (2003)
 *   - Brambilla, "Kinetic Theory of Plasma Waves" (1998)
 *   - MIT 22.611 / Princeton PHY 521
 *
 * Knowledge Levels:
 *   L4 — Fundamental Laws: Cold plasma dispersion, CMA diagram
 *   L5 — Computational: Dispersion root finding
 *   L6 — Canonical Systems: Langmuir, ion-acoustic, Alfven, whistler
 *   L8 — Advanced: Nonlinear wave coupling
 */

#ifndef PLASMA_WAVES_H
#define PLASMA_WAVES_H

#include "plasma_constants.h"

/* ================================================================
 * L6.1: Electrostatic Waves in Unmagnetized Plasma
 * ================================================================ */

/**
 * Langmuir wave (electron plasma wave) dispersion:
 *
 * omega^2 = omega_pe^2 + 3 * v_th,e^2 * k^2
 *
 * This is the Bohm-Gross dispersion relation.
 * For k -> 0: omega -> omega_pe
 * For k -> inf: omega ~ sqrt(3) * v_th,e * k (thermal correction)
 *
 * Complexity: O(1)
 */
double langmuir_dispersion(double k, double omega_pe, double v_the);

/**
 * Ion acoustic wave dispersion:
 *
 * omega^2 = k^2 * c_s^2 / (1 + k^2 * lambda_De^2)
 *
 * For k*lambda_De << 1: omega = k * c_s (sound-like)
 * For k*lambda_De >> 1: omega -> omega_pi (ion plasma oscillation)
 *
 * Complexity: O(1)
 */
double ion_acoustic_dispersion(double k, double cs, double lambda_De);

/**
 * Ion acoustic wave dispersion with Landau damping.
 *
 * Returns both real frequency and damping rate.
 *
 * omega_r = k * c_s / sqrt(1 + k^2 * lambda_De^2)
 * gamma = -omega_r * sqrt(pi/8) * (omega_r/(k*v_th,e))^3
 *         * [sqrt(m_e/m_i) + (T_e/T_i)^(3/2)*exp(-T_e/(2*T_i))]
 *
 * Complexity: O(1)
 */
void ion_acoustic_damped(double k, double cs, double lambda_De,
                         double Te, double Ti, double mi,
                         double *omega_r, double *gamma);

/* ================================================================
 * L6.2: Cold Plasma Dispersion (Stix)
 * ================================================================ */

/** Cold plasma dielectric tensor components (Stix notation) */
typedef struct {
    double S;  /**< Sum term: (R+L)/2 */
    double D;  /**< Difference term: (R-L)/2 */
    double P;  /**< Parallel term: 1 - sum_s (omega_ps^2/omega^2) */
    double R;  /**< Right-hand circular: S + D */
    double L;  /**< Left-hand circular: S - D */
} ColdPlasmaDielectric;

/**
 * Compute Stix cold plasma dielectric tensor elements.
 *
 * S = 1 - sum_s omega_ps^2 / (omega^2 - omega_cs^2)
 * D = sum_s (omega_cs/omega) * omega_ps^2 / (omega^2 - omega_cs^2)
 * P = 1 - sum_s omega_ps^2 / omega^2
 *
 * where s runs over all plasma species.
 *
 * Reference: Stix (1992), "Waves in Plasmas", Ch. 1
 *
 * Complexity: O(n_species)
 */
void cold_plasma_dielectric(double omega,
                            const double *wp_species,
                            const double *wc_species,
                            const int *signs,
                            int n_species,
                            ColdPlasmaDielectric *diel);

/**
 * Cold plasma dispersion relation for angle theta:
 *
 * A*n^4 - B*n^2 + C = 0
 *
 * where n = ck/omega is the refractive index,
 * A = S*sin^2(theta) + P*cos^2(theta)
 * B = R*L*sin^2(theta) + P*S*(1+cos^2(theta))
 * C = P*R*L
 *
 * Returns the two solutions n^2 for the fast/slow modes.
 *
 * Reference: Stix (1992), Ch. 2
 */
void cold_plasma_refractive_index(const ColdPlasmaDielectric *diel,
                                  double cos_theta,
                                  double *n2_fast, double *n2_slow);

/**
 * Check for resonance (n^2 -> infinity): solve A(omega) = 0
 *
 * This gives the resonance cones where waves are absorbed.
 */
int cold_plasma_resonance_condition(const ColdPlasmaDielectric *diel,
                                    double cos_theta);

/**
 * Check for cutoff (n^2 -> 0): solve C(omega) = 0 = P*R*L
 *
 * Cutoffs are: P=0 (plasma frequency cutoff),
 *              R=0 (R-cutoff, lower hybrid),
 *              L=0 (L-cutoff, ion cyclotron)
 */
void cold_plasma_cutoffs(const ColdPlasmaDielectric *diel,
                         int *has_P_cut, int *has_R_cut, int *has_L_cut);

/* ================================================================
 * L6.3: MHD Waves (revisited from full dispersion)
 * ================================================================ */

/**
 * Alfven wave dispersion (shear):
 * omega = k_parallel * v_A
 *
 * This is the incompressible shear Alfven wave.
 * Complexity: O(1)
 */
double alfven_wave_dispersion(double k_parallel, double vA);

/**
 * Fast magnetosonic wave dispersion:
 * omega = k * v_fast
 *
 * v_fast^2 = (1/2)*(c_s^2+v_A^2+sqrt((c_s^2+v_A^2)^2-4*c_s^2*v_A^2*cos^2(theta)))
 */
double fast_wave_dispersion(double k, double cs, double vA, double cos_theta);

/**
 * Slow magnetosonic wave dispersion:
 * omega = k * v_slow
 */
double slow_wave_dispersion(double k, double cs, double vA, double cos_theta);

/**
 * Whistler wave (R-wave, helicon) dispersion:
 *
 * n^2 = 1 - omega_pe^2 / (omega * (omega - omega_ce))
 *
 * For omega_ci << omega << omega_ce:
 *   omega = omega_ce * (k*c/omega_pe)^2  (whistler branch)
 *
 * Complexity: O(1)
 */
double whistler_dispersion(double k, double omega_pe, double omega_ce);

/**
 * Whistler group velocity:
 * v_g = d omega / dk = 2 * omega_ce * c^2 * k / omega_pe^2
 *
 * This is the speed at which whistler wave packets propagate.
 * Used in ionospheric/space plasma diagnostics.
 */
double whistler_group_velocity(double k, double omega_pe, double omega_ce);

/**
 * Lower hybrid wave dispersion:
 *
 * omega_LH^2 = (omega_ci*omega_ce + omega_pi^2) / (1 + omega_pe^2/omega_ce^2)
 *            ~ sqrt(omega_ci * omega_ce) for dense plasmas
 */
double lower_hybrid_dispersion(double omega_pi, double omega_ci, double omega_ce);

/* ================================================================
 * L6.4: Instabilities
 * ================================================================ */

/**
 * Two-stream instability growth rate.
 *
 * When an electron beam passes through a stationary plasma,
 * the dispersion relation:
 *   1 = omega_pe^2/omega^2 + omega_b^2/(omega - k*v_b)^2
 *
 * has complex roots for k < omega_pe/v_b.
 *
 * Maximum growth rate:
 *   gamma_max = (sqrt(3)/2^(4/3)) * (n_b/n_0)^(1/3) * omega_pe
 *
 * Reference: Buneman (1958), Phys. Rev. 115, 503
 *
 * Complexity: O(1)
 */
double two_stream_growth_rate(double n_beam, double n_plasma,
                              double omega_pe);

/**
 * Weibel instability growth rate.
 *
 * In a plasma with anisotropic temperature (T_perp > T_par),
 * the Weibel instability spontaneously generates magnetic fields.
 *
 * gamma_max = omega_pe * v_th/c * sqrt((T_perp/T_par - 1) * (beta/2))
 *
 * Reference: Weibel (1959), Phys. Rev. Lett. 2, 83
 */
double weibel_growth_rate(double omega_pe, double v_th, double T_perp,
                          double T_par, double beta);

/**
 * Firehose instability condition.
 *
 * Occurs when p_parallel > p_perp + B^2/mu_0.
 *
 * Returns 1 if firehose unstable, 0 otherwise.
 */
int firehose_unstable(double p_parallel, double p_perp, double B);

/**
 * Mirror instability condition.
 *
 * Occurs when p_perp/p_parallel > 1 + B^2/(2*mu_0*p_parallel).
 *
 * Returns 1 if mirror unstable, 0 otherwise.
 */
int mirror_unstable(double p_parallel, double p_perp, double B);

/**
 * Kelvin-Helmholtz instability growth rate (velocity shear).
 *
 * For two plasmas with relative velocity v0:
 *   gamma_KH = (k * v0 / 2) * sqrt((rho1*rho2)/((rho1+rho2)^2))
 */
double kelvin_helmholtz_growth(double k, double v0, double rho1, double rho2);

/**
 * Rayleigh-Taylor instability growth rate (density gradient + acceleration).
 *
 * gamma_RT = sqrt(g_eff * k * (rho_heavy - rho_light)/(rho_heavy + rho_light))
 */
double rayleigh_taylor_growth(double g_eff, double k,
                              double rho_heavy, double rho_light);

/* ================================================================
 * L5: Dispersion Solver (Numerical Root Finding)
 * ================================================================ */

/**
 * Double root finding callback type for dispersion relations.
 * @param omega Trial frequency [rad/s]
 * @param params User parameters array
 * @return D(omega, k) — the dispersion function (root when D=0)
 */
typedef double (*DispersionFunc)(double omega, const double *params);

/**
 * Find all real roots of D(omega) = 0 in [omega_min, omega_max]
 * using a sign-change scan with Brent refinement.
 *
 * @param D The dispersion function
 * @param params Parameters for D
 * @param omega_min Lower bound [rad/s]
 * @param omega_max Upper bound [rad/s]
 * @param n_scan Number of scanning intervals
 * @param roots Output array (caller pre-allocated)
 * @param max_roots Maximum number of roots to find
 * @return Number of roots found
 */
int find_dispersion_roots(DispersionFunc D, const double *params,
                          double omega_min, double omega_max,
                          int n_scan, double *roots, int max_roots);

/* ================================================================
 * L8: Nonlinear Three-Wave Coupling
 * ================================================================ */

/**
 * Three-wave resonant coupling condition.
 *
 * Two waves (omega1,k1) and (omega2,k2) can parametrically couple
 * to a third wave (omega3,k3) if:
 *   omega1 + omega2 = omega3
 *   k1 + k2 = k3
 *
 * Check if three given waves satisfy the resonance conditions
 * within tolerance.
 *
 * @return 1 if resonant, 0 otherwise
 */
int three_wave_resonance(double w1, double k1, double w2, double k2,
                         double w3, double k3, double tolerance);

/**
 * Parametric decay instability growth rate.
 *
 * A large-amplitude pump wave (omega0,k0) decays into two daughter
 * waves (omega1,k1) and (omega2,k2). The growth rate is:
 *
 * gamma = sqrt(V^2 * |k0|^2 / (4 * omega1 * omega2))
 *
 * where V is the coupling coefficient (oscillation velocity amplitude).
 *
 * Reference: Sagdeev & Galeev (1969), "Nonlinear Plasma Theory"
 */
double parametric_decay_growth(double V_coupling, double k0,
                                double omega1, double omega2);

#endif /* PLASMA_WAVES_H */
