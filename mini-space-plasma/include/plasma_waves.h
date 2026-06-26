/**
 * plasma_waves.h — Plasma Waves and Instabilities
 *
 * References:
 *   Stix "Waves in Plasmas" (1992)
 *   Swanson "Plasma Waves" (2003)
 *   Kivelson & Russell §5 (Plasma Waves in Space)
 *   MIT 22.611 Lectures 15-18
 *
 * Knowledge: L1 (definitions), L2 (wave concepts), L3 (dispersion),
 *   L4 (Vlasov-Maxwell), L5 (dispersion solvers), L6 (wave modes)
 */

#ifndef PLASMA_WAVES_H
#define PLASMA_WAVES_H

#include "space_plasma.h"

/*──────────────────────────────────────────────────────────────────────
 * L1: Wave Mode Identifiers
 *──────────────────────────────────────────────────────────────────────*/

typedef enum {
    WAVE_LANGMUIR,        /* Electron plasma oscillations           */
    WAVE_ION_ACOUSTIC,    /* Ion acoustic waves                     */
    WAVE_ALFVEN_SHEAR,    /* Shear (torsional) Alfven wave          */
    WAVE_FAST_MAGNETOSONIC,/* Fast (compressional) magnetosonic      */
    WAVE_SLOW_MAGNETOSONIC,/* Slow magnetosonic wave                 */
    WAVE_WHISTLER,        /* Whistler mode (R-wave, ω_ci < ω < ω_ce) */
    WAVE_LOWER_HYBRID,    /* Lower hybrid waves                     */
    WAVE_UPPER_HYBRID,    /* Upper hybrid waves                     */
    WAVE_BERNSTEIN,       /* Electron Bernstein modes               */
    WAVE_KINETIC_ALFVEN,  /* Kinetic Alfven wave                    */
    WAVE_MAGNETOACOUSTIC_GRAVITY, /* Gravity-modified MHD waves     */
    WAVE_NONE
} wave_mode_t;

/** @brief Cold plasma dielectric tensor (Stix parameters) */
typedef struct {
    double S;   /* Sum: (R + L)/2                     */
    double D;   /* Difference: (R - L)/2              */
    double P;   /* Plasma: 1 - ω_p²/ω²               */
    double R;   /* Right-hand: 1 - Σ ω_p²/(ω(ω+ω_c)) */
    double L;   /* Left-hand: 1 - Σ ω_p²/(ω(ω-ω_c))  */
} stix_tensor_t;

/** @brief Wave polarization (for given k, ω) */
typedef struct {
    double E_vec[3];        /* Electric field eigenvector            */
    double pol_ratio;       /* |E_x/E_y|                            */
    double helicity;        /* Degree of circular polarization       */
    double ellipticity;     /* b/a for polarization ellipse          */
} wave_polarization_t;

/*──────────────────────────────────────────────────────────────────────
 * L1: Cold Plasma Dispersion Relation
 *──────────────────────────────────────────────────────────────────────*/

/**
 * @brief Compute Stix dielectric tensor elements
 *
 * For cold, magnetized plasma with B along z:
 *   R = 1 - Σ_s ω_ps² / (ω (ω + ω_cs))
 *   L = 1 - Σ_s ω_ps² / (ω (ω - ω_cs))
 *   P = 1 - Σ_s ω_ps² / ω²
 *   S = (R + L)/2,  D = (R - L)/2
 *
 * @param omega      Wave frequency [rad/s]
 * @param omega_pe   Electron plasma frequency [rad/s]
 * @param omega_ce   Electron gyrofrequency [rad/s] (positive)
 * @param omega_pi   Ion plasma frequency [rad/s]
 * @param omega_ci   Ion gyrofrequency [rad/s] (positive)
 * @param eps        Output: Stix tensor
 */
void stix_dielectric_tensor(double omega, double omega_pe, double omega_ce,
                            double omega_pi, double omega_ci, stix_tensor_t *eps);

/**
 * @brief Cold plasma dispersion relation (Appleton-Hartree)
 *
 * For wave with wave-vector at angle θ to B₀:
 *   A n⁴ - B n² + C = 0
 *
 * where n = ck/ω is the refractive index.
 *
 *   A = S sin²θ + P cos²θ
 *   B = R L sin²θ + P S (1 + cos²θ)
 *   C = P R L
 *
 * @param omega   Angular frequency [rad/s]
 * @param k       Wave number [rad/m]
 * @param theta   Angle between k and B [rad]
 * @param eps     Stix tensor
 * @return        0 if dispersion satisfied (within tol), else residual
 */
double cold_plasma_dispersion_residual(double omega, double k, double theta,
                                       const stix_tensor_t *eps);

/**
 * @brief Solve for refractive index n = ck/ω given ω and θ
 *
 * Returns both solutions (fast and slow branches).
 *
 * @param omega    Frequency [rad/s]
 * @param theta    Propagation angle [rad]
 * @param eps      Stix tensor
 * @param n_fast   Output: fast mode refractive index
 * @param n_slow   Output: slow mode refractive index
 * @return         0 on success, -1 if no propagating solution
 */
int cold_plasma_refractive_index(double omega, double theta,
                                 const stix_tensor_t *eps,
                                 double *n_fast, double *n_slow);

/*──────────────────────────────────────────────────────────────────────
 * L2: Specific Wave Mode Dispersion Relations
 *──────────────────────────────────────────────────────────────────────*/

/**
 * @brief Langmuir wave dispersion (Bohm-Gross)
 *
 * ω² = ω_pe² + 3 k² v_the²
 *
 * Thermal correction for warm electrons. In cold plasma limit:
 * ω = ω_pe (electron plasma oscillations).
 *
 * @param k        Wavenumber [rad/m]
 * @param omega_pe Electron plasma frequency [rad/s]
 * @param v_the    Electron thermal velocity [m/s]
 * @return         Angular frequency [rad/s]
 */
double langmuir_frequency(double k, double omega_pe, double v_the);

/**
 * @brief Ion acoustic wave frequency (in long-wavelength limit)
 *
 * ω = k c_s / sqrt(1 + k²λ_De²)
 *
 * For kλ_De ≪ 1: ω ≈ k c_s (sound wave)
 * For kλ_De ≫ 1: ω ≈ ω_pi (ion plasma oscillation)
 *
 * @param k         Wavenumber [rad/m]
 * @param c_s       Ion sound speed [m/s]
 * @param lambda_De Electron Debye length [m]
 * @return          Angular frequency [rad/s]
 */
double ion_acoustic_frequency(double k, double c_s, double lambda_De);

/**
 * @brief Ion acoustic wave damping rate (Landau damping)
 *
 * γ/ω ≈ -√(π/8) * (c_s/v_the)³ * exp(-c_s²/(2v_the²))  for T_e ≫ T_i
 *
 * @param omega    Wave frequency [rad/s]
 * @param c_s      Sound speed [m/s]
 * @param v_the    Electron thermal velocity [m/s]
 * @return         Damping rate γ [s⁻¹] (negative = damped)
 */
double ion_acoustic_landau_damping(double omega, double c_s, double v_the);

/**
 * @brief Alfven wave frequency (shear Alfven)
 *
 * ω = k_par * v_A  where k_par = k cos(θ)
 *
 * Non-dispersive in ideal MHD.
 *
 * @param k      Wavenumber [rad/m]
 * @param theta  Angle between k and B [rad]
 * @param v_A    Alfven speed [m/s]
 * @return       Angular frequency [rad/s]
 */
double alfven_wave_frequency(double k, double theta, double v_A);

/**
 * @brief Fast magnetosonic wave frequency
 *
 * ω² = k² v_f²  where v_f² = (v_A² + c_s² + sqrt((v_A²+c_s²)² - 4v_A²c_s²cos²θ))/2
 *
 * @param k      Wavenumber [rad/m]
 * @param v_A    Alfven speed [m/s]
 * @param c_s    Sound speed [m/s]
 * @param theta  Angle between k and B [rad]
 * @return       Angular frequency [rad/s]
 */
double fast_magnetosonic_frequency(double k, double v_A, double c_s, double theta);

/**
 * @brief Slow magnetosonic wave frequency
 *
 * ω² = k² v_s²  where v_s² = (v_A² + c_s² - sqrt((v_A²+c_s²)² - 4v_A²c_s²cos²θ))/2
 */
double slow_magnetosonic_frequency(double k, double v_A, double c_s, double theta);

/**
 * @brief Whistler wave dispersion (R-mode, ω_ci ≪ ω ≪ ω_ce)
 *
 * ω = ω_ce * (k c/ω_pe)² * cos(θ)
 *
 * @param k         Wavenumber [rad/m]
 * @param omega_pe  Electron plasma frequency [rad/s]
 * @param omega_ce  Electron gyrofrequency [rad/s]
 * @param theta     Propagation angle [rad]
 * @return          Angular frequency [rad/s]
 */
double whistler_frequency(double k, double omega_pe, double omega_ce, double theta);

/**
 * @brief Lower hybrid wave frequency
 *
 * ω_LH ≈ sqrt(ω_ci ω_ce)  (when ω_pi ≫ ω_ci)
 *
 * @param omega_ci  Ion gyrofrequency [rad/s]
 * @param omega_ce  Electron gyrofrequency [rad/s]
 * @return          Lower hybrid frequency [rad/s]
 */
double lower_hybrid_wave_frequency(double omega_ci, double omega_ce);

/**
 * @brief Kinetic Alfven wave frequency (including finite k_perp ρ_i)
 *
 * ω² = k_par² v_A² * (1 + k_perp² ρ_s²) / (1 + k_perp² ρ_i²)
 *
 * where ρ_s = c_s/ω_ci is the ion sound Larmor radius.
 *
 * @param k_par    Parallel wavenumber [rad/m]
 * @param k_perp   Perpendicular wavenumber [rad/m]
 * @param v_A      Alfven speed [m/s]
 * @param rho_i    Ion Larmor radius [m]
 * @param rho_s    Ion sound Larmor radius [m]
 * @return         Angular frequency [rad/s]
 */
double kinetic_alfven_frequency(double k_par, double k_perp, double v_A,
                                double rho_i, double rho_s);

/*──────────────────────────────────────────────────────────────────────
 * L3: Wave Polarization
 *──────────────────────────────────────────────────────────────────────*/

/**
 * @brief Compute wave polarization for cold plasma at given (ω, k, θ)
 *
 * @param omega   Frequency [rad/s]
 * @param k       Wavenumber [rad/m]
 * @param theta   Propagation angle [rad]
 * @param eps     Stix dielectric tensor
 * @param pol     Output: polarization
 * @return        0 on success
 */
int wave_polarization(double omega, double k, double theta,
                      const stix_tensor_t *eps, wave_polarization_t *pol);

/*──────────────────────────────────────────────────────────────────────
 * L3: Resonance and Cutoff Conditions
 *──────────────────────────────────────────────────────────────────────*/

/**
 * @brief Find cutoff frequencies (n → 0)
 *
 * P = 0 cutoff: ω = ω_p   (plasma cutoff)
 * R = 0 cutoff: ω = (ω_ce + sqrt(ω_ce² + 4ω_pe²))/2  (R-cutoff)
 * L = 0 cutoff: ω = (-ω_ce + sqrt(ω_ce² + 4ω_pe²))/2 (L-cutoff)
 *
 * @param omega_pe  Electron plasma frequency [rad/s]
 * @param omega_ce  Electron gyrofrequency [rad/s]
 * @param cutoffs   Output: [0]=P, [1]=R, [2]=L cutoff frequencies [rad/s]
 */
void wave_cutoff_frequencies(double omega_pe, double omega_ce, double cutoffs[3]);

/**
 * @brief Find resonance frequencies (n → ∞)
 *
 * R → ∞: ω = ω_ce (electron cyclotron resonance)
 * S → ∞: ω = ω_LH (lower hybrid resonance)
 * cos²θ/R + sin²θ/S → ∞: ω = ω_UH (upper hybrid, for perp θ=π/2)
 *
 * @param omega_pe  Electron plasma frequency [rad/s]
 * @param omega_ce  Electron gyrofrequency [rad/s]
 * @param omega_pi  Ion plasma frequency [rad/s]
 * @param omega_ci  Ion gyrofrequency [rad/s]
 * @param theta     Propagation angle [rad]
 * @param res       Output: resonance frequencies [rad/s]
 */
void wave_resonance_frequencies(double omega_pe, double omega_ce,
                                double omega_pi, double omega_ci,
                                double theta, double res[3]);

/*──────────────────────────────────────────────────────────────────────
 * L5: CMA (Clemmow-Mullaly-Allis) Diagram Classifier
 *──────────────────────────────────────────────────────────────────────*/

/**
 * @brief Determine wave propagation bands from CMA parameters
 *
 * Parameters: α = ω_pe²/ω², β = ω_ce/ω
 *
 * @param alpha   CMA alpha parameter
 * @param beta    CMA beta parameter
 * @return        Bitmask of propagating wave modes
 */
int cma_propagation_bands(double alpha, double beta);

/**
 * @brief Compute group velocity for given wave mode
 *
 * v_g = dω/dk (numerical derivative via central difference)
 *
 * @param omega       Frequency function f(k,theta)
 * @param k           Wavenumber [rad/m]
 * @param theta       Angle [rad]
 * @param v_g_par     Output: parallel group velocity [m/s]
 * @param v_g_perp    Output: perpendicular group velocity [m/s]
 * @param dk          Finite difference step [rad/m]
 * @param dtheta      Finite difference step [rad]
 */
void group_velocity(double (*omega_fn)(double, double), double k, double theta,
                    double *v_g_par, double *v_g_perp, double dk, double dtheta);

/*──────────────────────────────────────────────────────────────────────
 * L5: General dispersion solver (Newton-Raphson in ω)
 *──────────────────────────────────────────────────────────────────────*/

/**
 * @brief Find frequency for given k using Newton-Raphson
 *
 * Solves D(ω,k,θ) = 0 starting from initial guess ω0.
 *
 * @param D        Dispersion function D(ω,k,θ)
 * @param k        Wavenumber [rad/m]
 * @param theta    Angle [rad]
 * @param omega0   Initial guess [rad/s]
 * @param tol      Convergence tolerance [rad/s]
 * @param max_iter Maximum iterations
 * @return         Converged ω [rad/s], or -1 on failure
 */
double solve_dispersion_omega(double (*D)(double,double,double),
                              double k, double theta, double omega0,
                              double tol, int max_iter);

/**
 * @brief Compute dielectric tensor for warm plasma (electrostatic limit)
 *
 * For Langmuir waves with thermal corrections.
 *
 * @param omega   Frequency [rad/s]
 * @param k       Wavenumber [rad/m]
 * @param omega_p Plasma frequency [rad/s]
 * @param v_th    Thermal velocity [m/s]
 * @return        Dielectric function ε(ω,k)
 */
double warm_plasma_dielectric_function(double omega, double k,
                                       double omega_p, double v_th);

#endif /* PLASMA_WAVES_H */
