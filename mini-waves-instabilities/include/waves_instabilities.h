/**
 * waves_instabilities.h — Plasma Waves and Instabilities
 *
 * Comprehensive treatment of linear and nonlinear wave phenomena,
 * kinetic and fluid instabilities in magnetized and unmagnetized plasmas.
 *
 * References:
 *   - Stix, "Waves in Plasmas" (1992)
 *   - Swanson, "Plasma Waves" (2003, 2nd ed)
 *   - Krall & Trivelpiece, "Principles of Plasma Physics" (1973)
 *   - Goldston & Rutherford, "Introduction to Plasma Physics" (1995)
 *   - Mikhailovskii, "Theory of Plasma Instabilities" (1974)
 *   - Sagdeev & Galeev, "Nonlinear Plasma Theory" (1969)
 *   - MIT 22.611 · Princeton PHY 521 · Berkeley PHYS 242
 *
 * Knowledge Levels:
 *   L1 — Definitions: plasma frequency, Debye length, growth/damping rates
 *   L2 — Core Concepts: wave-particle resonance, normal modes, instability
 *   L3 — Math Structures: dispersion relation D(w,k)=0, dielectric tensor
 *   L4 — Fundamental Laws: Vlasov-Maxwell, MHD energy principle
 *   L5 — Computational: Brent root finding, shooting method, eigenvalue
 *   L6 — Canonical Systems: Langmuir, IAW, whistler, two-stream, Weibel
 *   L7 — Applications: fusion, space weather, laser-plasma
 *   L8 — Advanced: parametric decay, quasilinear, weak turbulence
 */

#ifndef WAVES_INSTABILITIES_H
#define WAVES_INSTABILITIES_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================
 * Physical Constants (NIST CODATA 2018)
 * ================================================================ */
#define PLASMA_PI          3.14159265358979323846
#define EPSILON0           8.8541878128e-12
#define MU0                1.25663706212e-6
#define E_CHARGE           1.602176634e-19
#define M_ELECTRON         9.1093837015e-31
#define M_PROTON           1.67262192369e-27
#define K_BOLTZMANN        1.380649e-23
#define C_LIGHT            2.99792458e8
#define AMU                1.66053906660e-27

/* ================================================================
 * L1: Core Data Structures
 * ================================================================ */

/** Plasma species descriptor with charge, mass, temperature */
typedef struct {
    double n0;          /**< Equilibrium density [m^-3] */
    double T;           /**< Temperature [J] */
    double m;           /**< Mass [kg] */
    double q;           /**< Charge [C] (signed) */
    double B0;          /**< Background magnetic field [T] */
} PlasmaSpecies;

/** Derived plasma parameters for a single species */
typedef struct {
    double omega_p;     /**< Plasma frequency: sqrt(n q^2/(eps0 m)) [rad/s] */
    double omega_c;     /**< Cyclotron frequency: |q|B/m [rad/s] */
    double v_th;        /**< Thermal velocity: sqrt(2 k_B T/m) [m/s] */
    double lambda_D;    /**< Debye length: sqrt(eps0 k_B T/(n e^2)) [m] */
    double r_L;         /**< Larmor radius: m v_perp/(|q| B) [m] */
    double beta;        /**< Plasma beta: 2 mu0 n k_B T/B^2 */
    double c_s;         /**< Ion sound speed: sqrt(k_B T_e/m_i) [m/s] */
    double v_A;         /**< Alfven speed: B/sqrt(mu0 n m_i) [m/s] */
} PlasmaParams;

/** Complex frequency: omega = omega_r + i*gamma */
typedef struct {
    double omega_r;     /**< Real frequency [rad/s] */
    double gamma;       /**< Growth rate (>0) or damping rate (<0) [1/s] */
} ComplexOmega;

/** Full dispersion solution at one wavenumber */
typedef struct {
    int n_roots;                    /**< Number of roots found */
    ComplexOmega roots[8];          /**< Complex frequency roots */
    double k;                       /**< Wavenumber magnitude [rad/m] */
    double theta;                   /**< Propagation angle from B [rad] */
    int mode_type[8];               /**< Mode classification per root */
} DispersionResult;

/** Wave mode type identification */
enum {
    WAVE_NONE = 0,
    WAVE_LANGMUIR = 1,
    WAVE_ION_ACOUSTIC = 2,
    WAVE_WHISTLER = 3,
    WAVE_ALFVEN = 4,
    WAVE_FAST_MAGNETOSONIC = 5,
    WAVE_SLOW_MAGNETOSONIC = 6,
    WAVE_LOWER_HYBRID = 7,
    WAVE_UPPER_HYBRID = 8,
    WAVE_BERNSTEIN = 9,
    WAVE_ION_CYCLOTRON = 10,
    WAVE_O_MODE = 11,
    WAVE_X_MODE = 12,
    WAVE_DRIFT = 13,
    WAVE_ITG = 14,
    WAVE_TEM = 15,
    WAVE_KINETIC_ALFVEN = 16,
    WAVE_R_CUTOFF = 17,
    WAVE_L_CUTOFF = 18
};

/* ================================================================
 * L1: Core Parameter Computation
 * ================================================================ */

/** Compute all derived plasma parameters from species data */
void compute_plasma_params(const PlasmaSpecies *species,
                           PlasmaParams *params);

/** Generic plasma frequency: omega_p = sqrt(n q^2/(eps0 m)) */
double plasma_frequency(double n, double q, double m);

/** Electron plasma frequency [rad/s]: sqrt(n_e e^2/(eps0 m_e)) */
double electron_plasma_frequency(double n_e);

/** Ion plasma frequency [rad/s] */
double ion_plasma_frequency(double n_i, double Z_i, double m_i);

/** Cyclotron (gyro) frequency [rad/s]: |q| B / m */
double cyclotron_frequency(double q, double B, double m);

/** Debye length [m]: sqrt(eps0 k_B T/(n e^2)) */
double debye_length(double n, double T);

/** Thermal velocity [m/s]: sqrt(2 k_B T/m) */
double thermal_velocity(double T, double m);

/** Ion acoustic speed [m/s]: sqrt(k_B T_e/m_i) */
double ion_sound_speed(double T_e, double m_i);

/** Alfven speed [m/s]: B/sqrt(mu0 n m_i) */
double alfven_speed(double B, double n, double m_i);

/** Plasma beta: 2 mu0 n k_B T/B^2 */
double plasma_beta(double n, double T, double B);

/** Larmor (gyro) radius [m]: m v_perp/(|q| B) */
double larmor_radius(double m, double v_perp, double q, double B);

/** Coulomb logarithm: ln(12 pi n lambda_D^3) */
double coulomb_logarithm(double n, double lambda_D);

/** Electron-ion collision frequency [1/s], Spitzer formula */
double ei_collision_frequency(double n, double T_e, double ln_Lambda);

/* ================================================================
 * L4: Electrostatic Waves — Unmagnetized Plasma
 * ================================================================ */

/**
 * Langmuir wave dispersion — Bohm-Gross relation:
 *   omega^2 = omega_pe^2 + 3 * v_th_e^2 * k^2
 *
 * For k->0: omega->omega_pe. For large k: omega~sqrt(3) v_th_e k.
 *
 * Reference: Bohm & Gross, Phys. Rev. 75, 1851 (1949)
 * Complexity: O(1)
 */
double langmuir_wave_dispersion(double k, double omega_pe, double v_th_e);

/**
 * Langmuir wave with weak Landau damping.
 * Returns complex omega; gamma < 0 means damping.
 *
 * Weak-damping formula (k lambda_De << 1):
 *   omega_r = omega_pe * (1 + 1.5 * k^2 * lambda_De^2)
 *   gamma = -omega_pe * sqrt(pi/8) * (omega_pe/(k v_th_e))^3
 *           * exp(-omega_pe^2/(2 k^2 v_th_e^2) - 1.5)
 *
 * Reference: Landau, J. Phys. USSR 10, 25 (1946)
 * Complexity: O(1)
 */
ComplexOmega langmuir_wave_damped(double k, double omega_pe, double v_th_e);

/**
 * Ion acoustic wave dispersion (fluid):
 *   omega = k * c_s / sqrt(1 + k^2 * lambda_De^2)
 *
 * Sound-like for k*lambda_De << 1: omega ~ k c_s
 * Ion plasma oscillation for k*lambda_De >> 1: omega ~ omega_pi
 * Only propagates when T_e >> T_i (otherwise heavy Landau damping).
 *
 * Complexity: O(1)
 */
double ion_acoustic_wave_dispersion(double k, double c_s, double lambda_De);

/**
 * Ion acoustic wave with electron Landau damping.
 *
 * gamma = -|omega_r| * sqrt(pi/8) *
 *         [sqrt(m_e/m_i) + (T_e/T_i)^(3/2) * exp(-T_e/(2 T_i) - 1.5)]
 *
 * Reference: Chen, "Introduction to Plasma Physics" (2016), Ch. 4
 * Complexity: O(1)
 */
ComplexOmega ion_acoustic_wave_damped(double k, double c_s, double lambda_De,
                                       double T_e, double T_i, double m_i);

/* ================================================================
 * L4: Magnetized Cold Plasma — Stix Framework
 * ================================================================ */

/** Cold plasma dielectric tensor (Stix notation) */
typedef struct {
    double S;  /**< (R+L)/2 — perpendicular dielectric */
    double D;  /**< (R-L)/2 — Hall/conduction term */
    double P;  /**< Parallel dielectric: 1 - sum omega_ps^2/omega^2 */
    double R;  /**< S+D — right-hand circular polarization */
    double L;  /**< S-D — left-hand circular polarization */
} StixDielectric;

/**
 * Compute cold plasma dielectric tensor.
 *
 * For each species s:
 *   eps_s = omega_ps^2/(omega^2 - omega_cs^2)
 *   S = 1 - sum_s eps_s
 *   D = sum_s (omega_cs/omega) * eps_s
 *   P = 1 - sum_s omega_ps^2/omega^2
 *
 * Reference: Stix, "Waves in Plasmas" (1992), Ch. 1-2
 * Complexity: O(n_species)
 */
void stix_dielectric_tensor(double omega,
                            const double wp_species[],
                            const double wc_species[],
                            int n_species,
                            StixDielectric *diel);

/**
 * Cold plasma dispersion relation: A n^4 - B n^2 + C = 0
 *
 * With n = c k/omega:
 *   A = S sin^2(theta) + P cos^2(theta)
 *   B = R L sin^2(theta) + P S (1 + cos^2(theta))
 *   C = P R L
 *
 * Returns n^2 for the two modes. n^2 < 0 => evanescent (set to -1).
 *
 * Reference: Stix (1992), Ch. 2; Swanson (2003), Ch. 3
 * Complexity: O(1)
 */
void cold_plasma_n2(const StixDielectric *diel, double cos_theta,
                    double *n2_mode1, double *n2_mode2);

/**
 * Identify which wave mode propagates at (omega, theta).
 * @return integer wave mode constant
 */
int identify_wave_mode(const StixDielectric *diel,
                       double cos_theta, double omega, double omega_ce);

/** CMA diagram classification */
typedef struct {
    double alpha;   /**< omega_pe^2/omega^2 */
    double beta;    /**< omega_ce/omega */
    int region;     /**< CMA region 0-13 */
    char label[64]; /**< Region description */
} CMARegion;

/** Classify (omega_pe, omega_ce, omega) into CMA diagram region */
CMARegion cma_classify(double omega_pe, double omega_ce, double omega);

/**
 * Check if at cold plasma resonance: A(omega,theta) = 0.
 * @return 1 if resonant, 0 otherwise.
 */
int cold_plasma_is_resonance(double omega, double omega_pe,
                              double omega_ce, double cos_theta);

/**
 * Find cutoff frequencies for cold magnetized plasma.
 * Cutoffs: P=0 => omega=omega_pe; R=0 => omega_R; L=0 => omega_L.
 * Stores up to 3 cutoff values in cutoffs[]; returns count.
 */
int cold_plasma_cutoffs(double omega_pe, double omega_ce,
                        double cutoffs[3]);

/* ================================================================
 * L6: Canonical Wave Modes — Magnetized Plasma
 * ================================================================ */

/**
 * Whistler wave (R-mode, helicon) dispersion.
 *
 * n^2 = 1 + omega_pe^2/(omega*(omega_ce - omega))
 *
 * For omega_ci << omega << omega_ce:
 *   omega = omega_ce * (k*c/omega_pe)^2
 *
 * The quadratic dispersion causes higher frequencies
 * to have larger group velocity -> descending whistler tone.
 *
 * Reference: Storey, Phil. Trans. A 246, 113 (1953)
 * Complexity: O(1)
 */
double whistler_dispersion(double k, double omega_pe, double omega_ce);

/** Whistler group velocity v_g = d(omega)/dk */
double whistler_group_velocity(double k, double omega_pe, double omega_ce);

/**
 * Lower hybrid frequency:
 *   1/omega_LH^2 = 1/(omega_ci*omega_ce) + 1/omega_pi^2
 *
 * For dense plasma with omega_pi^2 >> omega_ci*omega_ce:
 *   omega_LH ~ sqrt(omega_ci*omega_ce)
 *
 * Important for lower hybrid current drive (LHCD) in tokamaks
 * (L7 application: ITER/Fusion).
 */
double lower_hybrid_frequency(double omega_ci, double omega_ce, double omega_pi);

/** Upper hybrid: omega_UH = sqrt(omega_pe^2 + omega_ce^2) */
double upper_hybrid_frequency(double omega_pe, double omega_ce);

/** R-cutoff: omega_R = omega_ce/2 + sqrt(omega_pe^2 + omega_ce^2/4) */
double r_cutoff_frequency(double omega_pe, double omega_ce);

/** L-cutoff: omega_L = -omega_ce/2 + sqrt(omega_pe^2 + omega_ce^2/4) */
double l_cutoff_frequency(double omega_pe, double omega_ce);

/**
 * O-mode (ordinary mode) perpendicular propagation.
 * E parallel to B0. Dispersion same as unmagnetized:
 *   n^2 = 1 - omega_pe^2/omega^2
 */
double omode_dispersion(double k, double omega_pe);

/**
 * X-mode (extraordinary mode) perpendicular propagation.
 * E perpendicular to B0.
 *   n^2 = 1 - (omega_pe^2/omega^2)*(omega^2 - omega_pe^2)/
 *              (omega^2 - omega_UH^2)
 *
 * @param is_propagating output: 1 if n^2 > 0, 0 otherwise
 * @return n^2 value
 */
double xmode_dispersion(double k, double omega_pe, double omega_ce,
                         int *is_propagating);

/**
 * Ion cyclotron wave (L-wave) near omega_ci.
 * Approximate dispersion for omega < omega_ci:
 *   n_parallel^2 = 1 + omega_pi^2/(omega_ci*(omega_ci - omega))
 */
double ion_cyclotron_dispersion(double k_parallel, double omega_pi,
                                 double omega_ci);

/* ================================================================
 * L4: MHD Waves
 * ================================================================ */

/** Shear Alfven wave: omega = k_parallel * v_A */
double alfven_wave_dispersion(double k_parallel, double v_A);

/** Fast magnetosonic speed */
double fast_magnetosonic_speed(double c_s, double v_A, double cos_theta);

/** Slow magnetosonic speed */
double slow_magnetosonic_speed(double c_s, double v_A, double cos_theta);

/** Fast magnetosonic wave: omega = k * v_fast */
double fast_magnetosonic_dispersion(double k, double c_s, double v_A,
                                     double cos_theta);

/** Slow magnetosonic wave: omega = k * v_slow */
double slow_magnetosonic_dispersion(double k, double c_s, double v_A,
                                     double cos_theta);

/** All three MHD wave speeds at given propagation angle */
void mhd_wave_speeds(double c_s, double v_A, double cos_theta,
                     double *v_alfven, double *v_fast, double *v_slow);

/**
 * Kinetic Alfven wave (includes finite Larmor radius effects).
 *
 * omega^2 = k_parallel^2 v_A^2 * [1 + k_perp^2*rho_eff^2]
 *
 * with rho_eff^2 = (3/4)*rho_i^2 + (T_e/T_i)*rho_i^2
 *
 * Reference: Hasegawa, JGR 81, 5083 (1976)
 */
double kinetic_alfven_dispersion(double k_parallel, double k_perp,
                                  double v_A, double rho_i,
                                  double Te_over_Ti);

#ifdef __cplusplus
}
#endif

#endif /* WAVES_INSTABILITIES_H */
