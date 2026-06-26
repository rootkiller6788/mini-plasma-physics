/**
 * plasma_params.c — Fundamental Plasma Parameter Computations
 *
 * Implements all L1 plasma definitions from Goldston & Rutherford,
 * Chen, and standard plasma physics references.
 *
 * References:
 *   - Goldston & Rutherford, "Introduction to Plasma Physics" (1995)
 *   - Chen, "Introduction to Plasma Physics and Controlled Fusion" (2016)
 *   - Spitzer, "Physics of Fully Ionized Gases" (1956)
 *   - NRL Plasma Formulary (2019)
 *
 * Knowledge Coverage:
 *   L1 Definitions: Complete (Debye length, plasma frequency, thermal
 *                    velocity, Alfven speed, plasma beta, gyrofrequency,
 *                    gyroradius, collision frequency, Coulomb logarithm,
 *                    Saha equation)
 *
 * Each function implements one independent knowledge point.
 */
#include "plasma_params.h"
#include "plasma_constants.h"
#include <math.h>
#include <stdlib.h>

/* ================================================================
 * Helper: Maxwellian distribution
 * ================================================================ */

/**
 * Maxwell-Boltzmann 3D velocity distribution
 *
 * f_M(v) = (m / (2*pi*kB*T))^(3/2) * exp(-m*v^2 / (2*kB*T))
 *
 * Returns the probability density at speed v [s^3/m^3].
 *
 * Theorem: In thermal equilibrium, the velocity distribution of a
 * classical ideal gas follows f_M(v). Derived by Maxwell (1860) and
 * Boltzmann (1872) from kinetic theory.
 *
 * Complexity: O(1)
 */
double maxwellian_3d(double v, double T, double m) {
    double a = m / (2.0 * K_B * T);
    double coeff = a / M_PI;
    double coeff_sqrt = sqrt(coeff);
    /* f(v) in 3D includes the 4*pi*v^2 spherical shell factor:
     *   f_M(v) dv = 4*pi * (m/(2*pi*kB*T))^(3/2) * v^2 * exp(-m*v^2/(2*kB*T)) dv
     */
    return 4.0 * M_PI * coeff_sqrt * coeff * v * v * exp(-a * v * v);
}

/**
 * Maxwellian distribution in energy
 *
 * f_E(E) dE = 2/sqrt(pi) * (1/(kB*T))^(3/2) * sqrt(E) * exp(-E/(kB*T)) dE
 *
 * Complexity: O(1)
 */
double maxwellian_energy(double E, double T) {
    if (E < 0.0) return 0.0;
    double kT = K_B * T;
    return (2.0 / sqrt(M_PI)) * pow(1.0 / kT, 1.5)
           * sqrt(E) * exp(-E / kT);
}

/* ================================================================
 * L1.2: Debye Length
 * ================================================================ */

double debye_length(double Te, double ne) {
    if (ne <= 0.0 || Te <= 0.0) return INFINITY;
    return sqrt(EPSILON_0 * K_B * Te / (ne * E_CHARGE * E_CHARGE));
}

double debye_length_electron(double Te, double ne) {
    return debye_length(Te, ne);
}

double debye_length_ion(double Ti, double ni, int Z) {
    if (ni <= 0.0 || Ti <= 0.0 || Z <= 0) return INFINITY;
    double Z_d = (double)Z;
    return sqrt(EPSILON_0 * K_B * Ti / (ni * Z_d * Z_d * E_CHARGE * E_CHARGE));
}

double debye_length_total(double Te, double Ti, double ne, double ni, int Z) {
    double lDe = debye_length_electron(Te, ne);
    double lDi = debye_length_ion(Ti, ni, Z);
    if (lDe == INFINITY || lDi == INFINITY) {
        return (lDe < lDi) ? lDe : lDi;
    }
    /* 1/lambda_D^2 = 1/lambda_De^2 + 1/lambda_Di^2 */
    return 1.0 / sqrt(1.0 / (lDe * lDe) + 1.0 / (lDi * lDi));
}

/* ================================================================
 * L1.3: Plasma Frequency
 * ================================================================ */

double plasma_frequency(double n, double m) {
    if (n <= 0.0 || m <= 0.0) return 0.0;
    return sqrt(n * E_CHARGE * E_CHARGE / (EPSILON_0 * m));
}

double electron_plasma_frequency(double ne) {
    return plasma_frequency(ne, M_ELECTRON);
}

double ion_plasma_frequency(double ni, double mi, int Z) {
    if (ni <= 0.0 || mi <= 0.0 || Z <= 0) return 0.0;
    double Z_d = (double)Z;
    return sqrt(ni * Z_d * Z_d * E_CHARGE * E_CHARGE / (EPSILON_0 * mi));
}

double plasma_frequency_hz(double n, double m) {
    return plasma_frequency(n, m) / (2.0 * M_PI);
}

double lower_hybrid_frequency(double ni, double mi, int Z, double B) {
    double w_pi = ion_plasma_frequency(ni, mi, Z);
    double w_ci = ion_gyrofrequency(B, mi, Z);
    double w_ce = electron_gyrofrequency(B);
    if (w_pi <= 0.0 || w_ci <= 0.0 || w_ce <= 0.0) return 0.0;
    /* 1/w_LH^2 = 1/w_pi^2 + 1/(w_ci * w_ce) */
    double inv_sq = 1.0 / (w_pi * w_pi) + 1.0 / (w_ci * w_ce);
    return 1.0 / sqrt(inv_sq);
}

double upper_hybrid_frequency(double ne, double B) {
    double w_pe = electron_plasma_frequency(ne);
    double w_ce = electron_gyrofrequency(B);
    /* w_UH^2 = w_pe^2 + w_ce^2 */
    return sqrt(w_pe * w_pe + w_ce * w_ce);
}

/* ================================================================
 * L1.4: Thermal Velocity and Sound Speed
 * ================================================================ */

double thermal_velocity(double T, double m) {
    if (T <= 0.0 || m <= 0.0) return 0.0;
    return sqrt(2.0 * K_B * T / m);
}

double electron_thermal_velocity(double Te) {
    return thermal_velocity(Te, M_ELECTRON);
}

double ion_thermal_velocity(double Ti, double mi) {
    return thermal_velocity(Ti, mi);
}

double thermal_velocity_rms(double T, double m) {
    if (T <= 0.0 || m <= 0.0) return 0.0;
    return sqrt(3.0 * K_B * T / m);
}

double thermal_velocity_mean(double T, double m) {
    if (T <= 0.0 || m <= 0.0) return 0.0;
    return sqrt(8.0 * K_B * T / (M_PI * m));
}

double ion_sound_speed(double Te, double Ti, double mi) {
    if (mi <= 0.0) return 0.0;
    /* c_s = sqrt((gamma_e kB Te + gamma_i kB Ti) / mi)
     * Using adiabatic indices: gamma_e = 1, gamma_i = 3 */
    double num = K_B * Te + 3.0 * K_B * Ti;
    if (num < 0.0) num = K_B * Te; /* safeguard if Ti < 0 */
    return sqrt(num / mi);
}

double ion_sound_speed_cold(double Te, double mi) {
    if (Te <= 0.0 || mi <= 0.0) return 0.0;
    return sqrt(K_B * Te / mi);
}

double electron_sound_speed(double Te) {
    if (Te <= 0.0) return 0.0;
    return sqrt(K_B * Te / M_ELECTRON);
}

/* ================================================================
 * L1.5: Alfven Speed and Plasma Beta
 * ================================================================ */

double alfven_speed(double B, double rho) {
    if (rho <= 0.0) return INFINITY;
    return B / sqrt(MU_0 * rho);
}

double plasma_beta(double pressure, double B) {
    if (B == 0.0) return INFINITY;
    return 2.0 * MU_0 * pressure / (B * B);
}

double species_beta(double n, double T, double B) {
    if (B == 0.0) return INFINITY;
    return 2.0 * MU_0 * n * K_B * T / (B * B);
}

/* ================================================================
 * L1.6: Gyromotion Parameters
 * ================================================================ */

double gyrofrequency(double B, double q_abs, double m) {
    if (m <= 0.0) return 0.0;
    return q_abs * B / m;
}

double electron_gyrofrequency(double B) {
    return E_CHARGE * B / M_ELECTRON;
}

double ion_gyrofrequency(double B, double mi, int Z) {
    if (mi <= 0.0 || Z <= 0) return 0.0;
    return (double)Z * E_CHARGE * B / mi;
}

double gyroradius(double v_perp, double B, double q_abs, double m) {
    if (B == 0.0 || q_abs <= 0.0 || m <= 0.0) return INFINITY;
    return m * v_perp / (q_abs * B);
}

double larmor_radius_thermal(double T, double B, double q_abs, double m) {
    if (B == 0.0 || q_abs <= 0.0 || m <= 0.0 || T <= 0.0) return INFINITY;
    return sqrt(m * K_B * T) / (q_abs * B);
}

double electron_larmor_radius_thermal(double Te, double B) {
    return larmor_radius_thermal(Te, B, E_CHARGE, M_ELECTRON);
}

double ion_larmor_radius_thermal(double Ti, double B, double mi, int Z) {
    if (Z <= 0) return INFINITY;
    return larmor_radius_thermal(Ti, B, (double)Z * E_CHARGE, mi);
}

/* ================================================================
 * L1.7: Collision Parameters
 * ================================================================ */

double coulomb_logarithm(double n, double Te) {
    if (n <= 0.0 || Te <= 0.0) return 0.0;
    double lD = debye_length(Te, n);
    /* Number of particles in a Debye sphere */
    double ND = (4.0 * M_PI / 3.0) * n * lD * lD * lD;
    if (ND <= 1.0) return 0.0;
    /* ln Lambda = ln(12*pi*n*lambda_D^3) = ln(9*ND) */
    return log(9.0 * ND);
}

double coulomb_logarithm_from_g(double g) {
    if (g <= 0.0 || g >= 1.0) return 0.0;
    return -log(g);
}

double electron_ion_collision_frequency(double ne, double ni, double Te, int Z) {
    if (ne <= 0.0 || ni <= 0.0 || Te <= 0.0 || Z <= 0) return 0.0;
    double lnL = coulomb_logarithm(ne, Te);
    if (lnL <= 0.0) return 0.0;
    double Z_d = (double)Z;
    double kTe = K_B * Te;
    /* Spitzer electron-ion collision frequency:
     * nu_ei = (4*sqrt(2*pi)/3) * (n_i Z^2 e^4 ln Lambda)
     *         / ((4*pi*epsilon_0)^2 * sqrt(m_e) * (kB*Te)^(3/2)) */
    double num = (4.0 * sqrt(2.0 * M_PI) / 3.0) * ni * Z_d * Z_d
                 * pow(E_CHARGE, 4.0) * lnL;
    double denom = pow(4.0 * M_PI * EPSILON_0, 2.0) * sqrt(M_ELECTRON)
                   * pow(kTe, 1.5);
    return num / denom;
}

double ion_electron_collision_frequency(double ne, double ni, double Te,
                                        double mi, int Z) {
    if (mi <= 0.0) return 0.0;
    double nu_ei = electron_ion_collision_frequency(ne, ni, Te, Z);
    return nu_ei * M_ELECTRON / mi;
}

double electron_mean_free_path(double ne, double ni, double Te, int Z) {
    double v_the = electron_thermal_velocity(Te);
    double nu_ei = electron_ion_collision_frequency(ne, ni, Te, Z);
    if (nu_ei <= 0.0) return INFINITY;
    return v_the / nu_ei;
}

double spitzer_resistivity(double ne, double Te, int Z, double ln_Lambda) {
    if (ne <= 0.0 || Te <= 0.0 || Z <= 0) return INFINITY;
    double kTe = K_B * Te;
    double Z_d = (double)Z;
    /* eta_parallel = (pi * Z * e^2 * sqrt(m_e) * ln Lambda)
     *                / ((4*pi*epsilon_0)^2 * (kB*Te)^(3/2)) */
    double num = M_PI * Z_d * E_CHARGE * E_CHARGE * sqrt(M_ELECTRON) * ln_Lambda;
    double denom = pow(4.0 * M_PI * EPSILON_0, 2.0) * pow(kTe, 1.5);
    return num / denom;
}

/* ================================================================
 * L1.8: Saha Equation and Ionization
 * ================================================================ */

double saha_ionization_fraction(double T, double n_total, double E_ion,
                                 int g_i, int g_n) {
    if (T <= 0.0 || n_total <= 0.0 || g_i <= 0 || g_n <= 0) return 0.0;
    double kT = K_B * T;
    /* Thermal de Broglie wavelength factor:
     * (2*pi*m_e*kB*T / h^2)^(3/2) */
    double lambda_db_factor = pow(2.0 * M_PI * M_ELECTRON * kT
                                  / (H_PLANCK * H_PLANCK), 1.5);
    double g_factor = (double)(g_i * 2) / (double)g_n; /* g_e = 2 */
    double S = g_factor * lambda_db_factor * exp(-E_ion / kT);
    /* S = n_i * n_e / n_n, with n_i = n_e = alpha * n_total, n_n = (1-alpha) * n_total
     * => S = alpha^2 * n_total / (1 - alpha)
     * => alpha^2 * n_total + S * alpha - S = 0
     * => alpha = (-S + sqrt(S^2 + 4*S*n_total)) / (2*n_total)  */
    double a = n_total;
    double b = S;
    double c = -S;
    double disc = b * b - 4.0 * a * c;
    if (disc < 0.0) return 0.0;
    double alpha = (-b + sqrt(disc)) / (2.0 * a);
    if (alpha < 0.0) alpha = 0.0;
    if (alpha > 1.0) alpha = 1.0;
    return alpha;
}

double ionization_percent(double T, double n_total, double E_ion,
                          int g_i, int g_n) {
    return 100.0 * saha_ionization_fraction(T, n_total, E_ion, g_i, g_n);
}

/* ================================================================
 * L2: Quasi-neutrality
 * ================================================================ */

/**
 * Check if a plasma satisfies quasi-neutrality condition.
 *
 * Quasi-neutrality requires: |n_i - n_e| / n_e << 1
 *
 * The condition is valid on length scales L >> lambda_D.
 *
 * Returns 1 if quasi-neutral (relative difference < 1%), 0 otherwise.
 */
int is_quasi_neutral(double ne, double ni, int Z) {
    if (ne <= 0.0) return 0;
    double ni_eff = ni * (double)Z;
    return (fabs(ni_eff - ne) / ne) < 0.01;
}

/**
 * Compute the plasma parameter g = 1/(n lambda_D^3).
 *
 * g << 1: ideal plasma (weakly coupled)
 * g ~ 1: strongly coupled plasma
 * g > 1: non-ideal (liquid-like or solid-like behavior)
 *
 * Complexity: O(1)
 */
double plasma_parameter(double n, double Te) {
    double lD = debye_length(Te, n);
    if (lD == INFINITY) return INFINITY;
    return 1.0 / (n * lD * lD * lD);
}

/**
 * Fill a PlasmaRegime structure for given electron parameters.
 */
void classify_plasma(double ne, double Te, PlasmaRegime *regime) {
    regime->debye_length = debye_length(Te, ne);
    regime->plasma_parameter = plasma_parameter(ne, Te);
    regime->coulomb_logarithm = coulomb_logarithm(ne, Te);
    regime->plasma_frequency = electron_plasma_frequency(ne);
    regime->is_ideal = (regime->plasma_parameter < 0.01) ? 1 : 0;
    regime->is_collisionless = (regime->coulomb_logarithm > 10.0) ? 1 : 0;
}

/* ================================================================
 * L2: Debye Shielding Potential
 * ================================================================ */

/**
 * Screened Coulomb potential (Yukawa potential)
 *
 * phi(r) = (q / (4*pi*epsilon_0 * r)) * exp(-r / lambda_D)
 *
 * This is the potential around a test charge q in a plasma,
 * accounting for Debye shielding by the surrounding charges.
 *
 * For r << lambda_D: phi(r) ~ q / (4*pi*epsilon_0*r) (bare Coulomb)
 * For r >> lambda_D: phi(r) ~ 0 (fully screened)
 *
 * Theorem (Debye-Huckel, 1923): The linearized Poisson-Boltzmann
 * equation for a plasma gives the Yukawa potential.
 */
double debye_shielding_potential(double r, double q, double lambda_D) {
    if (r <= 0.0 || lambda_D <= 0.0) return 0.0;
    return (q / (4.0 * M_PI * EPSILON_0 * r)) * exp(-r / lambda_D);
}

/**
 * Number of particles in a Debye sphere.
 *
 * N_D = (4*pi/3) * n * lambda_D^3
 *
 * For collective behavior to dominate over individual particle
 * interactions, we need N_D >> 1. Typical values:
 *   - Fusion plasma: N_D ~ 10^6 - 10^8
 *   - Space plasma:  N_D ~ 10^10
 *   - Glow discharge: N_D ~ 10^4
 */
double debye_sphere_count(double n, double Te) {
    double lD = debye_length(Te, n);
    if (lD == INFINITY) return 0.0;
    return (4.0 * M_PI / 3.0) * n * lD * lD * lD;
}

/* ================================================================
 * L2: Plasma Parameter Space (Hiroshima diagram)
 * ================================================================ */

/**
 * Compute the plasma parameter space coordinates.
 *
 * This identifies where a given plasma sits in the density-temperature
 * diagram (also known as the "Hiroshima diagram" in plasma physics).
 *
 * @param n Density [m^-3]
 * @param Te Electron temperature [K]
 * @param out_lambda_D Output: Debye length [m]
 * @param out_ND Output: particles per Debye sphere
 * @param out_regime Output string: "fusion", "space", "industrial", etc.
 */
const char* plasma_regime_name(double n, double Te) {
    double lD = debye_length(Te, n);
    double ND = (4.0 * M_PI / 3.0) * n * lD * lD * lD;

    if (Te > 1.0e7 && n > 1.0e19) return "Magnetic Fusion (MCF)";
    if (Te > 1.0e6 && n > 1.0e25) return "Inertial Fusion (ICF)";
    if (Te > 1.0e5 && n < 1.0e10)  return "Space/Solar Wind";
    if (n > 1.0e10 && n < 1.0e16 && Te < 5.0e4) return "Ionosphere";
    if (n > 1.0e15 && Te < 1.0e4) return "Glow Discharge";
    if (n > 1.0e20 && Te < 1.0e4) return "Arc Discharge";
    if (n > 1.0e26 && Te < 1.0e5) return "Warm Dense Matter";
    if (ND < 1.0) return "Strongly Coupled (Non-ideal)";
    if (Te > 1.0e8) return "Relativistic Plasma";
    return "Laboratory Plasma";
}
