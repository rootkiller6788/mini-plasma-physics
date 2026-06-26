/**
 * plasma_params.h — Fundamental Plasma Parameters
 *
 * Defines the core plasma parameters that characterize any plasma system.
 * These are the L1 definitions: every plasma physicist must know these.
 *
 * References:
 *   - Goldston & Rutherford, "Introduction to Plasma Physics" (1995), Ch. 1-2
 *   - Chen, "Introduction to Plasma Physics and Controlled Fusion" (2016), Ch. 1-4
 *   - MIT 22.611 / Princeton PHY 521 / Berkeley PHYS 242
 *
 * Knowledge Level: L1 — Definitions (Complete)
 *
 * Convention: All formulas use SI units.
 *   T  — temperature [K] (multiply by K_B for [J])
 *   n  — number density [m^-3]
 *   B  — magnetic field [T]
 *   m  — particle mass [kg]
 *   q  — particle charge [C]
 *   Z  — ion charge state
 */

#ifndef PLASMA_PARAMS_H
#define PLASMA_PARAMS_H

#include "plasma_constants.h"
#include <math.h>

/* ================================================================
 * L1.1: Plasma Parameter Structures
 * ================================================================ */

/** Plasma species descriptor */
typedef struct {
    double mass;        /**< Particle mass [kg] */
    double charge;      /**< Particle charge [C] (signed) */
    double density;     /**< Number density [m^-3] */
    double temperature; /**< Temperature [K] */
    int    Z;           /**< Charge state (|charge|/e, unsigned) */
    double A;           /**< Mass number (mass / proton_mass) */
} PlasmaSpecies;

/** Collective plasma state */
typedef struct {
    PlasmaSpecies *species;  /**< Array of plasma species */
    int            n_species;/**< Number of species */
    double         B;        /**< Magnetic field magnitude [T] */
    double         pressure; /**< Total plasma pressure [Pa] */
    double         volume;   /**< Plasma volume [m^3] */
} PlasmaState;

/** Dimensionless plasma regime classification */
typedef struct {
    double debye_length;      /**< lambda_D, Debye length [m] */
    double plasma_parameter;  /**< g = 1/(n lambda_D^3) */
    double coulomb_logarithm; /**< ln Lambda */
    double plasma_frequency;  /**< omega_p [rad/s] */
    int    is_ideal;          /**< 1 if g << 1 */
    int    is_collisionless;  /**< 1 if collisionless */
} PlasmaRegime;

/* ================================================================
 * L1.2: Debye Shielding Parameters
 * ================================================================ */

/**
 * Debye length lambda_D [m]
 * lambda_D = sqrt(epsilon_0 k_B T_e / (n_e e^2))
 */
double debye_length(double Te, double ne);

/** Electron Debye length [m] */
double debye_length_electron(double Te, double ne);

/** Ion Debye length [m] */
double debye_length_ion(double Ti, double ni, int Z);

/** Total Debye length: 1/lambda_D^2 = 1/lambda_De^2 + 1/lambda_Di^2 */
double debye_length_total(double Te, double Ti, double ne, double ni, int Z);

/* ================================================================
 * L1.3: Plasma Frequency
 * ================================================================ */

/** Plasma frequency omega_p = sqrt(n e^2/(epsilon_0 m)) [rad/s] */
double plasma_frequency(double n, double m);

/** Electron plasma frequency [rad/s] */
double electron_plasma_frequency(double ne);

/** Ion plasma frequency [rad/s] */
double ion_plasma_frequency(double ni, double mi, int Z);

/** Plasma frequency f_p [Hz] */
double plasma_frequency_hz(double n, double m);

/** Lower hybrid frequency [rad/s] */
double lower_hybrid_frequency(double ni, double mi, int Z, double B);

/** Upper hybrid frequency [rad/s] */
double upper_hybrid_frequency(double ne, double B);

/* ================================================================
 * L1.4: Thermal Velocity and Sound Speed
 * ================================================================ */

/** Thermal velocity v_th = sqrt(2 kB T / m) [m/s] */
double thermal_velocity(double T, double m);

/** Electron thermal velocity [m/s] */
double electron_thermal_velocity(double Te);

/** Ion thermal velocity [m/s] */
double ion_thermal_velocity(double Ti, double mi);

/** RMS thermal speed v_rms = sqrt(3 kB T / m) [m/s] */
double thermal_velocity_rms(double T, double m);

/** Mean thermal speed v_mean = sqrt(8 kB T/(pi m)) [m/s] */
double thermal_velocity_mean(double T, double m);

/** Ion sound speed (Bohm speed) c_s [m/s] */
double ion_sound_speed(double Te, double Ti, double mi);

/** Ion sound speed, cold ion limit [m/s] */
double ion_sound_speed_cold(double Te, double mi);

/** Electron sound speed [m/s] */
double electron_sound_speed(double Te);

/* ================================================================
 * L1.5: Alfven Speed and Plasma Beta
 * ================================================================ */

/** Alfven speed v_A = B/sqrt(mu_0 rho) [m/s] */
double alfven_speed(double B, double rho);

/** Plasma beta = 2 mu_0 p / B^2 */
double plasma_beta(double pressure, double B);

/** Species beta = 2 mu_0 n kB T / B^2 */
double species_beta(double n, double T, double B);

/* ================================================================
 * L1.6: Gyromotion Parameters
 * ================================================================ */

/** Cyclotron frequency omega_c = |q| B / m [rad/s] */
double gyrofrequency(double B, double q_abs, double m);

/** Electron cyclotron frequency [rad/s] */
double electron_gyrofrequency(double B);

/** Ion cyclotron frequency [rad/s] */
double ion_gyrofrequency(double B, double mi, int Z);

/** Gyroradius r_L = m v_perp / (|q| B) [m] */
double gyroradius(double v_perp, double B, double q_abs, double m);

/** Thermal Larmor radius rho_L = sqrt(m kB T)/(|q| B) [m] */
double larmor_radius_thermal(double T, double B, double q_abs, double m);

/** Electron thermal Larmor radius [m] */
double electron_larmor_radius_thermal(double Te, double B);

/** Ion thermal Larmor radius [m] */
double ion_larmor_radius_thermal(double Ti, double B, double mi, int Z);

/* ================================================================
 * L1.7: Collision Parameters
 * ================================================================ */

/** Coulomb logarithm ln Lambda = ln(12 pi n lambda_D^3) */
double coulomb_logarithm(double n, double Te);

/** Coulomb logarithm from plasma parameter g */
double coulomb_logarithm_from_g(double g);

/** Electron-ion collision frequency (Spitzer) [1/s] */
double electron_ion_collision_frequency(double ne, double ni, double Te, int Z);

/** Ion-electron collision frequency [1/s] */
double ion_electron_collision_frequency(double ne, double ni, double Te,
                                        double mi, int Z);

/** Electron mean free path [m] */
double electron_mean_free_path(double ne, double ni, double Te, int Z);

/** Spitzer resistivity [Ohm*m] */
double spitzer_resistivity(double ne, double Te, int Z, double ln_Lambda);

/* ================================================================
 * L1.8: Saha Equation
 * ================================================================ */

/** Saha ionization fraction alpha in [0, 1] */
double saha_ionization_fraction(double T, double n_total, double E_ion,
                                 int g_i, int g_n);

/** Ionization percent [0, 100] */
double ionization_percent(double T, double n_total, double E_ion,
                          int g_i, int g_n);

/* ================================================================
 * L2: Collective Plasma Properties
 * ================================================================ */

/** Check quasi-neutrality: |n_i - n_e|/n_e < 1% */
int is_quasi_neutral(double ne, double ni, int Z);

/** Plasma parameter g = 1/(n lambda_D^3) */
double plasma_parameter(double n, double Te);

/** Classify plasma into regime */
void classify_plasma(double ne, double Te, PlasmaRegime *regime);

/** Debye shielded potential (Yukawa) */
double debye_shielding_potential(double r, double q, double lambda_D);

/** Number of particles in Debye sphere */
double debye_sphere_count(double n, double Te);

/** Plasma regime name string */
const char* plasma_regime_name(double n, double Te);

/* Maxwellian 3D distribution f_M(v) */
double maxwellian_3d(double v, double T, double m);

/* Maxwellian energy distribution f_E(E) */
double maxwellian_energy(double E, double T);

#endif /* PLASMA_PARAMS_H */
