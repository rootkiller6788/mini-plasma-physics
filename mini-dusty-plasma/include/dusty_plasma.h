#ifndef DUSTY_PLASMA_H
#define DUSTY_PLASMA_H
/**
 * @file  dusty_plasma.h
 * @brief Main header for dusty plasma (complex plasma) physics library.
 *
 * Implements the theoretical framework of Shukla & Mamun (2002),
 * "Introduction to Dusty Plasma Physics", with foundations from
 * Goldston & Rutherford (1995) "Introduction to Plasma Physics".
 *
 * Dusty plasma = electron-ion plasma + charged dust grains (nm to um).
 * Key dimensionless parameters: Havnes P, Coulomb coupling Gamma.
 *
 * MIT 22.611 / 22.612 · Stanford PHYSICS 370 · Princeton PHY 535
 */
#include "dusty_constants.h"
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================
 * L1 — Core Type Definitions
 * ================================================================ */

/** Material composition of a dust grain. */
typedef enum {
    DUST_MATERIAL_SILICA = 0,   /**< SiO2 — most common in experiments */
    DUST_MATERIAL_CARBON = 1,   /**< Graphite / amorphous carbon */
    DUST_MATERIAL_IRON   = 2,   /**< Metallic iron */
    DUST_MATERIAL_MF     = 3,   /**< Melamine-formaldehyde */
    DUST_MATERIAL_ICE    = 4,   /**< Water ice — astrophysical */
    DUST_MATERIAL_CUSTOM = 5    /**< User-defined material properties */
} DustMaterial;

/** Charging model selection. */
typedef enum {
    CHARGE_MODEL_OML       = 0, /**< Orbital Motion Limited (Mott-Smith & Langmuir) */
    CHARGE_MODEL_OML_PLUS  = 1, /**< OML + photoelectric emission */
    CHARGE_MODEL_FULL      = 2  /**< Full: OML + photo + secondary + thermionic */
} ChargeModel;

/** Regime classification for a dusty plasma. */
typedef enum {
    DUST_REGIME_ISOLATED    = 0, /**< P << 1, isolated grain limit */
    DUST_REGIME_COLLECTIVE  = 1, /**< P ~ 1, collective effects onset */
    DUST_REGIME_DOMINATED   = 2  /**< P >> 1, dust-dominated plasma */
} DustRegime;

/**
 * @brief Physical description of a single dust grain.
 *
 * Grain radius a determines the regime: sub-micron (quantum effects),
 * micron-sized (classical, OML valid), or macroscopic.
 */
typedef struct {
    double radius;          /**< Grain radius [m] */
    double mass;            /**< Grain mass [kg] */
    double charge;          /**< Instantaneous charge Qd [C] */
    double charge_number;   /**< Charge number Zd = |Qd|/e */
    double material_density;/**< Mass density [kg/m^3] */
    double work_function;   /**< Work function [eV] */
    double surface_temp;    /**< Grain surface temperature [K] */
    DustMaterial material;  /**< Material type */
    double se_yield_max;    /**< Maximum secondary electron yield */
} DustGrain;

/**
 * @brief State variables of a dusty plasma (local fluid element).
 *
 * The plasma is described as a three-component fluid: electrons, ions, dust.
 * Quasineutrality: n_i = n_e + Z_d n_d (with sign convention Z_d > 0).
 */
typedef struct {
    double n_e;             /**< Electron number density [m^-3] */
    double n_i;             /**< Ion number density [m^-3] */
    double n_d;             /**< Dust number density [m^-3] */
    double n_n;             /**< Neutral gas number density [m^-3] */
    double T_e;             /**< Electron temperature [K] (or energy units) */
    double T_i;             /**< Ion temperature [K] */
    double T_d;             /**< Dust kinetic temperature [K] */
    double T_n;             /**< Neutral gas temperature [K] */
    double B;               /**< Magnetic field strength [T] */
    double Z_d;             /**< Average dust charge number */
    double m_i;             /**< Ion mass [kg] */
    double m_n;             /**< Neutral particle mass [kg] */
    double phi_float;       /**< Floating potential of dust [V] */
    double lambda_D;        /**< Total Debye length [m] */
    double lambda_De;       /**< Electron Debye length [m] */
    double lambda_Di;       /**< Ion Debye length [m] */
} DustPlasmaState;

/**
 * @brief Dispersion relation result for a wave mode.
 */
typedef struct {
    double omega_r;         /**< Real frequency [rad/s] */
    double omega_i;         /**< Growth/damping rate [rad/s] (positive = unstable) */
    double k;               /**< Wavenumber [rad/m] */
    double phase_velocity;  /**< v_phi = omega_r / k [m/s] */
    double group_velocity;  /**< v_g = d omega_r / d k [m/s] */
} WaveMode;

/**
 * @brief Yukawa (screened Coulomb) potential between two dust grains.
 *
 * phi(r) = (Q / (4 pi eps0 r)) * exp(-r / lambda_D)
 *
 * This is the fundamental interaction potential in dusty plasmas.
 * For kappa = r / lambda_D -> 0, reduces to bare Coulomb.
 */
typedef struct {
    double Q1;              /**< Charge of first grain [C] */
    double Q2;              /**< Charge of second grain [C] */
    double screening_length;/**< Debye screening length lambda_D [m] */
    double r;               /**< Separation distance [m] */
} YukawaPair;

/**
 * @brief Components of the total current to a dust grain.
 *
 * Sign convention: positive = net current TO the grain.
 * At floating potential: I_total = 0.
 */
typedef struct {
    double I_e;             /**< Electron collection current [A] */
    double I_i;             /**< Ion collection current [A] */
    double I_photo;         /**< Photoelectric emission current [A] */
    double I_secondary;     /**< Secondary electron emission current [A] */
    double I_thermionic;    /**< Thermionic emission current [A] */
    double I_total;         /**< Net current [A] */
} DustChargingCurrent;

/**
 * @brief Force components acting on a single dust grain.
 *
 * Force balance determines grain equilibrium position and dynamics.
 * Key forces: gravity, electrostatic, ion drag, neutral drag.
 */
typedef struct {
    double F_gravity;       /**< Gravitational force [N] (downward positive) */
    double F_electric;      /**< Electrostatic force [N] */
    double F_ion_drag;      /**< Ion drag force [N] (collection + orbital) */
    double F_neutral_drag;  /**< Neutral drag force [N] */
    double F_thermophoretic;/**< Thermophoretic force [N] */
    double F_yukawa;        /**< Net Yukawa force from neighbor grains [N] */
    double F_total;         /**< Net force [N] */
    double F_radiation;     /**< Radiation pressure force [N] */
} DustForceResult;

/**
 * @brief Pair correlation function g(r) for strongly coupled dusty plasma.
 *
 * In the strongly coupled regime (Gamma > 1), the liquid/solid structure
 * is characterized by g(r) and its Fourier transform S(k).
 */
typedef struct {
    size_t n_bins;          /**< Number of radial bins */
    double dr;              /**< Bin width [m] */
    double *g;              /**< Radial distribution function values */
    double *r;              /**< Radial bin centers [m] */
    double Gamma;           /**< Coulomb coupling parameter */
    double kappa;           /**< Screening parameter a/lambda_D */
} PairCorrelation;

/**
 * @brief Dust crystal lattice description (2D hexagonal by default).
 */
typedef struct {
    size_t n_grains;         /**< Number of grains in crystal */
    double lattice_constant; /**< Inter-particle spacing [m] */
    double coupling;         /**< Coulomb coupling parameter Gamma */
    double screening;        /**< Screening parameter kappa */
    double **positions;      /**< Grain positions [n_grains][3] */
    double eq_temp;          /**< Equilibrium kinetic temperature [K] */
} DustCrystal;

/**
 * @brief Particle trajectory state (phase space point).
 */
typedef struct {
    double x, y, z;         /**< Position [m] */
    double vx, vy, vz;      /**< Velocity [m/s] */
    double t;               /**< Time [s] */
} DustTrajectory;

/* ================================================================
 * L2 — Core Plasma Parameter Functions
 * ================================================================ */

/**
 * @brief Compute electron Debye length.
 * lambda_De = sqrt(eps0 * k_B * T_e / (n_e * e^2))
 *
 * Theorem: Debye shielding — the characteristic length over which
 * a test charge is screened in a plasma.
 * Ref: Goldston & Rutherford §2.3
 */
double dust_debye_electron(double n_e, double T_e);

/**
 * @brief Compute ion Debye length.
 * lambda_Di = sqrt(eps0 * k_B * T_i / (n_i * e^2))
 */
double dust_debye_ion(double n_i, double T_i);

/**
 * @brief Compute total Debye length for a dusty plasma.
 * lambda_D^(-2) = lambda_De^(-2) + lambda_Di^(-2)
 *
 * The dust component contribution is often neglected because
 * T_d << T_e, T_i for typical laboratory conditions.
 */
double dust_debye_total(double lambda_De, double lambda_Di);

/**
 * @brief Compute electron plasma frequency.
 * omega_pe = sqrt(n_e * e^2 / (eps0 * m_e))
 */
double dust_plasma_freq_electron(double n_e);

/**
 * @brief Compute ion plasma frequency.
 * omega_pi = sqrt(n_i * e^2 / (eps0 * m_i))
 */
double dust_plasma_freq_ion(double n_i, double m_i);

/**
 * @brief Compute dust plasma frequency.
 * omega_pd = sqrt(n_d * Z_d^2 * e^2 / (eps0 * m_d))
 *
 * This is the fundamental timescale for dust dynamics.
 * For micron-sized grains, omega_pd ~ 10-100 Hz (vs GHz for electrons).
 */
double dust_plasma_freq(double n_d, double Z_d, double m_d);

/**
 * @brief Compute Havnes parameter P.
 * P = Z_d * n_d / n_e
 *
 * P characterizes the relative importance of dust charge:
 *   P < 0.1: isolated grain regime
 *   P ~ 1:   collective effects
 *   P > 10:  dust-dominated plasma
 * Ref: Havnes et al. (1987), J. Geophys. Res.
 */
double dust_havnes_parameter(double Z_d, double n_d, double n_e);

/**
 * @brief Classify the dusty plasma regime based on Havnes P.
 */
DustRegime dust_classify_regime(double P);

/**
 * @brief Compute Coulomb coupling parameter for dust grains.
 * Gamma = Q_d^2 / (4 pi eps0 a k_B T_d)
 *
 * This is the ratio of inter-particle potential energy to thermal energy.
 * Gamma > 1:  strongly coupled (liquid-like)
 * Gamma > 170: crystallization (solid-like, for Yukawa systems)
 * Ref: Ikezi (1986), Phys. Fluids
 */
double dust_coulomb_coupling(double Q_d, double a, double T_d);

/**
 * @brief Compute Yukawa coupling parameter including screening.
 * Gamma* = Gamma * exp(-kappa),  where kappa = a / lambda_D
 *
 * The Yukawa coupling is the relevant parameter for screened systems.
 * Phase diagram: Gamma* vs kappa determines solid/liquid/gas phase.
 * Ref: Hamaguchi et al. (1997), Phys. Rev. E
 */
double dust_yukawa_coupling(double Gamma, double kappa);

/**
 * @brief Check if dust should crystallize based on coupling.
 *
 * For Yukawa systems: crystallization at Gamma* > Gamma_crit(kappa).
 * For kappa < 1: Gamma_crit ~ 170 (pure Coulomb limit).
 * For kappa > 5: Gamma_crit ~ 106 (strong screening limit).
 * Ref: Vaulina & Khrapak (2000), JETP
 */
int dust_crystal_condition(double Gamma, double kappa);

/**
 * @brief Compute critical coupling for Yukawa crystallization.
 *
 * Empirical fit from molecular dynamics simulations:
 * Gamma_crit = 170 * exp(-kappa) + 106 * (1 - exp(-kappa))
 * valid for 0 < kappa < 10.
 */
double dust_critical_coupling(double kappa);

/* ================================================================
 * L2 — Grain Initialization Utilities
 * ================================================================ */

/**
 * @brief Initialize a spherical dust grain.
 *
 * Computes mass from radius and material density assuming spherical shape:
 * m_d = (4/3) * pi * a^3 * rho
 */
DustGrain dust_grain_init(double radius, double density, DustMaterial material);

/**
 * @brief Compute grain mass for a given radius and density.
 * m = (4*pi/3) * rho * a^3
 */
double dust_grain_mass(double radius, double density);

/**
 * @brief Initialize plasma state from key parameters.
 *
 * Automatically computes Debye lengths and enforces quasineutrality
 * by adjusting n_i = n_e + Z_d * n_d.
 */
DustPlasmaState dust_plasma_state_init(
    double n_e, double T_e, double T_i, double T_d,
    double n_d, double Z_d, double m_i, double B);

/**
 * @brief Compute dust acoustic speed.
 * c_da = sqrt(Z_d * k_B * T_e / m_d)
 *
 * In the DAW, electron pressure provides the restoring force
 * while dust provides the inertia. Valid when T_i << Z_d * T_e.
 * Ref: Rao, Shukla & Yu (1990), Planet. Space Sci.
 */
double dust_acoustic_speed(double Z_d, double T_e, double m_d);

/**
 * @brief Compute ion acoustic speed (DIAW phase velocity).
 * c_s = sqrt(k_B * T_e / m_i)
 */
double dust_ion_acoustic_speed(double T_e, double m_i);

/* ================================================================
 * L4-L7 — Sheath and Void Physics (declared here for convenience)
 * ================================================================ */

/**
 * @brief Sheath electric field — exponential model.
 * E(z) = E_0 * exp(-z / lambda_D)
 */
double dust_sheath_electric_field(double z, double E_0, double lambda_D);

/**
 * @brief Sheath potential — exponential model.
 */
double dust_sheath_potential(double z, double phi_wall, double lambda_D);

/**
 * @brief Modified Bohm criterion with dust.
 */
double dust_modified_bohm_velocity(
    double T_e, double T_i, double m_i, double n_d, double n_i);

/**
 * @brief Dust charge in sheath as function of height.
 */
double dust_charge_in_sheath(double z, double Z_d0, double lambda_D);

/**
 * @brief Dust void size estimation.
 */
double dust_void_radius_estimate(
    double lambda_D, double F_id0, double Z_d0, double E_0);

/**
 * @brief Dust density profile in void transition.
 */
double dust_density_void_profile(
    double r, double n_d0, double R_void, double delta);

/**
 * @brief Debye-Hückel potential around a dust grain.
 */
double dust_debye_huckel_potential(
    double r, double a, double phi_s, double lambda_D);

/**
 * @brief Electric field from Debye-Hückel potential.
 */
double dust_debye_huckel_field(
    double r, double a, double phi_s, double lambda_D);

#ifdef __cplusplus
}
#endif

#endif /* DUSTY_PLASMA_H */
