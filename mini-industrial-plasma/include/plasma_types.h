/**
 * plasma_types.h — Core Plasma Data Types for Industrial Plasma Physics
 *
 * Reference: Lieberman & Lichtenberg, "Principles of Plasma Discharges
 *   and Materials Processing" (2nd Ed, 2005), Ch.1-4
 * Course: MIT 22.611 / Berkeley EECS 245 / Stanford EE 414
 *
 * Defines the fundamental data structures for low-temperature,
 * industrial processing plasmas:
 *   - Plasma state (density, temperature, potentials)
 *   - Species (electrons, ions, neutrals, radicals)
 *   - Collision cross sections and rate coefficients
 *   - Discharge geometry and operating parameters
 */

#ifndef PLASMA_TYPES_H
#define PLASMA_TYPES_H

#include <stddef.h>
#include <stdint.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_E
#define M_E  2.71828182845904523536
#endif

/* ============================================================
 * Physical Constants (CODATA 2018)
 * ============================================================ */

#define PLASMA_E_CHARGE   1.602176634e-19
#define PLASMA_M_E        9.1093837015e-31
#define PLASMA_M_P        1.67262192369e-27
#define PLASMA_M_AR       6.6335209e-26
#define PLASMA_EPS0       8.8541878128e-12
#define PLASMA_K_B        1.380649e-23
#define PLASMA_HBAR       1.054571817e-34
#define PLASMA_A0         5.29177210903e-11

/* ============================================================
 * L1: Core Definitions — Plasma State
 * ============================================================ */

typedef enum {
    SPECIES_ELECTRON = 0,
    SPECIES_ION_POS   = 1,
    SPECIES_ION_NEG   = 2,
    SPECIES_NEUTRAL   = 3,
    SPECIES_RADICAL   = 4,
    SPECIES_PHOTON    = 5,
    SPECIES_METASTABLE = 6,
    SPECIES_DIMER     = 7,
    SPECIES_EXCITED   = 8
} PlasmaSpeciesType;

typedef struct {
    double n_e;
    double n_i;
    double n_n;
    double T_e;
    double T_i;
    double T_g;
    double V_plasma;
    double V_float;
    double lam_D;
    double omega_pe;
    double nu_en;
    double nu_ei;
    double nu_iz;
    double pressure;
    double B_field;
} PlasmaState;

typedef struct {
    char   name[32];
    double mass;
    int    charge_number;
    double ionization_energy;
    double excitation_energy;
    double polarizability;
    PlasmaSpeciesType type;
} PlasmaSpecies;

typedef enum {
    REACTION_IONIZATION      = 0,
    REACTION_EXCITATION      = 1,
    REACTION_DISSOCIATION    = 2,
    REACTION_RECOMBINATION   = 3,
    REACTION_ATTACHMENT      = 4,
    REACTION_DETACHMENT      = 5,
    REACTION_ELASTIC         = 6,
    REACTION_CHARGE_EXCHANGE = 7,
    REACTION_PENNING         = 8,
    REACTION_RADICAL_RECOMB  = 9
} ReactionType;

typedef struct {
    char         equation[128];
    ReactionType type;
    double       threshold;
    double       rate_constant_A;
    double       rate_constant_n;
    double       cross_section_max;
    double       energy_at_max;
    int          reactant_ids[4];
    int          product_ids[4];
    int          n_reactants;
    int          n_products;
} PlasmaReaction;

typedef struct {
    double sigma0;
    double eps0;
    double a;
    double sigma1;
    double eps1;
} CrossSectionModel;

typedef struct {
    double electrode_area;
    double ground_area;
    double gap_length;
    double chamber_radius;
    double chamber_height;
    double wafer_diameter;
    double pump_speed;
    double gas_flow_rate;
} DischargeGeometry;

typedef struct {
    double power;
    double frequency;
    double dc_bias;
    double voltage_amplitude;
    double match_efficiency;
    int    is_pulsed;
    double pulse_duty;
    double pulse_freq;
} DischargeParams;

typedef struct {
    double mobility;
    double diffusion;
    double ambipolar_diff;
    double thermal_cond;
    double viscosity;
    double elec_conductivity;
} TransportCoeffs;

typedef struct {
    int     n_species;
    int     n_reactions;
    double *concentrations;
    double *production_rates;
    double *loss_rates;
    double *rate_constants;
} RateMatrix;

typedef enum {
    EEDF_SOLVER_MAXWELLIAN = 0,
    EEDF_SOLVER_DRUYVESTEYN = 1,
    EEDF_SOLVER_BOLTZMANN_TWO_TERM = 2,
    EEDF_SOLVER_MONTE_CARLO = 3,
    EEDF_SOLVER_DISCRETE_ORDINATES = 4
} EEDFSolverType;

typedef struct {
    int      n_grid;
    double   emin, emax;
    double  *energy_grid;
    double  *f0;
    double  *f1;
    double   mean_energy;
    double   E_over_N;
    EEDFSolverType solver;
} EEDFState;

typedef struct {
    int      n_points;
    double  *x;
    double  *potential;
    double  *ion_density;
    double  *electron_density;
    double  *electric_field;
    double   sheath_width;
    double   ion_flux;
    double   ion_energy;
    double   V_wall;
} SheathSolution;

typedef struct {
    int     n_points;
    double *pd_values;
    double *V_breakdown;
    double  V_min;
    double  pd_min;
    double  gamma_se;
    double  A_coeff;
    double  B_coeff;
} PaschenCurve;

typedef struct {
    double etch_rate_chemical;
    double etch_rate_physical;
    double etch_rate_total;
    double anisotropy;
    double selectivity_mask;
    double selectivity_etchstop;
    double uniformity_3sigma;
    double aspect_ratio;
} EtchModel;

typedef struct {
    double dep_rate;
    double refractive_index;
    double stress;
    double hydrogen_content;
    double uniformity;
    double step_coverage;
} DepositionModel;

typedef struct {
    double n_neg;
    double alpha;
    double T_neg;
    double detachment_rate;
    double ion_ion_plasma_freq;
} ElectronegativeParams;

/* inline helpers */
static inline double eV_to_K(double eV) { return eV * 11604.518; }

static inline double K_to_eV(double K)  { return K / 11604.518; }

static inline double debye_length_electron(double n_e, double T_e_eV) {
    if (n_e <= 0.0 || T_e_eV <= 0.0) return INFINITY;
    return sqrt(PLASMA_EPS0 * T_e_eV / (n_e * PLASMA_E_CHARGE));
}

static inline double plasma_freq_electron(double n_e) {
    if (n_e <= 0.0) return 0.0;
    return sqrt(n_e * PLASMA_E_CHARGE * PLASMA_E_CHARGE
               / (PLASMA_EPS0 * PLASMA_M_E));
}

static inline double plasma_freq_ion(double n_i, double m_i) {
    if (n_i <= 0.0 || m_i <= 0.0) return 0.0;
    return sqrt(n_i * PLASMA_E_CHARGE * PLASMA_E_CHARGE / (PLASMA_EPS0 * m_i));
}

static inline double mean_free_path(double n_gas, double cross_section) {
    if (n_gas <= 0.0 || cross_section <= 0.0) return INFINITY;
    return 1.0 / (n_gas * cross_section);
}

#endif /* PLASMA_TYPES_H */
