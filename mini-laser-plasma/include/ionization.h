#ifndef IONIZATION_H
#define IONIZATION_H

/*
 * ionization.h -- Field and collisional ionization models
 *
 * Models for plasma formation through laser-driven ionization:
 * barrier suppression ionization (BSI), tunnel ionization (ADK
 * theory), and collisional (avalanche) ionization.
 *
 * References:
 *   - Ammosov, Delone & Krainov (1986), Sov. Phys. JETP 64, 1191
 *   - Augst et al. (1991), Phys. Rev. A 43, 5031
 *   - Keldysh (1965), Sov. Phys. JETP 20, 1307
 *   - Penetrante & Bardsley (1991), Phys. Rev. A 43, 3100
 *
 * Knowledge Layers:
 *   L1: Ionization potential, Keldysh parameter, ADK rate
 *   L2: BSI, tunnel ionization, multiphoton ionization regimes
 *   L4: Keldysh theory, ADK tunneling formula
 *   L5: Rate equation integration for ionization dynamics
 *
 * Courses: MIT 22.611, Princeton PHY 525, ETH 402-0841
 */

#include "plasma_constants.h"
#include "plasma_params.h"

/* ============================================================
 *  L1: Ionization State Definitions
 * ============================================================ */

/** Atomic data for ionization calculations */
typedef struct {
    int    Z;           /* atomic number                          */
    double Ip_eV;       /* ionization potential for charge Z [eV] */
    double Ip_J;        /* ionization potential [J]               */
    double n_eff;       /* effective principal quantum number     */
    int    l;           /* orbital angular momentum quantum num   */
    int    m;           /* magnetic quantum number                */
} AtomicIonState;

/** Ionization dynamics state for rate-equation integration */
typedef struct {
    double *populations;  /* population fraction for each charge state */
    int     Z_max;        /* maximum atomic number                     */
    double  ne;           /* current electron density [m^-3]           */
    double  Z_avg;        /* average ionization state                  */
    double  time;         /* current time [s]                          */
} IonizationState;

/* ============================================================
 *  L1/L2: Keldysh Parameter (Adiabaticity)
 * ============================================================ */

/**
 * keldysh_parameter -- Distinguishes ionization regimes
 *
 *   gamma_K = omega * sqrt(2 me Ip) / (e E0)
 *           = omega * sqrt(2 Ip / Up)
 *
 * gamma_K << 1 : tunnel ionization (ADK regime)
 * gamma_K >> 1 : multiphoton ionization (MPI regime)
 * gamma_K ~ 1  : transitional
 *
 * This is the single most important dimensionless parameter
 * for strong-field ionization physics.
 * [Keldysh 1965]
 *
 * Complexity: O(1)
 */
double keldysh_parameter(double intensity, double lambda_m, double Ip_eV);

/**
 * barrier_suppression_intensity -- BSI threshold intensity
 *
 *   I_BSI = (pi^2 eps0^3 c / 2) * Ip^4 / (Z^2 e^6)
 *         ~ 4e9 * Ip_eV^4 / Z^2  [W/cm^2]
 *
 * Intensity at which the Coulomb barrier is suppressed below
 * the ionization potential.  Above this, classical over-the-
 * barrier ionization occurs.
 * [Augst et al. 1991]
 */
double barrier_suppression_intensity(double Ip_eV, double Z);

/* ============================================================
 *  L4: ADK Tunneling Ionization Rate
 * ============================================================ */

/**
 * adk_tunneling_rate -- ADK cycle-averaged ionization rate
 *
 *   W_ADK = C_{n*l*} Ip (3E/pi E_atom)^{1/2}
 *           (2E_atom/E)^{2n* - |m| - 1}
 *           exp(-2 E_atom / (3E))
 *
 * where E_atom = (2 Ip)^{3/2} is the atomic field strength,
 * E is the laser electric field, and C_{n*l*} is the ADK
 * coefficient.
 * [Ammosov, Delone & Krainov 1986]
 *
 * Complexity: O(1)
 * Theorem: Quasi-static tunneling in alternating field
 *          (Keldysh-Faisal-Reiss theory, gamma_K << 1 limit)
 */
double adk_tunneling_rate(double intensity, double lambda_m,
                           double Ip_eV, int Z, int l, int m);

/**
 * adk_coefficient -- Compute C_{n*l*} factor
 *
 * Prefactor involving factorial expressions over effective
 * quantum numbers.  [ADK 1986, Eq. 11]
 */
double adk_coefficient(double n_eff, int l, int m);

/**
 * effective_principal_quantum_number -- n* from Ip
 *
 *   n* = Z / sqrt(2 Ip / E_h)   where E_h = Hartree energy
 *
 * Relates the ionization potential to an effective principal
 * quantum number for hydrogenic scaling of ionization rates.
 */
double effective_principal_quantum_number(int Z, double Ip_eV);

/* ============================================================
 *  L2: Collisional (Avalanche) Ionization
 * ============================================================ */

/**
 * avalanche_ionization_rate -- Impact ionization rate
 *
 *   W_av = (nu_ei / Ip) * Up
 *
 * Rate at which free electrons gain enough energy from the laser
 * field to collisionally ionize neutrals in a cascade process.
 * Dominant for long pulses and high densities.
 * [Penetrante & Bardsley 1991]
 */
double avalanche_ionization_rate(double intensity, double lambda_m,
                                  double Ip_eV, double nu_ei);

/**
 * ionization_balance -- Steady-state ionization degree
 *
 * Solves dZ/dt = W_field + (ne sigma v) - (recombination)
 * for the equilibrium Z at given intensity and density.
 */
double ionization_balance(double intensity, double lambda_m,
                           double ne, double Ip_eV, int Z_atom);

/**
 * critical_free_electrons -- Free electrons needed for avalanche
 *
 * Number of seed electrons required to initiate avalanche
 * ionization at a given intensity.
 */
double critical_free_electrons(double intensity, double lambda_m,
                                double Ip_eV);

/* ============================================================
 *  L5: Rate Equation Integration
 * ============================================================ */

/**
 * ionization_state_alloc -- Allocate ionization state
 *
 * Returns NULL on failure.
 */
IonizationState *ionization_state_alloc(int Z_max);

/**
 * ionization_state_free -- Free ionization state
 */
void ionization_state_free(IonizationState *is);

/**
 * integrate_ionization_rate -- Advance ionization state by dt
 *
 * Integrates the coupled rate equations for all charge states
 * over time step dt using a 4th-order Runge-Kutta scheme.
 *
 *   dN_Z/dt = W_{Z-1} N_{Z-1} - (W_Z + R_Z) N_Z + R_{Z+1} N_{Z+1}
 *
 * Complexity: O(Z_max) per step
 */
int integrate_ionization_rate(IonizationState *is, double dt,
                               double intensity, double lambda_m,
                               double Te_eV);

/**
 * ionization_fraction -- Fraction ionized to charge Z
 *
 * Returns N_Z / N_total for charge state Z.
 */
double ionization_fraction(const IonizationState *is, int Z);

#endif /* IONIZATION_H */
