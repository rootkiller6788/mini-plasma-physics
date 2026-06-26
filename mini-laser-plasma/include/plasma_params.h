#ifndef PLASMA_PARAMS_H
#define PLASMA_PARAMS_H

/*
 * plasma_params.h -- Fundamental plasma parameter calculations
 *
 * Defines all core plasma state parameters and provides function
 * signatures for computing quantities that characterize a laser-
 * irradiated plasma.
 *
 * References:
 *   - Goldston & Rutherford (1995), Ch. 1-4
 *   - Chen (2016) "Introduction to Plasma Physics and Controlled Fusion"
 *   - Kruer (1988) "The Physics of Laser Plasma Interactions", Ch. 1-2
 *
 * Knowledge Layers:
 *   L1: plasma frequency, Debye length, critical density, beta
 *   L3: plasma dielectric function, collision operator
 *   L4: cold plasma dispersion, Bohm-Gross, Appleton-Hartree
 *
 * Courses: MIT 22.611, Caltech Ph 106, Princeton PHY 525
 */

#include <math.h>

/* ============================================================
 *  Plasma State Structures (L1 Definitions)
 * ============================================================ */

/**
 * PlasmaState carries all independent parameters describing a
 * laser-irradiated plasma slab.
 */
typedef struct {
    double ne;         /* electron density [m^-3]                    */
    double Te;         /* electron temperature [eV]                  */
    double Ti;         /* ion temperature [eV]                       */
    double Z;          /* avg ionization state (dimensionless)       */
    double A;          /* avg ion mass number [amu]                  */
    double B_field;    /* magnetic field [T] (0 = unmagnetized)      */
    double lam_laser;  /* laser vacuum wavelength [m]                */
    double I_laser;    /* laser intensity [W/m^2]                    */
} PlasmaState;

/**
 * PlasmaDerived holds every quantity derivable from PlasmaState.
 * All values in SI units unless noted.
 */
typedef struct {
    double wp;            /* electron plasma frequency [rad/s]       */
    double wpi;           /* ion plasma frequency [rad/s]            */
    double nc;            /* critical density [m^-3]                 */
    double lambda_D;      /* Debye length [m]                        */
    double N_D;           /* particles in Debye sphere               */
    double v_the;         /* electron thermal velocity [m/s]         */
    double v_thi;         /* ion thermal velocity [m/s]              */
    double cs;            /* ion sound speed [m/s]                   */
    double nu_ei;         /* e-i collision frequency [s^-1]          */
    double ln_Lambda;     /* Coulomb logarithm [dimensionless]       */
    double beta;          /* plasma beta = p_kin / p_mag             */
    double a0;            /* normalized vector potential             */
    double ne_over_nc;    /* density ratio                           */
    double lambda_plasma; /* in-plasma wavelength [m]                */
    double skin_depth;    /* collisionless skin depth c/wp [m]       */
    double omega_laser;   /* laser angular frequency [rad/s]         */
} PlasmaDerived;

/* ============================================================
 *  L1: Core Definitions -- each function = one knowledge point
 * ============================================================ */

/**
 * plasma_frequency -- Electron plasma frequency (Langmuir freq.)
 *
 *   omega_p = sqrt(ne e^2 / (eps0 me))
 *
 * Fundamental oscillation frequency of electrons. Sets time scale
 * for all collective plasma phenomena.
 * [Goldston & Rutherford Sec 1.3]
 *
 * Complexity: O(1)
 * Theorem: linearised cold-fluid continuity + Poisson => harmonic eq.
 */
double plasma_frequency(double ne);

/**
 * ion_plasma_frequency -- Ion plasma frequency
 *
 *   omega_pi = sqrt(ni Z^2 e^2 / (eps0 mi))
 *           = omega_p sqrt(Z me / A mp)
 */
double ion_plasma_frequency(double ne, double Z, double A);

/**
 * critical_density -- Critical electron density for EM propagation
 *
 *   nc = eps0 me omega^2 / e^2
 *
 * EM waves propagate only where ne < nc. The single most
 * important threshold in laser-plasma physics.  [Kruer Sec 1.2]
 *
 * Complexity: O(1)
 * Theorem: EM dispersion omega^2 = wp^2 + c^2 k^2; cutoff k->0.
 */
double critical_density(double lambda_m);

/**
 * debye_length -- Electron Debye screening length
 *
 *   lambda_D = sqrt(eps0 kB Te / (ne e^2))
 *
 * Distance over which an external charge is screened by electrons.
 * [Goldston & Rutherford Sec 1.4]
 */
double debye_length(double ne, double Te_eV);

/**
 * debye_sphere_particles -- Particles in a Debye sphere
 *
 *   N_D = (4pi/3) ne lambda_D^3
 *
 * The "plasma parameter". Collective behaviour requires N_D >> 1.
 * [Chen Sec 1.3]
 */
double debye_sphere_particles(double ne, double lambda_D);

/**
 * electron_thermal_velocity -- Maxwellian electron thermal speed
 *
 *   v_the = sqrt(e Te / me)      [Te in eV]
 */
double electron_thermal_velocity(double Te_eV);

/**
 * ion_thermal_velocity -- Maxwellian ion thermal speed
 *
 *   v_thi = sqrt(e Ti / (A mp))   [Ti in eV]
 */
double ion_thermal_velocity(double Ti_eV, double A);

/**
 * ion_sound_speed -- Ion acoustic (Bohm) speed
 *
 *   cs = sqrt(Z kB Te / mi) = sqrt(Z e Te / (A mp))
 *
 * Governs SBS, ion-front dynamics, Bohm-sheath criterion.
 * [Chen Sec 4.6]
 */
double ion_sound_speed(double Te_eV, double Z, double A);

/**
 * skin_depth -- Collisionless EM skin depth
 *
 *   delta = c / omega_p
 *
 * Penetration depth for omega << omega_p into overdense plasma.
 * [Goldston & Rutherford Sec 4.2]
 */
double skin_depth(double wp);

/* ============================================================
 *  L2: Collision Parameters (Core Concepts)
 * ============================================================ */

/**
 * coulomb_logarithm -- Compute ln Lambda
 *
 *   ln Lambda = max(ln(lambda_D / b_min), 2)
 *   b_min = max(Ze^2/(4pi eps0 me v_the^2), hbar/(me v_the))
 *
 * Accounts for cumulative small-angle Coulomb scatterings.
 * [Goldston & Rutherford Sec 2.3]
 */
double coulomb_logarithm(double ne, double Te_eV, double Z);

/**
 * electron_ion_collision_frequency -- Spitzer e-i momentum exchange
 *
 *   nu_ei = [4 sqrt(2)/3] Z ni e^4 lnL
 *           / [sqrt(pi) me^{1/2} (4pi eps0)^2 (kB Te)^{3/2}]
 *
 * Rate for directed electron momentum loss to ions.
 * [Goldston & Rutherford Sec 2.7]
 */
double electron_ion_collision_frequency(double ne, double Te_eV,
                                        double Z, double ln_Lambda);

/**
 * collision_mean_free_path -- Electron mean free path
 *
 *   lambda_mfp = v_the / nu_ei
 */
double collision_mean_free_path(double v_the, double nu_ei);

/**
 * spitzer_resistivity -- Spitzer-Harm parallel resistivity
 *
 *   eta = pi Ze^2 sqrt(me) lnL / [(4pi eps0)^2 (kB Te)^{3/2}]
 *
 * DC resistivity in fully ionized plasma.
 * [Goldston & Rutherford Sec 2.9]
 */
double spitzer_resistivity(double Te_eV, double Z, double ln_Lambda);

/* ============================================================
 *  L3: Plasma Dielectric Response (Mathematical Structures)
 * ============================================================ */

/**
 * cold_plasma_permittivity -- Collisionless cold-plasma dielectric
 *
 *   eps(omega) = 1 - wp^2/omega^2
 *
 * Refractive index N = sqrt(eps). Real for omega > wp,
 * imaginary for omega < wp (evanescence).  [Chen Sec 4.3]
 */
double cold_plasma_permittivity(double omega, double wp);

/**
 * collisional_permittivity -- Drude-model dielectric
 *
 *   eps = 1 - wp^2/(omega^2 + nu^2) + i (nu/omega) wp^2/(omega^2 + nu^2)
 *
 * Returns Re(eps); imag part via *eps_imag.  [Kruer Sec 3.2]
 */
double collisional_permittivity(double omega, double wp, double nu,
                                double *eps_imag);

/**
 * warm_plasma_permittivity -- Bohm-Gross dielectric
 *
 *   eps(k,omega) = 1 - omega_p^2/(omega^2 - 3 k^2 v_the^2)
 *
 * Thermal correction for finite wavenumber k.
 * [Chen Sec 4.10]
 */
double warm_plasma_permittivity(double omega, double k, double wp,
                                double v_the);

/**
 * plasma_refractive_index -- Cold plasma refractive index
 *
 *   N = sqrt(1 - ne/nc)   for ne < nc
 */
double plasma_refractive_index(double ne, double nc);

/**
 * laser_wavelength_in_plasma -- Wavelength stretch inside plasma
 *
 *   lambda_plasma = lambda_vac / N
 */
double laser_wavelength_in_plasma(double lambda_vac, double ne, double nc);

/**
 * group_velocity -- EM group velocity in plasma
 *
 *   v_g / c = sqrt(1 - ne/nc)
 */
double group_velocity(double ne, double nc);

/* ============================================================
 *  L1/L4: Plasma Pressure and Beta
 * ============================================================ */

/**
 * plasma_beta -- Thermal / magnetic pressure ratio
 *
 *   beta = n kB T / (B^2/(2 mu0))
 *
 * beta < 1: magnetically confined; beta > 1: kinetically dominated.
 * [Chen Sec 6.4]
 */
double plasma_beta(double ne, double Te_eV, double Ti_eV, double B_field);

/**
 * electron_pressure -- Electron ideal-gas pressure
 *
 *   pe = ne kB Te
 */
double electron_pressure(double ne, double Te_eV);

/**
 * radiation_pressure -- Laser radiation pressure
 *
 *   P_rad = I/c (absorbing surface), 2I/c (reflecting)
 */
double radiation_pressure(double intensity, int reflecting);

/**
 * compute_all_derived -- Fill PlasmaDerived from PlasmaState
 *
 * Returns 0 on success, -1 for unphysical input.
 * Complexity: O(1)
 */
int compute_all_derived(const PlasmaState *ps, PlasmaDerived *pd);

#endif /* PLASMA_PARAMS_H */
