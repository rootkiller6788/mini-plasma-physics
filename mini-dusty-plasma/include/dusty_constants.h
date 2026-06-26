#ifndef DUSTY_CONSTANTS_H
#define DUSTY_CONSTANTS_H
/**
 * @file  dusty_constants.h
 * @brief Physical constants for dusty plasma physics.
 *
 * Reference: Goldston & Rutherford (1995), "Introduction to Plasma Physics"
 *            Shukla & Mamun (2002), "Introduction to Dusty Plasma Physics"
 *            CODATA 2018 recommended values.
 */

/* --- Fundamental Constants (SI units) --- */

/** Boltzmann constant [J/K] */
#define DUSTY_KB           1.380649e-23

/** Elementary charge [C] */
#define DUSTY_EC           1.602176634e-19

/** Electron mass [kg] */
#define DUSTY_ME           9.1093837015e-31

/** Proton mass [kg] */
#define DUSTY_MP           1.67262192369e-27

/** Vacuum permittivity [F/m] */
#define DUSTY_EPS0         8.8541878128e-12

/** Vacuum permeability [N/A^2] */
#define DUSTY_MU0          1.25663706212e-6

/** Speed of light [m/s] */
#define DUSTY_C            299792458.0

/** Planck constant [J s] */
#define DUSTY_H            6.62607015e-34

/** Reduced Planck constant [J s] */
#define DUSTY_HBAR         1.054571817e-34

/** Avogadro constant [1/mol] */
#define DUSTY_NA           6.02214076e23

/** Standard gravitational acceleration [m/s^2] */
#define DUSTY_G0           9.80665

/** Gravitational constant [N m^2/kg^2] */
#define DUSTY_G            6.67430e-11

/* --- Plasma Reference Values --- */

/** Electron temperature 1 eV in Kelvin [K] */
#define DUSTY_EV_IN_K      11604.518

/** Classical electron radius [m] */
#define DUSTY_RE_CLASSICAL 2.8179403262e-15

/** Bohr radius [m] */
#define DUSTY_A0           5.29177210903e-11

/** Ionization energy of hydrogen [J] */
#define DUSTY_EI_HYDROGEN  2.179e-18

/** Ionization energy of argon [J] */
#define DUSTY_EI_ARGON     2.523e-18

/* --- Dust Material Properties --- */

/** Mass density of silica (SiO2) dust [kg/m^3] */
#define DUSTY_RHO_SILICA   2200.0

/** Mass density of carbon dust [kg/m^3] */
#define DUSTY_RHO_CARBON   2000.0

/** Mass density of iron dust [kg/m^3] */
#define DUSTY_RHO_IRON     7874.0

/** Mass density of melamine-formaldehyde [kg/m^3] */
#define DUSTY_RHO_MF       1510.0

/** Work function of silica [eV] */
#define DUSTY_WF_SILICA    5.0

/** Work function of carbon (graphite) [eV] */
#define DUSTY_WF_CARBON    4.7

/** Secondary electron yield maximum for SiO2 */
#define DUSTY_SEY_SILICA   2.5

/* --- Mathematical Constants --- */

/** Pi */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/** Euler's number */
#ifndef M_E
#define M_E  2.71828182845904523536
#endif

#endif /* DUSTY_CONSTANTS_H */
