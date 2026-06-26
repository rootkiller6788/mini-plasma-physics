/**
 * plasma_constants.h — Physical Constants for Plasma Physics
 *
 * Reference: NIST CODATA 2018
 * Course: MIT 22.611 / Princeton PHY 521 / Oxford CMT / ETH 402-0841
 *
 * All values in SI units unless otherwise noted.
 */
#ifndef PLASMA_CONSTANTS_H
#define PLASMA_CONSTANTS_H

#include <math.h>

/* Pi (not always defined in strict C11) */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---------- Fundamental Constants ---------- */

/** Boltzmann constant [J/K] */
#define K_B         1.380649e-23

/** Reduced Planck constant [J·s] */
#define HBAR        1.054571817e-34

/** Speed of light in vacuum [m/s] */
#define C_LIGHT     2.99792458e8

/** Elementary charge [C] */
#define E_CHARGE    1.602176634e-19

/** Electron rest mass [kg] */
#define M_ELECTRON  9.1093837015e-31

/** Proton rest mass [kg] */
#define M_PROTON    1.67262192369e-27

/** Vacuum permittivity [F/m] */
#define EPSILON_0   8.8541878128e-12

/** Vacuum permeability [N/A²] */
#define MU_0        1.25663706212e-6

/** Avogadro constant [1/mol] */
#define N_AVOGADRO  6.02214076e23

/** Planck constant [J·s] */
#define H_PLANCK    6.62607015e-34

/** Stefan-Boltzmann constant [W/(m²·K⁴)] */
#define SIGMA_SB    5.670374419e-8

/** Universal gas constant [J/(mol·K)] */
#define R_GAS       8.314462618

/** Gravitational constant [N·m²/kg²] */
#define G_GRAV      6.67430e-11

/* ---------- Plasma Reference Values ---------- */

/** Tokamak core density (JET/ITER-like) [m⁻³] */
#define N_FUSION    1.0e20

/** Tokamak core temperature [K] ~10 keV */
#define T_FUSION    1.16e8

/** Solar wind density at 1 AU [m⁻³] */
#define N_SOLARWIND 5.0e6

/** Solar wind temperature [K] */
#define T_SOLARWIND 1.0e5

/** Ionosphere F-layer density [m⁻³] */
#define N_IONOSPHERE 1.0e12

/** Ionosphere temperature [K] */
#define T_IONOSPHERE 1.0e3

#endif /* PLASMA_CONSTANTS_H */
