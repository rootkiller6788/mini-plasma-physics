#ifndef PLASMA_CONSTANTS_H
#define PLASMA_CONSTANTS_H

/*
 * plasma_constants.h -- Physical constants for laser-plasma physics
 *
 * Source: CODATA 2018 recommended values (SI units).
 * References:
 *   - Goldston & Rutherford (1995) "Introduction to Plasma Physics", Ch. 1
 *   - Kruer (1988) "The Physics of Laser Plasma Interactions", App. A
 *   - Gibbon (2005) "Short Pulse Laser Interactions with Matter", Ch. 1
 *
 * Knowledge Layer: L1 (Definitions)
 * Courses: MIT 22.611, Princeton PHY 525, Stanford PHYSICS 370
 */

/* --- Fundamental Constants (CODATA 2018, SI units) --- */

/** Speed of light in vacuum [m/s] */
#define PLASMA_C        2.99792458e8

/** Elementary charge [C] */
#define PLASMA_E        1.602176634e-19

/** Electron rest mass [kg] */
#define PLASMA_ME       9.1093837015e-31

/** Proton rest mass [kg] */
#define PLASMA_MP       1.67262192369e-27

/** Vacuum permittivity [F/m] */
#define PLASMA_EPS0     8.8541878128e-12

/** Vacuum permeability [N/A^2] */
#define PLASMA_MU0      1.25663706212e-6

/** Boltzmann constant [J/K] */
#define PLASMA_KB       1.380649e-23

/** Planck constant [J.s] */
#define PLASMA_H        6.62607015e-34

/** Reduced Planck constant [J.s] */
#define PLASMA_HBAR     1.054571817e-34

/** Stefan-Boltzmann constant [W/(m^2.K^4)] */
#define PLASMA_SIGMA_SB 5.670374419e-8

/** Avogadro number [mol^-1] */
#define PLASMA_NA       6.02214076e23

/* --- Derived Plasma Scales --- */

/**
 * Classical electron radius [m]
 * r_e = e^2 / (4.pi.eps_0.m_e.c^2) ~ 2.81794e-15 m
 * Electrostatic self-energy equals rest mass energy scale.
 */
#define PLASMA_RE       ((PLASMA_E * PLASMA_E) / \
                        (4.0 * M_PI * PLASMA_EPS0 * PLASMA_ME * PLASMA_C * PLASMA_C))

/**
 * Thomson cross-section [m^2]
 * sigma_T = (8pi/3) r_e^2 ~ 6.6524e-29 m^2
 * Classical scattering cross-section for a free electron.
 */
#define PLASMA_SIGMA_T  ((8.0 * M_PI / 3.0) * PLASMA_RE * PLASMA_RE)

/**
 * Bohr radius [m]
 * a_0 = 4.pi.eps_0.hbar^2 / (m_e e^2) ~ 5.29177e-11 m
 * Atomic length unit, reference for tunneling ionization.
 */
#define PLASMA_A0       ((4.0 * M_PI * PLASMA_EPS0 * PLASMA_HBAR * PLASMA_HBAR) \
                        / (PLASMA_ME * PLASMA_E * PLASMA_E))

/**
 * Hartree energy [J]
 * E_h = m_e e^4 / (2 (4.pi.eps_0.hbar)^2) ~ 4.3597e-18 J = 27.2114 eV
 * Atomic energy unit for field ionization thresholds.
 */
#define PLASMA_EH       ((PLASMA_ME * PLASMA_E * PLASMA_E * PLASMA_E * PLASMA_E) \
                        / (2.0 * (4.0 * M_PI * PLASMA_EPS0 * PLASMA_HBAR) \
                           * (4.0 * M_PI * PLASMA_EPS0 * PLASMA_HBAR)))

/** Laser wavelength to angular frequency [rad/s]: omega = 2.pi.c / lambda */
#define LAMBDA_TO_OMEGA(lam)  (2.0 * M_PI * PLASMA_C / (lam))

/** Intensity from E0 [W/m^2]: I = eps_0 c |E|^2 / 2 */
#define INTENSITY_FROM_E0(E0)  (0.5 * PLASMA_EPS0 * PLASMA_C * (E0) * (E0))

/** E0 from intensity [V/m]: E = sqrt(2 I / (eps_0 c)) */
#define E0_FROM_INTENSITY(I)  (sqrt(2.0 * (I) / (PLASMA_EPS0 * PLASMA_C)))

/** Electronvolt to Joules */
#define EV_TO_J(eV)     ((eV) * PLASMA_E)

/** Joules to electronvolts */
#define J_TO_EV(J)      ((J) / PLASMA_E)

#endif /* PLASMA_CONSTANTS_H */
