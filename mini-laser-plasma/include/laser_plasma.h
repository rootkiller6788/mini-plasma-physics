#ifndef LASER_PLASMA_H
#define LASER_PLASMA_H

/*
 * laser_plasma.h -- Main laser-plasma interaction header
 *
 * Aggregates all sub-system headers and defines the top-level
 * structures for laser-plasma interaction physics.
 *
 * References:
 *   - Kruer (1988) "The Physics of Laser Plasma Interactions"
 *   - Gibbon (2005) "Short Pulse Laser Interactions with Matter"
 *   - Atzeni & Meyer-ter-Vehn (2004) "The Physics of Inertial Fusion"
 *
 * Knowledge Layers:
 *   L1: normalized vector potential, ponderomotive potential
 *   L2: ponderomotive force, relativistic transparency
 *   L4: laser-plasma coupling equations
 *
 * Courses: MIT 22.611, Stanford PHYSICS 370, ETH 402-0891
 */

#include "plasma_constants.h"
#include "plasma_params.h"

/* ============================================================
 *  L1: Laser-Plasma Key Parameters
 * ============================================================ */

/**
 * LaserPulse -- complete description of a laser pulse
 *
 * lambda0  : central vacuum wavelength [m]
 * intensity: peak intensity [W/m^2]
 * tau_fwhm : FWHM pulse duration [s]
 * w0       : 1/e^2 focal spot radius [m]
 * energy   : pulse energy [J]
 * contrast : intensity contrast ratio (ASE/pedestal)
 * pol      : polarization 0=linear, 1=circular
 */
typedef struct {
    double lambda0;
    double intensity;
    double tau_fwhm;
    double w0;
    double energy;
    double contrast;
    int    pol;
} LaserPulse;

/**
 * InteractionGeometry -- laser-target geometry parameters
 *
 * theta_inc : incidence angle from target normal [rad]
 * s_pol     : 0=s-polarized, 1=p-polarized
 * spot_area : focal spot area [m^2]
 * rayleigh  : Rayleigh range [m]
 * f_number  : focusing f-number
 */
typedef struct {
    double theta_inc;
    int    s_pol;
    double spot_area;
    double rayleigh;
    double f_number;
} InteractionGeometry;

/**
 * LaserPlasmaCoupling -- combined laser-plasma interaction state
 */
typedef struct {
    PlasmaState         plasma;
    PlasmaDerived       derived;
    LaserPulse          laser;
    InteractionGeometry geom;
} LaserPlasmaCoupling;

/* ============================================================
 *  L1/L2: Normalized Vector Potential and Ponderomotive Force
 * ============================================================ */

/**
 * normalized_vector_potential -- Dimensionless laser amplitude
 *
 *   a0 = e E0 / (me c omega)
 *      = sqrt(e^2 I lambda^2 / (2 pi^2 eps0 me^2 c^5))
 *      ~ 0.855 * lambda_mum * sqrt(I_18)
 *
 * a0 separates non-relativistic (<1) from relativistic (>1) regimes.
 * This is the single most important dimensionless parameter in
 * high-intensity laser-plasma physics.
 * [Gibbon Sec 2.1, Kruer Sec 6.2]
 *
 * Complexity: O(1)
 */
double normalized_vector_potential(double intensity, double lambda_m);

/**
 * normalized_vector_potential_E0 -- a0 from E-field amplitude
 *
 *   a0 = e E0 / (me c omega)
 *
 * Used when E-field is known directly (e.g., from PIC simulations).
 */
double normalized_vector_potential_E0(double E0, double lambda_m);

/**
 * relativistic_intensity -- Intensity where a0 = 1
 *
 *   I_rel = (pi^2/2) eps0 me^2 c^5 / (e^2 lambda^2)
 *         ~ 1.37e18 / lambda_mum^2  [W/cm^2]
 *
 * At this intensity, the electron quiver energy equals its rest mass
 * energy.  Above this threshold, relativistic effects dominate.
 * [Gibbon Sec 2.1]
 */
double relativistic_intensity(double lambda_m);

/**
 * ponderomotive_potential -- Non-relativistic ponderomotive potential
 *
 *   Phi_p = e^2 E0^2 / (4 me omega^2)
 *         = (e^2 / (2 eps0 me c)) * I / omega^2
 *
 * The time-averaged oscillation energy of an electron in the laser
 * field.  Acts as an effective potential that expels electrons from
 * regions of high intensity (ponderomotive force).
 * [Kruer Sec 6.2, Gibbon Sec 2.1]
 */
double ponderomotive_potential(double intensity, double lambda_m);

/**
 * ponderomotive_potential_relativistic -- Relativistic ponderomotive potential
 *
 *   Phi_p^rel = me c^2 (sqrt(1 + a0^2) - 1)
 *
 * For a0 >> 1, Phi_p^rel ~ me c^2 a0 = e E0 c / omega.
 * This is the correct potential for relativistic laser intensities.
 * [Gibbon Sec 2.1]
 */
double ponderomotive_potential_relativistic(double a0);

/**
 * ponderomotive_force -- Magnitude of ponderomotive force density
 *
 *   f_p = -ne d(Phi_p)/dx ~ -ne e^2/(4 me omega^2) d(E0^2)/dx
 *
 * The ponderomotive force pushes plasma away from regions of high
 * laser intensity, leading to density profile steepening, hole
 * boring, and wakefield excitation.
 */
double ponderomotive_force(double ne, double intensity, double lambda_m,
                           double gradient_scale_length);

/**
 * quiver_velocity -- Electron quiver velocity
 *
 *   v_osc = e E0 / (me omega) = a0 c
 *
 * Peak velocity of an electron oscillating in the laser electric
 * field.  For a0 >= 1, v_osc approaches c.
 */
double quiver_velocity(double a0);

/**
 * quiver_energy -- Electron quiver (ponderomotive) energy
 *
 *   U_p = e^2 E0^2 / (4 me omega^2) = Phi_p
 *
 * Equal to the ponderomotive potential.  Key energy scale for
 * above-threshold ionization (ATI) and HHG.
 * [Gibbon Sec 2.1]
 */
double quiver_energy(double intensity, double lambda_m);

/* ============================================================
 *  L2/L4: Laser Propagation Effects
 * ============================================================ */

/**
 * relativistic_self_focusing_threshold -- Critical power for self-focusing
 *
 *   P_c = 17.3 (nc/ne) GW
 *
 * When laser power exceeds P_c, relativistic mass increase of
 * electrons creates a positive refractive-index gradient on axis,
 * leading to self-focusing (relativistic filamentation).
 * [Gibbon Sec 5.1]
 *
 * Complexity: O(1)
 */
double relativistic_self_focusing_threshold(double ne, double nc);

/**
 * ponderomotive_self_focusing_threshold -- Ponderomotive channel power
 *
 *   P_c ~ 1.2e17 (nc/ne)^{3/2} lambda_mum^2  [W]
 *
 * Threshold above which ponderomotive expulsion of electrons creates
 * a density channel that guides the laser pulse.
 */
double ponderomotive_self_focusing_threshold(double ne, double nc,
                                              double lambda_m);

/**
 * relativistic_transparency_density -- Density for relativistic
 * transparency at given a0
 *
 *   n_c^rel = gamma * n_c = sqrt(1 + a0^2) * n_c
 *
 * Due to relativistic mass increase, the effective plasma frequency
 * decreases by factor 1/gamma, allowing the laser to propagate
 * through formally overdense plasma.
 * [Gibbon Sec 5.3]
 *
 * Complexity: O(1)
 */
double relativistic_transparency_density(double nc, double a0);

/**
 * filamentation_growth_rate -- Maximum growth rate for filamentation
 *
 *   Gamma_fil = (omega_p^2 / 8 omega) a0^2
 *
 * Growth rate of transverse intensity modulations due to
 * relativistic/ponderomotive filamentation instability.
 * [Kruer Sec 8.2]
 */
double filamentation_growth_rate(double wp, double omega, double a0);

/**
 * hole_boring_velocity -- Laser piston velocity
 *
 *   v_hb / c = sqrt[I / (2 mi ni c^3)] for non-relativistic
 *   v_hb / c ~ sqrt[I / (rho c^3)]
 *
 * Velocity at which the laser radiation pressure drills into
 * an overdense plasma (hole boring).
 * [Gibbon Sec 6.1]
 */
double hole_boring_velocity(double intensity, double mi, double ni);

/**
 * plasma_mirror_reflectivity -- Reflectivity of overdense plasma
 *
 *   R = (1 - sqrt(1 - nc/ne))^2 / (1 + sqrt(1 - nc/ne))^2
 *     for ne > nc (normal incidence)
 *
 * Fresnel-like formula for plasma mirror at ne > nc.
 */
double plasma_mirror_reflectivity(double ne, double nc);

#endif /* LASER_PLASMA_H */
