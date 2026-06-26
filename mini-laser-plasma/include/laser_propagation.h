#ifndef LASER_PROPAGATION_H
#define LASER_PROPAGATION_H

/*
 * laser_propagation.h -- Laser propagation in plasma
 *
 * Models for electromagnetic wave propagation in underdense and
 * near-critical plasmas: refraction, diffraction, self-focusing,
 * and relativistic transparency.
 *
 * References:
 *   - Gibbon (2005) "Short Pulse Laser Interactions with Matter", Ch. 5
 *   - Kruer (1988), Ch. 6-8
 *   - Max (1976), "Strong Self-Focusing of EM Waves in Plasmas"
 *
 * Knowledge Layers:
 *   L2: Self-focusing, relativistic transparency, hole boring
 *   L3: Paraxial wave equation, eikonal equation
 *   L4: Nonlinear Schrodinger equation for plasma
 *   L5: Split-step Fourier propagation, ray tracing
 *
 * Courses: MIT 22.611, Princeton PHY 535, ETH 402-0891
 */

#include "plasma_constants.h"
#include "plasma_params.h"

/* ============================================================
 *  L3: Beam Profile Structures
 * ============================================================ */

/**
 * GaussianBeam -- Paraxial Gaussian beam in vacuum
 *
 * w0      : waist radius at focus [m]
 * zR      : Rayleigh range [m]
 * lambda0 : vacuum wavelength [m]
 * power   : total beam power [W]
 */
typedef struct {
    double w0;
    double zR;
    double lambda0;
    double power;
} GaussianBeam;

/**
 * PlasmaWaveguide -- Plasma density channel for optical guiding
 *
 * depth   : density depression at axis [m^-3]
 * width   : channel half-width [m]
 * ne0     : background electron density [m^-3]
 */
typedef struct {
    double depth;
    double width;
    double ne0;
} PlasmaWaveguide;

/* ============================================================
 *  L2/L4: Refractive Effects
 * ============================================================ */

/**
 * ray_deflection_angle -- Deflection by transverse density gradient
 *
 *   dtheta/ds = (1/2nc) grad_perp(ne)
 *
 * Uses eikonal equation (ray optics) approximation. Valid when
 * density scale length >> wavelength and ne << nc.
 * [Kruer Sec 6.1]
 *
 * Complexity: O(1)
 */
double ray_deflection_angle(double ne_gradient, double nc, double path_length);

/**
 * density_scale_length -- Exponential density scale length
 *
 *   L_n = ne / |grad ne|
 *
 * Characteristic distance over which the density changes by a
 * factor of e. Key parameter for resonance absorption and
 * parametric instabilities.
 */
double density_scale_length(double ne, double ne_gradient);

/**
 * refraction_index_profile -- Refractive index along linear density ramp
 *
 *   N(z) = sqrt(1 - ne(z)/nc)
 *
 * Returns N at position z for a linear density profile from
 * ne_start to ne_end over length L.
 */
double refraction_index_profile(double z, double ne_start, double ne_end,
                                double nc, double L);

/* ============================================================
 *  L5: Ray Tracing and Propagation
 * ============================================================ */

/**
 * ray_trace_step -- Forward step in ray-tracing
 *
 * Advances ray position and direction by one step ds along the
 * ray path in a plasma with given density gradient.
 *
 *   dr/ds = k_hat
 *   dk_hat/ds = (1/2nc) grad_perp(ne)
 *
 * Complexity: O(1)
 */
void ray_trace_step(double *r, double *k_hat,
                    double ne_grad_x, double ne_grad_y,
                    double nc, double ds);

/**
 * snells_law_plasma -- Snell's law for stratified plasma
 *
 *   sin(theta)/sin(theta0) = sqrt(1 - ne(0)/nc) / sqrt(1 - ne(z)/nc)
 *
 * Refraction in a stratified plasma slab. Returns theta(z).
 */
double snells_law_plasma(double theta0, double ne0, double ne_z, double nc);

/**
 * turning_point_density -- Density at which a ray turns
 *
 *   ne_turn = nc cos^2(theta0)
 *
 * An obliquely incident ray with angle theta0 reaches its turning
 * point where the local refractive index vanishes for that angle.
 */
double turning_point_density(double theta0, double nc);

/* ============================================================
 *  L5: Paraxial Wave Equation (Numerical Propagation)
 * ============================================================ */

/**
 * ParaxialState -- State for split-step paraxial propagation
 *
 * E_real, E_imag : complex electric field envelope
 *                 on the transverse grid
 * Nx, Ny         : grid dimensions
 * dx, dy         : grid spacing [m]
 * z              : current propagation distance [m]
 */
typedef struct {
    double *E_real;
    double *E_imag;
    int     Nx, Ny;
    double  dx, dy;
    double  z;
    double  k0;
    double  lambda0;
} ParaxialState;

/**
 * paraxial_state_alloc -- Allocate paraxial propagation state
 *
 * Returns pointer to allocated state or NULL on failure.
 * Caller must free with paraxial_state_free.
 */
ParaxialState *paraxial_state_alloc(int Nx, int Ny, double dx, double dy,
                                    double lambda0);

/**
 * paraxial_state_free -- Release paraxial states
 */
void paraxial_state_free(ParaxialState *ps);

/**
 * nonlinear_phase_shift -- Relativistic nonlinear phase
 *
 *   dphi/dz = (omega_p^2 / (2 k0 c^2)) * (1 - 1/sqrt(1 + a^2))
 *
 * Phase accumulated per unit propagation distance due to
 * relativistic electron mass increase.
 * [Gibbon Sec 5.1]
 */
double nonlinear_phase_shift(double a_squared, double wp, double k0);

/**
 * paraxial_step_vacuum -- Free-space paraxial propagation step
 *
 * Advances the field by dz using the paraxial (Fresnel) approximation.
 * Implements split-step: advance phase in real space, then propagate
 * in Fourier space.
 *
 * Complexity: O(Nx Ny log(Nx Ny)) per step
 */
int paraxial_step_vacuum(ParaxialState *ps, double dz);

/**
 * paraxial_step_plasma -- Propagation step with plasma response
 *
 * Same as vacuum step but adds the nonlinear plasma phase term
 * from the relativistic mass increase and ponderomotive density
 * modification.
 */
int paraxial_step_plasma(ParaxialState *ps, double dz,
                         double ne, double nc);

/**
 * paraxial_run_to -- Propagate beam to target z
 *
 * Convenience wrapper that calls paraxial_step_vacuum or
 * paraxial_step_plasma repeatedly to reach distance z.
 *
 * Complexity: O(N_steps Nx Ny log(Nx Ny))
 */
int paraxial_run_to(ParaxialState *ps, double z_target,
                    double ne, double nc, double dz);

#endif /* LASER_PROPAGATION_H */
