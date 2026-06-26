#ifndef PARTICLE_MOTION_H
#define PARTICLE_MOTION_H

/*
 * particle_motion.h -- Single-particle motion in laser fields
 *
 * Models for single-electron dynamics in intense EM fields:
 * relativistic equations of motion, Boris pusher, ponderomotive
 * scattering, and vacuum laser acceleration.
 *
 * References:
 *   - Landau & Lifshitz (1975) "The Classical Theory of Fields"
 *   - Gibbon (2005), Ch. 2
 *   - Hartemann et al. (1995), Phys. Rev. E 51, 4833
 *   - Boris (1970), Proc. 4th Conf. Num. Sim. Plasmas
 *
 * Knowledge Layers:
 *   L1: Lorentz force, relativistic momentum, ponderomotive scattering
 *   L3: Canonical momentum, phase-space structures
 *   L4: Relativistic equations of motion, Larmor formula
 *   L5: Boris integrator, adaptive RK45
 *
 * Courses: MIT 8.07, Stanford PHYSICS 370, Princeton PHY 525
 */

#include "plasma_constants.h"
#include "plasma_params.h"

/* ============================================================
 *  L1: Particle State Definitions
 * ============================================================ */

/** Particle3D -- 6D phase-space coordinate for a single particle */
typedef struct {
    double x, y, z;     /* position [m]                            */
    double px, py, pz;  /* relativistic momentum [kg.m/s]          */
    double gamma;       /* Lorentz factor sqrt(1+p^2/(m^2 c^2))    */
    double mass_kg;     /* rest mass [kg]                          */
    double charge_C;    /* electric charge [C]                     */
} Particle3D;

/** EMField -- Electromagnetic field at a point */
typedef struct {
    double Ex, Ey, Ez;  /* electric field [V/m]                    */
    double Bx, By, Bz;  /* magnetic field [T]                      */
} EMField;

/** PlaneWaveField -- Linearly polarized plane wave parameters */
typedef struct {
    double lambda0;     /* wavelength [m]                          */
    double a0;          /* normalized vector potential             */
    double omega0;      /* angular frequency [rad/s]               */
    double k0;          /* wavenumber [rad/m]                      */
    double E0;          /* peak electric field [V/m]               */
    double B0;          /* peak magnetic field [T] = E0/c          */
    int    polarization; /* 0=linear, 1=circular                    */
    double phi_cep;     /* carrier-envelope phase [rad]            */
} PlaneWaveField;

/** ParticleTrack -- recorded trajectory (for diagnostics) */
typedef struct {
    double *t;          /* time samples [s]                        */
    double *x, *y, *z; /* position history [m]                     */
    double *px, *py, *pz; /* momentum history [kg.m/s]             */
    double *gamma;      /* Lorentz factor history                  */
    int     n_samples;
    int     capacity;
} ParticleTrack;

/* ============================================================
 *  L1/L4: Particle Dynamics
 * ============================================================ */

/**
 * particle_init -- Initialize particle at rest
 *
 * Sets position to (x0,y0,z0) and momentum to (px0,py0,pz0).
 * Computes gamma from the momentum.
 */
void particle_init(Particle3D *p, double x0, double y0, double z0,
                   double px0, double py0, double pz0,
                   double mass_kg, double charge_C);

/**
 * lorentz_force -- Compute Lorentz force on a particle
 *
 *   F = q (E + v x B)
 *
 * Returns force vector components through fx, fy, fz pointers.
 * Uses relativistic velocity v = p / (gamma m).
 *
 * Complexity: O(1)
 * Theorem: Lorentz force = rate of change of relativistic momentum
 */
void lorentz_force(const Particle3D *p, const EMField *field,
                   double *fx, double *fy, double *fz);

/**
 * update_gamma -- Recompute Lorentz factor from momentum
 *
 *   gamma = sqrt(1 + p^2/(m^2 c^2))
 */
void update_gamma(Particle3D *p);

/**
 * kinetic_energy -- Relativistic kinetic energy
 *
 *   K = (gamma - 1) m c^2
 */
double kinetic_energy(const Particle3D *p);

/**
 * rel_kinetic_energy_eV -- Relativistic kinetic energy in eV
 */
double rel_kinetic_energy_eV(const Particle3D *p);

/* ============================================================
 *  L4/L5: Numerical Integrators
 * ============================================================ */

/**
 * boris_push -- Boris relativistic particle pusher
 *
 * The industry-standard integrator for relativistic charged
 * particles in PIC codes.  Splits the motion:
 *   1. Half-step E-field acceleration
 *   2. Full-step B-field rotation
 *   3. Half-step E-field acceleration
 *
 * This scheme is symplectic for B=0 and volume-preserving for
 * the full Lorentz force, ensuring long-term stability.
 * [Boris 1970; Birdsall & Langdon (1991) Sec 4-4]
 *
 * Complexity: O(1)
 * Stability: Unconditionally stable for uniform fields
 */
void boris_push(Particle3D *p, const EMField *field, double dt);

/**
 * boris_push_N_steps -- Apply N Boris steps
 *
 * Advances particle through N time steps of size dt in
 * a time-independent EM field.
 *
 * Complexity: O(N)
 */
void boris_push_N_steps(Particle3D *p, const EMField *field,
                         double dt, int N);

/**
 * rk4_push -- 4th-order Runge-Kutta particle pusher
 *
 * Higher-order (but slower) integrator.  Useful for accuracy
 * benchmarks against the Boris scheme.
 *
 * Complexity: O(1) per step (4 field evaluations)
 */
void rk4_push(Particle3D *p, const EMField *field, double dt);

/* ============================================================
 *  L2: Ponderomotive Scattering
 * ============================================================ */

/**
 * ponderomotive_scattering_angle -- Angular deflection by laser pulse
 *
 * An electron initially at rest, struck by a finite-duration
 * laser pulse, acquires no net energy (in the plane-wave limit)
 * but is deflected by the ponderomotive potential gradient at
 * the pulse edges.
 *
 * Returns scattering angle theta [rad] for given pulse parameters
 * and impact parameter.
 * [Hartemann et al. 1995]
 */
double ponderomotive_scattering_angle(double a0, double w0,
                                       double impact_param,
                                       double tau_fwhm);

/**
 * ponderomotive_energy_after_passage -- Residual energy
 *
 * For a finite-width (Gaussian) laser beam, the ponderomotive
 * potential has a radial gradient that imparts a net energy
 * to electrons as they exit the focal region.
 *
 *   Delta_E ~ me c^2 * a0^2 * (lambda0/w0)^2
 *
 * Returns residual kinetic energy [eV].
 */
double ponderomotive_energy_after_passage(double a0, double w0,
                                           double lambda0);

/**
 * plane_wave_solution_8 -- Analytic orbit in plane wave
 *
 * An exact solution exists for a charged particle in a
 * plane EM wave of arbitrary intensity (Volkov solution).
 * Returns (px, py, pz) after time t for a particle initially
 * at rest.
 *
 * Complexity: O(1)
 * Theorem: Exact solution to relativistic Lorentz equation
 *          for plane wave field (Volkov 1935)
 */
void plane_wave_orbit_analytic(double t, double a0, double omega0,
                                int polarization, double mass_kg,
                                double *px, double *py, double *pz,
                                double *gamma_out);

/* ============================================================
 *  L3/L4: EM Wave Properties
 * ============================================================ */

/**
 * plane_wave_field_at -- EM field at (t, z) for plane wave
 *
 * Fills *field with E and B for a linearly or circularly
 * polarized plane wave propagating in +z direction.
 *
 *   Ex = E0 cos(kz - omega t + phi_cep)   (linear)
 *   By = B0 cos(kz - omega t + phi_cep)   (linear)
 */
void plane_wave_field_at(double t, double z,
                          const PlaneWaveField *pw,
                          EMField *field);

/**
 * laser_intensity_from_a0 -- Intensity from a0
 *
 *   I = (pi/2) eps0 me^2 c^5 a0^2 / (e^2 lambda^2)
 *
 * Inverse of the a0(intensity) formula.
 */
double laser_intensity_from_a0(double a0, double lambda_m);

/**
 * plane_wave_init -- Initialize PlaneWaveField from laser params
 */
void plane_wave_init(PlaneWaveField *pw, double lambda_m,
                     double a0, int polarization, double phi_cep);

/* ============================================================
 *  L5: Trajectory Diagnostics
 * ============================================================ */

/**
 * particle_track_alloc -- Allocate trajectory storage
 */
ParticleTrack *particle_track_alloc(int capacity);

/**
 * particle_track_free -- Free trajectory storage
 */
void particle_track_free(ParticleTrack *track);

/**
 * particle_track_record -- Record current particle state
 *
 * Returns 0 on success, -1 if capacity exceeded.
 */
int particle_track_record(ParticleTrack *track, double t,
                           const Particle3D *p);

/**
 * particle_track_max_energy -- Maximum kinetic energy in track [eV]
 */
double particle_track_max_energy(const ParticleTrack *track);

/**
 * particle_track_final_energy -- Final kinetic energy in track [eV]
 */
double particle_track_final_energy(const ParticleTrack *track);

#endif /* PARTICLE_MOTION_H */
