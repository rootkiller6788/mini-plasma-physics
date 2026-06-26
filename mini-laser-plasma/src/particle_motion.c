/*
 * particle_motion.c -- Single-particle motion in laser fields
 *
 * Implements relativistic particle dynamics in intense EM fields,
 * including the Boris pusher, ponderomotive scattering, and the
 * analytic Volkov solution.
 *
 * References:
 *   - Landau & Lifshitz (1975) "The Classical Theory of Fields"
 *   - Gibbon (2005), Ch. 2
 *   - Hartemann et al. (1995), Phys. Rev. E 51, 4833
 *   - Boris (1970), Proc. 4th Conf. Num. Sim. Plasmas
 *   - Birdsall & Langdon (1991) "Plasma Physics via Computer Simulation"
 *
 * Knowledge Layers: L1, L3, L4, L5
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "plasma_constants.h"
#include "plasma_params.h"
#include "particle_motion.h"

static const double c   = PLASMA_C;
static const double e   = PLASMA_E;
static const double me  = PLASMA_ME;
static const double mp  = PLASMA_MP;
static const double eps0 = PLASMA_EPS0;

/* ============================================================
 *  L1: particle_init -- Initialize particle state
 *
 * Sets the 6D phase-space position and computes the Lorentz
 * factor gamma from the initial momentum.
 *
 * Complexity: O(1)
 * ============================================================ */
void particle_init(Particle3D *p, double x0, double y0, double z0,
                   double px0, double py0, double pz0,
                   double mass_kg, double charge_C)
{
    if (!p) return;
    p->x = x0; p->y = y0; p->z = z0;
    p->px = px0; p->py = py0; p->pz = pz0;
    p->mass_kg = mass_kg;
    p->charge_C = charge_C;
    update_gamma(p);
}

/* ============================================================
 *  L1: lorentz_force -- Relativistic Lorentz force
 *
 *  F = q (E + v x B)
 *
 * where v = p / (gamma m) is the relativistic 3-velocity.
 *
 * This is the fundamental force law for charged particles
 * in electromagnetic fields.  The relativistic formulation
 * uses the proper time to maintain covariance:
 *
 *   dp^mu/dtau = q F^{mu nu} u_nu
 *
 * Source: Landau & Lifshitz, Sec 17
 *
 * Complexity: O(1)
 * ============================================================ */
void lorentz_force(const Particle3D *p, const EMField *field,
                   double *fx, double *fy, double *fz)
{
    if (!p || !field || !fx || !fy || !fz) return;

    double gamma = p->gamma;
    double mass = p->mass_kg;
    double q    = p->charge_C;

    /* Velocity: v = p / (gamma m) */
    double vx = p->px / (gamma * mass);
    double vy = p->py / (gamma * mass);
    double vz = p->pz / (gamma * mass);

    /* F = q (E + v x B) */
    *fx = q * (field->Ex + vy * field->Bz - vz * field->By);
    *fy = q * (field->Ey + vz * field->Bx - vx * field->Bz);
    *fz = q * (field->Ez + vx * field->By - vy * field->Bx);
}

/* ============================================================
 *  L3: update_gamma -- Recompute Lorentz factor
 *
 *  gamma = sqrt(1 + p^2/(m^2 c^2))
 *
 * This is the relativistic energy-momentum relation.
 * The total energy is E = gamma m c^2.
 *
 * Complexity: O(1)
 * ============================================================ */
void update_gamma(Particle3D *p)
{
    if (!p) return;
    double mass = p->mass_kg;
    double p_sq = p->px * p->px + p->py * p->py + p->pz * p->pz;
    p->gamma = sqrt(1.0 + p_sq / (mass * mass * c * c));
}

/* ============================================================
 *  L1: kinetic_energy -- Relativistic kinetic energy [J]
 *
 *  K = (gamma - 1) m c^2
 *
 * For v << c: K ~ (1/2) m v^2 (classical limit)
 * For v -> c: K >> m c^2 (ultra-relativistic)
 *
 * Complexity: O(1)
 * ============================================================ */
double kinetic_energy(const Particle3D *p)
{
    if (!p) return 0.0;
    return (p->gamma - 1.0) * p->mass_kg * c * c;
}

/* ============================================================
 *  L1: rel_kinetic_energy_eV -- Kinetic energy in eV
 *
 * Complexity: O(1)
 * ============================================================ */
double rel_kinetic_energy_eV(const Particle3D *p)
{
    if (!p) return 0.0;
    return kinetic_energy(p) / e;
}

/* ============================================================
 *  L5: boris_push -- Boris relativistic particle pusher
 *
 * Algorithm:
 *   1. Half-step E acceleration:
 *        p^- = p^n + q E dt/2
 *   2. Full-step B rotation:
 *        p^+ = p^- + (p^- + p^- x t) x s  (Boris rotation)
 *        where t = q B dt / (2 gamma m)
 *              s = 2 t / (1 + |t|^2)
 *   3. Half-step E acceleration:
 *        p^{n+1} = p^+ + q E dt/2
 *   4. Position update:
 *        x^{n+1} = x^n + p^{n+1} dt / (gamma^{n+1} m)
 *
 * The Boris scheme is volume-preserving in phase space for
 * the Lorentz force, which ensures long-term numerical
 * stability (no artificial heating or cooling).
 *
 * Source: Boris (1970); Birdsall & Langdon Sec 4-4
 *
 * Complexity: O(1)
 * Stability: Unconditionally stable for uniform B (omega_c*dt
 *            can be arbitrarily large, though accuracy degrades)
 * ============================================================ */
void boris_push(Particle3D *p, const EMField *field, double dt)
{
    if (!p || !field || dt <= 0.0) return;

    double m = p->mass_kg;
    double q = p->charge_C;

    /* Step 1: Half-step E-field acceleration */
    double px_minus = p->px + 0.5 * q * field->Ex * dt;
    double py_minus = p->py + 0.5 * q * field->Ey * dt;
    double pz_minus = p->pz + 0.5 * q * field->Ez * dt;

    /* Compute gamma from p^- */
    double p2_minus = px_minus*px_minus + py_minus*py_minus
                      + pz_minus*pz_minus;
    double gamma_minus = sqrt(1.0 + p2_minus / (m * m * c * c));

    /* Step 2: B-field rotation (Boris rotation) */
    double t_x = 0.5 * q * field->Bx * dt / (gamma_minus * m);
    double t_y = 0.5 * q * field->By * dt / (gamma_minus * m);
    double t_z = 0.5 * q * field->Bz * dt / (gamma_minus * m);

    double t2 = t_x*t_x + t_y*t_y + t_z*t_z;
    double s_x = 2.0 * t_x / (1.0 + t2);
    double s_y = 2.0 * t_y / (1.0 + t2);
    double s_z = 2.0 * t_z / (1.0 + t2);

    /* p' = p^- + p^- x t */
    double px_prime = px_minus + py_minus * t_z - pz_minus * t_y;
    double py_prime = py_minus + pz_minus * t_x - px_minus * t_z;
    double pz_prime = pz_minus + px_minus * t_y - py_minus * t_x;

    /* p^+ = p^- + p' x s */
    double px_plus = px_minus + py_prime * s_z - pz_prime * s_y;
    double py_plus = py_minus + pz_prime * s_x - px_prime * s_z;
    double pz_plus = pz_minus + px_prime * s_y - py_prime * s_x;

    /* Step 3: Half-step E-field acceleration */
    p->px = px_plus + 0.5 * q * field->Ex * dt;
    p->py = py_plus + 0.5 * q * field->Ey * dt;
    p->pz = pz_plus + 0.5 * q * field->Ez * dt;

    /* Update gamma */
    update_gamma(p);

    /* Step 4: Position update */
    p->x += p->px * dt / (p->gamma * m);
    p->y += p->py * dt / (p->gamma * m);
    p->z += p->pz * dt / (p->gamma * m);
}

/* ============================================================
 *  L5: boris_push_N_steps -- Multiple Boris steps
 *
 * Advances particle through N steps in a static EM field.
 * Useful for computing trajectories over many laser periods.
 *
 * Complexity: O(N)
 * ============================================================ */
void boris_push_N_steps(Particle3D *p, const EMField *field,
                         double dt, int N)
{
    for (int i = 0; i < N; i++) {
        boris_push(p, field, dt);
    }
}

/* ============================================================
 *  L5: rk4_push -- 4th-order Runge-Kutta integrator
 *
 * Does 4 field evaluations per step, giving O(dt^4) local
 * error vs O(dt^2) for Boris.  Used for accuracy benchmarks.
 *
 * dp/dt = F(p, x, field)
 * dx/dt = p / (gamma m)
 *
 * Complexity: O(1) per step (4 evaluations)
 * ============================================================ */
void rk4_push(Particle3D *p, const EMField *field, double dt)
{
    if (!p || !field || dt <= 0.0) return;

    double m = p->mass_kg;
    double q = p->charge_C;

    /* Helper: compute acceleration F/m */
    /* F/m = (q/m) (E + v x B) */
    #define ACCEL(px,py,pz,gm, ax,ay,az) do { \
        double vx_ = (px)/((gm)*(m)); \
        double vy_ = (py)/((gm)*(m)); \
        double vz_ = (pz)/((gm)*(m)); \
        double qom_ = (q)/(m); \
        ax = qom_ * (field->Ex + vy_*field->Bz - vz_*field->By); \
        ay = qom_ * (field->Ey + vz_*field->Bx - vx_*field->Bz); \
        az = qom_ * (field->Ez + vx_*field->By - vy_*field->Bx); \
    } while(0)

    double px0 = p->px, py0 = p->py, pz0 = p->pz;
    double x0 = p->x, y0 = p->y, z0 = p->z;
    double g0 = p->gamma;

    /* k1 */
    double k1_px, k1_py, k1_pz;
    double k1_vx, k1_vy, k1_vz;
    ACCEL(px0, py0, pz0, g0, k1_px, k1_py, k1_pz);
    k1_vx = px0 / (g0 * m); k1_vy = py0 / (g0 * m); k1_vz = pz0 / (g0 * m);

    /* k2 (midpoint) */
    double pm2_x = px0 + 0.5*dt*k1_px;
    double pm2_y = py0 + 0.5*dt*k1_py;
    double pm2_z = pz0 + 0.5*dt*k1_pz;
    double gm2 = sqrt(1.0 + (pm2_x*pm2_x+pm2_y*pm2_y+pm2_z*pm2_z)/(m*m*c*c));
    double k2_px, k2_py, k2_pz;
    double k2_vx, k2_vy, k2_vz;
    ACCEL(pm2_x, pm2_y, pm2_z, gm2, k2_px, k2_py, k2_pz);
    k2_vx = pm2_x / (gm2 * m); k2_vy = pm2_y / (gm2 * m); k2_vz = pm2_z / (gm2 * m);

    /* k3 */
    double pm3_x = px0 + 0.5*dt*k2_px;
    double pm3_y = py0 + 0.5*dt*k2_py;
    double pm3_z = pz0 + 0.5*dt*k2_pz;
    double gm3 = sqrt(1.0 + (pm3_x*pm3_x+pm3_y*pm3_y+pm3_z*pm3_z)/(m*m*c*c));
    double k3_px, k3_py, k3_pz;
    double k3_vx, k3_vy, k3_vz;
    ACCEL(pm3_x, pm3_y, pm3_z, gm3, k3_px, k3_py, k3_pz);
    k3_vx = pm3_x / (gm3 * m); k3_vy = pm3_y / (gm3 * m); k3_vz = pm3_z / (gm3 * m);

    /* k4 */
    double pm4_x = px0 + dt*k3_px;
    double pm4_y = py0 + dt*k3_py;
    double pm4_z = pz0 + dt*k3_pz;
    double gm4 = sqrt(1.0 + (pm4_x*pm4_x+pm4_y*pm4_y+pm4_z*pm4_z)/(m*m*c*c));
    double k4_px, k4_py, k4_pz;
    double k4_vx, k4_vy, k4_vz;
    ACCEL(pm4_x, pm4_y, pm4_z, gm4, k4_px, k4_py, k4_pz);
    k4_vx = pm4_x / (gm4 * m); k4_vy = pm4_y / (gm4 * m); k4_vz = pm4_z / (gm4 * m);

    /* Update */
    p->px = px0 + dt/6.0 * (k1_px + 2.0*k2_px + 2.0*k3_px + k4_px);
    p->py = py0 + dt/6.0 * (k1_py + 2.0*k2_py + 2.0*k3_py + k4_py);
    p->pz = pz0 + dt/6.0 * (k1_pz + 2.0*k2_pz + 2.0*k3_pz + k4_pz);

    p->x = x0 + dt/6.0 * (k1_vx + 2.0*k2_vx + 2.0*k3_vx + k4_vx);
    p->y = y0 + dt/6.0 * (k1_vy + 2.0*k2_vy + 2.0*k3_vy + k4_vy);
    p->z = z0 + dt/6.0 * (k1_vz + 2.0*k2_vz + 2.0*k3_vz + k4_vz);

    update_gamma(p);

    #undef ACCEL
}

/* ============================================================
 *  L2: ponderomotive_scattering_angle
 *
 * An electron initially at rest that is struck by a finite-
 * duration laser pulse is deflected by the radial ponderomotive
 * potential gradient (phi_p ~ a0^2 exp(-2r^2/w0^2)).
 *
 * The scattering angle is:
 *   theta ~ a0 * (lambda / w0) * exp(-2 b^2 / w0^2)
 *
 * where b is the impact parameter.
 *
 * Source: Hartemann et al. (1995), Eq. 25
 *
 * Complexity: O(1)
 * ============================================================ */
double ponderomotive_scattering_angle(double a0, double w0,
                                       double impact_param,
                                       double tau_fwhm)
{
    if (w0 <= 0.0) return 0.0;

    double b_norm = impact_param / w0;
    double radial_profile = exp(-2.0 * b_norm * b_norm);

    /* Theta ~ a0 * (c * tau_fwhm / w0) * radial_gradient */
    /* Simplified: theta ~ a0 * radial_profile (order of magnitude) */
    double theta = a0 * radial_profile;

    /* Clamp to physical range */
    if (theta > M_PI) theta = M_PI;
    return theta;
}

/* ============================================================
 *  L2: ponderomotive_energy_after_passage
 *
 * For a finite-width Gaussian beam, electrons acquire net
 * energy from the radial ponderomotive potential gradient
 * as they leave the focal region:
 *
 *   Delta_E ~ me c^2 * a0^2 * (lambda / w0)^2
 *
 * For tightly focused beams (w0 ~ lambda), this can be
 * significant.
 *
 * Source: Hartemann et al. (1995)
 *
 * Complexity: O(1)
 * ============================================================ */
double ponderomotive_energy_after_passage(double a0, double w0,
                                           double lambda0)
{
    if (w0 <= 0.0 || lambda0 <= 0.0) return 0.0;

    double delta_E_J = me * c * c * a0 * a0
                       * (lambda0 * lambda0) / (w0 * w0);
    return delta_E_J / e;  /* return in eV */
}

/* ============================================================
 *  L3/L4: plane_wave_orbit_analytic -- Volkov solution
 *
 * An exact solution to the relativistic Lorentz equation
 * for a charged particle in a plane EM wave of arbitrary
 * intensity.
 *
 * For linear polarization along x:
 *   px(t) = -(e E0 / omega) [cos(phi) - cos(phi0)]
 *   py(t) = 0
 *   pz(t) = (e^2 E0^2 / (2 m c omega^2))
 *            [sin(2phi) - sin(2phi0)] / 4
 *
 * where phi = omega t - k z(t) is the laser phase.
 *
 * Source: Volkov (1935); Landau & Lifshitz, Sec 47
 * Theorem: The Volkov states are exact solutions to the
 *          Dirac equation in a plane EM wave.  The classical
 *          limit gives the above orbit.
 *
 * Complexity: O(1)
 * ============================================================ */
void plane_wave_orbit_analytic(double t, double a0, double omega0,
                                int polarization, double mass_kg,
                                double *px, double *py, double *pz,
                                double *gamma_out)
{
    if (!px || !py || !pz || !gamma_out) return;

    double phi = omega0 * t;
    double cos_phi = cos(phi);
    double sin_phi = sin(phi);

    /* Normalized momentum components for linear polarization along x */
    double px_norm, py_norm, pz_norm;

    if (polarization == 0) {
        /* Linear polarization along x */
        px_norm = a0 * (1.0 - cos_phi);
        py_norm = 0.0;
        pz_norm = (a0 * a0 / 4.0) * (2.0 * phi - sin(2.0 * phi));
    } else {
        /* Circular polarization */
        px_norm = a0 * (1.0 - cos_phi);
        py_norm = a0 * sin_phi;
        pz_norm = a0 * a0 * (phi - sin_phi * cos_phi);
    }

    double gamma = 1.0 + (px_norm*px_norm + py_norm*py_norm
                           + pz_norm*pz_norm) / 2.0;

    double p_scale = mass_kg * c;
    *px = px_norm * p_scale;
    *py = py_norm * p_scale;
    *pz = pz_norm * p_scale;
    *gamma_out = gamma;
}

/* ============================================================
 *  L3: plane_wave_field_at -- EM field at (t, z) for plane wave
 *
 * For a plane wave propagating in +z:
 *
 * Linear (along x):
 *   Ex = E0 cos(kz - omega t + phi_cep)
 *   By = B0 cos(kz - omega t + phi_cep) = Ex / c
 *   Ey = Ez = Bx = Bz = 0
 *
 * Circular:
 *   Ex = E0 cos(kz - omega t + phi_cep)
 *   Ey = E0 sin(kz - omega t + phi_cep)
 *   Bx = -Ey / c, By = Ex / c
 *   Ez = Bz = 0
 *
 * This satisfies Maxwell's equations in vacuum:
 *   div E = 0, div B = 0
 *   curl E = -dB/dt, curl B = (1/c^2) dE/dt
 *
 * Complexity: O(1)
 * ============================================================ */
void plane_wave_field_at(double t, double z,
                          const PlaneWaveField *pw,
                          EMField *field)
{
    if (!pw || !field) return;

    double phase = pw->k0 * z - pw->omega0 * t + pw->phi_cep;
    double cos_ph = cos(phase);
    double sin_ph = sin(phase);

    if (pw->polarization == 0) {
        /* Linear polarization along x */
        field->Ex = pw->E0 * cos_ph;
        field->Ey = 0.0;
        field->Ez = 0.0;
        field->Bx = 0.0;
        field->By = pw->B0 * cos_ph;
        field->Bz = 0.0;
    } else {
        /* Circular polarization */
        field->Ex = pw->E0 * cos_ph;
        field->Ey = pw->E0 * sin_ph;
        field->Ez = 0.0;
        field->Bx = -pw->B0 * sin_ph;
        field->By =  pw->B0 * cos_ph;
        field->Bz = 0.0;
    }
}

/* ============================================================
 *  L1: laser_intensity_from_a0 -- Convert a0 to intensity
 *
 *  I = (pi/2) eps0 me^2 c^5 a0^2 / (e^2 lambda^2)
 *
 * Inverse of the normalized_vector_potential formula.
 *
 * For a0=1, lambda=1 um: I ~ 1.37e18 W/cm^2
 *
 * Complexity: O(1)
 * ============================================================ */
double laser_intensity_from_a0(double a0, double lambda_m)
{
    if (lambda_m <= 0.0) return 0.0;
    double prefactor = M_PI / 2.0 * eps0 * me * me
                       * c * c * c * c * c;
    return prefactor * a0 * a0 / (e * e * lambda_m * lambda_m);
}

/* ============================================================
 *  L1: plane_wave_init -- Initialize plane wave from laser params
 *
 * Complexity: O(1)
 * ============================================================ */
void plane_wave_init(PlaneWaveField *pw, double lambda_m,
                     double a0, int polarization, double phi_cep)
{
    if (!pw) return;

    pw->lambda0 = lambda_m;
    pw->a0 = a0;
    pw->omega0 = 2.0 * M_PI * c / lambda_m;
    pw->k0 = 2.0 * M_PI / lambda_m;
    pw->E0 = me * c * pw->omega0 * a0 / e;
    pw->B0 = pw->E0 / c;
    pw->polarization = polarization;
    pw->phi_cep = phi_cep;
}

/* ============================================================
 *  L5: particle_track_alloc -- Trajectory storage
 *
 * Complexity: O(1) allocation
 * ============================================================ */
ParticleTrack *particle_track_alloc(int capacity)
{
    if (capacity < 1) return NULL;

    ParticleTrack *track = (ParticleTrack *)malloc(sizeof(ParticleTrack));
    if (!track) return NULL;

    size_t sz = (size_t)capacity * sizeof(double);
    track->t     = (double *)malloc(sz);
    track->x     = (double *)malloc(sz);
    track->y     = (double *)malloc(sz);
    track->z     = (double *)malloc(sz);
    track->px    = (double *)malloc(sz);
    track->py    = (double *)malloc(sz);
    track->pz    = (double *)malloc(sz);
    track->gamma = (double *)malloc(sz);

    if (!track->t || !track->x || !track->y || !track->z
        || !track->px || !track->py || !track->pz || !track->gamma) {
        free(track->t); free(track->x); free(track->y); free(track->z);
        free(track->px); free(track->py); free(track->pz);
        free(track->gamma);
        free(track);
        return NULL;
    }

    track->n_samples = 0;
    track->capacity = capacity;
    return track;
}

/* ============================================================
 *  L5: particle_track_free
 * ============================================================ */
void particle_track_free(ParticleTrack *track)
{
    if (!track) return;
    free(track->t);
    free(track->x); free(track->y); free(track->z);
    free(track->px); free(track->py); free(track->pz);
    free(track->gamma);
    free(track);
}

/* ============================================================
 *  L5: particle_track_record -- Record current state
 *
 * Returns 0 on success, -1 if capacity exceeded.
 *
 * Complexity: O(1)
 * ============================================================ */
int particle_track_record(ParticleTrack *track, double t,
                           const Particle3D *p)
{
    if (!track || !p) return -1;
    if (track->n_samples >= track->capacity) return -1;

    int i = track->n_samples;
    track->t[i]     = t;
    track->x[i]     = p->x;
    track->y[i]     = p->y;
    track->z[i]     = p->z;
    track->px[i]    = p->px;
    track->py[i]    = p->py;
    track->pz[i]    = p->pz;
    track->gamma[i] = p->gamma;
    track->n_samples++;

    return 0;
}

/* ============================================================
 *  L5: particle_track_max_energy -- Maximum kinetic energy
 *
 * Complexity: O(n_samples)
 * ============================================================ */
double particle_track_max_energy(const ParticleTrack *track)
{
    if (!track || track->n_samples == 0) return 0.0;

    double max_gamma = 0.0;
    for (int i = 0; i < track->n_samples; i++) {
        if (track->gamma[i] > max_gamma)
            max_gamma = track->gamma[i];
    }

    return (max_gamma - 1.0) * me * c * c / e;  /* eV */
}

/* ============================================================
 *  L5: particle_track_final_energy -- Final kinetic energy
 *
 * Complexity: O(1)
 * ============================================================ */
double particle_track_final_energy(const ParticleTrack *track)
{
    if (!track || track->n_samples == 0) return 0.0;

    int last = track->n_samples - 1;
    double gamma_last = track->gamma[last];
    return (gamma_last - 1.0) * me * c * c / e;  /* eV */
}
