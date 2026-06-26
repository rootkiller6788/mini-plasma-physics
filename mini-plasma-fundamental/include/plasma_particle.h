/**
 * plasma_particle.h — Particle Methods for Plasma Simulation
 *
 * Single-particle motion, Boris particle pusher, PIC method basics,
 * and orbit integration in prescribed electromagnetic fields.
 *
 * References:
 *   - Birdsall & Langdon, "Plasma Physics via Computer Simulation" (1991)
 *   - Hockney & Eastwood, "Computer Simulation Using Particles" (1988)
 *   - Boris (1970), Proc. 4th Conf. Num. Sim. Plasmas
 *   - MIT 22.611 / Berkeley PHYS 242
 *
 * Knowledge Levels:
 *   L2 — Core Concepts: ExB drift, grad-B drift, curvature drift
 *   L4 — Fundamental: Lorentz force law, guiding center equations
 *   L5 — Computational: Boris pusher, PIC cycle, leapfrog
 *   L6 — Canonical: Particle orbits in tokamak fields
 */

#ifndef PLASMA_PARTICLE_H
#define PLASMA_PARTICLE_H

#include "plasma_constants.h"

/* ================================================================
 * L2: Single Particle State
 * ================================================================ */

/** Single particle in 3D electromagnetic fields */
typedef struct {
    double x, y, z;     /**< Position [m] */
    double vx, vy, vz;  /**< Velocity [m/s] */
    double mass;        /**< Mass [kg] */
    double charge;      /**< Charge [C] (signed, e.g. -e for electron) */
    int    alive;       /**< 1 if particle is active, 0 if lost */
} Particle3D;

/** Guiding center variables (reduced phase space) */
typedef struct {
    double X, Y, Z;     /**< Guiding center position [m] */
    double v_parallel;  /**< Parallel velocity [m/s] */
    double mu;          /**< Magnetic moment [J/T] (adiabatic invariant) */
    double mass;
    double charge;
} GuidingCenter;

/** Collection of particles for ensemble methods */
typedef struct {
    Particle3D *particles;  /**< Particle array */
    int         n_particles;/**< Number of particles */
    double      weight;     /**< Super-particle weight (real/particle ratio) */
} ParticleEnsemble;

/* ================================================================
 * L5: Boris Particle Pusher
 * ================================================================ */

/**
 * Boris particle push — one full time step.
 *
 * The Boris scheme splits the Lorentz force advance into:
 *   v^- = v^n + (q*dt/(2m)) * E^n                       [half E-accel]
 *   v^+ - v^- = (q*dt/(2m)) * (v^+ + v^-) x B^n         [rotation]
 *   v^{n+1} = v^+ + (q*dt/(2m)) * E^n                   [half E-accel]
 *   x^{n+1} = x^n + dt * v^{n+1}                        [position update]
 *
 * This is explicit, second-order accurate, symplectic for B=0,
 * and conserves energy exactly for E=0.
 *
 * Reference: Boris (1970)
 *
 * Complexity: O(1) per particle
 */
void boris_push(Particle3D *p, double Ex, double Ey, double Ez,
                double Bx, double By, double Bz, double dt);

/**
 * Boris push with spatially-varying fields (E,B as function pointers).
 *
 * @param E_field Function returning (Ex,Ey,Ez) at (x,y,z,t)
 * @param B_field Function returning (Bx,By,Bz) at (x,y,z,t)
 */
typedef void (*FieldFunc3D)(double x, double y, double z, double t,
                            double *Ex, double *Ey, double *Ez,
                            double *Bx, double *By, double *Bz);

void boris_push_field(Particle3D *p, double t, double dt,
                      FieldFunc3D E_field, FieldFunc3D B_field);

/**
 * Push an ensemble of particles through one time step.
 *
 * Complexity: O(n_particles)
 */
void boris_push_ensemble(ParticleEnsemble *ensemble,
                         FieldFunc3D E_field, FieldFunc3D B_field,
                         double t, double dt);

/* ================================================================
 * L2: Particle Drifts
 * ================================================================ */

/**
 * ExB drift velocity [m/s]
 *
 * v_ExB = (E x B) / B^2
 *
 * This drift is independent of charge, mass, and energy.
 * All particles drift together at this velocity.
 *
 * Complexity: O(1)
 */
void exb_drift_velocity(double Ex, double Ey, double Ez,
                        double Bx, double By, double Bz,
                        double *vdx, double *vdy, double *vdz);

/**
 * Grad-B drift velocity [m/s]
 *
 * v_gradB = (mu / (q*B^2)) * B x grad|B|
 *
 * where mu = (1/2)*m*v_perp^2 / B is the magnetic moment.
 *
 * Ions and electrons drift in opposite directions → current.
 *
 * Complexity: O(1)
 */
void gradb_drift_velocity(double mu, double q, double B_mag,
                          double gradBx, double gradBy, double gradBz,
                          double Bx, double By, double Bz,
                          double *vdx, double *vdy, double *vdz);

/**
 * Curvature drift velocity [m/s]
 *
 * v_curv = (m*v_par^2 / (q*B^2)) * (B x (B·grad)B) / B^2
 *
 * Arises from the centrifugal force on particles following
 * curved magnetic field lines.
 */
void curvature_drift_velocity(double m, double v_par, double q,
                              double Bx, double By, double Bz,
                              double curv_x, double curv_y, double curv_z,
                              double *vdx, double *vdy, double *vdz);

/**
 * Compute the total guiding center drift velocity (ExB + gradB + curvature).
 *
 * Complexity: O(1)
 */
void total_guiding_center_velocity(double Ex, double Ey, double Ez,
                                   double Bx, double By, double Bz,
                                   double gradB_mag, double v_par,
                                   double mu, double m, double q,
                                   double *vdx, double *vdy, double *vdz);

/**
 * Polarization drift velocity [m/s]
 *
 * v_pol = (m/(q*B^2)) * dE/dt
 *
 * This is important for time-varying electric fields.
 * Mass-dependent → ions dominate.
 */
void polarization_drift_velocity(double m, double q, double B_mag,
                                 double dEx_dt, double dEy_dt, double dEz_dt,
                                 double *vdx, double *vdy, double *vdz);

/* ================================================================
 * L6: Particle Orbits in Tokamak
 * ================================================================ */

/**
 * Banana orbit parameters for trapped particles in tokamak.
 *
 * In a tokamak, particles with v_par/v < sqrt(2*r/R) are trapped
 * in the magnetic mirror on the outboard side.
 *
 * Banana width:
 *   Delta_b = q * rho_L / sqrt(epsilon)
 *
 * where epsilon = r/R is the inverse aspect ratio.
 *
 * @param rho_L Larmor radius [m]
 * @param epsilon r/R inverse aspect ratio
 * @param q_safety Safety factor
 * @return Banana width [m]
 */
double banana_width(double rho_L, double epsilon, double q_safety);

/**
 * Trapped particle fraction in a tokamak.
 *
 * f_trapped = sqrt(2 * epsilon)
 *
 * For typical tokamaks (epsilon ~ 0.3): f_trapped ~ 77%
 *
 * Complexity: O(1)
 */
double trapped_particle_fraction(double epsilon);

/**
 * Bounce frequency of trapped particles [rad/s]
 *
 * omega_b = sqrt(epsilon) * v_th / (q * R0)
 */
double bounce_frequency(double epsilon, double v_th, double q_safety, double R0);

/**
 * Toroidal precession frequency [rad/s]
 *
 * omega_phi = q * v_th^2 / (omega_c * R0 * r) * (s + 1)
 *
 * where s = (r/q)*(dq/dr) is the magnetic shear.
 */
double toroidal_precession_frequency(double v_th, double omega_c,
                                     double R0, double r, double q,
                                     double shear);

/* ================================================================
 * L5: PIC Method Basics
 * ================================================================ */

/** 1D electrostatic PIC grid */
typedef struct {
    int     nx;         /**< Number of grid cells */
    double  dx;         /**< Cell size [m] */
    double *rho;        /**< Charge density on grid [C/m^3] */
    double *phi;        /**< Electrostatic potential on grid [V] */
    double *E;          /**< Electric field on grid [V/m] */
} PICGrid1D;

/**
 * First-order (NGP — Nearest Grid Point) charge deposition.
 *
 * Deposit particle charge to nearest grid point.
 *
 * Complexity: O(1) per particle
 */
void pic_deposit_ngp(PICGrid1D *grid, double x, double q);

/**
 * First-order (CIC — Cloud-in-Cell) charge deposition.
 *
 * Linear weighting (tent function) from particle to two nearest grid points.
 *
 * Complexity: O(1) per particle
 */
void pic_deposit_cic(PICGrid1D *grid, double x, double q);

/**
 * CIC field interpolation: get E-field at particle position.
 *
 * Complexity: O(1) per particle
 */
double pic_interpolate_cic(const PICGrid1D *grid, double x);

/**
 * Solve Poisson equation on PIC grid: d^2 phi/dx^2 = -rho/epsilon_0.
 *
 * Uses tridiagonal solver with Dirichlet (phi=0) boundary conditions.
 *
 * Complexity: O(nx)
 */
void pic_solve_poisson(PICGrid1D *grid);

/**
 * Compute E-field from potential: E = -d phi/dx.
 *
 * Uses central differences interior, one-sided at boundaries.
 *
 * Complexity: O(nx)
 */
void pic_compute_e_field(PICGrid1D *grid);

/**
 * One full PIC cycle:
 *   1. Deposit charge (particles -> grid)
 *   2. Solve Poisson (density -> potential)
 *   3. Compute E-field (potential -> field)
 *   4. Push particles (field -> particles)
 *   5. Move particles (velocity -> position)
 *
 * Complexity: O(n_particles + nx) per time step
 */
void pic_cycle_1d(PICGrid1D *grid, ParticleEnsemble *ensemble,
                  double dt, double t);

#endif /* PLASMA_PARTICLE_H */
