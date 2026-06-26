/**
 * plasma_kinetic.h — Kinetic Theory of Plasmas
 *
 * Structures and APIs for Vlasov-Fokker-Planck kinetic description
 * of plasmas. Covers the velocity-space distribution function approach
 * required when MHD approximations break down.
 *
 * References:
 *   - Vlasov (1938), JETP 8, 291
 *   - Landau (1946), JETP 16, 574 — Landau damping
 *   - Rosenbluth, MacDonald & Judd (1957), Phys. Rev. 107, 1
 *   - Balescu (1960), Phys. Fluids 3, 52
 *
 * Knowledge Levels:
 *   L2 — Core Concepts: velocity distribution, Landau damping
 *   L3 — Math Structures: Vlasov equation, Fokker-Planck operator
 *   L4 — Fundamental Laws: Boltzmann equation with Coulomb collisions
 */

#ifndef PLASMA_KINETIC_H
#define PLASMA_KINETIC_H

#include "plasma_constants.h"
#include <stddef.h>

/* ================================================================
 * L2/L3: Distribution Function Grid
 * ================================================================ */

/** 1D velocity-space grid for distribution function */
typedef struct {
    int     nv;         /**< Number of velocity grid points */
    double  v_min;      /**< Minimum velocity [m/s] */
    double  v_max;      /**< Maximum velocity [m/s] */
    double  dv;         /**< Grid spacing dv = (vmax - vmin) / (nv - 1) */
    double *v;          /**< Velocity grid array [nv] */
    double *f;          /**< Distribution function f(v) [nv] */
} Distribution1D;

/** 2D velocity-space grid (v_parallel, v_perp) for gyrokinetics */
typedef struct {
    int     nvpar;      /**< Number of parallel velocity points */
    int     nvperp;     /**< Number of perpendicular velocity points */
    double  vpar_min, vpar_max, dvpar;
    double  vperp_min, vperp_max, dvperp;
    double **f;         /**< f(v_par, v_perp) [nvpar][nvperp] */
} Distribution2D;

/** 1D spatial + 1D velocity phase-space grid (1D1V) */
typedef struct {
    int     nx;         /**< Number of spatial points */
    int     nv;         /**< Number of velocity points */
    double  x_min, x_max, dx;
    double  v_min, v_max, dv;
    double **f;         /**< f(x,v) [nx][nv] */
    double *E_field;    /**< Electric field E(x) [nx] */
    double  dt;         /**< Time step */
    double  t;          /**< Current time */
} PhaseSpace1D1V;

/* ================================================================
 * L2: Reduced Distribution Functions (Moments)
 * ================================================================ */

/**
 * Compute number density n(x,t) = integral f(x,v,t) dv
 *
 * Complexity: O(nv)
 */
double density_from_f1d1v(const PhaseSpace1D1V *ps, int ix);

/**
 * Compute mean velocity u(x,t) = (1/n) * integral v * f(x,v,t) dv
 *
 * Complexity: O(nv)
 */
double mean_velocity_from_f1d1v(const PhaseSpace1D1V *ps, int ix);

/**
 * Compute kinetic temperature (2nd moment)
 * T_kin = (m/kB) * (1/n) * integral (v-u)^2 f(x,v,t) dv
 *
 * Complexity: O(nv)
 */
double temperature_from_f1d1v(const PhaseSpace1D1V *ps, int ix, double mass);

/**
 * Compute heat flux q (3rd centered moment)
 * q = (m/2) * integral (v-u)^3 f(x,v,t) dv
 *
 * Complexity: O(nv)
 */
double heat_flux_from_f1d1v(const PhaseSpace1D1V *ps, int ix, double mass);

/**
 * Compute pressure tensor component P_xx
 * P_xx = m * integral (v-u)^2 f(x,v,t) dv
 *
 * Complexity: O(nv)
 */
double pressure_from_f1d1v(const PhaseSpace1D1V *ps, int ix, double mass);

/* ================================================================
 * L3: Vlasov Equation
 * ================================================================ */

/**
 * Initialize a Maxwellian distribution on a 1D1V phase space grid.
 *
 * f(x,v,t=0) = n0 * (m/(2*pi*kB*T))^(1/2) * exp(-m*v^2/(2*kB*T))
 *
 * Complexity: O(nx * nv)
 */
void init_maxwellian_1d1v(PhaseSpace1D1V *ps, double n0, double T0,
                          double mass);

/**
 * Compute the spatial advection term for the Vlasov equation:
 *   df/dt = -v * df/dx - (q/m) E(x) * df/dv
 *
 * Uses 2nd-order central differences for df/dx.
 *
 * Complexity: O(nx * nv)
 */
void vlasov_advection_x(const PhaseSpace1D1V *ps, double **dfdt);

/**
 * Compute the velocity advection term for Vlasov equation:
 *   df/dt = ... - (q/m) * E(x) * df/dv
 *
 * Uses 2nd-order central differences for df/dv.
 *
 * Complexity: O(nx * nv)
 */
void vlasov_advection_v(const PhaseSpace1D1V *ps, double charge, double mass,
                        double **dfdt);

/**
 * One full Vlasov time step using operator splitting (Strang splitting).
 *
 * Step 1: Half-step velocity advection (dt/2)
 * Step 2: Full-step spatial advection (dt)
 * Step 3: Half-step velocity advection (dt/2)
 *
 * Reference: Cheng & Knorr (1976), J. Comp. Phys. 22, 330.
 *
 * Complexity: O(nx * nv)
 */
void vlasov_step_strang(PhaseSpace1D1V *ps, double charge, double mass);

/**
 * Solve Poisson's equation in 1D to update E-field from density:
 *   d^2 phi / dx^2 = - (e / epsilon_0) * (n_i - n_e)
 *   E = -d phi / dx
 *
 * Uses tridiagonal solver for the discrete system.
 *
 * Complexity: O(nx)
 */
void poisson_solve_1d(PhaseSpace1D1V *ps, double n0_ion);

/* ================================================================
 * L4: Landau Damping Analysis
 * ================================================================ */

/**
 * Compute the plasma dispersion function Z(zeta)
 *
 * Z(zeta) = (1/sqrt(pi)) * integral_{-inf}^{inf}
 *           exp(-t^2) / (t - zeta) dt   (Im(zeta) > 0)
 *
 * For Im(zeta) <= 0, the analytic continuation is used.
 *
 * This is the Fried-Conte function, central to linear kinetic
 * theory of plasma waves.
 *
 * Reference: Fried & Conte (1961), "The Plasma Dispersion Function"
 *
 * Complexity: O(1) using rational approximation
 */
double plasma_dispersion_function(double zeta_real, double zeta_imag,
                                  double *Z_real, double *Z_imag);

/**
 * Compute Landau damping rate gamma_L for Langmuir waves.
 *
 * gamma_L = -omega_pe * sqrt(pi/8) * (1/(k*lambda_D)^3)
 *           * exp(-1/(2*(k*lambda_D)^2) - 3/2)
 *
 * For a Maxwellian plasma, this is the collisionless damping rate
 * of electron plasma waves due to resonant wave-particle interaction.
 *
 * Theorem (Landau, 1946): Even in a collisionless plasma,
 * electrostatic waves are damped by phase mixing of particles
 * traveling near the wave phase velocity.
 */
double landau_damping_rate(double k, double omega_pe, double lambda_D);

/**
 * Compute the phase velocity v_phi = omega/k for Langmuir waves
 * using the Bohm-Gross dispersion relation:
 *   omega^2 = omega_pe^2 + 3 * (kB*Te/m_e) * k^2
 */
double langmuir_phase_velocity(double k, double omega_pe, double Te);

/**
 * Check if a wave is strongly Landau-damped.
 * Condition: |gamma_L / omega_r| > 0.1
 *
 * Returns 1 if strongly damped, 0 otherwise.
 */
int is_strongly_landau_damped(double k, double omega_pe, double lambda_D);

/* ================================================================
 * L3: Fokker-Planck Collision Operator
 * ================================================================ */

/**
 * Compute the Rosenbluth potentials for Coulomb collisions.
 *
 * Given a distribution function f_b(v) of background particles "b",
 * compute the Rosenbluth-Trubnikov potentials:
 *
 * h_b(v) = integral f_b(v') / |v - v'| dv'   (Rosenbluth potential)
 * g_b(v) = integral f_b(v') * |v - v'| dv'   (Trubnikov potential)
 *
 * These are used to construct the Fokker-Planck collision operator.
 *
 * Reference: Rosenbluth et al. (1957), Phys. Rev. 107, 1
 *
 * @param v Velocity magnitude at which to evaluate
 * @param dist 1D isotropic distribution function
 * @param h Output: Rosenbluth potential
 * @param g Output: Trubnikov potential
 * @param dh Output: dh/dv
 * @param dg Output: dg/dv
 * @param d2g Output: d^2g/dv^2
 */
void rosenbluth_potentials(double v, const Distribution1D *dist,
                           double *h, double *g,
                           double *dh, double *dg, double *d2g);

/**
 * Fokker-Planck collision operator for electron-electron collisions.
 *
 * df/dt|_c = Gamma * (1/v^2) * d/dv [ ... ]
 *
 * where Gamma = e^4 ln Lambda / (4*pi*epsilon_0^2 * m^2)
 *
 * Complexity: O(nv) per evaluation
 */
void fokker_planck_ee(double *dfdt_coll, const Distribution1D *f,
                      int nv, double dv, double ln_Lambda);

/**
 * Pitch-angle scattering operator for electron-ion collisions.
 *
 * df/dt = nu_ei(v) * (1/2) * d/dmu [(1-mu^2) df/dmu]
 *
 * where mu = v_parallel / v = cos(theta) is the pitch angle.
 *
 * Complexity: O(nvperp * nvpar)
 */
void pitch_angle_scattering(Distribution2D *f, double nu_ei);

#endif /* PLASMA_KINETIC_H */
