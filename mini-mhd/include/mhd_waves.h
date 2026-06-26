/*
 * mhd_waves.h -- MHD Wave Dispersion and Propagation
 *
 * Reference: Freidberg "Ideal MHD" (2014) Ch.6
 *            Goedbloed & Poedts "Principles of MHD" (2004) Ch.5
 *            Boyd & Sanderson "Physics of Plasmas" (2003)
 *
 * Knowledge: L2 -- MHD wave modes (Alfven, fast, slow magnetosonic)
 *             L3 -- Eigenvalue analysis of MHD flux Jacobian
 *             L6 -- Alfven wave as canonical MHD system
 *
 * In MHD there are 3 distinct wave modes:
 *   1. Alfven wave -- incompressible, transverse, propagates along B
 *   2. Fast magnetosonic -- compressible, propagates fastest perpendicular to B
 *   3. Slow magnetosonic -- compressible, propagates slowly along B
 *
 * The phase speeds satisfy the MHD dispersion relation:
 *   (omega^2/k^2 - v_A^2 cos^2(theta)) *
 *   [omega^4/k^4 - (c_s^2 + v_A^2) omega^2/k^2 + c_s^2 v_A^2 cos^2(theta)] = 0
 */

#ifndef MHD_WAVES_H
#define MHD_WAVES_H

#include "mhd_defs.h"

/* ================================================================
 * L2 -- MHD Wave Speeds
 *
 * For a given direction n (unit vector) and angle theta between k and B:
 *
 * Alfven speed:              c_A = v_A |cos(theta)|
 * Fast magnetosonic:         c_f^2 = 0.5*(c_s^2 + v_A^2 + sqrt((c_s^2+v_A^2)^2 - 4 c_s^2 v_A^2 cos^2(theta)))
 * Slow magnetosonic:         c_sl^2 = 0.5*(c_s^2 + v_A^2 - sqrt((c_s^2+v_A^2)^2 - 4 c_s^2 v_A^2 cos^2(theta)))
 *
 * Ordering: c_sl <= c_A <= c_f always.
 * ================================================================ */

/*
 * mhd_wavespeeds_compute -- Compute all 4 characteristic speeds
 *
 * Given primitive state, propagation direction n, and gamma,
 * fills the MHDWaveSpeeds struct.
 *
 * The speeds are the eigenvalues of the MHD flux Jacobian A_n = dF_n/dU.
 * Complexity: O(1)
 */
void mhd_wavespeeds_compute(const MHDState *state,
                             double nx, double ny, double nz,
                             double gamma, MHDWaveSpeeds *ws);

/*
 * mhd_fast_magnetosonic_speed -- Fast magnetosonic speed
 * c_f^2 = 0.5 * [c_s^2 + v_A^2 + sqrt((c_s^2+v_A^2)^2 - 4 c_s^2 v_A^2 cos^2(theta))]
 *
 * This is the fastest compressional mode. Perpendicular to B (theta = pi/2):
 *   c_f = sqrt(c_s^2 + v_A^2) -- combines thermal + magnetic pressure.
 * Parallel to B (theta = 0):
 *   c_f = max(c_s, v_A) -- faster of sound and Alfven speed.
 *
 * Complexity: O(1)
 */
double mhd_fast_magnetosonic_speed(double cs, double va, double cos_theta);

/*
 * mhd_slow_magnetosonic_speed -- Slow magnetosonic speed
 * c_sl^2 = 0.5 * [c_s^2 + v_A^2 - sqrt((c_s^2+v_A^2)^2 - 4 c_s^2 v_A^2 cos^2(theta))]
 *
 * This is the slowest compressional mode. Perpendicular to B: c_sl = 0.
 * Parallel to B: c_sl = min(c_s, v_A).
 *
 * Complexity: O(1)
 */
double mhd_slow_magnetosonic_speed(double cs, double va, double cos_theta);

/*
 * mhd_alfven_wave_speed_projected -- Alfven speed projected to direction n
 * c_A = |v_A cos(theta)| = |(B . n)| / sqrt(mu_0 * rho)
 *
 * The Alfven wave propagates only along the magnetic field.
 * It is purely transverse (incompressible): div(v) = 0, no density perturbation.
 *
 * Complexity: O(1)
 */
double mhd_alfven_wave_speed_projected(const MHDState *state,
                                        double nx, double ny, double nz);

/*
 * mhd_friedrichs_diagram -- Compute phase speed as function of angle
 *
 * Returns the phase speed v_ph(theta) for the 3 wave modes.
 * The Friedrichs diagram is a polar plot of phase speeds:
 *   - Alfven: figure-8 (two tangent circles)
 *   - Fast: quasi-ellipse
 *   - Slow: figure-8 inside Alfven (for v_A > c_s) or vice versa
 *
 * Output: v_ph_alfven, v_ph_fast, v_ph_slow for given angle theta.
 *
 * Complexity: O(1)
 */
void mhd_friedrichs_diagram(double cs, double va, double theta,
                              double *v_alfven, double *v_fast, double *v_slow);

/*
 * mhd_group_velocity -- Compute group velocity for MHD waves
 *
 * The group velocity v_g = d(omega)/dk determines energy propagation.
 * For MHD waves, v_g || B for Alfven waves.
 * For fast mode, v_g is nearly isotropic.
 * For slow mode, v_g is nearly along B.
 *
 * Returns components of group velocity in (B_parallel, B_perp) frame.
 * Complexity: O(1)
 */
void mhd_group_velocity(double cs, double va, double theta,
                         int mode, /* 0=Alfven, 1=fast, 2=slow */
                         double *vg_par, double *vg_perp);

/*
 * mhd_wave_decomposition -- Decompose perturbation into MHD wave modes
 *
 * Given a perturbation (drho, dv, dB, dp), decomposes into
 * amplitudes of Alfven, fast, slow, and entropy modes.
 *
 * Based on the eigenvector decomposition of the MHD equations.
 * Reference: Roe & Balsara (1996), "Notes on the Eigensystem of MHD"
 *
 * Complexity: O(1)
 */
void mhd_wave_decomposition(const MHDState *state0,
                             const MHDState *perturbation,
                             double nx, double ny, double nz,
                             double gamma,
                             double *amp_alfven,
                             double *amp_fast,
                             double *amp_slow,
                             double *amp_entropy);

/*
 * mhd_wave_phase_velocity -- General phase velocity for given k vector
 * v_ph = (omega/k) k_hat
 *
 * Solves the full MHD dispersion relation for arbitrary k direction.
 *
 * Complexity: O(1)
 */
void mhd_wave_phase_velocity(const MHDState *state,
                              double kx, double ky, double kz,
                              double gamma,
                              double *vph_x, double *vph_y, double *vph_z,
                              int *num_modes);

/*
 * mhd_alfven_wave_dispersion -- Alfven wave dispersion relation with resistivity
 * omega = k_par * v_A  (ideal)
 * omega = k_par * v_A - i * eta * k^2 / 2  (resistive, leading order damping)
 *
 * Resistivity introduces damping of Alfven waves at rate ~ eta*k^2/2.
 *
 * Complexity: O(1)
 */
void mhd_alfven_dispersion(double k_par, double k_perp, double va, double eta,
                            double *omega_real, double *omega_imag);

/*
 * mhd_magnetosonic_dispersion -- Fast/slow magnetosonic dispersion with resistivity
 *
 * Returns complex frequency omega = omega_r + i*gamma for fast and slow modes.
 * Resistive damping scales as eta*k^2 for both compressional modes.
 *
 * Complexity: O(1)
 */
void mhd_magnetosonic_dispersion(double k_par, double k_perp,
                                  double cs, double va, double eta,
                                  double *omega_fast_r, double *omega_fast_i,
                                  double *omega_slow_r, double *omega_slow_i);

/*
 * mhd_wave_cfl -- Compute CFL-stable timestep
 * dt <= CFL * min(dx, dy, dz) / max_wavespeed
 *
 * Returns the maximum allowed timestep for explicit time stepping.
 *
 * Complexity: O(n) for n states, O(1) per state
 */
double mhd_wave_cfl(const MHDState *states, int n,
                     double dx, double dy, double dz,
                     double gamma, double cfl_number);

#endif /* MHD_WAVES_H */
