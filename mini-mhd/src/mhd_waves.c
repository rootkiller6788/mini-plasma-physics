/*
 * mhd_waves.c -- MHD Wave Dispersion and Propagation
 *
 * Implements MHD wave speed calculations, dispersion relations,
 * wave decomposition, and group velocity computations.
 *
 * Reference: Freidberg "Ideal MHD" (2014) Ch.6
 *            Goedbloed & Poedts "Principles of MHD" (2004) Ch.5
 *            Roe & Balsara (1996), "Notes on the Eigensystem of MHD"
 */

#include "mhd_waves.h"
#include "mhd_defs.h"
#include <math.h>
#include <stdio.h>

/* ================================================================
 * L2 -- Characteristic Wave Speeds
 *
 * The MHD system has 7 distinct eigenvalues (in 3D, one degenerate):
 *   lambda_1,7 = v_n +/- c_f   (fast magnetosonic)
 *   lambda_2,6 = v_n +/- c_s   (slow magnetosonic)
 *   lambda_3,5 = v_n +/- c_a   (Alfven)
 *   lambda_4   = v_n           (entropy / div B)
 *
 * where c_f >= c_a >= c_s always.
 *
 * The speeds are computed from the dispersion relation of the
 * linearized MHD equations. For propagation direction n
 * making angle theta with B:
 *
 * c_f^2 = 0.5*(cs^2 + va^2 + sqrt((cs^2+va^2)^2 - 4*cs^2*va^2*cos^2(theta)))
 * c_s^2 = 0.5*(cs^2 + va^2 - sqrt((cs^2+va^2)^2 - 4*cs^2*va^2*cos^2(theta)))
 * c_a   = va * |cos(theta)|
 * ================================================================ */

double mhd_fast_magnetosonic_speed(double cs, double va, double cos_theta) {
    double sum = cs*cs + va*va;
    double disc = sum*sum - 4.0*cs*cs*va*va*cos_theta*cos_theta;
    if (disc < 0.0) disc = 0.0;
    return sqrt(0.5 * (sum + sqrt(disc)));
}

double mhd_slow_magnetosonic_speed(double cs, double va, double cos_theta) {
    double sum = cs*cs + va*va;
    double disc = sum*sum - 4.0*cs*cs*va*va*cos_theta*cos_theta;
    if (disc < 0.0) disc = 0.0;
    double val = 0.5 * (sum - sqrt(disc));
    if (val < 0.0) val = 0.0;
    return sqrt(val);
}

double mhd_alfven_wave_speed_projected(const MHDState *state,
                                        double nx, double ny, double nz) {
    if (!state) return 0.0;

    double Bn = state->Bx*nx + state->By*ny + state->Bz*nz;
    double B_mag = mhd_vector_magnitude(state->Bx, state->By, state->Bz);

    /* cos(theta) = |B.n| / (|B| * |n|) = |B.n| / |B| */
    double norm_n = sqrt(nx*nx + ny*ny + nz*nz);
    if (norm_n < 1e-40 || B_mag < 1e-40) return 0.0;
    double cos_theta = fabs(Bn) / (B_mag * norm_n);

    double va = mhd_alfven_speed(B_mag, state->rho);
    return va * cos_theta;
}

/* ================================================================
 * L2 -- Complete Wave Speed Computation
 *
 * Fills the MHDWaveSpeeds struct with all 4 characteristic speeds
 * for a given propagation direction.
 * ================================================================ */

void mhd_wavespeeds_compute(const MHDState *state,
                             double nx, double ny, double nz,
                             double gamma, MHDWaveSpeeds *ws) {
    if (!state || !ws) return;

    double B_mag = mhd_vector_magnitude(state->Bx, state->By, state->Bz);
    double norm_n = sqrt(nx*nx + ny*ny + nz*nz);
    if (norm_n < 1e-40) {
        /* Degenerate direction: all speeds = 0 */
        ws->c_fast = ws->c_slow = ws->c_alfven = ws->c_sound = 0.0;
        return;
    }

    double Bn = state->Bx*nx + state->By*ny + state->Bz*nz;
    double cos_theta = 0.0;
    if (B_mag > 1e-40) {
        cos_theta = fabs(Bn) / (B_mag * norm_n);
    }

    double cs = mhd_sound_speed(state->p, state->rho, gamma);
    double va = mhd_alfven_speed(B_mag, state->rho);

    ws->c_sound  = cs;
    ws->c_alfven = va * cos_theta;
    ws->c_fast   = mhd_fast_magnetosonic_speed(cs, va, cos_theta);
    ws->c_slow   = mhd_slow_magnetosonic_speed(cs, va, cos_theta);
}

/* ================================================================
 * L6 -- Friedrichs Diagram
 *
 * Polar plot of phase speeds v_ph(theta) for the 3 MHD wave modes.
 * The Friedrichs diagram shows the anisotropy of wave propagation:
 *
 * - Alfven mode: two tangent circles (figure-8) -- v_ph = v_A |cos(theta)|
 * - Fast mode: quasi-ellipse, fastest perpendicular to B
 * - Slow mode: figure-8 inside Alfven circle (for v_A > c_s)
 *
 * The topology depends on the ratio v_A/c_s:
 *   v_A > c_s: fast mode is isotropic, slow mode highly anisotropic
 *   v_A < c_s: fast mode anisotropic, slow mode nearly isotropic
 * ================================================================ */

void mhd_friedrichs_diagram(double cs, double va, double theta,
                              double *v_alfven, double *v_fast, double *v_slow) {
    double ct = cos(theta);
    *v_alfven = va * fabs(ct);
    *v_fast   = mhd_fast_magnetosonic_speed(cs, va, ct);
    *v_slow   = mhd_slow_magnetosonic_speed(cs, va, ct);
}

/* ================================================================
 * L3 -- Group Velocity
 *
 * The group velocity v_g = d(omega)/dk determines the direction
 * and speed of energy propagation. For MHD waves:
 *
 * Alfven: v_g = v_A * sign(cos(theta)) * B_hat
 *   Energy propagates exactly along B.
 *
 * Fast mode: v_g is between k_hat and B_hat
 * Slow mode: v_g is between k_hat and -B_hat (for v_A > c_s)
 *
 * We compute v_g in the frame where B is along the z-axis
 * and k is in the xz-plane: k = (k*sin(theta), 0, k*cos(theta)).
 * ================================================================ */

void mhd_group_velocity(double cs, double va, double theta,
                         int mode,
                         double *vg_par, double *vg_perp) {
    double ct = cos(theta), st = sin(theta);

    if (mode == 0) {
        /* Alfven wave: v_g || B, v_g = v_A * sign(cos(theta)) */
        double sign = (ct >= 0) ? 1.0 : -1.0;
        *vg_par  = sign * va;
        *vg_perp = 0.0;
    } else {
        /* Fast (mode==1) or Slow (mode==2) magnetosonic */
        double v_ph = (mode == 1)
            ? mhd_fast_magnetosonic_speed(cs, va, ct)
            : mhd_slow_magnetosonic_speed(cs, va, ct);

        /* d(v_ph)/d(theta) via finite difference or analytic formula */
        /* Analytic: d(v_ph^2)/dtheta = +/- (2*cs^2*va^2*sin(2*theta)) / (2*sqrt(...)) */
        double sum = cs*cs + va*va;
        double disc = sum*sum - 4.0*cs*cs*va*va*ct*ct;
        if (disc < 1e-40) disc = 1e-40;
        double sqrt_disc = sqrt(disc);

        double sign = (mode == 1) ? 1.0 : -1.0;
        double dvph2_dtheta = sign * 4.0*cs*cs*va*va*ct*st / sqrt_disc;
        double dvph_dtheta = 0.0;
        if (v_ph > 1e-40) dvph_dtheta = 0.5 * dvph2_dtheta / v_ph;

        /* v_g = v_ph * k_hat + d(v_ph)/d(theta) * theta_hat */
        *vg_par  = v_ph * ct - dvph_dtheta * st;
        *vg_perp = v_ph * st + dvph_dtheta * ct;
    }
}

/* ================================================================
 * L3 -- Wave Decomposition
 *
 * Decomposes an MHD perturbation into Alfven, fast, slow, and entropy
 * wave amplitudes. Based on the eigenvector decomposition of the
 * linearized MHD flux Jacobian (Roe & Balsara 1996).
 *
 * The 8 degrees of freedom decompose as:
 *   1 entropy wave    (density/temperature perturbation at constant p+B^2/2mu_0)
 *   2 Alfven waves    (+/- propagating, incompressible transverse)
 *   2 fast waves      (+/- propagating, compressible)
 *   2 slow waves      (+/- propagating, compressible)
 *   1 div(B) mode     (should be zero for physical states)
 *
 * We simplify by returning the amplitude of each mode type
 * (summing +/- pairs into a single amplitude).
 * ================================================================ */

void mhd_wave_decomposition(const MHDState *state0,
                             const MHDState *pert,
                             double nx, double ny, double nz,
                             double gamma,
                             double *amp_alfven,
                             double *amp_fast,
                             double *amp_slow,
                             double *amp_entropy) {
    if (!state0 || !pert) {
        *amp_alfven = *amp_fast = *amp_slow = *amp_entropy = 0.0;
        return;
    }

    double B_mag = mhd_vector_magnitude(state0->Bx, state0->By, state0->Bz);
    double norm_n = sqrt(nx*nx + ny*ny + nz*nz);
    if (B_mag < 1e-40 || norm_n < 1e-40) {
        *amp_alfven = *amp_fast = *amp_slow = *amp_entropy = 0.0;
        return;
    }

    /* Unit vectors: propagation direction n_hat, B direction b_hat */
    double nx_u = nx/norm_n, ny_u = ny/norm_n, nz_u = nz/norm_n;
    double bx = state0->Bx/B_mag, by = state0->By/B_mag, bz = state0->Bz/B_mag;

    /* Perturbation in B parallel and perpendicular to B */
    double dB_par = pert->Bx*bx + pert->By*by + pert->Bz*bz;
    double dB_perp_x = pert->Bx - dB_par*bx;
    double dB_perp_y = pert->By - dB_par*by;
    double dB_perp_z = pert->Bz - dB_par*bz;
    double dB_perp = sqrt(dB_perp_x*dB_perp_x + dB_perp_y*dB_perp_y
                         + dB_perp_z*dB_perp_z);

    /* Perturbation in velocity */
    double dv_par = pert->vx*bx + pert->vy*by + pert->vz*bz;
    double dv_perp_x = pert->vx - dv_par*bx;
    double dv_perp_y = pert->vy - dv_par*by;
    double dv_perp_z = pert->vz - dv_par*bz;
    double dv_perp = sqrt(dv_perp_x*dv_perp_x + dv_perp_y*dv_perp_y
                         + dv_perp_z*dv_perp_z);

    double va = mhd_alfven_speed(B_mag, state0->rho);

    /* Alfven wave: transverse velocity = +/- transverse B/sqrt(mu_0*rho) */
    *amp_alfven = 0.5 * (fabs(dv_perp - dB_perp/sqrt(MHD_MU0*state0->rho))
                      + fabs(dv_perp + dB_perp/sqrt(MHD_MU0*state0->rho)));

    /* Fast and slow: compressible modes (simplified amplitude estimate) */
    double dp = pert->p;
    double drho = pert->rho;
    double cs = mhd_sound_speed(state0->p, state0->rho, gamma);

    double v_f = mhd_fast_magnetosonic_speed(cs, va,
                   fabs(nx_u*bx + ny_u*by + nz_u*bz));
    double v_s = mhd_slow_magnetosonic_speed(cs, va,
                   fabs(nx_u*bx + ny_u*by + nz_u*bz));

    *amp_fast = 0.5 * fabs(dp + state0->rho*v_f*dv_par)
              / (state0->rho*v_f*v_f + 1e-40);
    *amp_slow = 0.5 * fabs(dp + state0->rho*v_s*dv_par)
              / (state0->rho*v_s*v_s + 1e-40);

    /* Entropy: density perturbation at constant total pressure */
    double dp_tot = dp + (B_mag*dB_par)*MHD_MU0_INV;
    *amp_entropy = fabs(drho - dp_tot/(cs*cs + 1e-40));
}

/* ================================================================
 * L3 -- Phase Velocity for Arbitrary k
 *
 * For a given k-vector, computes the 3 MHD wave phase velocities:
 * v_ph = (omega/k) * k_hat where omega/k satisfies the dispersion relation.
 * ================================================================ */

void mhd_wave_phase_velocity(const MHDState *state,
                              double kx, double ky, double kz,
                              double gamma,
                              double *vph_x, double *vph_y, double *vph_z,
                              int *num_modes) {
    double k_mag = sqrt(kx*kx + ky*ky + kz*kz);
    if (k_mag < 1e-40 || !state) {
        vph_x[0]=vph_y[0]=vph_z[0]=vph_x[1]=vph_y[1]=vph_z[1]=0.0;
        vph_x[2]=vph_y[2]=vph_z[2]=0.0;
        *num_modes = 0;
        return;
    }

    double B_mag = mhd_vector_magnitude(state->Bx, state->By, state->Bz);
    double k_dot_B = kx*state->Bx + ky*state->By + kz*state->Bz;
    double cos_theta = fabs(k_dot_B) / (k_mag * B_mag + 1e-40);

    double cs = mhd_sound_speed(state->p, state->rho, gamma);
    double va = mhd_alfven_speed(B_mag, state->rho);

    double vA_proj = va * cos_theta;
    double vf = mhd_fast_magnetosonic_speed(cs, va, cos_theta);
    double vs = mhd_slow_magnetosonic_speed(cs, va, cos_theta);

    /* Phase velocities: v_ph = v_phase * k_hat */
    double kx_u = kx/k_mag, ky_u = ky/k_mag, kz_u = kz/k_mag;

    /* Alfven */
    vph_x[0] = vA_proj * kx_u; vph_y[0] = vA_proj * ky_u; vph_z[0] = vA_proj * kz_u;
    /* Fast */
    vph_x[1] = vf * kx_u; vph_y[1] = vf * ky_u; vph_z[1] = vf * kz_u;
    /* Slow */
    vph_x[2] = vs * kx_u; vph_y[2] = vs * ky_u; vph_z[2] = vs * kz_u;

    *num_modes = 3;
}

/* ================================================================
 * L5 -- Alfven Wave Dispersion with Resistivity
 *
 * Ideal: omega = k_par * v_A
 * Resistive: omega ~ k_par * v_A - i * eta * k^2 / 2 (damping)
 *
 * The imaginary part gives the damping rate.
 * ================================================================ */

void mhd_alfven_dispersion(double k_par, double k_perp, double va, double eta,
                            double *omega_real, double *omega_imag) {
    double k2 = k_par*k_par + k_perp*k_perp;
    *omega_real = k_par * va;
    *omega_imag = -0.5 * eta * k2;
}

/* ================================================================
 * L5 -- Magnetosonic Dispersion with Resistivity
 *
 * Fast and slow modes also experience resistive damping.
 *
 * omega^2 ~ k^2 * v_{f,s}^2 - i * eta * k^2 * omega
 * => omega ~ k*v_{f,s} - i * eta*k^2/2 (same leading-order damping as Alfven)
 * ================================================================ */

void mhd_magnetosonic_dispersion(double k_par, double k_perp,
                                  double cs, double va, double eta,
                                  double *omega_fast_r, double *omega_fast_i,
                                  double *omega_slow_r, double *omega_slow_i) {
    double k_mag = sqrt(k_par*k_par + k_perp*k_perp);
    double cos_theta = fabs(k_par) / (k_mag + 1e-40);
    double k2 = k_mag * k_mag;

    double vf = mhd_fast_magnetosonic_speed(cs, va, cos_theta);
    double vs = mhd_slow_magnetosonic_speed(cs, va, cos_theta);

    *omega_fast_r = k_mag * vf;
    *omega_fast_i = -0.5 * eta * k2;
    *omega_slow_r = k_mag * vs;
    *omega_slow_i = -0.5 * eta * k2;
}

/* ================================================================
 * L5 -- CFL Timestep for MHD Waves
 *
 * dt = CFL * min(dx, dy, dz) / max_i(|v_n| + c_f_i)
 *
 * This ensures numerical stability for explicit time integration.
 * ================================================================ */

double mhd_wave_cfl(const MHDState *states, int n,
                     double dx, double dy, double dz,
                     double gamma, double cfl_number) {
    if (!states || n <= 0) return 0.0;

    double max_speed = 0.0;
    double dl_min = dx;
    if (dy < dl_min) dl_min = dy;
    if (dz < dl_min) dl_min = dz;

    for (int i = 0; i < n; i++) {
        double B_mag = mhd_vector_magnitude(states[i].Bx, states[i].By, states[i].Bz);
        double cs = mhd_sound_speed(states[i].p, states[i].rho, gamma);
        double va = mhd_alfven_speed(B_mag, states[i].rho);
        double cf = mhd_fast_magnetosonic_speed(cs, va, 0.0);
        /* Use worst-case (perpendicular) fast speed */
        cf = sqrt(cs*cs + va*va);

        double v_mag = mhd_vector_magnitude(states[i].vx, states[i].vy, states[i].vz);
        double speed = v_mag + cf;
        if (speed > max_speed) max_speed = speed;
    }

    if (max_speed < 1e-40) return INFINITY;
    return cfl_number * dl_min / max_speed;
}
