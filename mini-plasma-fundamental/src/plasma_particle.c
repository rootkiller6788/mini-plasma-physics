/**
 * plasma_particle.c -- Particle Methods for Plasma Simulation
 *
 * Implements Boris particle pusher, guiding center drifts,
 * PIC method basics, and tokamak orbit parameters.
 *
 * References:
 *   - Birdsall & Langdon, "Plasma Physics via Computer Simulation" (1991)
 *   - Boris (1970), Proc. 4th Conf. Num. Sim. Plasmas
 *   - Hockney & Eastwood, "Computer Simulation Using Particles" (1988)
 *
 * Knowledge Coverage:
 *   L2: Particle drifts (ExB, grad-B, curvature, polarization)
 *   L4: Lorentz force, guiding center equations
 *   L5: Boris pusher, PIC deposit/interpolate, leapfrog
 *   L6: Tokamak particle orbits (banana, bounce, precession)
 */
#include "plasma_particle.h"
#include "plasma_params.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ================================================================
 * L5: Boris Particle Pusher
 * ================================================================ */

void boris_push(Particle3D *p, double Ex, double Ey, double Ez,
                double Bx, double By, double Bz, double dt) {
    if (!p) return;
    double q_over_m = p->charge / p->mass;
    double h = 0.5 * q_over_m * dt;

    /* Half electric acceleration */
    double vx_minus = p->vx + h * Ex;
    double vy_minus = p->vy + h * Ey;
    double vz_minus = p->vz + h * Ez;

    /* Magnetic rotation via Boris scheme */
    double B2 = Bx*Bx + By*By + Bz*Bz;
    double t_mag = h / sqrt(1.0 + h*h * B2);
    double bx = Bx * t_mag, by = By * t_mag, bz = Bz * t_mag;
    double b2 = bx*bx + by*by + bz*bz;

    /* v' = v^- + v^- x b */
    double vpx = vx_minus + (vy_minus * bz - vz_minus * by);
    double vpy = vy_minus + (vz_minus * bx - vx_minus * bz);
    double vpz = vz_minus + (vx_minus * by - vy_minus * bx);

    /* v^+ = v^- + (2/(1+b^2)) * (v' x b) */
    double scale = 2.0 / (1.0 + b2);
    double vx_plus = vx_minus + scale * (vpy * bz - vpz * by);
    double vy_plus = vy_minus + scale * (vpz * bx - vpx * bz);
    double vz_plus = vz_minus + scale * (vpx * by - vpy * bx);

    /* Half electric acceleration */
    p->vx = vx_plus + h * Ex;
    p->vy = vy_plus + h * Ey;
    p->vz = vz_plus + h * Ez;

    /* Position update */
    p->x += p->vx * dt;
    p->y += p->vy * dt;
    p->z += p->vz * dt;
}

void boris_push_field(Particle3D *p, double t, double dt,
                      FieldFunc3D E_field, FieldFunc3D B_field) {
    if (!p || !E_field || !B_field) return;
    double Ex, Ey, Ez, Bx, By, Bz;
    E_field(p->x, p->y, p->z, t, &Ex, &Ey, &Ez, &Bx, &By, &Bz);
    double Bx_copy = Bx, By_copy = By, Bz_copy = Bz;
    B_field(p->x, p->y, p->z, t, &Ex, &Ey, &Ez, &Bx, &By, &Bz);
    boris_push(p, Ex, Ey, Ez, Bx_copy, By_copy, Bz_copy, dt);
}

void boris_push_ensemble(ParticleEnsemble *ensemble,
                         FieldFunc3D E_field, FieldFunc3D B_field,
                         double t, double dt) {
    if (!ensemble || !E_field || !B_field) return;
    for (int i = 0; i < ensemble->n_particles; i++) {
        if (ensemble->particles[i].alive) {
            boris_push_field(&ensemble->particles[i], t, dt,
                             E_field, B_field);
        }
    }
}

/* ================================================================
 * L2: Particle Drifts
 * ================================================================ */

void exb_drift_velocity(double Ex, double Ey, double Ez,
                        double Bx, double By, double Bz,
                        double *vdx, double *vdy, double *vdz) {
    double B2 = Bx*Bx + By*By + Bz*Bz;
    if (B2 < 1e-30) { *vdx = *vdy = *vdz = 0.0; return; }
    double inv_B2 = 1.0 / B2;
    *vdx = (Ey*Bz - Ez*By) * inv_B2;
    *vdy = (Ez*Bx - Ex*Bz) * inv_B2;
    *vdz = (Ex*By - Ey*Bx) * inv_B2;
}

void gradb_drift_velocity(double mu, double q, double B_mag,
                          double gradBx, double gradBy, double gradBz,
                          double Bx, double By, double Bz,
                          double *vdx, double *vdy, double *vdz) {
    if (q == 0.0 || B_mag < 1e-30) {
        *vdx = *vdy = *vdz = 0.0; return;
    }
    double coeff = mu / (q * B_mag * B_mag * B_mag);
    *vdx = coeff * (By*gradBz - Bz*gradBy);
    *vdy = coeff * (Bz*gradBx - Bx*gradBz);
    *vdz = coeff * (Bx*gradBy - By*gradBx);
}

void curvature_drift_velocity(double m, double v_par, double q,
                              double Bx, double By, double Bz,
                              double curv_x, double curv_y, double curv_z,
                              double *vdx, double *vdy, double *vdz) {
    double B2 = Bx*Bx + By*By + Bz*Bz;
    if (q == 0.0 || B2 < 1e-30) {
        *vdx = *vdy = *vdz = 0.0; return;
    }
    double coeff = m * v_par * v_par / (q * B2 * B2);
    *vdx = coeff * (By*curv_z - Bz*curv_y);
    *vdy = coeff * (Bz*curv_x - Bx*curv_z);
    *vdz = coeff * (Bx*curv_y - By*curv_x);
}

void total_guiding_center_velocity(double Ex, double Ey, double Ez,
                                   double Bx, double By, double Bz,
                                   double gradB_mag, double v_par,
                                   double mu, double m, double q,
                                   double *vdx, double *vdy, double *vdz) {
    (void)v_par; (void)m; /* reserved for curvature drift in full implementation */
    double v_ExB[3], v_gradB[3];
    exb_drift_velocity(Ex, Ey, Ez, Bx, By, Bz,
                       &v_ExB[0], &v_ExB[1], &v_ExB[2]);
    double B_mag = sqrt(Bx*Bx + By*By + Bz*Bz);
    if (B_mag < 1e-30) {
        *vdx = v_ExB[0]; *vdy = v_ExB[1]; *vdz = v_ExB[2]; return;
    }
    double gradBx = gradB_mag * Bx / B_mag;
    double gradBy = gradB_mag * By / B_mag;
    double gradBz = gradB_mag * Bz / B_mag;
    gradb_drift_velocity(mu, q, B_mag, gradBx, gradBy, gradBz,
                         Bx, By, Bz, &v_gradB[0], &v_gradB[1], &v_gradB[2]);
    *vdx = v_ExB[0] + v_gradB[0];
    *vdy = v_ExB[1] + v_gradB[1];
    *vdz = v_ExB[2] + v_gradB[2];
}

void polarization_drift_velocity(double m, double q, double B_mag,
                                 double dEx_dt, double dEy_dt, double dEz_dt,
                                 double *vdx, double *vdy, double *vdz) {
    if (q == 0.0 || B_mag < 1e-30) {
        *vdx = *vdy = *vdz = 0.0; return;
    }
    double coeff = m / (q * B_mag * B_mag);
    *vdx = coeff * dEx_dt;
    *vdy = coeff * dEy_dt;
    *vdz = coeff * dEz_dt;
}

/* ================================================================
 * L6: Tokamak Orbit Parameters
 * ================================================================ */

double banana_width(double rho_L, double epsilon, double q_safety) {
    if (epsilon <= 0.0) return 0.0;
    return q_safety * rho_L / sqrt(epsilon);
}

double trapped_particle_fraction(double epsilon) {
    if (epsilon < 0.0) return 0.0;
    if (epsilon > 0.5) epsilon = 0.5;
    return sqrt(2.0 * epsilon);
}

double bounce_frequency(double epsilon, double v_th, double q_safety,
                        double R0) {
    if (q_safety <= 0.0 || R0 <= 0.0) return 0.0;
    return sqrt(epsilon) * v_th / (q_safety * R0);
}

double toroidal_precession_frequency(double v_th, double omega_c,
                                     double R0, double r, double q,
                                     double shear) {
    if (omega_c <= 0.0 || R0 <= 0.0 || r <= 0.0) return 0.0;
    return q * v_th * v_th / (omega_c * R0 * r) * (shear + 1.0);
}

/* ================================================================
 * L5: PIC Method Basics
 * ================================================================ */

void pic_deposit_ngp(PICGrid1D *grid, double x, double q) {
    if (!grid || !grid->rho) return;
    int ix = (int)(x / grid->dx);
    if (ix >= 0 && ix < grid->nx) {
        grid->rho[ix] += q / grid->dx;
    }
}

void pic_deposit_cic(PICGrid1D *grid, double x, double q) {
    if (!grid || !grid->rho) return;
    double x_over_dx = x / grid->dx;
    int ix = (int)x_over_dx;
    double wx = x_over_dx - ix;
    double q_over_dx = q / grid->dx;
    if (ix >= 0 && ix < grid->nx) {
        grid->rho[ix] += (1.0 - wx) * q_over_dx;
    }
    if (ix + 1 >= 0 && ix + 1 < grid->nx) {
        grid->rho[ix + 1] += wx * q_over_dx;
    }
}

double pic_interpolate_cic(const PICGrid1D *grid, double x) {
    if (!grid || !grid->E) return 0.0;
    double x_over_dx = x / grid->dx;
    int ix = (int)x_over_dx;
    double wx = x_over_dx - ix;
    double E = 0.0;
    if (ix >= 0 && ix < grid->nx) {
        E += (1.0 - wx) * grid->E[ix];
    }
    if (ix + 1 >= 0 && ix + 1 < grid->nx) {
        E += wx * grid->E[ix + 1];
    }
    return E;
}

void pic_solve_poisson(PICGrid1D *grid) {
    if (!grid || !grid->rho || !grid->phi) return;
    int nx = grid->nx;
    double dx = grid->dx;
    double *a = (double*)malloc(nx * sizeof(double));
    double *b = (double*)malloc(nx * sizeof(double));
    double *c = (double*)malloc(nx * sizeof(double));
    double *rhs = (double*)malloc(nx * sizeof(double));
    if (!a || !b || !c || !rhs) {
        free(a); free(b); free(c); free(rhs); return;
    }
    double dx2 = dx * dx;
    for (int i = 0; i < nx; i++) {
        a[i] = 1.0; b[i] = -2.0; c[i] = 1.0;
        rhs[i] = -grid->rho[i] * dx2 / EPSILON_0;
    }
    b[0] = 1.0; c[0] = 0.0; rhs[0] = 0.0;
    a[nx-1] = 0.0; b[nx-1] = 1.0; rhs[nx-1] = 0.0;
    for (int i = 1; i < nx; i++) {
        double w = a[i] / b[i-1];
        b[i] -= w * c[i-1];
        rhs[i] -= w * rhs[i-1];
    }
    grid->phi[nx-1] = rhs[nx-1] / b[nx-1];
    for (int i = nx - 2; i >= 0; i--) {
        grid->phi[i] = (rhs[i] - c[i] * grid->phi[i+1]) / b[i];
    }
    free(a); free(b); free(c); free(rhs);
}

void pic_compute_e_field(PICGrid1D *grid) {
    if (!grid || !grid->phi || !grid->E) return;
    double dx = grid->dx;
    for (int i = 1; i < grid->nx - 1; i++) {
        grid->E[i] = -(grid->phi[i+1] - grid->phi[i-1]) / (2.0 * dx);
    }
    grid->E[0] = -(grid->phi[1] - grid->phi[0]) / dx;
    grid->E[grid->nx - 1] =
        -(grid->phi[grid->nx - 1] - grid->phi[grid->nx - 2]) / dx;
}

void pic_cycle_1d(PICGrid1D *grid, ParticleEnsemble *ensemble,
                  double dt, double t) {
    if (!grid || !ensemble) return;
    (void)t; /* t is available for time-dependent fields */
    memset(grid->rho, 0, grid->nx * sizeof(double));
    for (int i = 0; i < ensemble->n_particles; i++) {
        if (ensemble->particles[i].alive) {
            pic_deposit_cic(grid, ensemble->particles[i].x,
                            ensemble->particles[i].charge
                            * ensemble->weight);
        }
    }
    pic_solve_poisson(grid);
    pic_compute_e_field(grid);
    for (int i = 0; i < ensemble->n_particles; i++) {
        if (ensemble->particles[i].alive) {
            double E_at_p = pic_interpolate_cic(grid,
                                                ensemble->particles[i].x);
            boris_push(&ensemble->particles[i], E_at_p, 0.0, 0.0,
                       0.0, 0.0, 0.0, dt);
        }
    }
}

/* ================================================================
 * L6: Neoclassical Transport Quantities
 * ================================================================ */

double bootstrap_current_density(double epsilon, double B_theta,
                                  double dp_dr) {
    if (B_theta == 0.0) return 0.0;
    return -C_LIGHT * sqrt(epsilon) * dp_dr / B_theta;
}

double ware_pinch_velocity(double E_phi, double B_theta) {
    if (B_theta == 0.0) return 0.0;
    return E_phi / B_theta;
}

double pfirsch_schlueter_factor(double epsilon, double r, double R0) {
    if (R0 <= 0.0) return 0.0;
    return 2.0 * epsilon * r / R0;
}
