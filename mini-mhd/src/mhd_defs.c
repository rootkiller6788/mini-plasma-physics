/*
 * mhd_defs.c -- Core MHD Definitions Implementation
 *
 * Implements the non-inline functions declared in mhd_defs.h:
 * state conversions, flux computation, wave speeds, plasma parameters.
 *
 * Reference: Goldston & Rutherford (1995), Freidberg (2014)
 */

#include "mhd_defs.h"
#include <stdio.h>
#include <stdlib.h>

/* ================================================================
 * L1/L3 -- State Conversion Functions
 * ================================================================ */

/*
 * mhd_primitive_to_conserved -- Convert primitive to conserved variables
 *
 * Mapping:
 *   D  = rho
 *   M  = rho * v
 *   B  = B         (unchanged)
 *   E  = p/(gamma-1) + rho*v^2/2 + B^2/(2*mu_0)
 *
 * This is a direct (non-iterative) conversion.
 * The inverse requires solving a quadratic for pressure.
 */
void mhd_primitive_to_conserved(const MHDState *prim, double gamma,
                                 MHDConserved *cons) {
    if (!prim || !cons) return;

    cons->D       = prim->rho;
    cons->Mx      = prim->rho * prim->vx;
    cons->My      = prim->rho * prim->vy;
    cons->Mz      = prim->rho * prim->vz;
    cons->Bx      = prim->Bx;
    cons->By      = prim->By;
    cons->Bz      = prim->Bz;

    double v_mag   = mhd_vector_magnitude(prim->vx, prim->vy, prim->vz);
    double B_mag   = mhd_vector_magnitude(prim->Bx, prim->By, prim->Bz);
    cons->E        = mhd_total_energy_density(prim->rho, v_mag, B_mag,
                                               prim->p, gamma);
}

/*
 * mhd_conserved_to_primitive -- Convert conserved to primitive variables
 *
 * Inverse mapping. Requires solving for pressure from total energy:
 *   p = (gamma-1) * (E - M^2/(2*D) - B^2/(2*mu_0))
 *
 * Guard: if D <= 0 or p < 0, state is unphysical (set to small positive).
 */
void mhd_conserved_to_primitive(const MHDConserved *cons, double gamma,
                                 MHDState *prim) {
    if (!cons || !prim) return;

    prim->rho = cons->D;
    if (prim->rho < 1e-40) {
        /* Unphysical: zero density. Clamp to small value. */
        prim->rho = 1e-40;
        prim->vx = prim->vy = prim->vz = 0.0;
        prim->Bx = cons->Bx; prim->By = cons->By; prim->Bz = cons->Bz;
        prim->p = 1e-40;
        return;
    }

    double inv_rho = 1.0 / prim->rho;
    prim->vx = cons->Mx * inv_rho;
    prim->vy = cons->My * inv_rho;
    prim->vz = cons->Mz * inv_rho;
    prim->Bx = cons->Bx;
    prim->By = cons->By;
    prim->Bz = cons->Bz;

    double v_mag  = mhd_vector_magnitude(prim->vx, prim->vy, prim->vz);
    double B_mag  = mhd_vector_magnitude(prim->Bx, prim->By, prim->Bz);
    double ke     = mhd_kinetic_energy_density(prim->rho, v_mag);
    double me     = mhd_magnetic_energy_density(B_mag);

    prim->p = (gamma - 1.0) * (cons->E - ke - me);
    if (prim->p < 0.0) {
        /* Negative pressure: vacuum or numerical error floor */
        prim->p = 1e-40;
    }
}

/* ================================================================
 * L2/L4 -- MHD Flux Computation
 *
 * The flux vector F(U) in direction n = (nx, ny, nz):
 *
 * F_rho        = rho * v_n
 * F_rho_v      = rho * v * v_n + p_tot * n - (B_n/mu_0) * B
 * F_B          = v_n * B - B_n * v
 * F_E          = (E + p_tot) * v_n - (B_n/mu_0) * (v . B)
 *
 * where v_n = v.n, B_n = B.n, p_tot = p + B^2/(2*mu_0).
 *
 * Note: F_Bn = 0 (B_n is transported as a zero flux in the normal direction).
 * The Powell 8-wave formulation adds a source term for div(B) errors.
 * ================================================================ */

void mhd_flux_compute(const MHDState *state, double nx, double ny, double nz,
                       double gamma, MHDFlux *flux) {
    if (!state || !flux) return;

    double vn  = state->vx * nx + state->vy * ny + state->vz * nz;
    double Bn  = state->Bx * nx + state->By * ny + state->Bz * nz;
    double v_dot_B = state->vx * state->Bx + state->vy * state->By
                   + state->vz * state->Bz;
    double B2  = state->Bx * state->Bx + state->By * state->By
               + state->Bz * state->Bz;
    double p_tot = state->p + 0.5 * B2 * MHD_MU0_INV;
    double inv_mu0 = MHD_MU0_INV;

    double v2 = state->vx * state->vx + state->vy * state->vy
              + state->vz * state->vz;
    double E_total = state->p / (gamma - 1.0)
                   + 0.5 * state->rho * v2 + 0.5 * B2 * inv_mu0;

    /* Mass flux */
    flux->D_flux = state->rho * vn;

    /* Momentum flux: rho*v*v_n + p_tot*n - (B_n/mu_0)*B */
    flux->Mx_flux = state->rho * state->vx * vn + p_tot * nx
                  - (Bn * inv_mu0) * state->Bx;
    flux->My_flux = state->rho * state->vy * vn + p_tot * ny
                  - (Bn * inv_mu0) * state->By;
    flux->Mz_flux = state->rho * state->vz * vn + p_tot * nz
                  - (Bn * inv_mu0) * state->Bz;

    /* Induction flux: v_n*B - B_n*v  (F_Bn = 0 identically) */
    flux->Bx_flux = vn * state->Bx - Bn * state->vx;
    flux->By_flux = vn * state->By - Bn * state->vy;
    flux->Bz_flux = vn * state->Bz - Bn * state->vz;

    /* Energy flux: (E + p_tot)*v_n - (B_n/mu_0)*(v.B) */
    flux->E_flux = (E_total + p_tot) * vn - (Bn * inv_mu0) * v_dot_B;
}

/* ================================================================
 * L2 -- Maximum Wave Speed (for CFL condition)
 *
 * The maximum signal speed in direction n is:
 *   lambda_max = max(|v_n - c_f|, |v_n - c_s|, |v_n - c_a|,
 *                     |v_n + c_a|, |v_n + c_s|, |v_n + c_f|)
 *
 * but we can simply take: lambda_max = |v_n| + c_f
 * since c_f >= c_a, c_s always.
 * ================================================================ */

double mhd_max_wavespeed(const MHDState *state, double nx, double ny, double nz,
                          double gamma) {
    if (!state) return 0.0;

    double vn    = state->vx * nx + state->vy * ny + state->vz * nz;
    double Bn    = state->Bx * nx + state->By * ny + state->Bz * nz;
    double B_mag = mhd_vector_magnitude(state->Bx, state->By, state->Bz);

    double cs = mhd_sound_speed(state->p, state->rho, gamma);
    double va = mhd_alfven_speed(B_mag, state->rho);
    double ca = Bn / sqrt(MHD_MU0 * state->rho);

    /* c_f^2 = 0.5*(cs^2 + va^2 + sqrt((cs^2+va^2)^2 - 4*cs^2*ca^2)) */
    double sum = cs*cs + va*va;
    double disc = sum*sum - 4.0*cs*cs*ca*ca;
    if (disc < 0.0) disc = 0.0;
    double cf = sqrt(0.5 * (sum + sqrt(disc)));

    return fabs(vn) + cf;
}

/* ================================================================
 * L2 -- PlasmaParameters Computation
 *
 * Given a state and physical inputs, fills all dimensionless parameters.
 * Uses the inline functions from mhd_defs.h to compute each parameter.
 * ================================================================ */

void mhd_compute_all_parameters(const MHDState *state, double gamma,
                                 double L_char, double eta, double nu,
                                 double m_i, double n, double T_e, double T_i,
                                 double Omega, PlasmaParameters *params) {
    if (!state || !params) return;

    double B_mag = mhd_vector_magnitude(state->Bx, state->By, state->Bz);
    double v_mag = mhd_vector_magnitude(state->vx, state->vy, state->vz);

    params->beta           = mhd_plasma_beta(state->p, B_mag);
    params->mach_alfven    = mhd_mach_alfven(v_mag, B_mag, state->rho);
    params->mach_sonic     = mhd_mach_sonic(v_mag, state->p, state->rho, gamma);
    params->magnetic_reynolds = mhd_magnetic_reynolds(v_mag, L_char, eta);
    params->lundquist      = mhd_lundquist(B_mag, state->rho, L_char, eta);
    params->magnetic_prandtl = mhd_magnetic_prandtl(nu, eta);
    params->hartmann       = mhd_hartmann(B_mag, L_char, state->rho, nu, eta);
    params->ion_inertial_length = mhd_ion_inertial_length(n, m_i);
    params->ion_larmor_radius   = mhd_ion_larmor_radius(T_i, B_mag, m_i);
    params->debye_length   = mhd_debye_length(T_e, n);
    params->elsasser       = mhd_elsasser(B_mag, state->rho, eta, Omega);

    /* Electron and ion plasma beta */
    double factor = 2.0 * MHD_MU0 * n * MHD_KB / (B_mag * B_mag + 1e-60);
    params->electron_plasma_beta = factor * T_e;
    params->ion_plasma_beta      = factor * T_i;
}

/* ================================================================
 * L1 -- MHDGeometry initialization
 * ================================================================ */

void mhd_geometry_init(MHDGeometry *geom,
                        MHDCoordinateSystem coord,
                        double Lx, double Ly, double Lz,
                        int nx, int ny, int nz) {
    if (!geom) return;
    geom->coord = coord;
    geom->Lx = Lx; geom->Ly = Ly; geom->Lz = Lz;
    geom->nx = nx; geom->ny = ny; geom->nz = nz;
    geom->dx = (nx > 1) ? Lx / (double)nx : 0.0;
    geom->dy = (ny > 1) ? Ly / (double)ny : 0.0;
    geom->dz = (nz > 1) ? Lz / (double)nz : 0.0;
    geom->x0 = 0.0; geom->y0 = 0.0; geom->z0 = 0.0;
}

/*
 * mhd_state_init_uniform -- Initialize a uniform MHD state
 * Useful for setting initial conditions in simulations.
 */
void mhd_state_init_uniform(MHDState *state,
                             double rho, double vx, double vy, double vz,
                             double Bx, double By, double Bz, double p) {
    if (!state) return;
    state->rho = rho;
    state->vx  = vx;  state->vy = vy;  state->vz = vz;
    state->Bx  = Bx;  state->By = By;  state->Bz = Bz;
    state->p   = p;
}

/*
 * mhd_state_print -- Print MHD state for debugging
 */
void mhd_state_print(const MHDState *state) {
    if (!state) { printf("MHDState: NULL\n"); return; }
    printf("MHDState: rho=%.6e v=(%.6e,%.6e,%.6e) "
           "B=(%.6e,%.6e,%.6e) p=%.6e\n",
           state->rho, state->vx, state->vy, state->vz,
           state->Bx, state->By, state->Bz, state->p);
}

/*
 * mhd_params_print -- Print plasma parameters
 */
void mhd_params_print(const PlasmaParameters *p) {
    if (!p) { printf("PlasmaParameters: NULL\n"); return; }
    printf("Plasma Parameters:\n");
    printf("  beta              = %.6e\n", p->beta);
    printf("  M_A (Alfven Mach) = %.6e\n", p->mach_alfven);
    printf("  M_s (Sonic Mach)  = %.6e\n", p->mach_sonic);
    printf("  R_m (Mag Reynolds)= %.6e\n", p->magnetic_reynolds);
    printf("  S   (Lundquist)   = %.6e\n", p->lundquist);
    printf("  P_m (Mag Prandtl) = %.6e\n", p->magnetic_prandtl);
    printf("  Ha  (Hartmann)    = %.6e\n", p->hartmann);
    printf("  d_i (Ion skin)    = %.6e m\n", p->ion_inertial_length);
    printf("  rho_i (Larmor)    = %.6e m\n", p->ion_larmor_radius);
    printf("  lambda_D (Debye)  = %.6e m\n", p->debye_length);
    printf("  Lambda (Elsasser) = %.6e\n", p->elsasser);
}
