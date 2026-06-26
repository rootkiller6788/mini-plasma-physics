/*
 * laser_propagation.c -- Laser propagation in plasma
 *
 * Implements models for electromagnetic wave propagation in
 * underdense plasmas: normalized vector potential, ponderomotive
 * force, relativistic transparency, self-focusing, and beam
 * propagation.
 *
 * References:
 *   - Gibbon (2005) "Short Pulse Laser Interactions with Matter"
 *   - Kruer (1988) "The Physics of Laser Plasma Interactions"
 *   - Max (1976), Phys. Fluids 19, 74
 *
 * Knowledge Layers: L1, L2, L4, L5
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "plasma_constants.h"
#include "plasma_params.h"
#include "laser_plasma.h"
#include "laser_propagation.h"

static const double c   = PLASMA_C;
static const double e   = PLASMA_E;
static const double me  = PLASMA_ME;
static const double eps0 = PLASMA_EPS0;

/* ============================================================
 *  L1: normalized_vector_potential
 *
 *  a0 = e E0 / (me c omega) = sqrt(e^2 I lambda^2 / (2 pi^2 eps0 me^2 c^5))
 *     ~ 0.855 * lambda_mum * sqrt(I_18)
 *
 * where I_18 = I / 10^18 W/cm^2.
 *
 * a0 is the single most important dimensionless parameter in
 * high-intensity laser-plasma physics.
 *
 * a0 << 1 : non-relativistic (ponderomotive force proportional to I)
 * a0 ~ 1  : mildly relativistic (1.37e18 W/cm^2 at 1 um)
 * a0 >> 1 : ultra-relativistic (v_osc ~ c, radiation reaction matters)
 *
 * Source: Gibbon Sec 2.1, Eq. 2.6
 *
 * Complexity: O(1)
 * ============================================================ */
double normalized_vector_potential(double intensity, double lambda_m)
{
    if (intensity <= 0.0 || lambda_m <= 0.0) return 0.0;

    double omega = 2.0 * M_PI * c / lambda_m;
    double E0 = sqrt(2.0 * intensity / (eps0 * c));
    return e * E0 / (me * c * omega);
}

/* ============================================================
 *  L1: normalized_vector_potential_E0 -- a0 from E-field directly
 *
 *  a0 = e E0 / (me c omega)
 *
 * Used when the electric field amplitude is known directly,
 * as in PIC simulations.
 *
 * Complexity: O(1)
 * ============================================================ */
double normalized_vector_potential_E0(double E0, double lambda_m)
{
    if (E0 <= 0.0 || lambda_m <= 0.0) return 0.0;
    double omega = 2.0 * M_PI * c / lambda_m;
    return e * E0 / (me * c * omega);
}

/* ============================================================
 *  L1: relativistic_intensity
 *
 *  I_rel = (pi^2/2) eps0 me^2 c^5 / (e^2 lambda^2)
 *        ~ 1.37e18 / lambda_mum^2  [W/cm^2]
 *
 * This is the intensity at which a0 = 1. It provides the
 * boundary between non-relativistic and relativistic regimes.
 *
 * For Ti:Sapphire (0.8 um): I_rel ~ 2.14e18 W/cm^2
 * For Nd:glass (1.053 um): I_rel ~ 1.24e18 W/cm^2
 * For CO2 (10.6 um):      I_rel ~ 1.22e16 W/cm^2
 *
 * Source: Gibbon Sec 2.1, Eq. 2.7
 *
 * Complexity: O(1)
 * ============================================================ */
double relativistic_intensity(double lambda_m)
{
    if (lambda_m <= 0.0) return 0.0;

    double prefactor = 2.0 * M_PI * M_PI * eps0 * me * me
                       * c * c * c * c * c;
    double denom = e * e * lambda_m * lambda_m;
    return prefactor / denom;
}

/* ============================================================
 *  L1: ponderomotive_potential -- Non-relativistic
 *
 *  Phi_p = e^2 E0^2 / (4 me omega^2)
 *        = (e^2 / (2 eps0 me c)) * I / omega^2
 *
 * This is the time-averaged oscillation (quiver) energy of
 * a free electron in the laser electric field.  It acts as
 * an effective repulsive potential:
 *
 *   f_p = -grad Phi_p   (ponderomotive force)
 *
 * For I=10^18 W/cm^2, lambda=1 um:
 *   Phi_p ~ 93 keV  (already relativistic!)
 *
 * Source: Kruer Sec 6.2; Gibbon Sec 2.1
 *
 * Complexity: O(1)
 * ============================================================ */
double ponderomotive_potential(double intensity, double lambda_m)
{
    if (intensity <= 0.0 || lambda_m <= 0.0) return 0.0;

    double omega = 2.0 * M_PI * c / lambda_m;
    double E0_sq = 2.0 * intensity / (eps0 * c);
    return e * e * E0_sq / (4.0 * me * omega * omega);
}

/* ============================================================
 *  L1: ponderomotive_potential_relativistic -- Relativistic
 *
 *  Phi_p^rel = me c^2 (sqrt(1 + a0^2) - 1)
 *
 * Non-relativistic limit (a0 << 1): Phi_p^rel ~ (1/2) me c^2 a0^2
 *                                  = e^2 E0^2 / (4 me omega^2)
 * which recovers the classical ponderomotive potential.
 *
 * Ultra-relativistic limit (a0 >> 1):
 *   Phi_p^rel ~ me c^2 a0 ~ me c^2 * e E0 / (me c omega) = e E0 c / omega
 *
 * Source: Gibbon Sec 2.1, Eq. 2.4
 *
 * Complexity: O(1)
 * ============================================================ */
double ponderomotive_potential_relativistic(double a0)
{
    if (a0 < 0.0) a0 = 0.0;
    return me * c * c * (sqrt(1.0 + a0 * a0) - 1.0);
}

/* ============================================================
 *  L2: ponderomotive_force -- Ponderomotive force density
 *
 *  f_p = -ne d Phi_p / dx
 *      ~ -ne * (e^2 / (4 me omega^2)) * d(E0^2)/dx
 *
 * The ponderomotive force pushes plasma from high-intensity to
 * low-intensity regions.  Key effects include:
 *   - Density profile steepening at the critical surface
 *   - Hole boring into overdense targets
 *   - Self-focusing channel formation
 *   - Wakefield generation
 *
 * Source: Kruer Sec 6.2
 *
 * Complexity: O(1)
 * ============================================================ */
double ponderomotive_force(double ne, double intensity,
                           double lambda_m, double gradient_scale_length)
{
    if (ne <= 0.0 || intensity <= 0.0 || lambda_m <= 0.0
        || gradient_scale_length <= 0.0)
        return 0.0;

    double phi_p = ponderomotive_potential(intensity, lambda_m);
    return ne * phi_p / gradient_scale_length;
}

/* ============================================================
 *  L1: quiver_velocity
 *
 *  v_osc = e E0 / (me omega) = a0 c
 *
 * The peak velocity of an electron oscillating in the laser
 * electric field.  For a0 = 1, v_osc = c on the forward
 * half-cycle.
 *
 * This is the key velocity parameter for analyzing parametric
 * instabilities where the oscillatory motion beats with
 * density perturbations.
 *
 * Complexity: O(1)
 * ============================================================ */
double quiver_velocity(double a0)
{
    return a0 * c;
}

/* ============================================================
 *  L1: quiver_energy
 *
 *  U_p = e^2 E0^2 / (4 me omega^2) = Phi_p
 *
 * The ponderomotive (quiver) energy equals the ponderomotive
 * potential.  It determines:
 *   - Above-Threshold Ionization (ATI) energy cutoff: N omega ~ 10 U_p
 *   - High-Harmonic Generation (HHG) cutoff: N omega ~ Ip + 3.17 U_p
 *   - Transition from MPI to tunnel ionization (gamma_K ~ 1)
 *
 * Source: Gibbon Sec 2.1
 *
 * Complexity: O(1)
 * ============================================================ */
double quiver_energy(double intensity, double lambda_m)
{
    return ponderomotive_potential(intensity, lambda_m);
}

/* ============================================================
 *  L2: relativistic_self_focusing_threshold
 *
 *  P_c = 17.3 (nc / ne) GW
 *
 * When the laser power exceeds P_c, relativistic mass increase
 * of electrons on axis produces a positive refractive-index
 * step that acts as a converging lens, leading to self-focusing.
 *
 * This is a whole-beam instability.  The beam collapses until
 * some saturation mechanism intervenes (e.g., ponderomotive
 * density depletion).
 *
 * Source: Gibbon Sec 5.1, Eq. 5.1
 * Theorem: From nonlinear Schrodinger equation with Kerr-type
 *          nonlinearity in plasma: P_c ~ (c me^2 c^3 / e^2) (omega^2/wp^2)
 *
 * Complexity: O(1)
 * ============================================================ */
double relativistic_self_focusing_threshold(double ne, double nc)
{
    if (ne <= 0.0 || nc <= 0.0) return DBL_MAX;
    double P_c_W = 17.3e9 * (nc / ne);
    return P_c_W;
}

/* ============================================================
 *  L2: ponderomotive_self_focusing_threshold
 *
 *  P_c ~ 1.2e17 (nc/ne)^{3/2} lambda_mum^2  [W]
 *
 * Ponderomotive expulsion of electrons from the laser axis
 * creates a density channel.  When the laser power exceeds
 * P_c, the channel depth is sufficient to guide the beam.
 *
 * Source: Gibbon Sec 5.1
 *
 * Complexity: O(1)
 * ============================================================ */
double ponderomotive_self_focusing_threshold(double ne, double nc,
                                              double lambda_m)
{
    if (ne <= 0.0 || nc <= 0.0 || lambda_m <= 0.0) return DBL_MAX;
    double lambda_um = lambda_m * 1e6;
    double ratio = nc / ne;
    return 1.2e17 * pow(ratio, 1.5) * lambda_um * lambda_um;
}

/* ============================================================
 *  L2: relativistic_transparency_density
 *
 *  n_c^rel = gamma n_c = sqrt(1 + a0^2) n_c
 *
 * Electrons driven to relativistic velocities have increased
 * effective mass m_eff = gamma m_e, which reduces the plasma
 * frequency to omega_p^eff = omega_p / sqrt(gamma).
 *
 * Consequently, a laser with a0 > 1 can propagate through
 * plasma that is formally overdense (ne > n_c).  This is
 * the physical basis for:
 *   - Relativistic self-induced transparency (RSIT)
 *   - Hole boring with relativistic correction
 *   - Laser-driven ion acceleration (RPA, LS)
 *
 * Source: Gibbon Sec 5.3, Eq. 5.22
 *
 * Complexity: O(1)
 * ============================================================ */
double relativistic_transparency_density(double nc, double a0)
{
    if (nc <= 0.0 || a0 < 0.0) return 0.0;
    double gamma = sqrt(1.0 + a0 * a0);
    return gamma * nc;
}

/* ============================================================
 *  L2: filamentation_growth_rate
 *
 *  Gamma_fil = (omega_p^2 / (8 omega)) * a0^2
 *
 * Transverse intensity modulations grow exponentially at this
 * rate due to the relativistic/ponderomotive nonlinearity.
 * Filamentation breaks up the laser beam into small-scale
 * "hot spots" that can enhance local absorption.
 *
 * Source: Kruer Sec 8.2; Max (1976)
 *
 * Complexity: O(1)
 * ============================================================ */
double filamentation_growth_rate(double wp, double omega, double a0)
{
    if (omega <= 0.0 || wp <= 0.0) return 0.0;
    return (wp * wp) / (8.0 * omega) * a0 * a0;
}

/* ============================================================
 *  L2: hole_boring_velocity -- Laser piston in overdense plasma
 *
 *  v_hb / c = sqrt(I / (2 mi ni c^3))
 *
 * The laser radiation pressure acts as a piston, pushing the
 * critical surface into the overdense plasma at velocity v_hb.
 *
 * For a solid-density target (ni ~ 10^23 cm^-3, A=12):
 *   I = 10^20 W/cm^2 => v_hb/c ~ 0.015 (~ 4500 km/s)
 *
 * This is fundamental to:
 *   - Laser-driven shock ignition
 *   - Hole-boring RPA ion acceleration
 *   - Fast ignition cone compression
 *
 * Source: Gibbon Sec 6.1, Eq. 6.3
 *
 * Complexity: O(1)
 * ============================================================ */
double hole_boring_velocity(double intensity, double mi, double ni)
{
    if (intensity <= 0.0 || mi <= 0.0 || ni <= 0.0) return 0.0;
    double v_over_c_sq = intensity / (2.0 * mi * ni * c * c * c);
    if (v_over_c_sq >= 1.0) return c;
    return c * sqrt(v_over_c_sq);
}

/* ============================================================
 *  L2: plasma_mirror_reflectivity
 *
 *  R = (1 - N)^2 / (1 + N)^2  where N = sqrt(1 - nc/ne)
 *
 * Fresnel formula for reflection at a sharp plasma-vacuum
 * interface.  For ne >> nc:
 *   R ~ 1 - 2 sqrt(nc/ne)  -> approaches 1
 *
 * For ne = nc: N = 0, R = 1 (perfect reflection)
 * For ne = 2 nc: N ~ 0.707, R ~ 0.03 (weak reflection)
 *
 * This formula ignores profile effects (resonance absorption)
 * and collisions.
 *
 * Complexity: O(1)
 * ============================================================ */
double plasma_mirror_reflectivity(double ne, double nc)
{
    if (nc <= 0.0) return 0.0;
    if (ne <= nc) return 1.0; /* totally reflecting if overdense */

    /* For ne > nc: underdense side sees step-down in density */
    double N = sqrt(1.0 - nc / ne);
    double temp = (1.0 - N) / (1.0 + N);
    return temp * temp;
}

/* ============================================================
 *  L5: ray_deflection_angle
 *
 *  dtheta/ds ~ (1/2nc) * grad_perp(ne)
 *
 * Using the eikonal (ray-optics) approximation for an EM wave
 * in a slowly varying plasma density.
 *
 * Source: Kruer Sec 6.1
 *
 * Complexity: O(1)
 * ============================================================ */
double ray_deflection_angle(double ne_gradient, double nc,
                            double path_length)
{
    if (nc <= 0.0) return 0.0;
    return fabs(ne_gradient) * path_length / (2.0 * nc);
}

/* ============================================================
 *  L2: density_scale_length
 *
 *  L_n = ne / |grad ne|
 *
 * The exponential density gradient scale length is a critical
 * parameter for:
 *   - Resonance absorption efficiency (phi depends on L_n)
 *   - Parametric instability thresholds (gain ~ L_n)
 *   - Laser propagation (validity of WKB approximation)
 *
 * Large L_n => gentle gradient => weak interaction
 * Small L_n => steep gradient => strong interaction
 *
 * Complexity: O(1)
 * ============================================================ */
double density_scale_length(double ne, double ne_gradient)
{
    if (ne <= 0.0 || fabs(ne_gradient) < 1e-30) return DBL_MAX;
    return ne / fabs(ne_gradient);
}

/* ============================================================
 *  L5: refraction_index_profile -- N along linear density ramp
 *
 *  ne(z) = ne_start + (ne_end - ne_start) * z / L
 *  N(z)  = sqrt(1 - ne(z)/nc)
 *
 * Used for computing the phase accumulated by a ray traversing
 * a density gradient.
 *
 * Complexity: O(1)
 * ============================================================ */
double refraction_index_profile(double z, double ne_start, double ne_end,
                                double nc, double L)
{
    if (L <= 0.0 || nc <= 0.0) return 1.0;
    double frac, ne_z;

    if (z <= 0.0) {
        ne_z = ne_start;
    } else if (z >= L) {
        ne_z = ne_end;
    } else {
        frac = z / L;
        ne_z = ne_start + (ne_end - ne_start) * frac;
    }

    if (ne_z < 0.0) ne_z = 0.0;
    if (ne_z >= nc) return 0.0;

    return sqrt(1.0 - ne_z / nc);
}

/* ============================================================
 *  L5: snells_law_plasma -- Snell's law for stratified plasma
 *
 *  sin(theta(z)) / sin(theta_0) = N_0 / N(z)
 *
 * A ray incident at angle theta_0 on a plasma with surface
 * density ne0 is refracted as it travels to regions of
 * different density.
 *
 * This is critical for understanding:
 *   - Laser ray trajectories in ICF hohlraums
 *   - Resonance absorption angle optimization
 *   - Density profile reconstruction from refractometry
 *
 * Complexity: O(1)
 * ============================================================ */
double snells_law_plasma(double theta0, double ne0, double ne_z, double nc)
{
    if (nc <= 0.0) return theta0;
    if (ne0 >= nc) return M_PI / 2.0;  /* incident from overdense */

    double N0 = sqrt(1.0 - ne0 / nc);
    double Nz = sqrt(1.0 - (ne_z < nc ? ne_z / nc : 0.999));

    if (Nz <= 0.0) return M_PI / 2.0;  /* turning point */

    double sin_theta_z = N0 / Nz * sin(theta0);

    if (sin_theta_z >= 1.0) return M_PI / 2.0;

    return asin(sin_theta_z);
}

/* ============================================================
 *  L4: turning_point_density
 *
 *  ne_turn = nc cos^2(theta_0)
 *
 * An obliquely incident ray with vacuum angle theta_0 on
 * a plasma with ne=0 at the boundary will turn around when
 * the local density reaches ne_turn.
 *
 * For theta_0 = 45 degrees: ne_turn = nc/2
 * For theta_0 = 60 degrees: ne_turn = nc/4 (quarter-critical!)
 *
 * The quarter-critical surface (nc/4) is special because
 * it supports Two-Plasmon Decay (omega0 = 2 omega_p).
 *
 * Complexity: O(1)
 * ============================================================ */
double turning_point_density(double theta0, double nc)
{
    if (nc <= 0.0) return 0.0;
    double cos_theta = cos(theta0);
    return nc * cos_theta * cos_theta;
}

/* ============================================================
 *  L5: ray_trace_step -- Forward integration step
 *
 * Advances (r, k_hat) by ds along the ray trajectory.
 *
 *   r_new     = r + k_hat * ds
 *   k_hat_new = k_hat + dk_hat/ds * ds
 *   where dk_hat/ds = grad_perp(ne) / (2 nc)
 *
 * This implements the ray equation derived from the eikonal
 * approximation to the Helmholtz equation.
 *
 * Complexity: O(1)
 * ============================================================ */
void ray_trace_step(double *r, double *k_hat,
                    double ne_grad_x, double ne_grad_y,
                    double nc, double ds)
{
    if (!r || !k_hat || nc <= 0.0 || ds <= 0.0) return;

    /* Position update */
    r[0] += k_hat[0] * ds;
    r[1] += k_hat[1] * ds;
    r[2] += k_hat[2] * ds;

    /* Direction update from density gradient */
    double dk_ds_x = ne_grad_x / (2.0 * nc);
    double dk_ds_y = ne_grad_y / (2.0 * nc);

    k_hat[0] += dk_ds_x * ds;
    k_hat[1] += dk_ds_y * ds;

    /* Normalize k_hat (maintain unit vector) */
    double norm = sqrt(k_hat[0]*k_hat[0] + k_hat[1]*k_hat[1]
                       + k_hat[2]*k_hat[2]);
    if (norm > 0.0) {
        k_hat[0] /= norm;
        k_hat[1] /= norm;
        k_hat[2] /= norm;
    }
}

/* ============================================================
 *  L5: nonlinear_phase_shift -- Relativistic nonlinear phase
 *
 *  dphi/dz = (wp^2 / (2 k0 c^2)) * (1 - 1/sqrt(1 + a^2))
 *
 * Phase accumulated per propagation distance due to the
 * relativistic electron mass increase.  This is the "Kerr"
 * nonlinearity of a plasma.
 *
 * Source: Gibbon Sec 5.1
 *
 * Complexity: O(1)
 * ============================================================ */
double nonlinear_phase_shift(double a_sq, double wp, double k0)
{
    if (k0 <= 0.0 || wp <= 0.0) return 0.0;
    if (a_sq < 0.0) a_sq = 0.0;

    double gamma_inv = 1.0 / sqrt(1.0 + a_sq);
    return (wp * wp) / (2.0 * k0 * c * c) * (1.0 - gamma_inv);
}

/* ============================================================
 *  L5: paraxial_state_alloc -- Allocate propagation grid
 *
 * Allocates a ParaxialState structure and its 2D complex
 * field arrays on a grid of Nx x Ny points.
 *
 * Complexity: O(Nx * Ny)
 * ============================================================ */
ParaxialState *paraxial_state_alloc(int Nx, int Ny, double dx, double dy,
                                    double lambda0)
{
    if (Nx < 1 || Ny < 1 || dx <= 0.0 || dy <= 0.0 || lambda0 <= 0.0)
        return NULL;

    ParaxialState *ps = (ParaxialState *)malloc(sizeof(ParaxialState));
    if (!ps) return NULL;

    ps->Nx = Nx;
    ps->Ny = Ny;
    ps->dx = dx;
    ps->dy = dy;
    ps->z  = 0.0;
    ps->lambda0 = lambda0;
    ps->k0 = 2.0 * M_PI / lambda0;

    size_t size = (size_t)Nx * Ny * sizeof(double);
    ps->E_real = (double *)malloc(size);
    ps->E_imag = (double *)malloc(size);

    if (!ps->E_real || !ps->E_imag) {
        free(ps->E_real);
        free(ps->E_imag);
        free(ps);
        return NULL;
    }

    memset(ps->E_real, 0, size);
    memset(ps->E_imag, 0, size);

    return ps;
}

/* ============================================================
 *  L5: paraxial_state_free
 * ============================================================ */
void paraxial_state_free(ParaxialState *ps)
{
    if (!ps) return;
    free(ps->E_real);
    free(ps->E_imag);
    free(ps);
}

/* ============================================================
 *  L5: paraxial_step_vacuum -- Vacuum propagation (linear)
 *
 * Implements the split-step Fourier method:
 *   1. Half-step of free-space diffraction in Fourier domain
 *   2. Full-step of phase accumulation in real space
 *
 * Without FFTW, this implementation uses a simple finite-
 * difference Laplacian for the diffraction step.
 *
 *   dE/dz = i/(2k0) Laplacian_perp(E)
 *
 * Complexity: O(Nx * Ny) per step (FD form)
 * ============================================================ */
int paraxial_step_vacuum(ParaxialState *ps, double dz)
{
    if (!ps || dz <= 0.0) return -1;
    if (ps->Nx < 3 || ps->Ny < 3) return -1;

    int Nx = ps->Nx, Ny = ps->Ny;
    double dx = ps->dx, dy = ps->dy;
    double coeff = dz / (2.0 * ps->k0);

    /* Compute Laplacian of E in real space (5-point stencil) */
    double *lap_real = (double *)malloc((size_t)Nx * Ny * sizeof(double));
    double *lap_imag = (double *)malloc((size_t)Nx * Ny * sizeof(double));
    if (!lap_real || !lap_imag) {
        free(lap_real); free(lap_imag);
        return -1;
    }

    for (int i = 1; i < Nx - 1; i++) {
        for (int j = 1; j < Ny - 1; j++) {
            int idx = i * Ny + j;
            /* d^2/dx^2 */
            double d2r_dx2 = (ps->E_real[(i+1)*Ny+j] - 2.0*ps->E_real[idx]
                              + ps->E_real[(i-1)*Ny+j]) / (dx * dx);
            double d2i_dx2 = (ps->E_imag[(i+1)*Ny+j] - 2.0*ps->E_imag[idx]
                              + ps->E_imag[(i-1)*Ny+j]) / (dx * dx);
            /* d^2/dy^2 */
            double d2r_dy2 = (ps->E_real[i*Ny+j+1] - 2.0*ps->E_real[idx]
                              + ps->E_real[i*Ny+j-1]) / (dy * dy);
            double d2i_dy2 = (ps->E_imag[i*Ny+j+1] - 2.0*ps->E_imag[idx]
                              + ps->E_imag[i*Ny+j-1]) / (dy * dy);

            lap_real[idx] = d2r_dx2 + d2r_dy2;
            lap_imag[idx] = d2i_dx2 + d2i_dy2;
        }
    }

    /* Update E: dEr/dz = -dEi_laplacian/(2k0), dEi/dz = dEr_laplacian/(2k0) */
    for (int i = 1; i < Nx - 1; i++) {
        for (int j = 1; j < Ny - 1; j++) {
            int idx = i * Ny + j;
            ps->E_real[idx] -= coeff * lap_imag[idx];
            ps->E_imag[idx] += coeff * lap_real[idx];
        }
    }

    ps->z += dz;
    free(lap_real);
    free(lap_imag);
    return 0;
}

/* ============================================================
 *  L5: paraxial_step_plasma -- Propagation with plasma nonlinearity
 *
 * Adds the nonlinear phase shift from the relativistic plasma
 * response to the vacuum propagation step.
 *
 *   dE/dz = i/(2k0) Laplacian_perp(E) + i dphi/dz(E) E
 *
 * where dphi/dz is the relativistic nonlinear phase shift.
 *
 * Complexity: O(Nx * Ny) per step
 * ============================================================ */
int paraxial_step_plasma(ParaxialState *ps, double dz,
                          double ne, double nc)
{
    (void)nc;   /* nc unused in current implementation, reserved for future */
    if (!ps || dz <= 0.0) return -1;

    /* First apply vacuum diffraction step */
    int ret = paraxial_step_vacuum(ps, dz);
    if (ret != 0) return ret;

    /* Then apply nonlinear phase from plasma response */
    double wp = plasma_frequency(ne);
    for (int i = 0; i < ps->Nx; i++) {
        for (int j = 0; j < ps->Ny; j++) {
            int idx = i * ps->Ny + j;
            double Er = ps->E_real[idx];
            double Ei = ps->E_imag[idx];
            double a_sq = (Er*Er + Ei*Ei) * (e*e) / (me*me * c*c * c*c);
            double dphi = nonlinear_phase_shift(a_sq, wp, ps->k0) * dz;

            /* Rotate: E -> E * exp(i * dphi) */
            double cos_phi = cos(dphi);
            double sin_phi = sin(dphi);
            ps->E_real[idx] = Er * cos_phi - Ei * sin_phi;
            ps->E_imag[idx] = Er * sin_phi + Ei * cos_phi;
        }
    }

    return 0;
}

/* ============================================================
 *  L5: paraxial_run_to -- Propagate to target distance
 *
 * Repeatedly steps the beam forward until z >= z_target
 * (or maximum steps reached).
 *
 * Complexity: O(N_steps * Nx * Ny)
 * ============================================================ */
int paraxial_run_to(ParaxialState *ps, double z_target,
                    double ne, double nc, double dz)
{
    if (!ps || z_target <= ps->z || dz <= 0.0) return -1;

    int max_steps = 10000;
    int step_count = 0;

    while (ps->z < z_target && step_count < max_steps) {
        double remaining = z_target - ps->z;
        double actual_dz = (dz < remaining) ? dz : remaining;

        int ret;
        if (ne > 0.0 && nc > 0.0) {
            ret = paraxial_step_plasma(ps, actual_dz, ne, nc);
        } else {
            ret = paraxial_step_vacuum(ps, actual_dz);
        }

        if (ret != 0) return -1;
        step_count++;
    }

    return (step_count < max_steps) ? 0 : -1;
}
