/*
 * plasma_params.c -- Fundamental plasma parameter implementations
 *
 * Each function implements exactly one knowledge point from the
 * laser-plasma physics curriculum.
 *
 * References:
 *   - Goldston & Rutherford (1995) "Introduction to Plasma Physics"
 *   - Chen (2016) "Introduction to Plasma Physics and Controlled Fusion"
 *   - Kruer (1988) "The Physics of Laser Plasma Interactions"
 *
 * Knowledge Layers: L1, L2, L3, L4
 */

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include "plasma_constants.h"
#include "plasma_params.h"

/* ============================================================
 *  Physical constants (aliased for readability)
 * ============================================================ */
static const double c   = PLASMA_C;
static const double e   = PLASMA_E;
static const double me  = PLASMA_ME;
static const double mp  = PLASMA_MP;
static const double eps0 = PLASMA_EPS0;
static const double kB  = PLASMA_KB;
static const double hbar = PLASMA_HBAR;
static const double mu0 = PLASMA_MU0;

/* ============================================================
 *  L1: plasma_frequency -- Electron plasma (Langmuir) frequency
 *
 *  omega_p = sqrt(ne e^2 / (eps0 me))
 *
 * This is the most fundamental frequency in plasma physics.
 * All collective electron phenomena occur on the time scale
 * omega_p^{-1}.
 *
 * For ne = 10^21 cm^-3 (typical laser-produced plasma):
 *   omega_p ~ 1.78e15 rad/s  (period ~ 3.5 fs)
 *
 * Source: Goldston & Rutherford Sec 1.3, Eq. 1.16
 *
 * Complexity: O(1)   Theorem: linearised cold-fluid oscillation
 * ============================================================ */
double plasma_frequency(double ne)
{
    if (ne <= 0.0) return 0.0;
    return sqrt(ne * e * e / (eps0 * me));
}

/* ============================================================
 *  L1: ion_plasma_frequency -- Ion plasma frequency
 *
 *  omega_pi = sqrt(ni Z^2 e^2 / (eps0 mi))
 *           = omega_pe * sqrt(Z me / (A mp))
 *
 * Because mp/me ~ 1836, ion motions are ~43x slower than
 * electron motions even for Z ~ 1.  This decoupling is key
 * to SBS (ion timescale instability vs SRS electron timescale).
 *
 * Source: Chen Sec 4.12
 *
 * Complexity: O(1)
 * ============================================================ */
double ion_plasma_frequency(double ne, double Z, double A)
{
    if (ne <= 0.0 || Z <= 0.0 || A <= 0.0) return 0.0;
    /* ni = ne/Z for quasi-neutrality */
    double ni = ne / Z;
    double mi = A * mp;
    return sqrt(ni * Z * Z * e * e / (eps0 * mi));
}

/* ============================================================
 *  L1: critical_density -- Critical electron density
 *
 *  nc = eps0 me omega^2 / e^2
 *     = 1.115e15 / (lambda_mum)^2  [cm^-3]
 *
 * For a 1 um laser: nc ~ 1.1e21 cm^-3
 * For a 0.35 um (3omega Nd:glass): nc ~ 9.0e21 cm^-3
 * For CO2 (10.6 um): nc ~ 1.0e19 cm^-3
 *
 * EM wave propagation requires ne < nc (underdense plasma).
 *
 * Source: Kruer Sec 1.2, Eq. 1.5
 * Theorem: EM dispersion omega^2 = wp^2 + c^2 k^2,
 *          cutoff occurs at k -> 0 => omega = wp
 *
 * Complexity: O(1)
 * ============================================================ */
double critical_density(double lambda_m)
{
    if (lambda_m <= 0.0) return 0.0;
    double omega = 2.0 * M_PI * c / lambda_m;
    return eps0 * me * omega * omega / (e * e);
}

/* ============================================================
 *  L1: debye_length -- Electron Debye screening length
 *
 *  lambda_D = sqrt(eps0 kB Te / (ne e^2))
 *
 * For a laser-produced plasma with ne = 10^21 cm^-3,
 * Te = 1 keV: lambda_D ~ 7.4 nm
 *
 * The Debye length sets the minimum spatial scale for
 * collective (fluid) plasma descriptions.  Structures
 * smaller than lambda_D require kinetic treatment.
 *
 * Source: Goldston & Rutherford Sec 1.4, Eq. 1.22
 * Theorem: From linearised Poisson-Boltzmann equation
 *          with Maxwell-Boltzmann electron distribution
 *
 * Complexity: O(1)
 * ============================================================ */
double debye_length(double ne, double Te_eV)
{
    if (ne <= 0.0 || Te_eV <= 0.0) return 0.0;
    /* lambda_D = sqrt(eps0 * k_B * T[K] / (ne * e^2))
     * with T[K] = Te_eV * e / k_B:
     * lambda_D = sqrt(eps0 * Te_eV * e / (ne * e^2))
     *          = sqrt(eps0 * Te_eV / (ne * e))
     */
    return sqrt(eps0 * Te_eV / (ne * e));
}

/* ============================================================
 *  L1: debye_sphere_particles -- Plasma parameter N_D
 *
 *  N_D = (4pi/3) ne lambda_D^3
 *
 * This is the "plasma parameter" that determines whether
 * the system exhibits collective behaviour.
 *
 * N_D >> 1    : collective plasma (weakly coupled)
 * N_D ~ 1     : strongly coupled plasma
 * N_D < 1     : not a plasma (gas of independent particles)
 *
 * For ICF plasmas: N_D ~ 10^2 - 10^5
 *
 * Source: Chen Sec 1.3, Eq. 1.6
 *
 * Complexity: O(1)
 * ============================================================ */
double debye_sphere_particles(double ne, double lambda_D)
{
    if (ne <= 0.0 || lambda_D <= 0.0) return 0.0;
    return (4.0 * M_PI / 3.0) * ne * lambda_D * lambda_D * lambda_D;
}

/* ============================================================
 *  L1: electron_thermal_velocity
 *
 *  v_the = sqrt(kB Te / me) = sqrt(e Te_eV / me)
 *
 * For Te = 1 keV: v_the ~ 1.33e7 m/s ~ 0.044 c
 * For Te = 100 eV: v_the ~ 4.2e6 m/s
 *
 * This velocity appears in:
 *   - Landau damping rate
 *   - Collision frequency scaling
 *   - Bohm-Gross dispersion
 *
 * Complexity: O(1)
 * ============================================================ */
double electron_thermal_velocity(double Te_eV)
{
    if (Te_eV <= 0.0) return 0.0;
    /* v_the = sqrt(kB T/me) with kB T = e * Te_eV */
    return sqrt(e * Te_eV / me);
}

/* ============================================================
 *  L1: ion_thermal_velocity
 *
 *  v_thi = sqrt(kB Ti / mi) = sqrt(e Ti_eV / (A mp))
 *
 * For A=12 (carbon), Ti=1 keV: v_thi ~ 9.0e4 m/s
 * Much slower than electron thermal velocity by ~ sqrt(me/mi).
 *
 * Complexity: O(1)
 * ============================================================ */
double ion_thermal_velocity(double Ti_eV, double A)
{
    if (Ti_eV <= 0.0 || A <= 0.0) return 0.0;
    /* v_thi = sqrt(kB T/mi) with kB T = e * Ti_eV */
    double mi = A * mp;
    return sqrt(e * Ti_eV / mi);
}

/* ============================================================
 *  L1: ion_sound_speed -- Bohm speed
 *
 *  cs = sqrt(Z kB Te / mi) = sqrt(Z e Te_eV / (A mp))
 *
 * The speed of ion-acoustic perturbations in a plasma with
 * hot electrons and cold ions.  For Te=1 keV, Z=6 (carbon),
 * A=12: cs ~ 2.2e5 m/s.
 *
 * This speed appears in:
 *   - Bohm sheath criterion: v_i >= cs at sheath edge
 *   - SBS frequency shift: delta_omega ~ 2 omega0 cs/c
 *   - Ion front expansion velocity
 *
 * Source: Chen Sec 4.6, Eq. 4.49
 *
 * Complexity: O(1)
 * ============================================================ */
double ion_sound_speed(double Te_eV, double Z, double A)
{
    if (Te_eV <= 0.0 || Z <= 0.0 || A <= 0.0) return 0.0;
    /* cs = sqrt(Z kB Te/mi) with kB Te = e * Te_eV */
    double mi = A * mp;
    return sqrt(Z * e * Te_eV / mi);
}

/* ============================================================
 *  L1: skin_depth -- Collisionless EM skin depth
 *
 *  delta = c / omega_p
 *
 * For ne = 10^21 cm^-3: delta ~ 0.17 um
 * For ne = 10^19 cm^-3: delta ~ 1.7 um
 *
 * An EM field of frequency omega << omega_p penetrates
 * into an overdense plasma only to a depth ~ delta.
 * This sets the absorption layer thickness.
 *
 * Source: Goldston & Rutherford Sec 4.2
 *
 * Complexity: O(1)
 * ============================================================ */
double skin_depth(double wp)
{
    if (wp <= 0.0) return DBL_MAX;
    return c / wp;
}

/* ============================================================
 *  L2: coulomb_logarithm -- ln Lambda
 *
 *  ln Lambda = max(ln(lambda_D / b_min), 2)
 *
 * where lambda_D is the Debye length (maximum impact parameter
 * for a single collision) and b_min is the minimum impact
 * parameter (90-degree deflection distance).
 *
 * For weakly coupled plasmas (N_D >> 1): ln Lambda ~ 5-20.
 * The Coulomb logarithm accounts for the collective shielding
 * of individual charges in a plasma.
 *
 * Source: Goldston & Rutherford Sec 2.3, Eq. 2.40
 *
 * Complexity: O(1)
 * ============================================================ */
double coulomb_logarithm(double ne, double Te_eV, double Z)
{
    if (ne <= 0.0 || Te_eV <= 0.0 || Z <= 0.0) return 2.0;

    double ld = debye_length(ne, Te_eV);
    double vth = electron_thermal_velocity(Te_eV);

    /* Minimum impact parameter: classical 90-degree deflection */
    double b_min_cl = (Z * e * e) / (4.0 * M_PI * eps0 * me * vth * vth);

    /* Quantum correction: de Broglie wavelength */
    double b_min_qm = hbar / (2.0 * me * vth);

    /* Use the larger of classical and quantum b_min */
    double b_min = (b_min_cl > b_min_qm) ? b_min_cl : b_min_qm;

    /* Guard against b_min >= ld (strongly coupled plasma) */
    if (b_min >= ld || b_min <= 0.0) return 2.0;

    double lnL = log(ld / b_min);
    /* The Coulomb logarithm is never less than 2 for a plasma */
    return (lnL > 2.0) ? lnL : 2.0;
}

/* ============================================================
 *  L2: electron_ion_collision_frequency -- Spitzer formula
 *
 *  nu_ei = [4 sqrt(2) / 3] Z ni e^4 ln Lambda
 *          / [sqrt(pi) me^{1/2} (4 pi eps0)^2 (kB Te)^{3/2}]
 *
 * This is the rate at which a drifting electron population
 * transfers momentum to ions via Coulomb collisions.
 *
 * For ne = 10^21 cm^-3, Te = 1 keV, Z = 6, lnL = 10:
 *   nu_ei ~ 3.5e12 s^-1  (mean free time ~ 0.3 ps)
 *
 * Source: Goldston & Rutherford Sec 2.7, Eq. 2.73
 * Also known as: Braginskii electron-ion collision frequency
 *
 * Complexity: O(1)
 * ============================================================ */
double electron_ion_collision_frequency(double ne, double Te_eV,
                                         double Z, double ln_Lambda)
{
    if (ne <= 0.0 || Te_eV <= 0.0 || Z <= 0.0 || ln_Lambda < 2.0)
        return 0.0;

    double Te_J = Te_eV * e;
    double ni = ne / Z;  /* quasi-neutrality */

    double prefactor = (4.0 * sqrt(2.0) / 3.0) * Z * ni
                       * pow(e, 4.0) * ln_Lambda;

    double denominator = sqrt(M_PI) * sqrt(me)
                         * pow(4.0 * M_PI * eps0, 2.0)
                         * pow(kB * Te_J, 1.5);

    return prefactor / denominator;
}

/* ============================================================
 *  L2: collision_mean_free_path
 *
 *  lambda_mfp = v_the / nu_ei
 *
 * Average distance an electron travels between significant
 * (90-degree cumulative) momentum-redirection events.
 *
 * For the Spitzer example above: lambda_mfp ~ 4 um
 * In ICF plasmas, this can be comparable to or larger than
 * the laser wavelength, making IB absorption weak for
 * high temperatures.
 *
 * Complexity: O(1)
 * ============================================================ */
double collision_mean_free_path(double v_the, double nu_ei)
{
    if (v_the <= 0.0 || nu_ei <= 0.0) return DBL_MAX;
    return v_the / nu_ei;
}

/* ============================================================
 *  L2: spitzer_resistivity -- DC resistivity
 *
 *  eta_parallel = pi Z e^2 sqrt(me) ln Lambda
 *                 / [(4 pi eps0)^2 (kB Te)^{3/2}]
 *
 * Parallel DC resistivity in a fully ionized, unmagnetized
 * plasma.  Scales as Te^{-3/2}, so hotter plasmas are
 * better conductors (important for fast ignition concepts).
 *
 * Source: Goldston & Rutherford Sec 2.9, Eq. 2.82
 *
 * Complexity: O(1)
 * ============================================================ */
double spitzer_resistivity(double Te_eV, double Z, double ln_Lambda)
{
    if (Te_eV <= 0.0 || Z <= 0.0 || ln_Lambda < 2.0) return 0.0;

    double Te_J = Te_eV * e;
    double numerator = M_PI * Z * e * e * sqrt(me) * ln_Lambda;
    double denominator = pow(4.0 * M_PI * eps0, 2.0)
                         * pow(kB * Te_J, 1.5);

    return numerator / denominator;
}

/* ============================================================
 *  L3: cold_plasma_permittivity
 *
 *  epsilon(omega) = 1 - (omega_p / omega)^2
 *
 * This is the real part of the unmagnetized cold plasma
 * dielectric function.  It gives the refractive index
 * N = sqrt(epsilon).
 *
 * omega > omega_p : epsilon in (0, 1) -> propagation (N real)
 * omega < omega_p : epsilon < 0      -> evanescence (N imag)
 * omega = omega_p : epsilon = 0      -> cutoff (N=0)
 *
 * Source: Chen Sec 4.3, Eq. 4.24
 * Theorem: From cold fluid + Maxwell: nabla x B = mu0 J + c^{-2} dE/dt
 *
 * Complexity: O(1)
 * ============================================================ */
double cold_plasma_permittivity(double omega, double wp)
{
    if (omega <= 0.0) return -DBL_MAX;
    double ratio = wp / omega;
    return 1.0 - ratio * ratio;
}

/* ============================================================
 *  L3: collisional_permittivity -- Drude model
 *
 *  epsilon(omega) = 1 - wp^2/(omega^2 + nu^2)
 *                        + i (nu/omega) wp^2/(omega^2 + nu^2)
 *
 * Returns Re(epsilon); imaginary part via *eps_imag.
 *
 * This is the key function for computing inverse bremsstrahlung
 * absorption coefficient:
 *   k_ib = (omega/c) Im(sqrt(epsilon))
 *        ~ (nu_ei/c) (ne/nc) / sqrt(1 - ne/nc)
 *
 * Source: Kruer Sec 3.2
 *
 * Complexity: O(1)
 * ============================================================ */
double collisional_permittivity(double omega, double wp, double nu,
                                double *eps_imag)
{
    if (omega <= 0.0) {
        if (eps_imag) *eps_imag = 0.0;
        return -DBL_MAX;
    }

    double denom = omega * omega + nu * nu;
    double delta = wp * wp / denom;

    double eps_real = 1.0 - delta;
    double eps_imag_val = (nu / omega) * delta;

    if (eps_imag) *eps_imag = eps_imag_val;
    return eps_real;
}

/* ============================================================
 *  L3: warm_plasma_permittivity -- Bohm-Gross dielectric
 *
 *  epsilon(k, omega) = 1 - wp^2 / (omega^2 - 3 k^2 v_the^2)
 *
 * Adding the thermal pressure term (3 k^2 v_the^2) transforms
 * the cold plasma oscillation into a propagating wave:
 *
 *   omega^2 = wp^2 + 3 k^2 v_the^2   (Bohm-Gross dispersion)
 *
 * The factor 3 comes from the adiabatic index gamma=3 for
 * 1D electron compression.  In 3D, gamma=5/3 would be used.
 *
 * Source: Chen Sec 4.10, Eq. 4.67
 * Theorem: Linearised warm fluid equations with adiabatic closure
 *
 * Complexity: O(1)
 * ============================================================ */
double warm_plasma_permittivity(double omega, double k, double wp,
                                double v_the)
{
    double thermal_term = 3.0 * k * k * v_the * v_the;
    double denom = omega * omega - thermal_term;

    /* Guard against resonance pole */
    if (fabs(denom) < 1e-30) {
        return (denom >= 0.0) ? DBL_MAX : -DBL_MAX;
    }

    return 1.0 - (wp * wp) / denom;
}

/* ============================================================
 *  L3: plasma_refractive_index
 *
 *  N = sqrt(1 - ne/nc)  for ne < nc
 *
 * Phase velocity v_phi = c / N (superluminal)
 * Group velocity  v_g = c * N  (subluminal)
 *
 * At the critical surface (ne = nc): N=0, causing:
 *   - Reflection of EM wave
 *   - Singular electric field (in cold theory)
 *   - Linear mode conversion to plasma waves (warm theory)
 *
 * Complexity: O(1)
 * ============================================================ */
double plasma_refractive_index(double ne, double nc)
{
    if (nc <= 0.0) return 1.0;
    if (ne >= nc) return 0.0;
    if (ne < 0.0) return 1.0;
    return sqrt(1.0 - ne / nc);
}

/* ============================================================
 *  L3: laser_wavelength_in_plasma
 *
 *  lambda_plasma = lambda_vac / N
 *
 * As the laser approaches the critical surface (ne -> nc),
 * N -> 0 and lambda_plasma -> infinity.  This wavelength
 * stretching is a key observable (e.g., in interferometry
 * diagnostics).
 *
 * Complexity: O(1)
 * ============================================================ */
double laser_wavelength_in_plasma(double lambda_vac, double ne, double nc)
{
    double N = plasma_refractive_index(ne, nc);
    if (N <= 0.0) return DBL_MAX;
    return lambda_vac / N;
}

/* ============================================================
 *  L3: group_velocity -- EM wave group velocity in plasma
 *
 *  v_g / c = sqrt(1 - ne/nc)
 *
 * The speed at which laser energy propagates.  This is
 * always less than c.  As ne -> nc, v_g -> 0, causing
 * "group velocity slowing" -- important for Raman
 * backscatter where slow scattered light has more time
 * to interact.
 *
 * Complexity: O(1)
 * ============================================================ */
double group_velocity(double ne, double nc)
{
    if (nc <= 0.0) return c;
    if (ne >= nc) return 0.0;
    if (ne < 0.0) return c;
    return c * sqrt(1.0 - ne / nc);
}

/* ============================================================
 *  L4: plasma_beta -- Ratio of thermal to magnetic pressure
 *
 *  beta = 2 mu0 n kB T / B^2
 *
 * beta << 1 : magnetically dominated (fusion confinement)
 * beta >> 1 : kinetically dominated (laser plasmas)
 * beta ~ 1  : equipartition
 *
 * For a laser-produced plasma (B=0 or weak self-generated B),
 * beta is formally infinite => kinetic regime.
 *
 * Source: Chen Sec 6.4
 *
 * Complexity: O(1)
 * ============================================================ */
double plasma_beta(double ne, double Te_eV, double Ti_eV, double B_field)
{
    if (B_field <= 0.0) return DBL_MAX;

    double Te_J = Te_eV * e;
    double Ti_J = Ti_eV * e;
    double p_thermal = ne * kB * (Te_J + Ti_J);
    double p_magnetic = B_field * B_field / (2.0 * mu0);

    if (p_magnetic <= 0.0) return DBL_MAX;
    return p_thermal / p_magnetic;
}

/* ============================================================
 *  L4: electron_pressure -- Ideal gas equation of state
 *
 *  p_e = n_e k_B T_e
 *
 * The simplest EOS for electrons in a plasma.  Valid when
 * degeneracy pressure is negligible (n_e << n_quantum where
 * n_quantum ~ (m_e kB Te/hbar^2)^{3/2}).
 *
 * Complexity: O(1)
 * ============================================================ */
double electron_pressure(double ne, double Te_eV)
{
    double Te_J = Te_eV * e;
    return ne * kB * Te_J;
}

/* ============================================================
 *  L2: radiation_pressure -- Laser radiation pressure
 *
 *  P_rad = I/c      for absorbing surface
 *  P_rad = 2 I/c    for perfectly reflecting surface
 *
 * For I = 10^19 W/cm^2: P_rad ~ 3.3e8 Pa = 3.3 Mbar
 * (absorbing case), sufficient to drive strong shocks.
 *
 * Complexity: O(1)
 * ============================================================ */
double radiation_pressure(double intensity, int reflecting)
{
    if (intensity <= 0.0) return 0.0;
    double p_absorbing = intensity / c;
    return reflecting ? 2.0 * p_absorbing : p_absorbing;
}

/* ============================================================
 *  L4: compute_all_derived -- Fill PlasmaDerived from PlasmaState
 *
 * One-shot function to compute all derived plasma parameters.
 * This is the "entry function" for most plasma calculations:
 * a user creates a PlasmaState, calls this, and gets back all
 * the physical quantities needed for further analysis.
 *
 * Returns 0 on success, -1 if input parameters are unphysical.
 *
 * Complexity: O(1)
 * ============================================================ */
int compute_all_derived(const PlasmaState *ps, PlasmaDerived *pd)
{
    if (!ps || !pd) return -1;
    if (ps->ne <= 0.0 || ps->Te <= 0.0 || ps->lam_laser <= 0.0)
        return -1;

    pd->omega_laser = 2.0 * M_PI * c / ps->lam_laser;
    pd->wp          = plasma_frequency(ps->ne);
    pd->wpi         = ion_plasma_frequency(ps->ne, ps->Z, ps->A);
    pd->nc          = critical_density(ps->lam_laser);
    pd->lambda_D    = debye_length(ps->ne, ps->Te);
    pd->N_D         = debye_sphere_particles(ps->ne, pd->lambda_D);
    pd->v_the       = electron_thermal_velocity(ps->Te);
    pd->v_thi       = ion_thermal_velocity(ps->Ti, ps->A);
    pd->cs          = ion_sound_speed(ps->Te, ps->Z, ps->A);
    pd->ln_Lambda   = coulomb_logarithm(ps->ne, ps->Te, ps->Z);
    pd->nu_ei       = electron_ion_collision_frequency(
                           ps->ne, ps->Te, ps->Z, pd->ln_Lambda);
    pd->skin_depth  = skin_depth(pd->wp);
    pd->beta        = plasma_beta(ps->ne, ps->Te, ps->Ti, ps->B_field);
    pd->ne_over_nc  = ps->ne / pd->nc;
    pd->a0          = 0.0;  /* computed by caller with laser intensity */
    pd->lambda_plasma = laser_wavelength_in_plasma(
                           ps->lam_laser, ps->ne, pd->nc);

    return 0;
}
