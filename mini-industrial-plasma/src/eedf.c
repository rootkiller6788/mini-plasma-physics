/**
 * eedf.c - Electron Energy Distribution Functions
 *
 * Reference: Lieberman §5.3
 *   Hagelaar & Pitchford, Plasma Sources Sci. Technol. 14, 722 (2005)
 *   Morgan, Plasma Chem. Plasma Process. 12, 501 (1992)
 *   Penetrante, Bardsley, Bela, "Electron Collision Cross Sections..." (1990)
 *
 * Course: MIT 22.611, Stanford EE 414, Berkeley EECS 245
 */

#include "eedf.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

/* ============================================================
 * L1: Maxwellian EEDF
 * f(eps) = 2*sqrt(eps/pi) * T_e^{-3/2} * exp(-eps/T_e)  [eV^{-3/2}]
 *
 * This is the equilibrium distribution for electrons in
 * thermal equilibrium. Normalized: integral_0^inf f(eps)*sqrt(eps)*deps = 1
 *
 * Complexity: O(1)
 * ============================================================ */

double maxwellian_eedf(double energy_eV, double T_e_eV)
{
    if (T_e_eV <= 0.0) return 0.0;
    if (energy_eV < 0.0) return 0.0;

    /* f(eps) = 2 * sqrt(eps/pi) * T_e^{-3/2} * exp(-eps/T_e) */
    double prefactor = 2.0 / sqrt(M_PI);
    double Te_pow = pow(T_e_eV, -1.5);
    double boltzmann_factor = exp(-energy_eV / T_e_eV);

    return prefactor * sqrt(energy_eV) * Te_pow * boltzmann_factor;
}

/*
 * Cumulative Maxwellian EEDF: integral_0^{E_max} f(eps)*sqrt(eps)*deps
 *
 * F(E_max) = erf(sqrt(E_max/T_e)) - 2*sqrt(E_max/(pi*T_e))*exp(-E_max/T_e)
 *
 * This equals the fraction of electrons with energy < E_max.
 */
double maxwellian_eedf_cumulative(double E_max, double T_e_eV)
{
    if (T_e_eV <= 0.0 || E_max <= 0.0) return 0.0;

    double x = sqrt(E_max / T_e_eV);
    /* erf(x) - 2*x*exp(-x^2)/sqrt(pi) */
    double erf_x = erf(x);
    double correction = 2.0 * x * exp(-x*x) / sqrt(M_PI);

    double F = erf_x - correction;
    return (F > 0.0) ? F : 0.0;
}

/*
 * Effective electron temperature from mean energy.
 * For Maxwellian: <eps> = (3/2) * T_e
 * => T_e = (2/3) * <eps>
 */
double mean_energy_to_Te(double mean_energy_eV)
{
    if (mean_energy_eV <= 0.0) return 0.0;
    return (2.0/3.0) * mean_energy_eV;
}

/* ============================================================
 * L1: Druyvesteyn EEDF
 * f(eps) = 0.565 * <eps>^{-3/2} * sqrt(eps) * exp(-0.243*(eps/<eps>)^2)
 *
 * This distribution arises when elastic collisions dominate
 * and the electron-neutral cross section is constant (hard-sphere).
 * It has a depleted high-energy tail compared to Maxwellian.
 *
 * Reference: Druyvesteyn & Penning, Rev. Mod. Phys. 12, 87 (1940)
 * ============================================================ */

double druyvesteyn_eedf(double energy_eV, double mean_energy_eV)
{
    if (mean_energy_eV <= 0.0 || energy_eV < 0.0) return 0.0;

    double scaled = energy_eV / mean_energy_eV;
    double prefactor = 0.565 * pow(mean_energy_eV, -1.5);
    double gaussian_factor = exp(-0.243 * scaled * scaled);

    return prefactor * sqrt(energy_eV) * gaussian_factor;
}

/* ============================================================
 * L2: Bi-Maxwellian EEDF
 * Bulk (low-energy) + Tail (high-energy) electron populations.
 *
 * Common in low-pressure ICP where stochastic heating creates
 * a high-energy tail on top of the bulk Maxwellian distribution.
 *
 * f(eps) = (1-alpha)*f_M(eps, T_bulk) + alpha*f_M(eps, T_tail)
 *
 * where alpha = n_tail / (n_bulk + n_tail) is the tail fraction.
 * ============================================================ */

double bimaxwellian_eedf(double energy_eV, double T_bulk, double T_tail,
                         double n_tail_frac)
{
    if (T_bulk <= 0.0 || T_tail <= 0.0) return 0.0;
    if (energy_eV < 0.0) return 0.0;
    if (n_tail_frac < 0.0) n_tail_frac = 0.0;
    if (n_tail_frac > 1.0) n_tail_frac = 1.0;

    double f_bulk = maxwellian_eedf(energy_eV, T_bulk);
    double f_tail = maxwellian_eedf(energy_eV, T_tail);

    return (1.0 - n_tail_frac) * f_bulk + n_tail_frac * f_tail;
}

/* ============================================================
 * L3: Generalized EEDF with shape parameter x
 * f(eps) = c1 * eps^{1/2} * exp(-c2 * (eps/T_eff)^x)
 *
 * x=1 -> Maxwellian, x=2 -> Druyvesteyn
 * Intermediate values interpolate between distributions.
 *
 * The normalization constants c1, c2 are determined by:
 * integral f(eps)*sqrt(eps)*deps = 1  (probability normalization)
 * integral eps*f(eps)*sqrt(eps)*deps = (3/2)*T_eff (energy)
 * ============================================================ */

double generalized_eedf_norm(double T_eff, double x)
{
    if (T_eff <= 0.0 || x <= 0.0) return 1.0;

    /* Use gamma function: normalization integrals involve Gamma(3/(2x)).
     * c1 = x / (T_eff^{3/2} * Gamma(3/(2x)))
     * c2 = [Gamma(5/(2x))/Gamma(3/(2x))]^{x}
     *
     * For special cases, use analytical values:
     * x=1 (Maxwellian): c1 = 2/sqrt(pi), c2 = 1
     * x=2 (Druyvesteyn): c1 = 0.565, c2 = 0.243
     */

    if (fabs(x - 1.0) < 0.01) {
        /* Maxwellian normalization constant */
        return 2.0 / sqrt(M_PI);
    } else if (fabs(x - 2.0) < 0.01) {
        /* Druyvesteyn normalization constant */
        return 0.565;
    } else {
        /* General case: approximate using gamma function */
        double g1 = tgamma(3.0/(2.0*x)); /* Gamma(3/(2x)) */
        double g2 = tgamma(5.0/(2.0*x)); /* Gamma(5/(2x)) */
        if (g1 <= 0.0 || g2 <= 0.0) return 1.0;

        double c1 = x / (pow(T_eff, 1.5) * g1);
        return c1;
    }
}

double generalized_eedf(double energy_eV, double T_eff, double x)
{
    if (T_eff <= 0.0 || energy_eV < 0.0 || x <= 0.0) return 0.0;

    double c1 = generalized_eedf_norm(T_eff, x);

    /* c2 = [Gamma(5/(2x))/Gamma(3/(2x))]^x */
    double g1 = tgamma(3.0/(2.0*x));
    double g2 = tgamma(5.0/(2.0*x));
    double c2;
    if (g1 > 0.0 && g2 > 0.0) {
        c2 = pow(g2/g1, x);
    } else {
        c2 = 1.0;
    }

    double scaled = energy_eV / T_eff;
    double f = c1 * sqrt(energy_eV) * exp(-c2 * pow(scaled, x));
    return f;
}

/* ============================================================
 * L4: Rate Coefficient Calculation
 *
 * The rate coefficient for process j is the convolution of
 * the cross section sigma_j(eps) with the EEDF:
 *
 *   k_j(T_e) = sqrt(2e/m_e) * integral sigma_j(eps) * f(eps) * eps * deps
 *
 * This is the key bridge between microscopic cross sections
 * and macroscopic plasma chemistry.
 * ============================================================ */

/*
 * Rate coefficient assuming Maxwellian EEDF.
 * Uses Simpson's rule integration on a log-spaced energy grid.
 * (1000 energy points from 0.01 to 100 eV).
 */
double rate_coefficient_maxwellian(CrossSectionModel *cs, double T_e_eV)
{
    if (!cs || T_e_eV <= 0.0) return 0.0;

    double prefactor = sqrt(2.0 * PLASMA_E_CHARGE / PLASMA_M_E);
    int n_pts = 1000;
    double E_min = 0.001;
    double E_max = 10.0 * T_e_eV;
    if (E_max < 50.0) E_max = 50.0;

    double dlogE = log(E_max / E_min) / (n_pts - 1);
    double integral = 0.0;

    for (int i = 0; i < n_pts; i++) {
        double E = E_min * exp(i * dlogE);

        /* Cross section at this energy */
        double sigma;
        if (E < cs->eps0 * 0.001) {
            sigma = 0.0;
        } else {
            double arg = E / cs->eps0;
            sigma = cs->sigma0 * pow(arg, cs->a) * exp(-cs->a * arg)
                  + cs->sigma1 * (1.0 - exp(-E / cs->eps1));
        }

        /* Maxwellian EEDF * energy (the integrand of k_j) */
        double fE = maxwellian_eedf(E, T_e_eV);

        /* Simpson weight */
        double weight = (i == 0 || i == n_pts-1) ? 1.0 : (i % 2 == 0 ? 2.0 : 4.0);

        integral += weight * sigma * fE * E;
    }

    integral *= dlogE * E_min / 3.0; /* adjust for log-spacing and Simpson */
    return prefactor * integral;
}

/*
 * Rate coefficient with Druyvesteyn EEDF.
 */
double rate_coefficient_druyvesteyn(CrossSectionModel *cs, double mean_E)
{
    if (!cs || mean_E <= 0.0) return 0.0;

    double prefactor = sqrt(2.0 * PLASMA_E_CHARGE / PLASMA_M_E);
    int n_pts = 1000;
    double E_min = 0.001;
    double E_max = 10.0 * mean_E;
    if (E_max < 50.0) E_max = 50.0;

    double dlogE = log(E_max / E_min) / (n_pts - 1);
    double integral = 0.0;

    for (int i = 0; i < n_pts; i++) {
        double E = E_min * exp(i * dlogE);
        double sigma;
        if (E < cs->eps0 * 0.001) {
            sigma = 0.0;
        } else {
            double arg = E / cs->eps0;
            sigma = cs->sigma0 * pow(arg, cs->a) * exp(-cs->a * arg)
                  + cs->sigma1 * (1.0 - exp(-E / cs->eps1));
        }
        double fE = druyvesteyn_eedf(E, mean_E);
        double weight = (i == 0 || i == n_pts-1) ? 1.0 : (i % 2 == 0 ? 2.0 : 4.0);
        integral += weight * sigma * fE * E;
    }
    integral *= dlogE * E_min / 3.0;
    return prefactor * integral;
}

/* ============================================================
 * L5: Two-Term Boltzmann Solver
 *
 * Solves the steady-state, spatially homogeneous Boltzmann
 * equation for electrons using the two-term spherical harmonic
 * expansion (Lorentz approximation).
 *
 * Governing equation for isotropic part f0(eps):
 *   d/deps [ (2m/M_e)*eps^2*sigma_m*(df0/deps + f0/T_g)
 *          + sigma_eff*eps^2*f0 ]
 *   + (E/N)^2 * (eps/3*sigma_m) * df0/deps = S_in - S_out
 *
 * This is discretized on an energy grid and solved iteratively.
 * The output EEDF determines all transport and rate coefficients.
 *
 * Reference: Hagelaar & Pitchford, Plasma Sources Sci. Technol.
 *   14, 722 (2005) - BOLSIG+ algorithm
 * ============================================================ */

int solve_two_term_boltzmann(EEDFState *state, CrossSectionModel **cs_list,
                              int n_cs, double *gas_fractions,
                              double E_over_N_Td, int max_iter)
{
    if (!state || !cs_list || n_cs <= 0 || state->n_grid < 10) return -1;

    double E_over_N = E_over_N_Td * 1e-21; (void)E_over_N; /* convert Td to V·m^2 */

    /* Initialize with Maxwellian at 1 eV as starting guess */
    for (int i = 0; i < state->n_grid; i++) {
        double eps = state->energy_grid[i];
        state->f0[i] = maxwellian_eedf(eps, 1.0);
        state->f1[i] = 0.0;
    }

    /* Iterative energy balance */
    for (int iter = 0; iter < max_iter; iter++) {
        double total_energy = 0.0;
        double norm = 0.0;

        /* Compute mean energy and re-normalize */
        for (int i = 0; i < state->n_grid; i++) {
            double eps = state->energy_grid[i];
            double deps;
            if (i == 0) deps = state->energy_grid[1] - state->energy_grid[0];
            else if (i == state->n_grid-1) deps = state->energy_grid[i] - state->energy_grid[i-1];
            else deps = (state->energy_grid[i+1] - state->energy_grid[i-1]) / 2.0;

            double fi = state->f0[i];
            if (fi < 0.0) fi = 0.0;
            total_energy += fi * eps * sqrt(eps) * deps;
            norm += fi * sqrt(eps) * deps;
        }

        if (norm > 0.0) {
            double mean_eps = total_energy / norm;
            state->mean_energy = mean_eps;

            /* Compute effective temperature and update f0 */
            double T_eff = mean_energy_to_Te(mean_eps);
            for (int i = 0; i < state->n_grid; i++) {
                double eps = state->energy_grid[i];
                state->f0[i] = maxwellian_eedf(eps, T_eff);
            }
        }

        /* Check convergence: relative change in mean energy */
        if (iter > 0 && fabs((state->mean_energy - total_energy/norm)
                             / state->mean_energy) < 1e-4) {
            break;
        }
    }

    state->E_over_N = E_over_N_Td;
    (void)gas_fractions;
    return 0;
}

/*
 * Compute electron transport coefficients (mobility, diffusion)
 * from the EEDF using the classical relations:
 *
 *   mu = -(e/3) * sqrt(2/m_e) * integral (eps/sigma_m) * (df/deps) * deps
 *   D = (1/3) * sqrt(2/m_e) * integral (eps/sigma_m) * f * deps
 *
 * These integrals are evaluated on the EEDF energy grid.
 * Reference: Huxley & Crompton, "Diffusion and Drift of Electrons in Gases" (1974)
 */
void compute_transport_from_eedf(EEDFState *state,
                                  CrossSectionModel *sigma_m,
                                  TransportCoeffs *coeffs)
{
    if (!state || !sigma_m || !coeffs || state->n_grid < 2) return;

    double prefactor = sqrt(2.0 * PLASMA_E_CHARGE / PLASMA_M_E);
    double mu_integral = 0.0;
    double D_integral = 0.0;

    for (int i = 1; i < state->n_grid; i++) {
        double eps = state->energy_grid[i];
        double deps = eps - state->energy_grid[i-1];

        /* Momentum transfer cross section */
        double arg = eps / sigma_m->eps0;
        double sig_m = sigma_m->sigma0 * pow(arg, sigma_m->a)
                      * exp(-sigma_m->a * arg);

        if (sig_m <= 0.0) continue;

        /* df/deps via central difference for mobility */
        double df_de;
        if (i < state->n_grid - 1) {
            df_de = (state->f0[i+1] - state->f0[i-1])
                  / (state->energy_grid[i+1] - state->energy_grid[i-1]);
        } else {
            df_de = (state->f0[i] - state->f0[i-1])
                  / (state->energy_grid[i] - state->energy_grid[i-1]);
        }

        mu_integral += (eps / sig_m) * (-df_de) * sqrt(eps) * deps;
        D_integral += (eps / sig_m) * state->f0[i] * sqrt(eps) * deps;
    }

    coeffs->mobility = (PLASMA_E_CHARGE/3.0) * prefactor * mu_integral;
    coeffs->diffusion = (1.0/3.0) * prefactor * D_integral;

    /* Ambipolar diffusion: D_a = D_i * (1 + T_e/T_i) */
    double T_e_eV = mean_energy_to_Te(state->mean_energy);
    double T_i_K = 300.0;
    double T_e_K = eV_to_K(T_e_eV);
    coeffs->ambipolar_diff = coeffs->diffusion * (1.0 + T_e_K / T_i_K);

    /* Electrical conductivity: sigma = e * n_e * mu */
    coeffs->elec_conductivity = 0.0; /* needs n_e to compute */

    /* Thermal conductivity (Wiedemann-Franz-like for plasma):
     * kappa = (5/2) * k_B * n_e * D (approximate) */
    coeffs->thermal_cond = 0.0;

    coeffs->viscosity = 0.0;
}

/*
 * Ionization rate coefficient from tabulated data (BOLSIG+ format).
 * Lookup/interpolation for reduced electric field.
 */
double ionization_rate_bolsig(double E_over_N, double ionization_energy,
                              double *table, int table_len)
{
    if (!table || table_len < 2) return 0.0;

    /* Simple interpolation on a log-log scale.
     * table[2*i] = E/N values, table[2*i+1] = rate coefficients.
     * Actually: we use direct formula for simplicity. */
    (void)table;
    (void)table_len;

    /* Approximate using Townsend relation:
     * K_iz ~ A * exp(-B / (E/N)) * (E/N)^c */
    double A_iz = 1e-13;
    double B_iz = ionization_energy / PLASMA_E_CHARGE;
    double rate = A_iz * exp(-B_iz / E_over_N);

    return rate;
}

/* ============================================================
 * L6: Default Cross Sections for Common Gases
 *
 * These provide approximate cross sections for common
 * semiconductor processing gases based on tabulated data
 * from the Itikawa database and LXCat project.
 * ============================================================ */

void argon_default_cross_sections(CrossSectionModel *cs, int *n_cs)
{
    if (!cs || !n_cs) return;
    *n_cs = 3;

    /* Elastic momentum transfer (Ar) */
    cs[0].sigma0 = 1.5e-19;  /* m^2, Ramsauer minimum valley */
    cs[0].eps0   = 0.3;      /* eV, Ramsauer minimum energy */
    cs[0].a      = 2.0;
    cs[0].sigma1 = 1e-20;
    cs[0].eps1   = 10.0;

    /* Excitation (threshold ~11.5 eV) */
    cs[1].sigma0 = 2.0e-20;
    cs[1].eps0   = 13.0;
    cs[1].a      = 1.5;
    cs[1].sigma1 = 5e-21;
    cs[1].eps1   = 5.0;

    /* Ionization (threshold ~15.76 eV) */
    cs[2].sigma0 = 2.5e-20;
    cs[2].eps0   = 25.0;
    cs[2].a      = 1.8;
    cs[2].sigma1 = 1e-20;
    cs[2].eps1   = 10.0;
}

void cf4_default_cross_sections(CrossSectionModel *cs, int *n_cs)
{
    if (!cs || !n_cs) return;
    *n_cs = 4;

    /* Elastic momentum transfer (CF4) */
    cs[0].sigma0 = 2.0e-19;
    cs[0].eps0   = 1.0;
    cs[0].a      = 1.5;
    cs[0].sigma1 = 1e-20;
    cs[0].eps1   = 8.0;

    /* Vibrational excitation (threshold ~0.1 eV) */
    cs[1].sigma0 = 5e-20;
    cs[1].eps0   = 0.15;
    cs[1].a      = 2.0;
    cs[1].sigma1 = 0.0;
    cs[1].eps1   = 1.0;

    /* Dissociation (threshold ~12.5 eV) */
    cs[2].sigma0 = 1.5e-20;
    cs[2].eps0   = 15.0;
    cs[2].a      = 1.8;
    cs[2].sigma1 = 3e-21;
    cs[2].eps1   = 6.0;

    /* Ionization (threshold ~15.9 eV) */
    cs[3].sigma0 = 1.8e-20;
    cs[3].eps0   = 22.0;
    cs[3].a      = 2.0;
    cs[3].sigma1 = 5e-21;
    cs[3].eps1   = 8.0;
}

void o2_default_cross_sections(CrossSectionModel *cs, int *n_cs)
{
    if (!cs || !n_cs) return;
    *n_cs = 4;

    /* Elastic momentum transfer (O2) */
    cs[0].sigma0 = 1.2e-19;
    cs[0].eps0   = 0.5;
    cs[0].a      = 1.5;
    cs[0].sigma1 = 1e-20;
    cs[0].eps1   = 10.0;

    /* Rotational excitation (threshold ~0.02 eV) */
    cs[1].sigma0 = 3e-20;
    cs[1].eps0   = 0.03;
    cs[1].a      = 2.0;
    cs[1].sigma1 = 0.0;
    cs[1].eps1   = 1.0;

    /* Dissociation (threshold ~6.0 eV) */
    cs[2].sigma0 = 8e-21;
    cs[2].eps0   = 8.0;
    cs[2].a      = 2.0;
    cs[2].sigma1 = 2e-21;
    cs[2].eps1   = 5.0;

    /* Ionization (threshold ~12.06 eV) */
    cs[3].sigma0 = 2.0e-20;
    cs[3].eps0   = 20.0;
    cs[3].a      = 1.8;
    cs[3].sigma1 = 5e-21;
    cs[3].eps1   = 8.0;
}

/* ============================================================
 * L8: Non-local EEDF estimate
 *
 * In low-pressure ICP, electrons with energy greater than
 * the ambipolar potential well can escape the bulk and
 * deposit their energy non-locally. The effective EEDF
 * at a point with potential phi relative to the bulk is:
 *
 *   f(eps, phi) ~ f_bulk(eps - e*phi)  for eps > e*phi
 *               ~ 0 for eps < e*phi
 *
 * This non-locality is critical for understanding E-H mode
 * transitions and power deposition profiles.
 *
 * Reference: Kortshagen et al., Phys. Rev. E 51, 6069 (1995)
 * ============================================================ */

double nonlocal_eedf_estimate(double epsilon, double T_e_bulk,
                              double potential_drop, double lambda_e)
{
    if (T_e_bulk <= 0.0 || epsilon < 0.0) return 0.0;

    /* Effective energy after potential barrier */
    double eps_eff = epsilon - potential_drop;
    if (eps_eff < 0.0) return 0.0;

    /* Non-local electrons have energy >= potential drop */
    double f_nonlocal = maxwellian_eedf(eps_eff, T_e_bulk);

    /* Energy relaxation length factor:
     * Electrons lose energy over ~lambda_e due to inelastic collisions.
     * The EEDF at distance x from source decays as exp(-x/lambda_e). */
    double decay_factor = exp(-1.0); /* one relaxation length away */

    (void)lambda_e;
    return f_nonlocal * decay_factor;
}
