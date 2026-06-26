/**
 * @file  dusty_crystal.c
 * @brief Dust crystal formation, structure and phase transitions.
 *
 * When dust grains in a plasma are strongly coupled (Gamma >> 1),
 * they can self-organize into crystalline structures — plasma
 * crystals. These are macroscopic (micron-sized grains, mm spacing)
 * and can be observed with visible-light cameras.
 *
 * L6-L8: Yukawa potential, crystal energy, pair correlations,
 *        phase diagram, nucleation.
 *
 * References:
 *   Thomas et al. (1994), Phys. Rev. Lett. 73, 652
 *   Hamaguchi et al. (1997), Phys. Rev. E 56, 4671
 *   Vaulina et al. (2002), Phys. Rev. Lett. 88, 035001
 */

#include "dusty_crystal.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* ================================================================
 * L6 — Phase Diagram and Crystal Condition
 * ================================================================ */

double dust_critical_coupling_yukawa(double kappa)
{
    /* Empirical phase boundary for Yukawa systems.
     * Gamma_crit(kappa) = 170/(1+0.18*kappa^2) + 106*(1-1/(1+0.06*kappa^2))
     * Vaulina et al. (2002), Phys. Rev. Lett. 88, 035001 */
    if (kappa < 0.0) kappa = 0.0;
    double k2 = kappa * kappa;
    double term1 = 170.0 / (1.0 + 0.18 * k2);
    double term2 = 106.0 * (1.0 - 1.0 / (1.0 + 0.06 * k2));
    return term1 + term2;
}

int dust_phase_determine(double Gamma, double kappa)
{
    /* 0=fluid, 1=strongly-coupled liquid, 2=solid */
    double Gamma_crit = dust_critical_coupling_yukawa(kappa);
    if (Gamma >= Gamma_crit) return 2;
    if (Gamma >= 1.0) return 1;
    return 0;
}

double dust_thermal_velocity(double T_d, double m_d)
{
    if (T_d <= 0.0 || m_d <= 0.0) return 0.0;
    return sqrt(DUSTY_KB * T_d / m_d);
}

double dust_einstein_frequency(double Q_d, double m_d, double a, double kappa)
{
    /* Omega_E^2 = (Q_d^2/(4*pi*eps0*m_d)) * (kappa^2/a^3) * (1+kappa)*exp(-kappa)
     * Nunomura et al. (2002), Phys. Rev. E 65, 066402 */
    if (a <= 0.0 || m_d <= 0.0 || Q_d <= 0.0) return 0.0;
    double prefactor = Q_d * Q_d / (4.0 * M_PI * DUSTY_EPS0 * m_d * a * a * a);
    double screening = kappa * kappa * (1.0 + kappa) * exp(-kappa);
    return sqrt(prefactor * screening);
}

/* ================================================================
 * L3 — Yukawa Potential and Force
 * ================================================================ */

double yukawa_potential(double Q1, double Q2, double r, double lambda_D)
{
    if (r <= 0.0 || lambda_D <= 0.0) return 0.0;
    return (Q1 * Q2 / (4.0 * M_PI * DUSTY_EPS0 * r)) * exp(-r / lambda_D);
}

double yukawa_force_magnitude(double Q1, double Q2, double r, double lambda_D)
{
    if (r <= 0.0 || lambda_D <= 0.0) return 0.0;
    double coulomb_force = Q1 * Q2 / (4.0 * M_PI * DUSTY_EPS0 * r * r);
    double yukawa_factor = (1.0 + r / lambda_D) * exp(-r / lambda_D);
    return coulomb_force * yukawa_factor;
}

void yukawa_force_vector(
    const double r1[3], const double r2[3],
    double Q1, double Q2, double lambda_D, double F_out[3])
{
    if (!r1 || !r2 || !F_out) return;
    double dx = r1[0] - r2[0];
    double dy = r1[1] - r2[1];
    double dz = r1[2] - r2[2];
    double r = sqrt(dx*dx + dy*dy + dz*dz);
    if (r < 1e-30) {
        F_out[0] = F_out[1] = F_out[2] = 0.0;
        return;
    }
    double F_mag = yukawa_force_magnitude(Q1, Q2, r, lambda_D);
    double inv_r = 1.0 / r;
    F_out[0] = F_mag * dx * inv_r;
    F_out[1] = F_mag * dy * inv_r;
    F_out[2] = F_mag * dz * inv_r;
}

double dust_crystal_total_energy(
    const double **positions, size_t n_grains,
    double Q_d, double lambda_D, double Lx, double Ly)
{
    if (!positions || n_grains < 2) return 0.0;
    double U_total = 0.0;
    for (size_t i = 0; i < n_grains; i++) {
        for (size_t j = i + 1; j < n_grains; j++) {
            double dx = positions[i][0] - positions[j][0];
            double dy = positions[i][1] - positions[j][1];
            double dz = positions[i][2] - positions[j][2];
            if (Lx > 0.0) dx -= Lx * round(dx / Lx);
            if (Ly > 0.0) dy -= Ly * round(dy / Ly);
            double r = sqrt(dx*dx + dy*dy + dz*dz);
            if (r < 1e-15) continue;
            U_total += yukawa_potential(Q_d, Q_d, r, lambda_D);
        }
    }
    return U_total;
}

double dust_yukawa_madelung(int lattice_type, double kappa)
{
    double alpha = 0.0;
    double a = 1.0;
    double lambda_D = (kappa > 0.0) ? a / kappa : 1e10;
    int n_shells = 20;

    if (lattice_type == 0) {
        for (int n = 1; n <= n_shells; n++) {
            double r = n * a;
            alpha += yukawa_potential(1.0, 1.0, r, lambda_D)
                     * (4.0 * M_PI * DUSTY_EPS0 * a);
        }
        alpha *= 2.0;
    } else if (lattice_type == 1) {
        for (int nx = -n_shells; nx <= n_shells; nx++) {
            for (int ny = -n_shells; ny <= n_shells; ny++) {
                if (nx == 0 && ny == 0) continue;
                double x = nx * a + ny * a * 0.5;
                double y = ny * a * 0.8660254037844386;
                double r = sqrt(x*x + y*y);
                alpha += yukawa_potential(1.0, 1.0, r, lambda_D)
                         * (4.0 * M_PI * DUSTY_EPS0 * a);
            }
        }
    } else {
        for (int nx = -n_shells; nx <= n_shells; nx++) {
            for (int ny = -n_shells; ny <= n_shells; ny++) {
                for (int nz = -n_shells; nz <= n_shells; nz++) {
                    if (nx == 0 && ny == 0 && nz == 0) continue;
                    double x = nx * a + 0.5 * a * (ny % 2);
                    double y = ny * a + 0.5 * a * (nz % 2);
                    double z = nz * a;
                    double r = sqrt(x*x + y*y + z*z);
                    alpha += yukawa_potential(1.0, 1.0, r, lambda_D)
                             * (4.0 * M_PI * DUSTY_EPS0 * a);
                }
            }
        }
    }
    return alpha;
}

/* ================================================================
 * L8 — Pair Correlation g(r) and Structure Factor S(k)
 * ================================================================ */

void dust_pair_correlation_hnc(
    double Gamma, double kappa, double r_max,
    int n_bins, double *r, double *g)
{
    if (!r || !g || n_bins <= 0) return;
    double dr = r_max / (double)n_bins;
    double a_ws = 1.0;
    double Gamma_eff = Gamma * exp(-kappa);

    for (int i = 0; i < n_bins; i++) {
        double ri = (i + 0.5) * dr;
        r[i] = ri;
        if (ri < 0.8 * a_ws) {
            g[i] = 0.0;
        } else {
            double r_norm = ri / a_ws;
            double r_peak = (Gamma_eff < 10.0) ? 1.1 : 1.05;
            double xi = 0.5 * pow(Gamma_eff, 0.25);
            double d_osc = a_ws * (0.9 + 0.1 / (1.0 + 0.01 * Gamma_eff));
            g[i] = 1.0 + 2.5 * exp(-fabs(r_norm - r_peak) / xi)
                   * cos(2.0 * M_PI * (r_norm - r_peak) * a_ws / d_osc);
            if (g[i] < 0.0) g[i] = 0.0;
            g[i] += (1.0 - g[i]) * (1.0 - exp(-ri / (3.0 * xi)));
        }
    }
}

void dust_structure_factor(
    const double *r, const double *g, int n_bins,
    double n_d, double k_max, int n_k,
    double *k_out, double *S_out)
{
    if (!r || !g || !k_out || !S_out || n_bins <= 0 || n_k <= 0) return;
    double dk = k_max / (double)(n_k - 1);
    for (int ik = 0; ik < n_k; ik++) {
        double k = (ik == 0) ? 1e-30 : ik * dk;
        k_out[ik] = k;
        double integral = 0.0;
        double dr_bin = r[1] - r[0];
        for (int i = 0; i < n_bins; i++) {
            double ri = r[i];
            double h = g[i] - 1.0;
            double kr = k * ri;
            double sinc = (kr < 1e-6) ? (1.0 - kr*kr/6.0) : (sin(kr) / kr);
            integral += ri * ri * h * sinc * dr_bin;
        }
        S_out[ik] = 1.0 + 4.0 * M_PI * n_d * integral;
    }
}

/* ================================================================
 * L6 — Melting Criteria
 * ================================================================ */

int dust_lindemann_melting(double rms_displacement, double a, double c_L)
{
    if (a <= 0.0) return 1;
    return (rms_displacement / a > c_L) ? 1 : 0;
}

double dust_rms_displacement_thermal(double T_d, double m_d, double Omega_E)
{
    if (Omega_E <= 0.0 || m_d <= 0.0) return 0.0;
    return sqrt(3.0 * DUSTY_KB * T_d / (m_d * Omega_E * Omega_E));
}

/* ================================================================
 * L8 — Nucleation and Grain Boundaries
 * ================================================================ */

double dust_nucleation_rate(
    double Gamma, double kappa, double T_d, double n_d)
{
    if (T_d <= 0.0 || n_d <= 0.0) return 0.0;
    double J0 = n_d * 50.0;
    double f_kappa = 1.0 / (1.0 + 0.1 * kappa * kappa);
    double barrier = Gamma * f_kappa;
    if (barrier > 50.0) return 0.0;
    return J0 * exp(-barrier);
}

double dust_grain_boundary_energy(double theta, double E_0, double A)
{
    if (theta <= 0.0) return 0.0;
    double theta_deg = theta * 180.0 / M_PI;
    if (theta_deg < 15.0) {
        return E_0 * theta * (A - log(theta));
    } else {
        double theta_sat = 15.0 * M_PI / 180.0;
        return E_0 * theta_sat * (A - log(theta_sat));
    }
}
