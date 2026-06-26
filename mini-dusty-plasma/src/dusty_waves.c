/**
 * @file  dusty_waves.c
 * @brief Wave dispersion relations in dusty plasmas.
 *
 * Dusty plasmas support a rich variety of wave modes that span
 * an enormous frequency range — from ~10 Hz (DAW) to ~GHz (electron
 * plasma waves). The key novelty is the dust-acoustic wave (DAW),
 * first predicted by Rao, Shukla & Yu (1990).
 *
 * L4-L6: DAW, DIAW, DLW dispersion relations and damping mechanisms.
 *
 * References:
 *   Rao, Shukla & Yu (1990), Planet. Space Sci. 38, 543
 *   Shukla & Silin (1992), Physica Scripta 45, 504
 *   Melandso (1996), Phys. Plasmas 3, 3890
 *   Piel (2010), "Plasma Physics", Ch. 11
 */

#include "dusty_waves.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* ================================================================
 * L4 — Dust-Acoustic Wave (DAW) Dispersion
 * ================================================================
 *
 * The DAW is the quintessential dusty plasma wave mode. It arises
 * because the massive dust grains provide inertia while the light
 * electrons provide pressure (restoring force).
 *
 * The dispersion relation is formally identical to the ion-acoustic
 * wave but with ions replaced by dust and shifted to much lower
 * frequencies (by factor sqrt(m_i/m_d) ~ 10^-6).
 */

WaveMode dust_acoustic_wave_dispersion(
    double k, double c_da, double lambda_D)
{
    /* Standard DAW dispersion (fluid model, cold ions).
     *
     * omega^2 = k^2 * c_da^2 / (1 + k^2 * lambda_D^2)
     *
     * Phase velocity: v_phi = c_da / sqrt(1 + k^2 * lambda_D^2)
     * Group velocity: v_g = c_da / (1 + k^2 * lambda_D^2)^(3/2)
     *
     * Long wavelength limit (k*lambda_D << 1):
     *   omega ≈ k * c_da  (sound-like, non-dispersive)
     *
     * Short wavelength limit (k*lambda_D >> 1):
     *   omega → omega_pd  (constant, dust plasma frequency)
     *
     * Rao, Shukla & Yu (1990), Eq. (5) */
    WaveMode wm;
    wm.k = k;

    if (k <= 0.0 || c_da <= 0.0) {
        wm.omega_r = 0.0;
        wm.omega_i = 0.0;
        wm.phase_velocity = 0.0;
        wm.group_velocity = 0.0;
        return wm;
    }

    double k2_lambda2 = k * k * lambda_D * lambda_D;
    double denom = 1.0 + k2_lambda2;

    wm.omega_r = sqrt(k * k * c_da * c_da / denom);
    wm.omega_i = 0.0; /* fluid model has no damping */
    wm.phase_velocity = wm.omega_r / k;
    wm.group_velocity = c_da / (denom * sqrt(denom));

    return wm;
}

WaveMode dust_acoustic_wave_with_ion_temp(
    double k, double Z_d, double T_e, double T_i, double m_d,
    double n_e, double n_i, double lambda_D)
{
    /* DAW dispersion including finite ion temperature.
     *
     * Finite T_i modifies the effective temperature:
     *
     * T_eff = T_i * T_e / (T_i + Z_d * T_e * n_e/n_i)
     *
     * This is a harmonic mean weighted by charge density ratios.
     * When T_i << Z_d * T_e * n_e/n_i, we recover T_eff ≈ T_i,
     * but typically T_eff lies between T_i and T_e*Z_d.
     *
     * The modified acoustic speed:
     * c_da_eff = sqrt(Z_d * k_B * T_eff / m_d)
     *
     * Finite T_i increases the phase velocity and introduces
     * ion Landau damping as an additional mechanism.
     *
     * Merlino et al. (1998), Phys. Plasmas 5, 1607 */
    WaveMode wm;
    wm.k = k;

    if (k <= 0.0 || Z_d <= 0.0 || T_e <= 0.0 || m_d <= 0.0
        || n_e <= 0.0 || n_i <= 0.0) {
        wm.omega_r = 0.0;
        wm.omega_i = 0.0;
        wm.phase_velocity = 0.0;
        wm.group_velocity = 0.0;
        return wm;
    }

    /* Effective temperature */
    double ZTe_ne_ni = Z_d * T_e * n_e / n_i;
    double T_eff = T_i * T_e / (T_i + ZTe_ne_ni);
    if (T_eff <= 0.0) T_eff = T_e;

    double c_da_eff = sqrt(Z_d * DUSTY_KB * T_eff / m_d);
    double k2_lambda2 = k * k * lambda_D * lambda_D;
    double denom = 1.0 + k2_lambda2;

    wm.omega_r = sqrt(k * k * c_da_eff * c_da_eff / denom);
    wm.omega_i = 0.0;
    wm.phase_velocity = wm.omega_r / k;
    wm.group_velocity = c_da_eff / (denom * sqrt(denom));

    return wm;
}

double dust_acoustic_wave_landau_damping(
    double k, double omega, double omega_pd, double v_thd)
{
    /* Landau damping of DAW by resonant dust particles.
     *
     * gamma_L / omega = -sqrt(pi/8) * (omega_pd/(k*v_thd))^3
     *                   * exp(-omega^2/(2*k^2*v_thd^2))
     *
     * This is the collisionless (Landau) damping from the kinetic
     * theory for a Maxwellian dust velocity distribution.
     *
     * Damping is exponentially small when:
     *   v_phi >> v_thd  i.e., omega/k >> sqrt(k_B*T_d/m_d)
     *
     * For DAW: v_phi ~ 30 mm/s, v_thd ~ 0.5 mm/s (T_d = 0.03 eV, m_d = 10^-15 kg)
     * → v_phi/v_thd ~ 60 → Landau damping negligible.
     *
     * For heavier damping (smaller Z_d or larger T_d), Landau
     * damping can suppress the DAW entirely.
     *
     * Rosenberg (1993), Planet. Space Sci. 41, 229 */
    if (k <= 0.0 || omega <= 0.0 || v_thd <= 0.0) return 0.0;

    double z = omega / (k * v_thd);
    double z2 = z * z;

    /* Guard against overflow in exp(z^2) */
    if (z2 > 50.0) return 0.0;

    double factor = sqrt(M_PI / 8.0)
                    * pow(omega_pd / (k * v_thd), 3.0)
                    * exp(-0.5 * z2);

    return -omega * factor;
}

double dust_acoustic_wave_collisional_damping(double nu_dn)
{
    /* Collisional damping from dust-neutral friction.
     *
     * gamma_coll = -nu_dn / 2
     *
     * This arises from the drag force term in the dust momentum
     * equation: m_d * n_d * nu_dn * v_d, which acts like a friction.
     *
     * The factor 1/2 comes from the energy dissipation rate:
     * dE/dt ∝ -nu_dn * E → E(t) = E_0 * exp(-nu_dn * t)
     * → amplitude ∝ exp(-nu_dn * t/2).
     *
     * For typical lab conditions: p = 10 Pa, nu_dn ~ 100 Hz,
     * so gamma_coll ~ -50 s^-1. This is often larger than
     * Landau damping and can completely suppress DAW propagation. */
    return -0.5 * nu_dn;
}

WaveMode dust_acoustic_wave_full(
    double k, double c_da, double lambda_D, double omega_pd,
    double v_thd, double nu_dn)
{
    /* Full DAW complex frequency combining both damping mechanisms.
     *
     * omega = omega_fluid + i * (gamma_Landau + gamma_coll)
     *
     * The wave is observable only if |omega_r| > |omega_i|,
     * i.e., the quality factor Q = |omega_r / (2*|omega_i|)| > 0.5. */
    WaveMode wm = dust_acoustic_wave_dispersion(k, c_da, lambda_D);
    double gamma_L = dust_acoustic_wave_landau_damping(
        k, wm.omega_r, omega_pd, v_thd);
    double gamma_c = dust_acoustic_wave_collisional_damping(nu_dn);
    wm.omega_i = gamma_L + gamma_c;
    return wm;
}

/* ================================================================
 * L4 — Dust-Ion-Acoustic Wave (DIAW)
 * ================================================================ */

WaveMode dust_ion_acoustic_wave_dispersion(
    double k, double c_s, double lambda_De)
{
    /* DIAW dispersion: ions provide inertia, electrons provide pressure.
     * Dust is stationary (frozen) on the DIAW timescale.
     *
     * omega^2 = k^2 * c_s^2 / (1 + k^2 * lambda_De^2)
     *
     * The dust modifies DIAW primarily by:
     * 1. Depleting free electrons (quasineutrality: n_e = n_i - Z_d*n_d)
     * 2. Providing an additional damping channel via charge fluctuations
     *
     * For typical lab parameters: f_DIAW ~ 100 kHz - 1 MHz,
     * vs f_DAW ~ 10-100 Hz — a 10^4 range!
     *
     * Shukla & Silin (1992), Physica Scripta 45, 504 */
    WaveMode wm;
    wm.k = k;

    if (k <= 0.0 || c_s <= 0.0) {
        wm.omega_r = 0.0;
        wm.omega_i = 0.0;
        wm.phase_velocity = 0.0;
        wm.group_velocity = 0.0;
        return wm;
    }

    double k2_lambda2 = k * k * lambda_De * lambda_De;
    double denom = 1.0 + k2_lambda2;

    wm.omega_r = sqrt(k * k * c_s * c_s / denom);
    wm.omega_i = 0.0;
    wm.phase_velocity = wm.omega_r / k;
    wm.group_velocity = c_s / (denom * sqrt(denom));

    return wm;
}

WaveMode dust_ion_acoustic_wave_charge_fluctuation(
    double k, double c_s, double lambda_De,
    double nd_over_ni, double tau_charge)
{
    /* DIAW modified by dust charge fluctuations.
     *
     * When dust charge fluctuates (due to discrete collection),
     * the waves can be damped or even become unstable.
     *
     * The growth rate from charge fluctuations:
     * gamma_cf ∝ nd_over_ni * omega * tau_charge / (1 + omega^2 * tau_charge^2)
     *
     * This is typically a damping mechanism (gamma < 0) because
     * charge fluctuations extract energy from the wave.
     *
     * Varma, Shukla & Krishan (1993), Phys. Rev. E 47, 3612 */
    WaveMode wm = dust_ion_acoustic_wave_dispersion(k, c_s, lambda_De);

    if (wm.omega_r <= 0.0 || tau_charge <= 0.0) return wm;

    double omega_tau = wm.omega_r * tau_charge;
    double gamma_cf = -nd_over_ni * wm.omega_r * omega_tau
                      / (1.0 + omega_tau * omega_tau);

    wm.omega_i = gamma_cf;
    return wm;
}

/* ================================================================
 * L8 — Dust Lattice Waves (DLW) - 1D Chain
 * ================================================================
 *
 * In strongly coupled (crystal) dusty plasmas, grains oscillate
 * about their equilibrium lattice positions. The normal modes are
 * dust lattice waves (DLW). Unlike DAW which is a fluid mode,
 * DLW is a discrete lattice mode analogous to phonons in solids.
 *
 * The dispersion relation involves a lattice sum over Yukawa
 * interactions between all particle pairs.
 */

WaveMode dust_lattice_wave_1d_longitudinal(
    double k, double a, double Q_d, double m_d,
    double kappa, int n_shells)
{
    /* 1D longitudinal DLW: grains oscillate along the chain.
     *
     * omega_L^2(k) = (Q_d^2 / (4*pi*eps0 * m_d * a^3))
     *   * sum_{n=1}^{N_shells} (2/n^3) * (1 + kappa*n)
     *   * exp(-kappa*n) * [1 - cos(k * n * a)]
     *
     * Each term in the sum represents the contribution from the
     * nth pair of neighbors (both directions).
     *
     * For nearest-neighbor only (N_shells=1, large kappa):
     * omega_L^2(k) ≈ 2*Omega_0^2 * [1 - cos(k*a)]
     *
     * where Omega_0 is the Einstein frequency.
     *
     * For long wavelength (k*a << 1):
     * omega ≈ c_DLW * k  (sound-like)
     * c_DLW = a * Omega_E is the dust-lattice sound speed.
     *
     * Melandso (1996), Phys. Plasmas 3, 3890 */
    WaveMode wm;
    wm.k = k;

    if (a <= 0.0 || Q_d <= 0.0 || m_d <= 0.0 || n_shells < 1) {
        wm.omega_r = 0.0;
        wm.omega_i = 0.0;
        wm.phase_velocity = 0.0;
        wm.group_velocity = 0.0;
        return wm;
    }

    double prefactor = Q_d * Q_d / (4.0 * M_PI * DUSTY_EPS0 * m_d * a * a * a);
    double sum = 0.0;

    for (int n = 1; n <= n_shells; n++) {
        double yukawa_factor = (1.0 + kappa * (double)n) * exp(-kappa * (double)n);
        double lattice_factor = 1.0 - cos(k * (double)n * a);
        sum += (2.0 / ((double)n * (double)n * (double)n))
               * yukawa_factor * lattice_factor;
    }

    if (sum < 0.0) sum = 0.0;
    wm.omega_r = sqrt(prefactor * sum);
    wm.omega_i = 0.0;

    if (k > 0.0 && wm.omega_r > 0.0) {
        wm.phase_velocity = wm.omega_r / k;

        /* Group velocity: analytical nearest-neighbor approximation.
         * For omega^2 = 2*Omega_0^2*(1 - cos(ka)):
         * v_g = d_omega/dk = Omega_0^2 * a * sin(ka) / omega.
         * This avoids recursive self-calls for the numerical derivative. */
        double Omega_0_sq = prefactor * (1.0 + kappa) * exp(-kappa) / (1.0 * 1.0 * 1.0);
        wm.group_velocity = Omega_0_sq * a * sin(k * a) / wm.omega_r;
    } else {
        wm.phase_velocity = 0.0;
        wm.group_velocity = 0.0;
    }

    return wm;
}

WaveMode dust_lattice_wave_1d_transverse(
    double k, double a, double Q_d, double m_d,
    double kappa, int n_shells)
{
    /* 1D transverse DLW: grains oscillate perpendicular to chain.
     *
     * omega_T^2(k) = (Q_d^2 / (4*pi*eps0 * m_d * a^3))
     *   * sum_{n=1}^{N_shells} (1/n^3) * (1 + kappa*n)
     *   * exp(-kappa*n) * [1 - cos(k * n * a)]
     *
     * The factor 1/n^3 vs 2/n^3 makes transverse modes have
     * LOWER frequency than longitudinal modes. The restoring
     * force for transverse displacement is weaker because the
     * projection of the Yukawa force is smaller.
     *
     * At k = 0: omega_T → 0 (Goldstone mode — broken translational
     * symmetry perpendicular to the chain).
     *
     * Melandso (1996), Phys. Plasmas 3, 3890 */
    WaveMode wm;
    wm.k = k;

    if (a <= 0.0 || Q_d <= 0.0 || m_d <= 0.0 || n_shells < 1) {
        wm.omega_r = 0.0;
        wm.omega_i = 0.0;
        wm.phase_velocity = 0.0;
        wm.group_velocity = 0.0;
        return wm;
    }

    double prefactor = Q_d * Q_d / (4.0 * M_PI * DUSTY_EPS0 * m_d * a * a * a);
    double sum = 0.0;

    for (int n = 1; n <= n_shells; n++) {
        double yukawa_factor = (1.0 + kappa * (double)n) * exp(-kappa * (double)n);
        double lattice_factor = 1.0 - cos(k * (double)n * a);
        sum += (1.0 / ((double)n * (double)n * (double)n))
               * yukawa_factor * lattice_factor;
    }

    if (sum < 0.0) sum = 0.0;
    wm.omega_r = sqrt(prefactor * sum);
    wm.omega_i = 0.0;

    if (k > 0.0 && wm.omega_r > 0.0) {
        wm.phase_velocity = wm.omega_r / k;
        /* Group velocity: analytical approximation (avoids recursion) */
        double Omega_0_sq_T = prefactor * 0.5 * (1.0 + kappa) * exp(-kappa);
        wm.group_velocity = Omega_0_sq_T * a * sin(k * a) / wm.omega_r;
    } else {
        wm.phase_velocity = 0.0;
        wm.group_velocity = 0.0;
    }

    return wm;
}

WaveMode dust_lattice_wave_2d_out_of_plane(
    double kx, double ky, double a, double Q_d, double m_d, double kappa)
{
    /* 2D hexagonal lattice DLW: out-of-plane (vertical) mode.
     *
     * omega_O^2(k) = Omega_0^2 * sum_{j=1}^6 (1 + kappa*a) * exp(-kappa*a)
     *                * (1/a^3) * [1 - cos(k · r_j)]
     *
     * where r_j are the 6 nearest-neighbor vectors of the hexagonal lattice:
     *   r_1 = a*(1, 0),       r_2 = a*(-1, 0)
     *   r_3 = a*(1/2, √3/2),  r_4 = a*(-1/2, -√3/2)
     *   r_5 = a*(-1/2, √3/2), r_6 = a*(1/2, -√3/2)
     *
     * The out-of-plane mode has the highest frequency at zone center
     * because all 6 neighbors contribute in phase.
     *
     * Nunomura et al. (2002), Phys. Rev. E 65, 066402 */
    WaveMode wm;
    wm.k = sqrt(kx*kx + ky*ky);

    if (a <= 0.0 || Q_d <= 0.0 || m_d <= 0.0) {
        wm.omega_r = 0.0;
        wm.omega_i = 0.0;
        wm.phase_velocity = 0.0;
        wm.group_velocity = 0.0;
        return wm;
    }

    /* Six nearest-neighbor vectors in hexagonal lattice */
    double nn_vectors[6][2] = {
        { 1.0,  0.0},
        {-1.0,  0.0},
        { 0.5,  0.8660254037844386},  /* sqrt(3)/2 */
        {-0.5, -0.8660254037844386},
        {-0.5,  0.8660254037844386},
        { 0.5, -0.8660254037844386}
    };

    double prefactor = Q_d * Q_d / (4.0 * M_PI * DUSTY_EPS0 * m_d);
    double yukawa_coupling = (1.0 + kappa) * exp(-kappa) / (a * a * a);
    double sum = 0.0;

    for (int j = 0; j < 6; j++) {
        double k_dot_r = kx * nn_vectors[j][0] * a + ky * nn_vectors[j][1] * a;
        sum += (1.0 - cos(k_dot_r));
    }

    wm.omega_r = sqrt(prefactor * yukawa_coupling * sum);
    wm.omega_i = 0.0;

    if (wm.k > 0.0 && wm.omega_r > 0.0) {
        wm.phase_velocity = wm.omega_r / wm.k;
    } else {
        wm.phase_velocity = 0.0;
    }
    wm.group_velocity = 0.0;

    return wm;
}

WaveMode dust_lattice_wave_2d_in_plane_longitudinal(
    double kx, double ky, double a, double Q_d, double m_d, double kappa)
{
    /* 2D in-plane longitudinal DLW.
     *
     * For in-plane motion, the projection of the force depends on
     * the angle between displacement and neighbor bond direction.
     * This gives a more complex tensor structure.
     *
     * Simplified: use effective spring constant from 2 neighbors
     * along the wavevector direction.
     *
     * Nunomura et al. (2002), Phys. Rev. E 65, 066402 */
    WaveMode wm;
    wm.k = sqrt(kx*kx + ky*ky);

    if (a <= 0.0 || Q_d <= 0.0 || m_d <= 0.0) {
        wm.omega_r = 0.0;
        wm.omega_i = 0.0;
        wm.phase_velocity = 0.0;
        wm.group_velocity = 0.0;
        return wm;
    }

    double prefactor = Q_d * Q_d / (4.0 * M_PI * DUSTY_EPS0 * m_d * a * a * a);
    double yukawa_spring = (1.0 + kappa) * exp(-kappa);

    /* Approximate: effective 1D longitudinal with renormalized spring */
    double sum = 2.0 * yukawa_spring * (1.0 - cos(wm.k * a));
    wm.omega_r = sqrt(prefactor * sum);
    wm.omega_i = 0.0;

    if (wm.k > 0.0 && wm.omega_r > 0.0) {
        wm.phase_velocity = wm.omega_r / wm.k;
    } else {
        wm.phase_velocity = 0.0;
    }
    wm.group_velocity = 0.0;

    return wm;
}

/* ================================================================
 * L4 — Magnetized Dust Modes
 * ================================================================ */

double dust_cyclotron_frequency(double Z_d, double B, double m_d)
{
    /* omega_cd = Z_d * e * B / m_d
     *
     * Dust cyclotron frequency is extremely small due to large m_d:
     * For Z_d=2000, B=0.1 T, m_d=10^-15 kg:
     * omega_cd ≈ 32 rad/s ≈ 5 Hz.
     *
     * Compare: omega_ce ≈ 1.8 GHz, omega_ci ≈ 2.4 MHz.
     *
     * This means dust is essentially unmagnetized in most
     * laboratory plasmas (omega_cd << nu_dn ~ 100 Hz). */
    if (m_d <= 0.0) return 0.0;
    return fabs(Z_d) * DUSTY_EC * B / m_d;
}

double dust_lower_hybrid_frequency(
    double Z_d, double B, double m_d, double m_i)
{
    /* Dust lower hybrid: geometric mean of dust and ion cyclotron.
     *
     * omega_LHd = sqrt(omega_cd * omega_ci)
     *
     * This mode is relevant when both ions and dust are magnetized,
     * which requires very strong B fields or very small grains.
     *
     * For dust to be magnetized: omega_cd > nu_dn
     * → B > m_d * nu_dn / (Z_d * e) ~ 0.3 T for micron grains. */
    double omega_cd = dust_cyclotron_frequency(Z_d, B, m_d);
    double omega_ci = DUSTY_EC * B / m_i;

    if (omega_cd <= 0.0 || omega_ci <= 0.0) return 0.0;
    return sqrt(omega_cd * omega_ci);
}

WaveMode dust_cyclotron_wave_dispersion(
    double k, double c_da, double lambda_D, double omega_cd)
{
    /* Electrostatic dust-cyclotron wave:
     *
     * omega^2 = omega_cd^2 + k^2 * c_da^2 / (1 + k^2 * lambda_D^2)
     *
     * This is the dust analog of the electrostatic ion cyclotron wave.
     * At k → 0: omega → omega_cd
     * At k → ∞: omega → omega_pd (if omega_cd < omega_pd)
     *
     * For most lab plasmas, omega_cd << omega_pd,
     * so the cyclotron branch is essentially flat. */
    WaveMode wm = dust_acoustic_wave_dispersion(k, c_da, lambda_D);
    wm.omega_r = sqrt(omega_cd * omega_cd + wm.omega_r * wm.omega_r);
    if (k > 0.0) wm.phase_velocity = wm.omega_r / k;

    return wm;
}
