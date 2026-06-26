#ifndef DUSTY_WAVES_H
#define DUSTY_WAVES_H
/**
 * @file  dusty_waves.h
 * @brief Wave modes in dusty plasmas — dispersion relations.
 *
 * Dusty plasmas support several novel wave modes beyond the usual
 * electron and ion plasma waves.
 *
 * Key modes:
 *   DAW  (Dust-Acoustic Wave)      — ultra-low frequency, driven by
 *                                     electron pressure + dust inertia.
 *   DIAW (Dust-Ion-Acoustic Wave)  — faster, ions provide inertia.
 *   DLW  (Dust-Lattice Wave)       — in strongly coupled crystals.
 *   DEDW (Dust-Electron-Drift Wave) — with E x B drifts.
 *
 * Refs:
 *   Rao, Shukla & Yu (1990), Planet. Space Sci. 38, 543 — DAW
 *   Shukla & Silin (1992), Phys. Scr. 45, 504 — DIAW
 *   Melandso (1996), Phys. Plasmas 3, 3890 — DLW
 *   Piel (2010), "Plasma Physics" Ch. 11 — pedagogical overview
 */
#include "dusty_plasma.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ================================================================
 * L4 — Dust-Acoustic Wave (DAW)
 * ================================================================ */

/**
 * @brief Dust-acoustic wave dispersion relation.
 *
 * omega^2 = (k^2 * c_da^2) / (1 + k^2 * lambda_D^2)
 *
 * where c_da = sqrt(Z_d * k_B * T_e / m_d) is the dust-acoustic speed
 * and lambda_D is the total Debye length.
 *
 * Long wavelength (k lambda_D << 1): omega ≈ k * c_da  (sound-like)
 * Short wavelength (k lambda_D >> 1): omega → omega_pd (constant)
 *
 * Complexity: O(1).
 */
WaveMode dust_acoustic_wave_dispersion(
    double k, double c_da, double lambda_D);

/**
 * @brief DAW dispersion relation including ion temperature effects.
 *
 * omega^2 = k^2 * (Z_d * k_B * T_eff / m_d) / (1 + k^2 * lambda_D^2)
 *
 * where T_eff = T_i * T_e / (T_i + Z_d * T_e * n_e/n_i)
 *
 * Ref: Merlino et al. (1998), Phys. Plasmas 5, 1607
 */
WaveMode dust_acoustic_wave_with_ion_temp(
    double k, double Z_d, double T_e, double T_i, double m_d,
    double n_e, double n_i, double lambda_D);

/**
 * @brief DAW Landau damping rate for Maxwellian dust.
 *
 * gamma / omega = -sqrt(pi/8) * (omega_pd / (k * v_thd))^3
 *                 * exp(-omega^2 / (2 * k^2 * v_thd^2))
 *
 * This is the collisionless damping due to resonant dust particles.
 * In practice, DAW is heavily damped unless T_d << Z_d * T_e.
 *
 * Ref: Rosenberg (1993), Planet. Space Sci. 41, 229
 */
double dust_acoustic_wave_landau_damping(
    double k, double omega, double omega_pd, double v_thd);

/**
 * @brief DAW collisional damping rate due to dust-neutral collisions.
 *
 * gamma_coll = nu_dn / 2
 *
 * where nu_dn is the dust-neutral collision frequency.
 * This is often the dominant damping mechanism in laboratory plasmas.
 */
double dust_acoustic_wave_collisional_damping(double nu_dn);

/**
 * @brief Full DAW complex frequency including both dampings.
 *
 * omega_complex = omega_r + i * gamma_total
 * gamma_total = gamma_Landau + gamma_coll
 */
WaveMode dust_acoustic_wave_full(
    double k, double c_da, double lambda_D, double omega_pd,
    double v_thd, double nu_dn);

/* ================================================================
 * L4 — Dust-Ion-Acoustic Wave (DIAW)
 * ================================================================ */

/**
 * @brief Dust-ion-acoustic wave dispersion relation.
 *
 * omega^2 = (k^2 * c_s^2) / (1 + k^2 * lambda_De^2)
 *
 * where c_s = sqrt(k_B * T_e / m_i) and lambda_De is the electron
 * Debye length. Dust is considered stationary on this timescale.
 *
 * DIAW is the usual ion-acoustic wave modified by dust.
 */
WaveMode dust_ion_acoustic_wave_dispersion(
    double k, double c_s, double lambda_De);

/**
 * @brief DIAW with dust charge fluctuation effects.
 *
 * When Z_d fluctuates, an additional damping/growth channel opens.
 * Ref: Varma, Shukla & Krishan (1993), Phys. Rev. E 47, 3612
 */
WaveMode dust_ion_acoustic_wave_charge_fluctuation(
    double k, double c_s, double lambda_De,
    double nd_over_ni, double tau_charge);

/* ================================================================
 * L6 — Dust Lattice Wave (DLW) in 1D Chain
 * ================================================================ */

/**
 * @brief 1D dust lattice wave dispersion (longitudinal mode).
 *
 * omega_L^2(k) = (Q_d^2 / (4 pi eps0 m_d a^3))
 *                * sum_{n=1}^{N} (2/n^3) * (1 + kappa * n)
 *                * exp(-kappa * n) * (1 - cos(k * n * a))
 *
 * where a is the inter-particle spacing and kappa = a / lambda_D.
 *
 * In the nearest-neighbor approximation (kappa >> 1):
 * omega_L^2 ≈ 2 * Omega_0^2 * (1 - cos(k*a))
 *
 * where Omega_0^2 = Q_d^2 / (4 * pi * eps0 * m_d * a^3) * (1 + kappa) * exp(-kappa)
 *
 * Complexity: O(N) per call, N ≈ 20 gives converged sum.
 */
WaveMode dust_lattice_wave_1d_longitudinal(
    double k, double a, double Q_d, double m_d,
    double kappa, int n_shells);

/**
 * @brief 1D dust lattice wave dispersion (transverse mode).
 *
 * omega_T^2(k) = (Q_d^2 / (4 pi eps0 m_d a^3))
 *                * sum_{n=1}^{N} (1/n^3) * (1 + kappa * n)
 *                * exp(-kappa * n) * (1 - cos(k * n * a))
 *
 * Transverse mode has lower frequency than longitudinal due to
 * weaker restoring force (1/n^3 vs 2/n^3 factor).
 */
WaveMode dust_lattice_wave_1d_transverse(
    double k, double a, double Q_d, double m_d,
    double kappa, int n_shells);

/**
 * @brief 2D hexagonal dust lattice wave dispersion (out-of-plane mode).
 *
 * omega_O^2(k) = Omega_0^2 * sum_{nn} (1 + kappa * r_nn)
 *                * exp(-kappa * r_nn) * (1 - cos(k · r_nn))
 *                * (1 / r_nn^3)
 *
 * where sum runs over 6 nearest neighbors at positions r_nn.
 * This is the out-of-plane (vertical) mode.
 */
WaveMode dust_lattice_wave_2d_out_of_plane(
    double kx, double ky, double a, double Q_d, double m_d, double kappa);

/**
 * @brief 2D in-plane longitudinal dust lattice wave.
 */
WaveMode dust_lattice_wave_2d_in_plane_longitudinal(
    double kx, double ky, double a, double Q_d, double m_d, double kappa);

/* ================================================================
 * L4 — Magnetized Dust Modes
 * ================================================================ */

/**
 * @brief Dust cyclotron frequency.
 * omega_cd = Z_d * e * B / m_d
 *
 * For micron-sized dust (m_d ~ 10^-15 kg), omega_cd ~ 10^-4 B,
 * extremely small compared to ion/electron cyclotron frequencies.
 */
double dust_cyclotron_frequency(double Z_d, double B, double m_d);

/**
 * @brief Dust lower-hybrid frequency.
 * omega_LHd ≈ sqrt(omega_cd * omega_ci)
 *
 * This mode couples dust and ion dynamics across B field.
 */
double dust_lower_hybrid_frequency(
    double Z_d, double B, double m_d, double m_i);

/**
 * @brief Electrostatic dust-cyclotron wave dispersion.
 *
 * omega^2 = omega_cd^2 + k^2 * c_da^2 / (1 + k^2 * lambda_D^2)
 *
 * Dust cyclotron motion modifies the DAW at very low frequencies.
 */
WaveMode dust_cyclotron_wave_dispersion(
    double k, double c_da, double lambda_D, double omega_cd);

#ifdef __cplusplus
}
#endif

#endif /* DUSTY_WAVES_H */
