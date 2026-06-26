/**
 * plasma_parameters.h — Plasma Parameter Calculations
 *
 * References:
 *   Goldston & Rutherford §2 (Debye shielding, plasma frequency)
 *   Kivelson & Russell §3 (Space plasma parameters)
 *   Chen "Introduction to Plasma Physics" §1-§4
 *   MIT 22.611 Lecture 2-3
 *
 * Knowledge coverage: L1 (definitions), L2 (core concepts),
 *   L4 (fundamental relations), L5 (computational)
 */

#ifndef PLASMA_PARAMETERS_H
#define PLASMA_PARAMETERS_H

#include "space_plasma.h"

/*──────────────────────────────────────────────────────────────────────
 * L1: Debye Length & Plasma Parameter
 *──────────────────────────────────────────────────────────────────────*/

/**
 * @brief Electron Debye length
 * λ_De = sqrt(ε₀ k_B T_e / (n_e e²))
 *
 * Debye length measures the distance over which charge screening occurs.
 * For λ_D ≪ L (system size), plasma exhibits collective behavior.
 *
 * @param n_e   Electron density [m⁻³]
 * @param T_e   Electron temperature [K]
 * @return      Debye length [m]
 *
 * Reference: Goldston & Rutherford Eq. (2.1)
 */
double sp_debye_length_e(double n_e, double T_e);

/**
 * @brief Ion Debye length
 * λ_Di = sqrt(ε₀ k_B T_i / (n_i Z² e²))
 *
 * @param n_i   Ion density [m⁻³]
 * @param T_i   Ion temperature [K]
 * @param Z     Ion charge state
 * @return      Ion Debye length [m]
 */
double sp_debye_length_i(double n_i, double T_i, double Z);

/**
 * @brief Total Debye length
 * λ_D⁻² = λ_De⁻² + λ_Di⁻²
 *
 * @param n_e   Electron density [m⁻³]
 * @param T_e   Electron temperature [K]
 * @param n_i   Ion density [m⁻³]
 * @param T_i   Ion temperature [K]
 * @param Z     Ion charge state
 * @return      Total Debye length [m]
 */
double sp_debye_length_total(double n_e, double T_e, double n_i, double T_i, double Z);

/**
 * @brief Number of particles in Debye sphere
 * N_D = (4π/3) n_e λ_D³
 *
 * Plasma parameter: N_D ≫ 1 indicates weakly coupled plasma.
 * N_D ≪ 1 indicates strongly coupled plasma.
 *
 * @param n_e       Electron density [m⁻³]
 * @param lambda_D  Debye length [m]
 * @return          N_D (dimensionless)
 */
double sp_debye_number(double n_e, double lambda_D);

/**
 * @brief Coulomb coupling parameter
 * Γ = e² / (4π ε₀ a k_B T) where a = (4πn/3)^(-1/3)
 *
 * Γ < 1: weakly coupled (ideal plasma)
 * Γ > 1: strongly coupled (non-ideal, liquid-like)
 *
 * @param n   Density [m⁻³]
 * @param T   Temperature [K]
 * @param Z   Charge state
 * @return    Γ (dimensionless)
 */
double sp_coulomb_coupling(double n, double T, double Z);

/*──────────────────────────────────────────────────────────────────────
 * L1: Plasma Frequency
 *──────────────────────────────────────────────────────────────────────*/

/**
 * @brief Electron plasma frequency
 * ω_pe = sqrt(n_e e² / (ε₀ m_e))
 *
 * Fundamental timescale of electron response. Electromagnetic waves
 * with ω < ω_pe cannot propagate through plasma.
 *
 * @param n_e   Electron density [m⁻³]
 * @return      Angular plasma frequency [rad/s]
 *
 * Reference: Goldston & Rutherford Eq. (2.4)
 */
double sp_plasma_freq_e(double n_e);

/**
 * @brief Ion plasma frequency
 * ω_pi = sqrt(n_i Z² e² / (ε₀ m_i))
 *
 * @param n_i   Ion density [m⁻³]
 * @param Z     Ion charge state
 * @param m_i   Ion mass [kg] (use SP_MP for protons)
 * @return      Angular ion plasma frequency [rad/s]
 */
double sp_plasma_freq_i(double n_i, double Z, double m_i);

/**
 * @brief Lower hybrid frequency
 * ω_LH⁻² = (ω_ci ω_ce)⁻¹ + ω_pi⁻²
 *
 * Important for lower hybrid waves and current drive.
 *
 * @param omega_ci  Ion gyrofrequency [rad/s]
 * @param omega_ce  Electron gyrofrequency [rad/s]
 * @param omega_pi  Ion plasma frequency [rad/s]
 * @return          Lower hybrid frequency [rad/s]
 */
double sp_lower_hybrid_freq(double omega_ci, double omega_ce, double omega_pi);

/**
 * @brief Upper hybrid frequency
 * ω_UH² = ω_pe² + ω_ce²
 *
 * @param omega_pe  Electron plasma frequency [rad/s]
 * @param omega_ce  Electron gyrofrequency [rad/s]
 * @return          Upper hybrid frequency [rad/s]
 */
double sp_upper_hybrid_freq(double omega_pe, double omega_ce);

/*──────────────────────────────────────────────────────────────────────
 * L1: Gyrofrequency & Larmor Radius
 *──────────────────────────────────────────────────────────────────────*/

/**
 * @brief Gyrofrequency (cyclotron frequency)
 * ω_c = |q| B / m
 *
 * @param B    Magnetic field magnitude [T]
 * @param m    Particle mass [kg]
 * @param q    Particle charge magnitude [C] (positive)
 * @return      Gyrofrequency [rad/s]
 *
 * Reference: Goldston & Rutherford Eq. (2.6)
 */
double sp_gyrofreq(double B, double m, double q);

/**
 * @brief Larmor (gyro) radius
 * r_L = m v_perp / (|q| B)
 *
 * @param v_perp  Perpendicular velocity [m/s]
 * @param B       Magnetic field magnitude [T]
 * @param m       Particle mass [kg]
 * @param q       Particle charge magnitude [C] (positive)
 * @return        Larmor radius [m]
 *
 * Reference: Goldston & Rutherford Eq. (2.7)
 */
double sp_larmor_radius(double v_perp, double B, double m, double q);

/**
 * @brief Compute gyroperiod
 * T_c = 2π / ω_c
 */
double sp_gyroperiod(double B, double m, double q);

/*──────────────────────────────────────────────────────────────────────
 * L2: Alfven Velocity & Sound Speed
 *──────────────────────────────────────────────────────────────────────*/

/**
 * @brief Alfven velocity
 * v_A = B / sqrt(μ₀ ρ)
 *
 * Speed of transverse MHD waves along magnetic field lines.
 * Energy propagates at v_A in the direction of B.
 *
 * @param B    Magnetic field magnitude [T]
 * @param rho  Mass density [kg/m³] = n_i m_i + n_e m_e ≈ n_i m_i
 * @return      Alfven velocity [m/s]
 *
 * Reference: Kivelson & Russell Eq. (4.23)
 */
double sp_alfven_speed(double B, double rho);

/**
 * @brief Ion sound speed
 * c_s = sqrt(γ Z k_B T_e / m_i)
 *
 * For ion acoustic waves, with cold ions (T_i ≪ T_e).
 * γ is the ratio of specific heats (typically 5/3 for adiabatic,
 * 1 for isothermal).
 *
 * @param T_e   Electron temperature [K]
 * @param m_i   Ion mass [kg]
 * @param Z     Ion charge state
 * @param gamma Polytropic index (1.0 = isothermal, 5/3 = adiabatic)
 * @return      Sound speed [m/s]
 */
double sp_sound_speed(double T_e, double m_i, double Z, double gamma);

/**
 * @brief Magnetosonic speed
 * v_ms = sqrt(v_A² + c_s²)
 *
 * @param v_A  Alfven speed [m/s]
 * @param c_s  Sound speed [m/s]
 * @return     Magnetosonic speed [m/s]
 */
double sp_magnetosonic_speed(double v_A, double c_s);

/*──────────────────────────────────────────────────────────────────────
 * L2: Plasma Beta
 *──────────────────────────────────────────────────────────────────────*/

/**
 * @brief Plasma beta
 * β = 2 μ₀ p / B²
 *
 * Ratio of thermal pressure to magnetic pressure.
 * β < 1: magnetically dominated (corona, magnetosphere)
 * β > 1: kinetically dominated (solar wind at 1 AU, fusion plasmas)
 *
 * @param p   Total pressure [Pa]
 * @param B   Magnetic field magnitude [T]
 * @return    β (dimensionless)
 *
 * Reference: Kivelson & Russell Eq. (4.34)
 */
double sp_plasma_beta(double p, double B);

/**
 * @brief Electron and ion beta separately
 * β_e = 2 μ₀ n_e k_B T_e / B²
 * β_i = 2 μ₀ n_i k_B T_i / B²
 *
 * @param n_e   Electron density [m⁻³]
 * @param T_e   Electron temperature [K]
 * @param n_i   Ion density [m⁻³]
 * @param T_i   Ion temperature [K]
 * @param B     Magnetic field [T]
 * @param beta_e Output: electron beta
 * @param beta_i Output: ion beta
 */
void sp_plasma_beta_components(double n_e, double T_e, double n_i, double T_i,
                               double B, double *beta_e, double *beta_i);

/*──────────────────────────────────────────────────────────────────────
 * L2: Collision Frequencies
 *──────────────────────────────────────────────────────────────────────*/

/**
 * @brief Electron-ion collision frequency (Spitzer)
 * ν_ei = (4√(2π)/(3)) n_i Z² e⁴ lnΛ / (√(m_e) (k_B T_e)^(3/2))
 *
 * Using simplified form: ν_ei [s⁻¹] ≈ 2.9e-12 * n[cm⁻³] * lnΛ * T_e[eV]^(-3/2)
 *
 * @param n_e      Electron density [m⁻³]
 * @param T_e      Electron temperature [K]
 * @param Z        Ion charge state
 * @param lnLambda Coulomb logarithm
 * @return         Collision frequency [s⁻¹]
 */
double sp_collision_freq_ei(double n_e, double T_e, double Z, double lnLambda);

/**
 * @brief Coulomb logarithm
 * lnΛ = ln(λ_D / b_min)
 * where b_min = Z e² / (4π ε₀ m_e v_th²) for classical plasmas
 *
 * In practice: lnΛ ≈ 30 - ½ ln(n_e[cm⁻³]) + ³⁄₂ ln(T_e[eV])
 * for T_e < 10 eV: add ½ ln(T_e)
 *
 * @param n_e  Electron density [m⁻³]
 * @param T_e  Electron temperature [K]
 * @return     Coulomb logarithm (dimensionless)
 */
double sp_coulomb_logarithm(double n_e, double T_e);

/**
 * @brief Mean free path for electron-ion collisions
 * λ_mfp = v_th / ν_ei
 *
 * @param v_th  Electron thermal velocity [m/s]
 * @param nu_ei Electron-ion collision frequency [s⁻¹]
 * @return      Mean free path [m]
 */
double sp_mean_free_path(double v_th, double nu_ei);

/**
 * @brief Spitzer resistivity
 * η = m_e ν_ei / (n_e e²)
 *
 * Plasma resistivity due to Coulomb collisions.
 * η ∝ T_e^(-3/2), independent of density to first order.
 *
 * @param n_e   Electron density [m⁻³]
 * @param T_e   Electron temperature [K]
 * @param Z     Ion charge state
 * @param lnL   Coulomb logarithm
 * @return      Resistivity [Ω·m]
 */
double sp_spitzer_resistivity(double n_e, double T_e, double Z, double lnL);

/*──────────────────────────────────────────────────────────────────────
 * L2: Dimensionless Plasma Numbers
 *──────────────────────────────────────────────────────────────────────*/

/**
 * @brief Lundquist number
 * S = μ₀ L v_A / η
 *
 * Measures ratio of advective to diffusive magnetic flux transport.
 * S ≫ 1: ideal MHD regime.
 * S ~ 10⁶-10¹² in space plasmas.
 *
 * @param L    Characteristic length scale [m]
 * @param v_A  Alfven speed [m/s]
 * @param eta  Resistivity [Ω·m]
 * @return     Lundquist number (dimensionless)
 */
double sp_lundquist_number(double L, double v_A, double eta);

/**
 * @brief Magnetic Reynolds number
 * Rm = μ₀ L v / η
 *
 * Determines if magnetic field is frozen into plasma (Rm ≫ 1).
 *
 * @param L    Characteristic length [m]
 * @param v    Characteristic velocity [m/s]
 * @param eta  Resistivity [Ω·m]
 * @return     Magnetic Reynolds number (dimensionless)
 */
double sp_magnetic_reynolds(double L, double v, double eta);

/**
 * @brief Ion inertial length (skin depth)
 * d_i = c / ω_pi
 *
 * Scale below which ions decouple from electrons.
 * Critical for Hall MHD and reconnection.
 *
 * @param n_i   Ion density [m⁻³]
 * @param m_i   Ion mass [kg]
 * @param Z     Ion charge state
 * @return      Ion inertial length [m]
 */
double sp_ion_inertial_length(double n_i, double m_i, double Z);

/**
 * @brief Electron inertial length
 * d_e = c / ω_pe
 */
double sp_electron_inertial_length(double n_e);

/**
 * @brief Compute plasma state from fundamental parameters
 *
 * Fills plasma_state_t struct with derived quantities:
 *   p_e = n_e k_B T_e, p_i = n_i k_B T_i
 *   is_magnetized = (ω_ci * τ_ii > 1)
 *
 * @param state  Pointer to plasma_state_t to be filled
 */
void sp_compute_plasma_state(plasma_state_t *state);

/**
 * @brief Electron plasma frequency in Hz: f_pe = ω_pe / (2π)
 */
double sp_plasma_freq_e_hz(double n_e);

/**
 * @brief Ion plasma frequency in Hz: f_pi = ω_pi / (2π)
 */
double sp_plasma_freq_i_hz(double n_i, double Z, double m_i);

/**
 * @brief Gyrofrequency in Hz: f_c = ω_c / (2π)
 */
double sp_gyrofreq_hz(double B, double m, double q);

#endif /* PLASMA_PARAMETERS_H */
