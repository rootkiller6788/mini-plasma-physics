/**
 * plasma_parameters.c — Fundamental Plasma Parameter Calculations
 *
 * Every function implements one well-defined plasma physics quantity.
 * References: Goldston & Rutherford (1995), Chen (1984), MIT 22.611.
 */
#include "../include/space_plasma.h"
#include "../include/plasma_parameters.h"
#include <math.h>
#include <stdio.h>

/*──────────────────────────────────────────────────────────────────────
 * Debye Length & Screening
 *──────────────────────────────────────────────────────────────────────*/

double sp_debye_length_e(double n_e, double T_e) {
    /* λ_De = sqrt(ε₀ k_B T_e / (n_e e²))
     * Goldston & Rutherford Eq. (2.1) */
    if (n_e <= 0.0 || T_e <= 0.0) return 0.0;
    return sqrt(SP_EPS0 * SP_KB * T_e / (n_e * SP_EC * SP_EC));
}

double sp_debye_length_i(double n_i, double T_i, double Z) {
    /* λ_Di = sqrt(ε₀ k_B T_i / (n_i Z² e²)) */
    if (n_i <= 0.0 || T_i <= 0.0 || Z <= 0.0) return 0.0;
    double Ze = Z * SP_EC;
    return sqrt(SP_EPS0 * SP_KB * T_i / (n_i * Ze * Ze));
}

double sp_debye_length_total(double n_e, double T_e, double n_i, double T_i, double Z) {
    /* λ_D⁻² = λ_De⁻² + λ_Di⁻² */
    double ld_e = sp_debye_length_e(n_e, T_e);
    double ld_i = sp_debye_length_i(n_i, T_i, Z);
    if (ld_e <= 0.0 && ld_i <= 0.0) return 0.0;
    if (ld_e <= 0.0) return ld_i;
    if (ld_i <= 0.0) return ld_e;
    return 1.0 / sqrt(1.0/(ld_e*ld_e) + 1.0/(ld_i*ld_i));
}

double sp_debye_number(double n_e, double lambda_D) {
    /* N_D = (4π/3) n_e λ_D³
     * Collective behavior requires N_D ≫ 1 */
    if (n_e <= 0.0 || lambda_D <= 0.0) return 0.0;
    return (4.0 * M_PI / 3.0) * n_e * lambda_D * lambda_D * lambda_D;
}

double sp_coulomb_coupling(double n, double T, double Z) {
    /* Γ = Z² e² / (4π ε₀ a k_B T)
     * a = (3/(4πn))^(1/3)  ← Wigner-Seitz radius
     * Γ > 1: strongly coupled, Γ < 1: weakly coupled */
    if (n <= 0.0 || T <= 0.0 || Z <= 0.0) return 0.0;
    double a = cbrt(3.0 / (4.0 * M_PI * n));
    double Ze = Z * SP_EC;
    return (Ze * Ze) / (4.0 * M_PI * SP_EPS0 * a * SP_KB * T);
}

/*──────────────────────────────────────────────────────────────────────
 * Plasma Frequency
 *──────────────────────────────────────────────────────────────────────*/

double sp_plasma_freq_e(double n_e) {
    /* ω_pe = sqrt(n_e e² / (ε₀ m_e))
     * Goldston & Rutherford Eq. (2.4) */
    if (n_e <= 0.0) return 0.0;
    return sqrt(n_e * SP_EC * SP_EC / (SP_EPS0 * SP_ME));
}

double sp_plasma_freq_i(double n_i, double Z, double m_i) {
    /* ω_pi = sqrt(n_i Z² e² / (ε₀ m_i)) */
    if (n_i <= 0.0 || Z <= 0.0 || m_i <= 0.0) return 0.0;
    double Ze = Z * SP_EC;
    return sqrt(n_i * Ze * Ze / (SP_EPS0 * m_i));
}

double sp_lower_hybrid_freq(double omega_ci, double omega_ce, double omega_pi) {
    /* ω_LH⁻² = (ω_ci ω_ce)⁻¹ + ω_pi⁻²
     * Approximation when ω_ci ≪ ω_LH ≪ ω_ce */
    if (omega_ci <= 0.0 || omega_ce <= 0.0 || omega_pi <= 0.0) return 0.0;
    double inv2 = 1.0/(omega_ci * omega_ce) + 1.0/(omega_pi * omega_pi);
    return 1.0 / sqrt(inv2);
}

double sp_upper_hybrid_freq(double omega_pe, double omega_ce) {
    /* ω_UH² = ω_pe² + ω_ce² */
    return sqrt(omega_pe*omega_pe + omega_ce*omega_ce);
}

/*──────────────────────────────────────────────────────────────────────
 * Gyrofrequency & Larmor Radius
 *──────────────────────────────────────────────────────────────────────*/

double sp_gyrofreq(double B, double m, double q) {
    /* ω_c = q B / m  (Goldston & Rutherford Eq. 2.6)
     * q is magnitude of charge (positive) */
    if (B <= 0.0 || m <= 0.0 || q <= 0.0) return 0.0;
    return q * B / m;
}

double sp_larmor_radius(double v_perp, double B, double m, double q) {
    /* r_L = m v_perp / (q B)  (Goldston & Rutherford Eq. 2.7) */
    if (B <= 0.0 || m <= 0.0 || q <= 0.0) return 0.0;
    return m * v_perp / (q * B);
}

double sp_gyroperiod(double B, double m, double q) {
    /* T_c = 2π / ω_c */
    double wc = sp_gyrofreq(B, m, q);
    if (wc <= 0.0) return 0.0;
    return 2.0 * M_PI / wc;
}

/*──────────────────────────────────────────────────────────────────────
 * Alfven Velocity & Sound Speed
 *──────────────────────────────────────────────────────────────────────*/

double sp_alfven_speed(double B, double rho) {
    /* v_A = B / sqrt(μ₀ ρ)
     * Kivelson & Russell Eq. (4.23) */
    if (rho <= 0.0) return 0.0;
    return B / sqrt(SP_MU0 * rho);
}

double sp_sound_speed(double T_e, double m_i, double Z, double gamma) {
    /* c_s = sqrt(γ Z k_B T_e / m_i)
     * For ion acoustic waves with cold ions.
     * γ = 1 for isothermal, 5/3 for adiabatic electrons */
    if (T_e <= 0.0 || m_i <= 0.0 || Z <= 0.0 || gamma <= 0.0) return 0.0;
    return sqrt(gamma * Z * SP_KB * T_e / m_i);
}

double sp_magnetosonic_speed(double v_A, double c_s) {
    /* v_ms = sqrt(v_A² + c_s²) */
    return sqrt(v_A*v_A + c_s*c_s);
}

/*──────────────────────────────────────────────────────────────────────
 * Plasma Beta
 *──────────────────────────────────────────────────────────────────────*/

double sp_plasma_beta(double p, double B) {
    /* β = 2 μ₀ p / B²
     * Kivelson & Russell Eq. (4.34) */
    if (B <= 0.0) return 0.0;
    return 2.0 * SP_MU0 * p / (B * B);
}

void sp_plasma_beta_components(double n_e, double T_e, double n_i, double T_i,
                               double B, double *beta_e, double *beta_i) {
    /* β_e = 2 μ₀ n_e k_B T_e / B²
     * β_i = 2 μ₀ n_i k_B T_i / B² */
    double B2_inv = 0.0;
    if (B > 0.0) B2_inv = 2.0 * SP_MU0 / (B * B);
    *beta_e = B2_inv * n_e * SP_KB * T_e;
    *beta_i = B2_inv * n_i * SP_KB * T_i;
}

/*──────────────────────────────────────────────────────────────────────
 * Collision Frequencies & Transport
 *──────────────────────────────────────────────────────────────────────*/

double sp_coulomb_logarithm(double n_e, double T_e) {
    /* Coulomb logarithm: lnΛ
     * For classical plasma: lnΛ ≈ ln(λ_D / b_min)
     * Practical formula (NRL formulary):
     *   T_e < 10 eV: lnΛ = 23 - ½ ln(n_e[cm⁻³]) + ³⁄₂ ln(T_e[eV])
     *   T_e > 10 eV: lnΛ = 24 - ½ ln(n_e[cm⁻³]) + ln(T_e[eV])
     *
     * Here n_e in [m⁻³], T_e in [K].
     * Convert: n[cm⁻³] = n[m⁻³] * 1e-6, T[eV] = k_B T / e
     */
    if (n_e <= 0.0 || T_e <= 0.0) return 0.0;
    double n_cm3 = n_e * 1.0e-6;
    double T_eV = SP_KB * T_e / SP_EC;
    double lnLambda;
    if (T_eV < 10.0) {
        lnLambda = 23.0 - 0.5 * log(n_cm3) + 1.5 * log(T_eV);
    } else {
        lnLambda = 24.0 - 0.5 * log(n_cm3) + log(T_eV);
    }
    /* Clamp to physically reasonable range */
    if (lnLambda < 2.0) lnLambda = 2.0;
    if (lnLambda > 40.0) lnLambda = 40.0;
    return lnLambda;
}

double sp_collision_freq_ei(double n_e, double T_e, double Z, double lnLambda) {
    /* ν_ei = (4√(2π)/3) n_i Z² e⁴ lnΛ / (√(m_e) (k_B T_e)^(3/2))
     * Using n_i ≈ n_e/Z for quasi-neutrality:
     *   ν_ei ≈ 2.9e-12 * n_e * Z * lnΛ * (T_e[eV])^(-3/2) [s⁻¹]
     * Input n_e in [m⁻³], T_e in [K].
     */
    if (n_e <= 0.0 || T_e <= 0.0 || Z <= 0.0) return 0.0;
    /* NRL Plasma Formulary (2018): nu[s^-1] ~ 2.91e-6 * n_cm3 * Z * lnL * T_eV^(-1.5)
     * Full MKS formula: nu = sqrt(2)*n*Z^2*e^4*lnL / (12*pi*eps0^2*sqrt(me)*(kTe)^(3/2)) */
    double kTe = SP_KB * T_e;
    double coeff = sqrt(2.0) * n_e * Z * Z * SP_EC * SP_EC * SP_EC * SP_EC * lnLambda;
    coeff /= (12.0 * M_PI * SP_EPS0 * SP_EPS0 * sqrt(SP_ME) * pow(kTe, 1.5));
    return coeff;
}

double sp_mean_free_path(double v_th, double nu_ei) {
    /* λ_mfp = v_th / ν_ei */
    if (nu_ei <= 0.0) return INFINITY;
    return v_th / nu_ei;
}

double sp_spitzer_resistivity(double n_e, double T_e, double Z, double lnL) {
    /* Spitzer resistivity: η = m_e ν_ei / (n_e e²)
     * η [Ω·m] ≈ 5.2e-5 * Z * lnΛ * (T_e[eV])^(-3/2)
     * Standard result: η_∥ (Spitzer) independent of density to first order */
    if (n_e <= 0.0 || T_e <= 0.0) return 0.0;
    double nu_ei = sp_collision_freq_ei(n_e, T_e, Z, lnL);
    return SP_ME * nu_ei / (n_e * SP_EC * SP_EC);
}

/*──────────────────────────────────────────────────────────────────────
 * Dimensionless Numbers
 *──────────────────────────────────────────────────────────────────────*/

double sp_lundquist_number(double L, double v_A, double eta) {
    /* S = μ₀ L v_A / η
     * For space plasmas S ~ 1e6-1e14 → nearly ideal MHD */
    if (eta <= 0.0) return INFINITY;
    return SP_MU0 * L * v_A / eta;
}

double sp_magnetic_reynolds(double L, double v, double eta) {
    /* Rm = μ₀ L v / η */
    if (eta <= 0.0) return INFINITY;
    return SP_MU0 * L * v / eta;
}

double sp_ion_inertial_length(double n_i, double m_i, double Z) {
    /* d_i = c / ω_pi = sqrt(ε₀ m_i / (n_i Z² e²))
     * Important for Hall MHD: ion and electron motions decouple below d_i */
    if (n_i <= 0.0 || m_i <= 0.0 || Z <= 0.0) return 0.0;
    double Ze = Z * SP_EC;
    return SP_C * sqrt(SP_EPS0) / sqrt(n_i * Ze * Ze / SP_EPS0);
    /* More directly: d_i = c / ω_pi = sqrt(m_i/(μ₀ n_i Z² e²)) */
}

double sp_electron_inertial_length(double n_e) {
    /* d_e = c / ω_pe = sqrt(m_e / (μ₀ n_e e²)) */
    if (n_e <= 0.0) return 0.0;
    return SP_C / sp_plasma_freq_e(n_e);
}

/*──────────────────────────────────────────────────────────────────────
 * Derived Plasma State
 *──────────────────────────────────────────────────────────────────────*/

void sp_compute_plasma_state(plasma_state_t *state) {
    /* Compute derived quantities:
     * p = n k_B T, magnetized flag */
    if (!state) return;
    state->p_e = state->n_e * SP_KB * state->T_e;
    state->p_i = state->n_i * SP_KB * state->T_i;

    /* Check if magnetized: need collision time and gyroperiod */
    double Bmag = SP_MAG3(state->B);
    if (Bmag > 0.0 && state->n_e > 0.0 && state->T_e > 0.0) {
        double omega_ci = sp_gyrofreq(Bmag, SP_MP, state->Z * SP_EC);
        double lnL = sp_coulomb_logarithm(state->n_e, state->T_e);
        double nu_ii = sp_collision_freq_ei(state->n_e, state->T_e, state->Z, lnL);
        /* Magnetized if gyroperiod ≪ collision time: ω_c / ν > 1 */
        state->is_magnetized = (nu_ii > 0.0 && omega_ci / nu_ii > 1.0) ? 1 : 0;
    } else {
        state->is_magnetized = 0;
    }
}

double sp_plasma_freq_e_hz(double n_e) {
    return sp_plasma_freq_e(n_e) / (2.0 * M_PI);
}

double sp_plasma_freq_i_hz(double n_i, double Z, double m_i) {
    return sp_plasma_freq_i(n_i, Z, m_i) / (2.0 * M_PI);
}

double sp_gyrofreq_hz(double B, double m, double q) {
    return sp_gyrofreq(B, m, q) / (2.0 * M_PI);
}
