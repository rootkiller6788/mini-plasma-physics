/**
 * fusion_heating.c — Plasma Heating and Current Drive
 *
 * Every function implements an independent physics knowledge point.
 * Refs: Stix "Waves in Plasmas" (1992), Wesson (2011),
 *       ITER Physics Basis (1999), Fisch RMP 59, 175 (1987)
 *
 * L2: NBI deposition, RF wave propagation, cyclotron resonance
 * L4: Alpha heating, ohmic heating, fusion gain
 * L5: NBI shine-through, orbit loss models
 * L7: ITER heating mix, DEMO heating requirements
 * L8: Current drive efficiency, helicity injection
 */

#include "include/fusion_plasma.h"
#include <math.h>
#include <assert.h>

/* ================================================================
 * L2: Neutral Beam Injection (NBI) — 4 functions
 * ================================================================ */

/**
 * nbi_stopping_cross_section — NBI stopping cross-section [m^2]
 *
 * Definition: sigma_stop = sum over plasma species of
 *          (Z_i * e^2 / (4*pi*epsilon0*E_beam))^2 * lnLambda * ...
 * Physics: cross-section for beam attenuation through
 *          ionization, charge exchange, and Coulomb collisions.
 * Simplified formula for D beam into D-T plasma.
 *
 * @param E_beam_keV  beam energy [keV]
 * @param ne          electron density [m^-3]
 * @param lnLambda    Coulomb logarithm
 * @return            effective stopping cross-section [m^2]
 */
double nbi_stopping_cross_section(double E_beam_keV, double ne, double lnLambda) {
    (void)ne;  /* density dependence negligible in simplified model */
    double E_J = E_beam_keV * 1e3 * E_CHARGE;
    double b90_sq = E_CHARGE * E_CHARGE / (4.0 * M_PI * EPSILON0 * E_J);
    b90_sq = b90_sq * b90_sq;
    return b90_sq * lnLambda * 4.0 * M_PI;
}

/**
 * nbi_shine_through_fraction — NBI shine-through fraction
 *
 * Definition: f_shine = exp(-a/lambda_stop)
 *          lambda_stop = 1/(n_e * sigma_stop)
 * Physics: fraction of neutral beam that passes through plasma
 *          without being ionized. Depends on beam energy and
 *          plasma line density.
 * Typical: f_shine < 0.05 for ITER (1 MeV D beams)
 *
 * @param a            plasma minor radius [m]
 * @param ne           electron density [m^-3]
 * @param sigma_stop   stopping cross-section [m^2]
 * @return             shine-through fraction
 */
double nbi_shine_through_fraction(double a, double ne, double sigma_stop) {
    double lambda_stop = 1.0 / (ne * sigma_stop);
    return exp(-a / lambda_stop);
}

/**
 * nbi_orbit_loss_fraction — NBI fast ion orbit loss fraction
 *
 * Definition: f_orbit ~ exp(-a/rho_fast) for co-injection
 * Physics: fast ions born near plasma edge may be on
 *          unconfined (loss) orbits. Depends on injection
 *          geometry and fast ion Larmor radius.
 * Typical: f_orbit ~ 0.05-0.15 for tangential co-injection
 *
 * @param a             minor radius [m]
 * @param rho_fast      fast ion Larmor radius [m]
 * @param tangency_radius  beam tangency radius / R
 * @return              orbit loss fraction
 */
double nbi_orbit_loss_fraction(double a, double rho_fast, double tangency_radius) {
    double rho_norm = rho_fast / a;
    double loss_factor = rho_norm / tangency_radius;
    if (loss_factor > 1.0) loss_factor = 1.0;
    return loss_factor;
}

/**
 * nbi_current_drive_efficiency — NBI current drive efficiency [A/W]
 *
 * Definition: eta_CD = I_CD / P_NBI
 *          ~ (n_e * R) scaling (Fisch, 1987)
 * Physics: fast ions from NBI can drive toroidal current.
 *          Efficiency peaks for off-axis injection with
 *          trapped particle effects.
 * ITER: eta_CD ~ 0.2-0.5 A/W (for 1 MeV negative ion NBI)
 *
 * @param ne_20       density [10^20 m^-3]
 * @param R           major radius [m]
 * @param Z_eff       effective charge
 * @return            current drive efficiency [A/W]
 */
double nbi_current_drive_efficiency(double ne_20, double R, double Z_eff) {
    double eta_0 = 0.3;  /* base efficiency factor */
    return eta_0 * R / (ne_20 * Z_eff);
}

/* ================================================================
 * L2: Radio Frequency (RF) Heating — 4 functions
 * ================================================================ */

/**
 * icrf_resonance_frequency — ICRF resonance frequency [Hz]
 *
 * Definition: f_ICRF = Omega_ci / (2*pi) = Z*e*B / (2*pi*mi)
 * Physics: Ion Cyclotron Range of Frequencies heating uses
 *          the fundamental ion cyclotron resonance.
 * For D at B=5T: f_ICRF ~ 38 MHz
 */
double icrf_resonance_frequency(double B, double mi_kg, double Z) {
    return ion_cyclotron_frequency(B, mi_kg, Z) / (2.0 * M_PI);
}

/**
 * ecrh_resonance_frequency — ECRH resonance frequency [Hz]
 *
 * Definition: f_ECRH = Omega_ce / (2*pi) = e*B / (2*pi*me)
 * Physics: Electron Cyclotron Resonance Heating uses the
 *          fundamental electron cyclotron resonance or
 *          2nd harmonic (X2) for core heating.
 * For B=5T: f_ECRH ~ 140 GHz (fundamental)
 */
double ecrh_resonance_frequency(double B) {
    return electron_cyclotron_frequency(B) / (2.0 * M_PI);
}

/**
 * lh_resonance_frequency — Lower Hybrid resonance frequency [Hz]
 *
 * Definition: f_LH ~ sqrt(Omega_ci * Omega_ce) / (2*pi)
 * Physics: Lower Hybrid frequency range. LH waves are
 *          efficient for current drive (Landau damping on
 *          electrons). Penetrate to plasma core.
 * For D at B=5T, ne=1e20: f_LH ~ 2-5 GHz
 */
double lh_resonance_frequency(double B, double mi_kg, double Z, double ne) {
    double omega_ci = ion_cyclotron_frequency(B, mi_kg, Z);
    double omega_ce = electron_cyclotron_frequency(B);
    double omega_pi = ion_plasma_frequency(ne, mi_kg, Z);
    double omega_lh_sq = (omega_ci * omega_ce * omega_pi * omega_pi) /
                          (omega_pi * omega_pi + omega_ce * omega_ce);
    return sqrt(omega_lh_sq) / (2.0 * M_PI);
}

/**
 * rf_power_deposition_profile — RF power deposition parametrization
 *
 * Simplified Gaussian deposition profile for RF heating.
 * P(r) = P_0 * exp(-(r - r_res)^2 / (2*Delta_r^2))
 *
 * @param r          radial position [m]
 * @param r_res      resonance location [m]
 * @param delta_r    deposition width [m]
 * @param P_total    total launched power [W]
 * @param a          plasma radius [m]
 * @return           local power density [W/m^3] (approximate)
 */
double rf_power_deposition(double r, double r_res, double delta_r,
                             double P_total, double a) {
    double arg = (r - r_res) / delta_r;
    double gaussian = exp(-0.5 * arg * arg);
    double norm = 1.0 / (sqrt(2.0 * M_PI) * delta_r * 2.0 * M_PI * M_PI);
    (void)a;
    return P_total * norm * gaussian;
}

/* ================================================================
 * L4: Alpha Particle Heating — 3 functions
 * ================================================================ */

/**
 * alpha_power_to_electrons — Alpha power to electrons [W]
 *
 * Definition: at E_alpha > E_crit, alpha energy goes to electrons.
 *          P_alpha_e = P_alpha * ln(E_alpha/E_crit)^2 /
 *                      (ln(E_alpha/E_crit)^2 + pi^2/4)
 * Physics: alpha slowing-down theory (Sivukhin, 1966).
 *          E_crit ~ 14.8*Te*A_alpha separates electron/ion heating.
 *
 * @param P_alpha   total alpha heating power [W]
 * @param Te_eV     electron temperature [eV]
 * @return          power to electron channel [W]
 */
double alpha_power_to_electrons(double P_alpha, double Te_eV) {
    double E_alpha_eV = E_ALPHA;
    double E_crit = alpha_critical_energy(Te_eV);
    if (E_crit >= E_alpha_eV) return 0.0;
    double ratio = E_alpha_eV / E_crit;
    double ln_ratio = log(ratio);
    double ln2 = ln_ratio * ln_ratio;
    double electron_fraction = ln2 / (ln2 + M_PI * M_PI / 4.0);
    return P_alpha * electron_fraction;
}

/**
 * alpha_power_to_ions — Alpha power to ions [W]
 *
 * Complement to alpha_power_to_electrons.
 * At E < E_crit, ion drag dominates.
 */
double alpha_power_to_ions(double P_alpha, double Te_eV) {
    double P_e = alpha_power_to_electrons(P_alpha, Te_eV);
    return P_alpha - P_e;
}

/**
 * alpha_ash_accumulation_time — Alpha ash (He) confinement time [s]
 *
 * Definition: tau_He* ~ tau_E / S_He
 *          S_He: helium recycling coefficient (>1)
 * Physics: helium ash from D-T reactions must be exhausted
 *          to prevent fuel dilution. Need tau_He-star / tau_E < 10-15
 *          for efficient fusion burn.
 *
 *  tau_E: energy confinement time [s]
 *  S_He:  helium recycling coefficient
 *  return: effective He particle confinement time [s]
 */
double alpha_ash_confinement_time(double tau_E, double S_He) {
    return tau_E / S_He;
}

/* ================================================================
 * L6: Heating System Integration — 2 functions
 * ================================================================ */

/**
 * heating_system_cost — Heating system cost estimate [arbitrary units]
 *
 * Simplified cost model for fusion heating systems.
 * NBI: most expensive (~4-5 USD/W)
 * ECRH: expensive (~3-4 USD/W, plus gyrotron development)
 * ICRF: moderate (~2-3 USD/W)
 * LH: cheapest (~1-2 USD/W)
 *
 * @param P_nbi    NBI power [W]
 * @param P_ecrh   ECRH power [W]
 * @param P_icrf   ICRF power [W]
 * @param P_lh     LH power [W]
 * @return         relative cost [arbitrary units]
 */
double heating_system_cost(double P_nbi, double P_ecrh, double P_icrf, double P_lh) {
    return 4.5 * P_nbi + 3.5 * P_ecrh + 2.5 * P_icrf + 1.5 * P_lh;
}

/**
 * optimum_heating_mix — Optimum heating mix for ITER-like parameters
 *
 * Based on ITER heating baseline: 33 MW NBI + 20 MW ECRH + 20 MW ICRF
 * This function computes the total and fills a HeatingSystem struct.
 *
 * @param hs  output: filled HeatingSystem structure
 */
void optimum_heating_mix(HeatingSystem *hs) {
    if (!hs) return;
    hs->P_nbi = 3.3e7;   /* 33 MW NBI (ITER baseline) */
    hs->P_icrf = 2.0e7;  /* 20 MW ICRF */
    hs->P_ecrh = 2.0e7;  /* 20 MW ECRH */
    hs->P_lh = 0.0;       /* LH not planned for ITER baseline */
    hs->P_ohmic = 1.0e6;  /* ~1 MW ohmic (negligible at high Te) */
    hs->P_alpha = 0.0;    /* filled by power balance calculation */
    hs->P_total = hs->P_nbi + hs->P_icrf + hs->P_ecrh + hs->P_lh + hs->P_ohmic;
}