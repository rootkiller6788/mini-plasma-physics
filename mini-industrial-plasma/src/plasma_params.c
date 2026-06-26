/**
 * plasma_params.c - Core Plasma Parameter Calculations
 *
 * Fundamental plasma physics computations:
 *   - Debye length, plasma frequency, collision frequencies
 *   - Quasineutrality, plasma beta, magnetization
 *   - Species initialization and reaction setup
 *   - Discharge state estimation from engineering parameters
 *
 * Reference: Lieberman §1-4, Goldston & Rutherford §1-3
 * Course: MIT 22.611, Caltech Ph 106, Stanford PHYSICS 370
 */

#include "plasma_types.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

/* ============================================================
 * L1: Fundamental Plasma Parameters
 * ============================================================ */

/*
 * Initialize a PlasmaState structure from basic parameters.
 *
 * Given n_e, T_e, and neutral gas pressure, compute all
 * dependent plasma parameters (Debye length, plasma frequency,
 * collision frequencies, etc.).
 *
 * This is the fundamental "plasma initialization" function
 * that bridges engineering inputs to plasma physics quantities.
 *
 * Complexity: O(1) — all analytical formulas
 */
void plasma_state_init(PlasmaState *ps, double n_e, double T_e_eV,
                       double T_i_K, double T_g_K, double pressure_pa,
                       double B_field_T)
{
    if (!ps) return;

    ps->n_e = n_e;
    ps->T_e = T_e_eV;
    ps->T_i = T_i_K;
    ps->T_g = T_g_K;
    ps->pressure = pressure_pa;
    ps->B_field = B_field_T;

    /* Neutral gas density from ideal gas law: n = p / (k_B * T_g) */
    if (T_g_K > 0.0) {
        ps->n_n = pressure_pa / (PLASMA_K_B * T_g_K);
    } else {
        ps->n_n = 0.0;
    }

    /* Ion density from quasineutrality: n_i = n_e (for singly charged ions) */
    ps->n_i = n_e;

    /* Debye length: lam_D = sqrt(eps0 * T_e[eV] / (n_e * e)) */
    ps->lam_D = debye_length_electron(n_e, T_e_eV);

    /* Electron plasma frequency: omega_pe = sqrt(n_e * e^2 / (eps0 * m_e)) */
    ps->omega_pe = plasma_freq_electron(n_e);

    /* Plasma potential (approximate from Bohm criterion and floating potential):
     * V_p ≈ V_f + T_e/2 * ln(m_i/(2*pi*m_e))
     * Using 0 as reference, V_p ~ T_e/2 * ln(m_i/(2*pi*m_e)) */
    double m_i = PLASMA_M_AR;
    double V_f = -0.5 * T_e_eV * log(m_i / (2.0 * M_PI * PLASMA_M_E));
    ps->V_float = V_f;
    ps->V_plasma = V_f + 0.5 * T_e_eV * log(m_i / (2.0 * M_PI * PLASMA_M_E));

    /* Collision frequencies:
     * Electron-neutral: nu_en = n_n * sigma_en * v_thermal_e */
    double v_th_e = sqrt(8.0 * PLASMA_E_CHARGE * fabs(T_e_eV)
                         / (M_PI * PLASMA_M_E));
    double sigma_en = 1e-19; /* typical elastic cross section [m^2] for Ar */
    ps->nu_en = ps->n_n * sigma_en * v_th_e;

    /* Electron-ion (Coulomb) collision frequency:
     * nu_ei = 2.9e-12 * n_e * Lambda * T_e^{-3/2} [s^-1] (Braginskii)
     * with Coulomb logarithm Lambda ~ 10-15 for low-T plasmas */
    double Lambda = 10.0;
    double n_e_20 = n_e * 1e-20; /* n_e in units of 10^20 m^-3 */
    double T_eV = fabs(T_e_eV);
    ps->nu_ei = 2.9e-12 * n_e_20 * Lambda * pow(T_eV, -1.5);

    /* Ionization frequency (approximate from Arrhenius-like rate):
     * nu_iz = n_n * K_iz(T_e) where K_iz ~ 1e-14 * exp(-E_iz/T_e) */
    double E_iz = 15.76; /* Ar ionization energy [eV] */
    double K_iz;
    if (T_eV > 0.0) {
        K_iz = 2.3e-14 * pow(T_eV, 0.59) * exp(-E_iz / T_eV);
    } else {
        K_iz = 0.0;
    }
    ps->nu_iz = ps->n_n * K_iz;
}

/*
 * Initialize a PlasmaSpecies structure.
 *
 * Each species in the plasma chemistry model requires its
 * fundamental properties: mass, charge, and energy thresholds.
 */
void plasma_species_init(PlasmaSpecies *sp, const char *name,
                         double mass_amu, int charge, double ioniz_eV,
                         double excit_eV, double polariz_m3,
                         PlasmaSpeciesType type)
{
    if (!sp) return;

    strncpy(sp->name, name, 31);
    sp->name[31] = '\0';
    sp->mass = mass_amu * 1.66053906660e-27; /* amu -> kg */
    sp->charge_number = charge;
    sp->ionization_energy = ioniz_eV;
    sp->excitation_energy = excit_eV;
    sp->polarizability = polariz_m3;
    sp->type = type;
}

/*
 * Initialize a PlasmaReaction structure.
 *
 * Sets up a single elementary reaction with Arrhenius rate constant:
 *   k(T_e) = A * T_e^n * exp(-E_threshold / T_e)
 */
void plasma_reaction_init(PlasmaReaction *rxn, const char *equation,
                          ReactionType type, double threshold_eV,
                          double A_factor, double n_exp,
                          double sigma_max, double E_at_max,
                          int *reactants, int n_r, int *products, int n_p)
{
    if (!rxn) return;

    strncpy(rxn->equation, equation, 127);
    rxn->equation[127] = '\0';
    rxn->type = type;
    rxn->threshold = threshold_eV;
    rxn->rate_constant_A = A_factor;
    rxn->rate_constant_n = n_exp;
    rxn->cross_section_max = sigma_max;
    rxn->energy_at_max = E_at_max;
    rxn->n_reactants = (n_r <= 4) ? n_r : 4;
    rxn->n_products = (n_p <= 4) ? n_p : 4;

    for (int i = 0; i < rxn->n_reactants && i < 4; i++)
        rxn->reactant_ids[i] = reactants[i];
    for (int i = 0; i < rxn->n_products && i < 4; i++)
        rxn->product_ids[i] = products[i];
}

/* ============================================================
 * L2: Plasma Quasineutrality and Screening
 * ============================================================ */

/*
 * Check if a plasma satisfies the quasineutrality condition.
 *
 * A plasma must satisfy three conditions simultaneously:
 * 1. Debye screening: lam_D << L (system size)
 * 2. Collective behavior: N_D = n * lam_D^3 >> 1 (many particles in Debye sphere)
 * 3. Plasma frequency >> collision frequency: omega_pe > nu_en
 *
 * Returns 1 if all conditions are met (true plasma), 0 otherwise.
 */
int plasma_quasineutrality_check(double n_e, double T_e_eV,
                                 double system_size, double nu_en)
{
    double lam_D = debye_length_electron(n_e, T_e_eV);
    double omega_pe = plasma_freq_electron(n_e);

    /* Condition 1: Debye length much less than system size */
    if (lam_D >= system_size * 0.1) return 0;

    /* Condition 2: Debye number N_D >> 1 (at least 10) */
    double N_D = n_e * lam_D * lam_D * lam_D;
    if (N_D < 10.0) return 0;

    /* Condition 3: Collective oscillations dominate collisions */
    if (omega_pe < nu_en) return 0;

    return 1;
}

/*
 * Compute the plasma parameter Lambda (Coulomb logarithm).
 *
 * Lambda = ln(12 * pi * n_e * lam_D^3) = ln(4 * pi * N_D)
 *
 * In industrial plasmas: Lambda ~ 10-15
 * In fusion plasmas: Lambda ~ 15-20
 * In space plasmas: Lambda ~ 20-30
 */
double plasma_coulomb_logarithm(double n_e, double T_e_eV)
{
    double lam_D = debye_length_electron(n_e, T_e_eV);
    if (lam_D <= 0.0) return 10.0;

    double N_D = n_e * lam_D * lam_D * lam_D;
    if (N_D <= 0.0) return 10.0;

    /* Lambda = ln(12 * pi * N_D) in the classical definition,
     * or more commonly Lambda = ln(4*pi*N_D) for Braginskii. */
    double Lambda = log(12.0 * M_PI * N_D);
    return (Lambda > 2.0) ? Lambda : 2.0;
}

/* ============================================================
 * L3: Magnetic Field Effects - Plasma Beta and Magnetization
 * ============================================================ */

/*
 * Plasma beta: ratio of thermal pressure to magnetic pressure.
 *
 * beta = n * k_B * T / (B^2 / (2*mu_0))
 *
 * beta << 1: magnetically confined (fusion, ICP at low power)
 * beta ~ 1: transitional
 * beta >> 1: flow-dominated (thermal plasma spray, arcs)
 *
 * For industrial ICP: B ~ 0.001 T, n ~ 1e17, T_e ~ 3 eV
 *   -> beta ~ 10^-3 (weakly magnetized)
 */
double plasma_beta(double n_e, double T_e_eV, double T_i_K, double B_field)
{
    if (B_field <= 0.0) return INFINITY;

    double p_e = n_e * PLASMA_E_CHARGE * T_e_eV; /* p_e = n_e * kT_e */
    double p_i = n_e * PLASMA_K_B * T_i_K;       /* p_i = n_i * kT_i */
    double p_mag = B_field * B_field / (2.0 * 1.25663706212e-6); /* B^2/(2*mu0) */

    return (p_e + p_i) / p_mag;
}

/*
 * Electron cyclotron frequency: omega_ce = e * B / m_e
 */
double electron_cyclotron_freq(double B_field)
{
    if (B_field <= 0.0) return 0.0;
    return PLASMA_E_CHARGE * B_field / PLASMA_M_E;
}

/*
 * Ion cyclotron frequency: omega_ci = e * B / m_i
 */
double ion_cyclotron_freq(double B_field, double m_i)
{
    if (B_field <= 0.0 || m_i <= 0.0) return 0.0;
    return PLASMA_E_CHARGE * B_field / m_i;
}

/*
 * Electron Larmor radius: r_L = m_e * v_perp / (e * B)
 * with v_perp from thermal speed.
 */
double electron_larmor_radius(double T_e_eV, double B_field)
{
    if (B_field <= 0.0) return INFINITY;

    double v_th = sqrt(2.0 * PLASMA_E_CHARGE * T_e_eV / PLASMA_M_E);
    return PLASMA_M_E * v_th / (PLASMA_E_CHARGE * B_field);
}

/*
 * Hall parameter: h = omega_ce * tau = omega_ce / nu_en
 *
 * h >> 1: magnetized electrons (E x B drift dominates)
 * h << 1: unmagnetized (collisions dominate)
 *
 * Critical for magnetron sputtering design.
 */
double hall_parameter(double B_field, double nu_en)
{
    if (nu_en <= 0.0 || B_field <= 0.0) return 0.0;
    double omega_ce = electron_cyclotron_freq(B_field);
    return omega_ce / nu_en;
}

/* ============================================================
 * L4: Energy and Power Flow
 * ============================================================ */

/*
 * Power balance in a steady-state discharge.
 *
 * At steady state, the RF/microwave power absorbed by the
 * plasma must balance all loss channels:
 *
 *   P_abs = P_ions + P_electrons + P_radiation + P_chemistry
 *
 * where:
 *   P_ions = Gamma_i * A_wall * (E_i + V_sheath)
 *   P_electrons = Gamma_e * A_wall * (2*T_e)
 *   P_radiation = n_e * n_n * K_exc * E_photon (volume)
 *   P_chemistry = dissociation/ionization energy sinks
 *
 * For low-T industrial plasmas, ion and electron wall losses
 * typically dominate (~70-90% of total power).
 */
double plasma_power_loss_to_walls(double n_e, double T_e_eV, double m_i,
                                   double wall_area, double V_sheath)
{
    if (n_e <= 0.0 || T_e_eV <= 0.0 || wall_area <= 0.0) return 0.0;

    double cs = sqrt(PLASMA_E_CHARGE * T_e_eV / m_i);
    double Gamma_i = n_e * cs; /* ion flux density */

    /* Ion energy loss: kinetic + recombination at wall */
    double E_i_loss = fabs(V_sheath) + 5.0 * T_e_eV; /* ion energy + ionization cost */

    /* Electron energy loss: kinetic energy carried to wall */
    double E_e_loss = 2.0 * T_e_eV;

    double P_loss = Gamma_i * wall_area * (E_i_loss + E_e_loss)
                   * PLASMA_E_CHARGE; /* convert eV to J */
    return P_loss;
}

/*
 * Estimate the ionization degree: alpha = n_i / (n_i + n_n)
 *
 * In low-temperature processing plasmas: alpha ~ 10^-6 to 10^-3
 * (weakly ionized). This is a key difference from fusion plasmas
 * where alpha ~ 1 (fully ionized).
 */
double ionization_degree(double n_e, double n_n)
{
    if (n_n <= 0.0) return 1.0;
    return n_e / (n_e + n_n);
}

/*
 * Collision mean free path for electrons at a given energy.
 *
 * lambda_mfp = 1 / (n_n * sigma_total)
 *
 * For 3 eV electrons in Ar at 1 Pa: lambda_mfp ~ 3 cm.
 * For sheaths (~1 mm) and features (~100 nm), this determines
 * whether transport is ballistic or diffusive.
 */
double electron_mean_free_path(double n_gas, double energy_eV,
                               CrossSectionModel *sigma_m)
{
    if (n_gas <= 0.0 || !sigma_m) return INFINITY;

    double sigma;
    if (energy_eV < sigma_m->eps0 * 0.001) {
        sigma = sigma_m->sigma0 * 1e-6;
    } else {
        double arg = energy_eV / sigma_m->eps0;
        sigma = sigma_m->sigma0 * pow(arg, sigma_m->a)
               * exp(-sigma_m->a * arg)
               + sigma_m->sigma1 * (1.0 - exp(-energy_eV / sigma_m->eps1));
    }

    if (sigma <= 0.0) return INFINITY;
    return 1.0 / (n_gas * sigma);
}

/* ============================================================
 * L5: Discharge Scaling Laws
 * ============================================================ */

/*
 * Similarity (scaling) laws for gas discharges.
 *
 * Discharges with the same p*d product (Paschen similarity)
 * and the same E/N (reduced electric field) will have similar
 * electron kinetics, even if dimensions differ by orders of
 * magnitude. This enables scaling from lab to industrial reactors.
 *
 * Scaling rules:
 * - Linear dimension L: all lengths scale with 1/p if pd conserved
 * - Voltage: V conserved if pd conserved (same E/p)
 * - Current density: j ~ p^2 (j/p^2 conserved)
 * - Power density: P/V ~ p (P/(V*p) conserved)
 *
 * Reference: von Engel, "Ionized Gases" (1965), Ch.4
 *            Francis, "Gas Discharge Physics" (1956)
 */
void discharge_similarity_scale(double pd_original, double p_original,
                                 double p_new, double *L_scale,
                                 double *V_scale, double *j_scale)
{
    if (p_original <= 0.0 || p_new <= 0.0) {
        *L_scale = 1.0; *V_scale = 1.0; *j_scale = 1.0;
        return;
    }

    /* Length scales inversely with pressure (to keep p*d const) */
    *L_scale = p_original / p_new;

    /* Voltage conserved (same E/p -> same V for same pd) */
    *V_scale = 1.0;

    /* Current density scales with p^2 */
    double p_ratio = p_new / p_original;
    *j_scale = p_ratio * p_ratio;

    (void)pd_original;
}

/*
 * Estimate RF plasma density from engineering parameters.
 *
 * Simplified global model for CCP/ICP density estimation.
 * Input: RF power, pressure, chamber dimensions, driving frequency.
 * Output: estimated n_e, T_e.
 *
 * Uses the global (volume-averaged) model scaling:
 *   n_e ~ P_rf / (e * u_B * A_eff * E_T)
 *   T_e determined by particle balance.
 */
void estimate_plasma_density_rf(double P_rf, double pressure_pa,
                                 double radius, double height,
                                 double freq_hz, double *n_e_out,
                                 double *T_e_out)
{
    /* Chamber volume and effective loss area */
    double V = M_PI * radius * radius * height;
    double A_eff = 2.0 * M_PI * radius * radius + 2.0 * M_PI * radius * height;

    /* Neutral density */
    double n_gas = pressure_pa / (PLASMA_K_B * 300.0);

    /* Particle balance determines T_e (for Argon):
     * K_iz(T_e) / u_B(T_e) = 1 / (n_gas * d_eff) */
    double d_eff = V / A_eff;
    double m_i = PLASMA_M_AR;
    double E_iz = 15.76;

    /* Iterative T_e solver */
    double T_e = 3.0;
    for (int iter = 0; iter < 10; iter++) {
        double u_B = sqrt(PLASMA_E_CHARGE * T_e / m_i);
        double K_iz = 2.3e-14 * pow(T_e, 0.59) * exp(-E_iz / T_e);
        double tau = d_eff / u_B;
        double balance = K_iz * n_gas * tau;

        if (balance > 0.0) {
            double T_e_new = E_iz / log(K_iz * n_gas * tau * 2.4);
            if (T_e_new > 0.1 && T_e_new < 20.0) {
                T_e = T_e_new;
            }
        }
        if (fabs(balance - 1.0) < 0.01) break;
    }

    /* Power balance determines n_e */
    double u_B = sqrt(PLASMA_E_CHARGE * T_e / m_i);
    double E_T = E_iz + 3.0 * T_e + 2.0 * T_e + 100.0; /* ~100V sheath */
    double n_e = P_rf / (PLASMA_E_CHARGE * u_B * A_eff * E_T);

    *n_e_out = (n_e > 0.0) ? n_e : 1e15;
    *T_e_out = (T_e > 0.0) ? T_e : 2.0;

    (void)freq_hz;
}

/* ============================================================
 * L7 Applications: Industrial Plasma Parameter Database
 * ============================================================ */

/*
 * Get typical plasma parameters for common industrial processes.
 *
 * These are representative values from literature and industry
 * (Lam Research, Applied Materials, Tokyo Electron tools).
 */
void industrial_plasma_defaults(const char *process_type,
                                 PlasmaState *ps,
                                 DischargeParams *dp,
                                 DischargeGeometry *geom)
{
    if (!ps || !dp || !geom) return;

    if (strstr(process_type, "RIE") || strstr(process_type, "etch")) {
        /* Reactive Ion Etching (CCP, 13.56 MHz, 1-100 mTorr) */
        plasma_state_init(ps, 1e16, 3.0, 300.0, 300.0, 5.0, 0.0);
        dp->power = 500.0;
        dp->frequency = 13.56e6;
        dp->dc_bias = -200.0;
        dp->voltage_amplitude = 500.0;
        dp->match_efficiency = 0.9;
        dp->is_pulsed = 0;
        dp->pulse_duty = 1.0;
        dp->pulse_freq = 0.0;
        geom->electrode_area = 0.05;
        geom->ground_area = 0.2;
        geom->gap_length = 0.03;
        geom->chamber_radius = 0.2;
        geom->chamber_height = 0.1;
        geom->wafer_diameter = 0.3;
        geom->pump_speed = 1.0;
        geom->gas_flow_rate = 100.0;

    } else if (strstr(process_type, "PECVD") || strstr(process_type, "deposition")) {
        /* PECVD (CCP, 13.56 MHz, 100-1000 mTorr) */
        plasma_state_init(ps, 5e15, 2.5, 300.0, 500.0, 100.0, 0.0);
        dp->power = 300.0;
        dp->frequency = 13.56e6;
        dp->dc_bias = -50.0;
        dp->voltage_amplitude = 200.0;
        dp->match_efficiency = 0.85;
        dp->is_pulsed = 0;
        dp->pulse_duty = 1.0;
        dp->pulse_freq = 0.0;
        geom->electrode_area = 0.07;
        geom->ground_area = 0.25;
        geom->gap_length = 0.025;
        geom->chamber_radius = 0.25;
        geom->chamber_height = 0.1;
        geom->wafer_diameter = 0.3;
        geom->pump_speed = 0.5;
        geom->gas_flow_rate = 500.0;

    } else if (strstr(process_type, "ICP") || strstr(process_type, "inductive")) {
        /* ICP (13.56 MHz, 1-50 mTorr, high density) */
        plasma_state_init(ps, 5e17, 3.5, 500.0, 500.0, 1.0, 0.001);
        dp->power = 1000.0;
        dp->frequency = 13.56e6;
        dp->dc_bias = -20.0;
        dp->voltage_amplitude = 100.0;
        dp->match_efficiency = 0.8;
        dp->is_pulsed = 0;
        dp->pulse_duty = 1.0;
        dp->pulse_freq = 0.0;
        geom->electrode_area = 0.03;
        geom->ground_area = 0.3;
        geom->gap_length = 0.1;
        geom->chamber_radius = 0.15;
        geom->chamber_height = 0.2;
        geom->wafer_diameter = 0.3;
        geom->pump_speed = 2.0;
        geom->gas_flow_rate = 50.0;

    } else if (strstr(process_type, "sputter") || strstr(process_type, "PVD")) {
        /* Magnetron sputtering (DC, 1-10 mTorr) */
        plasma_state_init(ps, 1e17, 3.0, 400.0, 400.0, 0.5, 0.05);
        dp->power = 2000.0;
        dp->frequency = 0.0; /* DC */
        dp->dc_bias = -400.0;
        dp->voltage_amplitude = 400.0;
        dp->match_efficiency = 0.9;
        dp->is_pulsed = 1;
        dp->pulse_duty = 0.8;
        dp->pulse_freq = 100e3;
        geom->electrode_area = 0.02;
        geom->ground_area = 0.5;
        geom->gap_length = 0.1;
        geom->chamber_radius = 0.3;
        geom->chamber_height = 0.3;
        geom->wafer_diameter = 0.3;
        geom->pump_speed = 1.0;
        geom->gas_flow_rate = 30.0;

    } else {
        /* Generic low-pressure plasma defaults */
        plasma_state_init(ps, 1e16, 3.0, 300.0, 300.0, 10.0, 0.0);
        dp->power = 300.0;
        dp->frequency = 13.56e6;
        dp->dc_bias = -100.0;
        dp->voltage_amplitude = 300.0;
        dp->match_efficiency = 0.85;
        dp->is_pulsed = 0;
        dp->pulse_duty = 1.0;
        dp->pulse_freq = 0.0;
        geom->electrode_area = 0.05;
        geom->ground_area = 0.2;
        geom->gap_length = 0.04;
        geom->chamber_radius = 0.2;
        geom->chamber_height = 0.15;
        geom->wafer_diameter = 0.3;
        geom->pump_speed = 1.0;
        geom->gas_flow_rate = 100.0;
    }
}
