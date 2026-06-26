/**
 * fusion_equilibrium.c — Equilibrium, Operating Points, and Power Balance
 *
 * Every function implements an independent physics knowledge point.
 * Refs: ITER Technical Basis (2002), Wesson (2011),
 *       Sorbom et al., FED 100, 378 (2015) — SPARC,
 *       JET Team, Nucl. Fusion 39, 1619 (1999),
 *       Maisonnier et al., FED 75-79, 1173 (2005) — DEMO
 *
 * L6: ITER, DEMO, SPARC, JET operating points
 * L7: Power balance, tritium breeding, neutron wall loading
 */

#include "include/fusion_plasma.h"
#include <math.h>
#include <assert.h>
#include <stdio.h>

/* ================================================================
 * L6: ITER Operating Point
 * ================================================================ */

/**
 * iter_operating_point — ITER baseline operating parameters
 *
 * ITER (International Thermonuclear Experimental Reactor)
 * Key parameters for 500 MW fusion power, Q=10 scenario:
 *   R = 6.2 m, a = 2.0 m, kappa = 1.7, delta = 0.33
 *   B = 5.3 T, Ip = 15 MA
 *   ne = 1.0e20 m^-3, Te = 10 keV, Ti = 10 keV
 *   P_fusion = 500 MW, Q = 10, tau_E ~ 3.7 s
 *
 * Source: ITER Technical Basis, IAEA (2002)
 */
void iter_operating_point(PlasmaParameters *p) {
    if (!p) return;
    p->R = 6.2;
    p->a = 2.0;
    p->elongation = 1.7;
    p->triangularity = 0.33;
    p->B = 5.3;
    p->Ip = 15.0e6;
    p->ne = 1.0e20;
    p->ni = 0.85e20;
    p->Te = 10.0e3;
    p->Ti = 10.0e3;
    p->Zeff = 1.6;
    p->V_p = plasma_volume(p->R, p->a, p->elongation);
}

/**
 * iter_power_balance — ITER power balance calculation
 *
 * Computes complete power balance for ITER parameters.
 *   P_fusion = 500 MW (design)
 *   P_alpha = 100 MW (20% of fusion power)
 *   P_aux = 50 MW (injected, Q=10 baseline)
 *   P_radiation = P_brem + P_cyc ~ 50 MW (core)
 *   P_transport = balance
 *
 * @param p   plasma parameters (ITER baseline)
 * @param eb  output: energy balance filled
 */
void iter_power_balance(const PlasmaParameters *p, EnergyBalance *eb) {
    if (!p || !eb) return;

    double ni = p->ni > 0 ? p->ni : p->ne;
    double nD = 0.5 * ni;
    double nT = 0.5 * ni;

    /* Compute D-T reactivity at Ti=10 keV */
    double Ti_keV = p->Ti / 1e3;
    double sigma_v = bosch_hale_sigma_v_dt(Ti_keV);

    /* Fusion power: P_fus = nD * nT * <sigma v> * E_fus */
    eb->P_fusion = nD * nT * sigma_v * E_FUSION_DT_J * p->V_p;

    /* Alpha power: 20% of fusion energy */
    double E_alpha_J = E_ALPHA * E_CHARGE;
    eb->P_alpha = nD * nT * sigma_v * E_alpha_J * p->V_p;

    /* Bremsstrahlung radiation */
    double Te_eV = p->Te;
    double P_brem_core = bremsstrahlung_power_density(p->ne, p->Zeff, Te_eV);
    double P_cyc_core = cyclotron_power_density(p->ne, Te_eV, p->B);
    eb->P_radiation = (P_brem_core + P_cyc_core) * p->V_p;

    /* Ohmic heating */
    double lnLambda = coulomb_logarithm(Te_eV, p->ne);
    double eta = spitzer_resistivity(p->Zeff, Te_eV, lnLambda);
    eb->P_ohmic = ohmic_heating_power(eta, p->Ip, p->a, p->elongation) * p->V_p;

    /* Auxiliary heating set to achieve Q=10 */
    eb->P_aux = eb->P_fusion / 10.0;

    /* Transport losses = everything else */
    eb->P_transport = eb->P_alpha + eb->P_ohmic + eb->P_aux - eb->P_radiation;

    /* Q */
    eb->Q = fusion_gain_Q(eb->P_fusion, eb->P_aux);

    /* Ignition check */
    double W = plasma_stored_energy(p);
    double tauE = energy_confinement_time(W, eb->P_radiation + eb->P_transport);
    double n_20 = p->ne / 1e20;
    double T_keV = p->Te / 1e3;
    eb->ignited = ignition_condition(n_20, T_keV, tauE);
}

/* ================================================================
 * L6: DEMO Operating Point
 * ================================================================ */

/**
 * demo_operating_point — DEMO (Demonstration Power Plant) parameters
 *
 * EU DEMO conceptual design:
 *   R = 9.0 m, a = 2.9 m, kappa = 1.65, delta = 0.33
 *   B = 5.7 T, Ip = 19.6 MA
 *   ne = 1.2e20 m^-3
 *   P_fusion ~ 2 GW, P_elec ~ 500 MW
 *   Q = 25-50, pulse length: 2 hours (steady-state goal)
 *
 * Source: Maisonnier et al., FED 75-79, 1173 (2005)
 */
void demo_operating_point(PlasmaParameters *p) {
    if (!p) return;
    p->R = 9.0;
    p->a = 2.9;
    p->elongation = 1.65;
    p->triangularity = 0.33;
    p->B = 5.7;
    p->Ip = 19.6e6;
    p->ne = 1.2e20;
    p->ni = 1.0e20;
    p->Te = 15.0e3;
    p->Ti = 15.0e3;
    p->Zeff = 1.8;
    p->V_p = plasma_volume(p->R, p->a, p->elongation);
}

/* ================================================================
 * L6: SPARC Operating Point
 * ================================================================ */

/**
 * sparc_operating_point — SPARC (CFS/MIT) compact high-field tokamak
 *
 * SPARC parameters:
 *   R = 1.85 m, a = 0.57 m, kappa = 1.7, delta = 0.3
 *   B = 12.2 T (HTS magnets), Ip = 8.7 MA
 *   ne = 4.0e20 m^-3, Te = 7 keV, Ti = 7 keV
 *   P_fusion ~ 50-140 MW, Q >= 2 (goal: Q~11)
 *
 * Key innovation: REBCO high-temperature superconductor
 * magnets enabling B > 10 T at R < 2 m.
 *
 * Source: Sorbom et al., FED 100, 378 (2015)
 *         Creely et al., JPP 86, 865860502 (2020)
 */
void sparc_operating_point(PlasmaParameters *p) {
    if (!p) return;
    p->R = 1.85;
    p->a = 0.57;
    p->elongation = 1.7;
    p->triangularity = 0.3;
    p->B = 12.2;
    p->Ip = 8.7e6;
    p->ne = 4.0e20;
    p->ni = 3.5e20;
    p->Te = 7.0e3;
    p->Ti = 7.0e3;
    p->Zeff = 1.5;
    p->V_p = plasma_volume(p->R, p->a, p->elongation);
}

/* ================================================================
 * L6: JET D-T Record
 * ================================================================ */

/**
 * jet_dt_record — JET D-T record discharge parameters
 *
 * JET (Joint European Torus) holds the world record for
 * magnetic fusion energy production:
 *   P_fusion = 16.1 MW (1997, D-T campaign)
 *   Q = 0.67 (near breakeven)
 *   W = 16 MJ stored energy
 *   triple product ~ 1.5e21 m^-3 keV s
 *
 * JET 2023 DTE3 campaign: sustained 5-10 s of D-T fusion,
 * demonstrating ITER-relevant scenarios.
 *
 * Source: JET Team, Nucl. Fusion 39, 1619 (1999)
 */
void jet_dt_record(PlasmaParameters *p) {
    if (!p) return;
    p->R = 2.96;
    p->a = 0.90;
    p->elongation = 1.6;
    p->triangularity = 0.24;
    p->B = 3.8;
    p->Ip = 4.0e6;
    p->ne = 4.0e19;
    p->ni = 3.5e19;
    p->Te = 10.0e3;
    p->Ti = 10.0e3;
    p->Zeff = 2.5;
    p->V_p = plasma_volume(p->R, p->a, p->elongation);
}

/* ================================================================
 * L7: Fusion Reactor Scaling — 4 functions
 * ================================================================ */

/**
 * fusion_power_scaling — Fusion power scaling with size [W]
 *
 * Definition: P_fus ~ n^2 * <sigma v> * V_p
 *          ~ (beta_N * B^2)^2 * V_p (at fixed T, shape)
 * For fixed beta_N, B, and T: P_fus ~ V_p ~ R * a^2
 * For fixed aspect ratio: P_fus ~ R^3
 *
 * @param p  plasma parameters
 * @return   estimated fusion power [W]
 */
double fusion_power_scaling(const PlasmaParameters *p) {
    if (!p) return 0.0;
    double Ti_keV = p->Ti / 1e3;
    double sigma_v = bosch_hale_sigma_v_dt(Ti_keV);
    double nD = 0.5 * p->ni;
    double nT = 0.5 * p->ni;
    return nD * nT * sigma_v * E_FUSION_DT_J * p->V_p;
}

/**
 * fusion_relevant_temperature — Optimal fusion temperature
 *
 * Definition: T_opt maximizes <sigma v>/T^2 (reactor efficiency)
 *          D-T: T_opt ~ 14 keV for <sigma v>/T^2 maximum
 *          But considering realistic profiles: T(0) ~ 15-25 keV
 * Physics: trade-off between cross-section (increases with T)
 *          and pressure (T appears in denominator of Lawson).
 */
double optimal_fusion_temperature_dt(void) {
    return 14.0;  /* keV — maximizes <sigma v>/T^2 */
}

/**
 * fusion_triple_product_target — Required triple product for Q=target
 *
 * Definition: n*T*tau_E >= (12/E_alpha) * (T/<sigma v>) * (1/Q) * const
 * Simplified: triple_prod_req = 3.0 * (1 + 5/Q)
 * Physics: scaling of required triple product with desired Q.
 * For Q=10: triple_prod_req ~ 4.5 (in 10^20 m^-3 keV s units)
 * For Q=inf (ignition): triple_prod_req = 3.0
 */
double required_triple_product_for_Q(double Q_target) {
    if (Q_target <= 0.0) return 1e10;
    double eta_abs = 0.2;  /* alpha absorption efficiency */
    return 3.0 * (1.0 + 1.0 / (Q_target * eta_abs));
}

/**
 * capital_cost_scaling — Fusion capital cost scaling [relative units]
 *
 * Definition: Cost ~ V_p^0.5 * B^1.5 (simplified)
 * Physics: larger volume and higher field increase magnet,
 *          vacuum vessel, and building costs.
 * Smaller reactors (SPARC, STEP) aim to reduce cost through
 * high-field HTS magnets despite higher B cost factor.
 *
 * @param p  plasma parameters
 * @return   relative capital cost
 */
double capital_cost_scaling(const PlasmaParameters *p) {
    if (!p) return 0.0;
    double volume_factor = pow(p->V_p, 0.5);
    double field_factor = pow(p->B, 1.5);
    return volume_factor * field_factor;
}

/* ================================================================
 * L7: Blanket and Tritium Breeding — 3 functions
 * ================================================================ */

/**
 * lithium_blanket_enrichment — Required Li-6 enrichment for TBR
 *
 * Definition: TBR ~ Li-6 enrichment * breeding capability
 * Natural Li: 7.5% Li-6
 * Enriched Li: 30-90% Li-6 (needed for TBR > 1.0)
 *
 * @param TBR_target       desired TBR (>1.0)
 * @param natural_capability  TBR with natural lithium
 * @return                 required Li-6 enrichment fraction
 */
double lithium_enrichment_for_tbr(double TBR_target, double natural_capability) {
    if (natural_capability <= 0.0) return 1.0;
    double enrichment = TBR_target / natural_capability;
    if (enrichment > 0.9) enrichment = 0.9;
    if (enrichment < 0.075) enrichment = 0.075;
    return enrichment;
}

/**
 * tritium_inventory_required — Required tritium inventory [kg]
 *
 * Definition: I_T = f_T * (P_fus * t_burn) / (E_fus * eta_burn)
 * Physics: tritium needed for a fusion burn campaign.
 *          Crucial for start-up inventory planning.
 * For ITER: ~0.5-1 kg T (total campaign)
 * For DEMO: ~5-10 kg T (annual consumption)
 *
 * @param P_fusion     fusion power [W]
 * @param t_burn       burn time [s]
 * @param f_burn       burn fraction
 * @return             tritium inventory [kg]
 */
double tritium_inventory_required(double P_fusion, double t_burn, double f_burn) {
    double E_DT_J = E_FUSION_DT * E_CHARGE;
    double reactions_needed = P_fusion * t_burn / E_DT_J;
    double T_kg = reactions_needed * M_T;
    if (f_burn > 0.0) T_kg /= f_burn;
    return T_kg;
}

/**
 * doubling_time_tritium — Tritium doubling time [years]
 *
 * Definition: t_2 = ln(2) * I_T / (S_production - S_consumption)
 * Physics: time to double the tritium inventory for
 *          starting up new reactors. t_2 < 5 years needed
 *          for reasonable fusion deployment rate.
 *
 * @param inventory   current T inventory [kg]
 * @param production  annual T production [kg/yr]
 * @param consumption annual T consumption [kg/yr]
 * @return            doubling time [years]
 */
double tritium_doubling_time(double inventory, double production, double consumption) {
    double net = production - consumption;
    if (net <= 0.0) return 1e10;  /* never doubles */
    return 0.693 * inventory / net;
}