/**
 * fusion_reactor.c — Fusion Reactor Engineering and Power Plant Design
 *
 * Every function implements an independent physics knowledge point.
 * Refs: Stacey "Fusion: An Introduction to the Physics and
 *       Technology of Magnetic Confinement Fusion" 2nd ed (2010),
 *       ITER Technical Basis (2002),
 *       Maisonnier et al., FED (2005) — DEMO,
 *       Sorbom et al., FED (2015) — SPARC/ARC
 *
 * L7: Power plant thermal efficiency, availability,
 *     levelized cost of electricity, tokamak reactor
 *     parameter optimization
 * L8: Liquid metal breeder blankets, advanced divertor,
 *     structural material limits (dpa, helium embrittlement)
 * L9: ARC-class compact reactors, stellarator reactor concepts,
 *     inertial fusion energy (IFE) reactor concepts
 */

#include "include/fusion_plasma.h"
#include "include/fusion_confinement.h"
#include "include/fusion_mhd.h"
#include <math.h>
#include <assert.h>

/* ================================================================
 * L7: Fusion Power Plant Engineering — 5 functions
 * ================================================================ */

/**
 * thermal_to_electric_efficiency — Thermal-electric conversion [fraction]
 *
 * Definition: eta_th = P_elec / P_thermal
 * Physics: Carnot efficiency limit modified by engineering factors.
 *          eta_th ~ 0.3-0.4 for steam Rankine cycle
 *          eta_th ~ 0.4-0.5 for supercritical CO2 Brayton cycle
 *          eta_th ~ 0.5-0.6 for advanced combined cycles
 *
 * @param T_hot     coolant hot leg temperature [K]
 * @param T_cold    coolant cold leg temperature [K]
 * @param eta_eng   engineering efficiency factor (0.5-0.7)
 * @return          thermal efficiency
 */
double thermal_to_electric_efficiency(double T_hot, double T_cold, double eta_eng) {
    if (T_hot <= T_cold) return 0.0;
    double eta_carnot = 1.0 - T_cold / T_hot;
    return eta_eng * eta_carnot;
}

/**
 * fusion_power_plant_net_power — Net electric power output [W]
 *
 * Definition: P_net = P_gross - P_recirc
 * Physics: fusion plant must subtract recirculating power
 *          for heating, magnets, cryogenics, pumps, etc.
 *          P_recirc/P_gross ~ 10-20% for tokamak.
 *
 * @param P_fusion       fusion power [W]
 * @param P_heating      auxiliary heating power [W]
 * @param eta_th         thermal conversion efficiency
 * @param P_housekeeping other plant power consumption [W]
 * @return               net electric power [W]
 */
double net_electric_power(double P_fusion, double P_heating, double eta_th,
                           double P_housekeeping) {
    double P_thermal = P_fusion + P_heating;
    double P_gross = eta_th * P_thermal;
    return P_gross - P_heating - P_housekeeping;
}

/**
 * plant_availability — Fusion plant availability
 *
 * Definition: A = t_operating / (t_operating + t_maintenance)
 * Physics: fusion plants require regular maintenance for
 *          in-vessel components (divertor, first wall, blankets).
 *          Availability crucial for economic viability.
 * ITER: A ~ 0.05-0.25 (experimental)
 * DEMO target: A ~ 0.3-0.5
 * Commercial target: A > 0.7
 *
 * @param t_operating     operating time per year [days]
 * @param t_maintenance   maintenance time per year [days]
 * @return                availability fraction
 */
double plant_availability(double t_operating, double t_maintenance) {
    double total = t_operating + t_maintenance;
    if (total <= 0.0) return 0.0;
    return t_operating / total;
}

/**
 * first_wall_neutron_damage — First wall neutron damage [dpa/yr]
 *
 * Definition: dpa_rate = NRT_factor * NWL * t_op * sigma_dpa / E_threshold
 * Physics: displacement per atom (dpa) measures radiation damage.
 *          EUROFER steel limit: ~50-70 dpa (end of life)
 *          Advanced ODS steels: ~150-200 dpa
 * For ITER: ~3 dpa over lifetime (low)
 * For DEMO: ~50-100 dpa over lifetime (significant)
 *
 * @param NWL            neutron wall loading [MW/m^2]
 * @param operation_years operating time at full power [years]
 * @param material_factor material-specific dpa rate [(dpa/yr)/(MW/m^2)]
 * @return               accumulated dpa
 */
double first_wall_dpa(double NWL, double operation_years, double material_factor) {
    return NWL * operation_years * material_factor;
}

/**
 * helium_production_rate — Helium production in structural materials [appm/yr]
 *
 * Definition: He_rate ~ NWL * cross_section_factor * density_factor
 * Physics: (n,alpha) reactions produce helium in structural
 *          materials, leading to embrittlement at grain boundaries.
 *          He limit for EUROFER: ~500-1000 appm
 *
 * @param NWL                neutron wall loading [MW/m^2]
 * @param He_production_coeff material coefficient [appm/(MW/m^2*yr)]
 * @param operation_years    years at full power
 * @return                   helium concentration [appm]
 */
double helium_production(double NWL, double He_production_coeff, double operation_years) {
    return NWL * He_production_coeff * operation_years;
}

/* ================================================================
 * L7: Levelized Cost of Electricity (LCOE) — 3 functions
 * ================================================================ */

/**
 * fusion_lcoe_simple — Simplified fusion LCOE [currency/MWh]
 *
 * Definition: LCOE = (C_capital * FCR + C_O&M + C_fuel) / E_annual
 *          FCR = fixed charge rate (amortization)
 *          E_annual = P_net * 8760 * A (MWh/year)
 *
 * @param capital_cost        overnight capital cost [currency]
 * @param fcr                 fixed charge rate [1/yr]
 * @param om_cost_annual      annual O&M cost [currency/yr]
 * @param fuel_cost_annual    annual fuel cost [currency/yr]
 * @param P_net_MW            net electric power [MW]
 * @param availability        plant availability [fraction]
 * @return                    LCOE [currency/MWh]
 */
double fusion_lcoe(double capital_cost, double fcr, double om_cost_annual,
                    double fuel_cost_annual, double P_net_MW, double availability) {
    double E_annual = P_net_MW * 8760.0 * availability;
    if (E_annual <= 0.0) return 1e10;
    double annual_cost = capital_cost * fcr + om_cost_annual + fuel_cost_annual;
    return annual_cost / E_annual;
}

/**
 * fusion_plant_capital_cost — Fusion plant capital cost estimate [USD]
 *
 * Simplified scaling based on ITER and DEMO studies:
 *   Cost ~ C0 * (P_fus/500)^0.6 * (R/6.2)^1.2
 * Physics: cost increases with size and power, but with
 *          economy of scale (exponent < 1).
 * ITER construction cost: ~20-25 BUSD (2020 estimate)
 * ARC (CFS) estimated cost: ~3-5 BUSD (compact HTS)
 *
 * @param P_fusion_MW  fusion power [MW]
 * @param R            major radius [m]
 * @return             estimated capital cost [BUSD]
 */
double fusion_capital_cost_busd(double P_fusion_MW, double R) {
    return 20.0 * pow(P_fusion_MW / 500.0, 0.6) * pow(R / 6.2, 1.2);
}

/**
 * fusion_competitiveness_index — Fusion economic competitiveness
 *
 * Definition: FCI = (LCOE_fusion / LCOE_reference)
 *          FCI < 1: fusion competitive
 *          FCI > 1: fusion not competitive without carbon pricing
 *
 * This simplified metric estimates whether fusion can compete
 * with other low-carbon energy sources.
 */
double fusion_competitiveness_index(double LCOE_fusion, double LCOE_reference) {
    if (LCOE_reference <= 0.0) return 1e10;
    return LCOE_fusion / LCOE_reference;
}

/* ================================================================
 * L8: Advanced Blanket Concepts — 4 functions
 * ================================================================ */

/**
 * liquid_metal_breeding_potential — Breeding potential of liquid metals
 *
 * Pb-Li eutectic (17% Li-6): natural TBR ~ 0.8-1.0
 * Pure Li: TBR ~ 1.2-1.5
 * FLiBe (LiF-BeF2): TBR ~ 0.9-1.1 (with Be neutron multiplier)
 *
 * @param Li6_enrichment   Li-6 enrichment fraction (0.075 natural, up to 0.9)
 * @param Pb_fraction      Pb fraction in Pb-Li (0 for pure Li, ~0.83 for PbLi)
 * @return                 estimated TBR
 */
double liquid_metal_tbr_estimate(double Li6_enrichment, double Pb_fraction) {
    double Li_density_factor = 1.0 - Pb_fraction;
    double base_TBR = 1.2 * Li_density_factor;
    return base_TBR * (Li6_enrichment / 0.075);
}

/**
 * beryllium_multiplier_factor — Beryllium neutron multiplication
 *
 * Be(n,2n) reaction multiplies neutrons for tritium breeding.
 * M_Be = (total neutrons out) / (source neutrons in)
 * Pure Be: M ~ 1.5-2.0
 * Be in Be12Ti: M ~ 1.3-1.5
 *
 * @param Be_thickness_cm   Be layer thickness [cm]
 * @return                  neutron multiplication factor
 */
double beryllium_neutron_multiplier(double Be_thickness_cm) {
    double M_sat = 1.85;
    double lambda = 5.0;
    return 1.0 + (M_sat - 1.0) * (1.0 - exp(-Be_thickness_cm / lambda));
}

/**
 * blanket_energy_multiplication — Total blanket energy multiplication
 *
 * Definition: M_E = E_total / E_neutron_DT
 * Physics: exothermic reactions (Li-6 + n -> T + He-4 + 4.78 MeV)
 *          in the blanket add to the fusion neutron energy.
 * For D-T fusion: E_neutron = 14.07 MeV
 * With breeding blanket: M_E ~ 1.1-1.3 (10-30% energy gain)
 *
 * @param TBR        tritium breeding ratio
 * @param Be_factor  Be neutron multiplication factor
 * @return           total energy multiplication factor
 */
double blanket_energy_multiplication(double TBR, double Be_factor) {
    double E_n = 14.07e6 * E_CHARGE;
    double E_exo = 4.78e6 * E_CHARGE;
    double E_total = E_n * Be_factor + E_exo * TBR;
    return E_total / E_n;
}

/**
 * structural_material_temperature_limit — Max operating temperature [K]
 *
 * Reduced Activation Ferritic-Martensitic (RAFM, e.g., EUROFER97):
 *   T_max ~ 550 C (823 K) limited by creep strength
 * Oxide Dispersion Strengthened (ODS) steels:
 *   T_max ~ 650-700 C (923-973 K)
 * SiC/SiC composites:
 *   T_max ~ 1000-1200 C (1273-1473 K)
 * Tungsten (divertor armor):
 *   T_max ~ 1200 C (1473 K) in plasma-facing, recrystallization limit
 *
 * Returns temperature-dependent allowable stress ratio.
 */
double structural_thermal_efficiency_factor(double T_K) {
    double T_ref = 823.0;
    if (T_K <= T_ref) return 1.0;
    double degradation = (T_K - T_ref) / (1273.0 - T_ref);
    if (degradation > 1.0) degradation = 1.0;
    return 1.0 - 0.5 * degradation * degradation;
}

/* ================================================================
 * L9: Advanced Reactor Concepts — 5 functions
 * ================================================================ */

/**
 * arc_reactor_magnet_cost — ARC-class HTS magnet cost [USD]
 *
 * ARC (Affordable, Robust, Compact) reactor concept uses
 * REBCO high-temperature superconductor (HTS) at ~20 K.
 * HTS enables B > 9 T at R ~ 3.3 m.
 * This function estimates magnet cost based on stored energy.
 *
 * @param B_max       peak field on coil [T]
 * @param R_coil      coil major radius [m]
 * @return            estimated magnet cost [USD]
 */
double arc_magnet_cost(double B_max, double R_coil) {
    double stored_energy_MJ = 0.5 * (B_max * B_max / MU0) *
                                (2.0 * M_PI * M_PI * R_coil * R_coil * R_coil * 0.3);
    double cost_per_MJ = 2.0e6;
    return stored_energy_MJ * cost_per_MJ;
}

/**
 * stellarator_coil_complexity — Stellarator coil complexity metric
 *
 * Definition: complexity ~ (1/N_fp) * sum |kappa_i|
 *          N_fp = number of field periods
 *          kappa_i = coil curvature
 * Physics: stellarator coils are non-planar 3D shapes.
 *          W7-X has 50 non-planar coils in 5 field periods.
 *          Complexity drives manufacturing cost.
 *
 * @param N_fp            number of field periods
 * @param avg_curvature   average coil curvature [1/m]
 * @param n_coils         number of coils
 * @return                complexity metric [arbitrary]
 */
double stellarator_coil_complexity(int N_fp, double avg_curvature, int n_coils) {
    return (double)n_coils * avg_curvature / (double)N_fp;
}

/**
 * stellarator_alpha_particle_confinement — Alpha confinement fraction
 *
 * Definition: f_conf = fraction of alpha orbits confined
 * Physics: stellarators must be optimized for alpha particle
 *          confinement. Ripple-induced losses are the main concern.
 *          W7-X optimized for low ripple (delta_B/B ~ 0.1%).
 *
 * @param ripple_amplitude   magnetic field ripple delta_B/B
 * @param epsilon_alpha      alpha Larmor radius / minor radius
 * @return                   confined fraction [0-1]
 */
double stellarator_alpha_confinement(double ripple_amplitude, double epsilon_alpha) {
    double loss_threshold = 0.01;
    double confinement = 1.0 / (1.0 + (ripple_amplitude / loss_threshold) *
                                       (epsilon_alpha / 0.01));
    if (confinement > 1.0) confinement = 1.0;
    return confinement;
}

/**
 * inertial_fusion_gain_curve — ICF gain curve (simplified)
 *
 * Definition: G = P_fus / P_driver
 * Physics: inertial confinement fusion gain as function of
 *          fuel rho*R (areal density, g/cm^2).
 *          G ~ (rho*R)^2 for rho*R < 3 g/cm^2 (central ignition)
 *          G ~ exp(rho*R) for rho*R > 3 g/cm^2 (propagating burn)
 *
 * @param rho_R     areal density [g/cm^2]
 * @return          target gain G
 */
double icf_gain_from_rhoR(double rho_R) {
    if (rho_R <= 1.0) return 0.01 * rho_R * rho_R;
    if (rho_R <= 3.0) return 0.01 * rho_R * rho_R;
    return 0.09 * exp(rho_R - 3.0);
}

/**
 * inertial_confinement_driver_efficiency — ICF driver efficiency
 *
 * Different ICF drivers have different wall-plug efficiencies:
 *   Flashlamp-pumped Nd:glass (NIF): eta ~ 0.5-1%
 *   Diode-pumped solid state (LIFE): eta ~ 10-15%
 *   Krypton Fluoride (KrF) laser: eta ~ 5-7%
 *   Heavy ion beam: eta ~ 25-35%
 *   Z-pinch (pulsed power): eta ~ 15-25%
 *
 * This function gives typical driver efficiency.
 */
double icf_driver_efficiency(int driver_type) {
    switch (driver_type) {
        case 0: return 0.008;   /* NIF-like flashlamp */
        case 1: return 0.12;    /* DPSSL (LIFE) */
        case 2: return 0.30;    /* Heavy ion */
        case 3: return 0.20;    /* Z-pinch */
        default: return 0.05;
    }
}

/* ================================================================
 * L9: Sustainability and Environment — 3 functions
 * ================================================================ */

/**
 * fusion_waste_disposal_rating — Waste disposal rating
 *
 * Definition: WDR = (contact dose rate after 50 yr) / (clearance limit)
 * Physics: reduced activation materials can achieve WDR < 1
 *          (shallow land burial) or WDR < 0.1 (clearance).
 *          Fusion has advantages over fission for waste.
 *
 * @param dose_rate_50yr      dose rate after 50 years cooling [microSv/hr]
 * @param clearance_limit     regulatory clearance limit [microSv/hr]
 * @return                    waste disposal rating
 */
double waste_disposal_rating(double dose_rate_50yr, double clearance_limit) {
    if (clearance_limit <= 0.0) return 1e10;
    return dose_rate_50yr / clearance_limit;
}

/**
 * fusion_safety_advantage_factor — Inherent fusion safety
 *
 * Physics: fusion has inherent safety advantages over fission:
 *   1. No nuclear chain reaction (no criticality accident)
 *   2. Limited fuel inventory in vessel (~1 g T at any time)
 *   3. Afterheat density ~10x lower than fission (easier cooling)
 *   4. No long-lived actinides (reduced long-term waste)
 *
 * This simplified metric compares accident source term.
 * Lower is safer.
 */
double fusion_safety_factor(double inventory_MJ, double afterheat_fraction) {
    return inventory_MJ * afterheat_fraction;
}

/**
 * carbon_emission_lifecycle — Fusion lifecycle carbon intensity [gCO2/kWh]
 *
 * Estimated lifecycle carbon emissions for fusion:
 *   ITER (experimental): ~20-50 gCO2/kWh (not optimized)
 *   DEMO: ~10-20 gCO2/kWh
 *   Commercial fusion: ~5-10 gCO2/kWh (target)
 *
 * Compare: Solar PV ~ 40 gCO2/kWh, Wind ~ 12 gCO2/kWh,
 *          Nuclear fission ~ 12 gCO2/kWh, Coal ~ 820 gCO2/kWh
 */
double fusion_carbon_intensity(double construction_CO2, double operation_CO2,
                                double lifetime_E) {
    if (lifetime_E <= 0.0) return 1e10;
    return (construction_CO2 + operation_CO2) / lifetime_E;
}