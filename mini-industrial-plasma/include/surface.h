#ifndef SURFACE_H
#define SURFACE_H

#include "plasma_types.h"
#include <stddef.h>

/*
 * surface.h — Plasma-Surface Interactions
 *
 * Reference: Lieberman §9 (Etching), §15 (PECVD), §16 (Sputtering)
 *   Coburn & Winters, J. Appl. Phys. 50, 3189 (1979)
 *   Steinbruchel, Appl. Phys. Lett. 55, 1960 (1989)
 * Course: MIT 22.611, Berkeley EECS 245
 */

/*
 * L1: Sputtering Yield — atoms removed per incident ion
 * Sigmund linear cascade theory (1969):
 * Y(E) = 0.042 * alpha * S_n(E) / U_s
 */

double sigmund_sputtering_yield(double E_ion_eV, double M_ion, double M_target,
                                double U_surface_eV, double alpha_efficiency);

/* Bohdansky formula (empirical fit for low-energy sputtering) */
double bohdansky_sputtering_yield(double E_ion_eV, double E_threshold,
                                  double Q_yield, double alpha_empirical);

/* Yamamura semi-empirical formula (most widely used for ion energies < 10 keV) */
double yamamura_sputtering_yield(double E_ion_eV, double M_ion, double M_target,
                                  double U_s, double Q_param, double s_param);

/*
 * L2: Ion-enhanced etching — the ion-neutral synergy
 * Coburn-Winters experiment: XeF2 + Ar+ on Si gives 10x faster
 * etching than either process alone.
 */

typedef struct {
    double ER_chemical;    /* chemical (spontaneous) rate [nm/min] */
    double ER_physical;    /* physical sputtering rate [nm/min] */
    double ER_ion_enhanced;/* synergistic ion-enhanced rate [nm/min] */
    double ER_total;
    double ion_flux;       /* ion flux [cm^-2 s^-1] */
    double neutral_flux;   /* reactive neutral flux [cm^-2 s^-1] */
} EtchRates;

/* Ion-enhanced etching yield model (Steinbruchel):
 * ER_total = ER_chem + ER_phys + beta * sqrt(J_ion) * ER_chem */
void compute_etch_rates(EtchRates *rates, double T_e_eV, double n_e,
                        double V_sheath, double n_reactive,
                        double sticking_coeff, const char *material,
                        const char *etchant_gas);

/*
 * L3: Feature profile evolution — aspect ratio dependent etching (ARDE)
 * RIE lag: narrower features etch slower
 */
double rie_lag_etch_rate(double open_area_rate, double aspect_ratio,
                         double sticking_coefficient);

/* Neutral shadowing in high aspect ratio features */
double neutral_flux_at_trench_bottom(double aspect_ratio,
                                     double neutral_mfp, double trench_width);

/* Ion angular distribution at feature bottom */
double ion_flux_at_trench_bottom(double aspect_ratio, double ion_energy_eV,
                                 double T_e_eV, double V_sheath);

/*
 * L4: Surface kinetics for PECVD
 */

typedef struct {
    double dep_rate;         /* nm/min */
    double SiH4_frac;        /* silane fraction in feed */
    double H2_dilution;      /* hydrogen dilution ratio */
    double T_substrate;      /* substrate temperature [K] */
    double P_rf;             /* RF power density [W/cm^2] */
    double film_stress;      /* MPa */
    double H_content_at_pct; /* at% hydrogen in film */
    double defect_density;   /* cm^-3 */
} PECVDState;

/* Arrhenius deposition rate: R = R0 * exp(-E_a/kT_sub) * [SiH3]^n */
double pecvd_deposition_rate(double T_substrate, double n_SiH3,
                             double E_activation, double R0);

/* Hydrogen content in a-Si:H films vs T_substrate */
double pecvd_hydrogen_content(double T_substrate, double deposition_rate);

/* Film stress model: thermal + intrinsic components */
double pecvd_film_stress(double T_substrate, double T_deposition,
                         double CTE_film, double CTE_substrate,
                         double E_modulus, double intrinsic_stress);

/*
 * L5: Surface reaction kinetics — Langmuir-Hinshelwood adsorption
 */
double surface_coverage_langmuir(double precursor_pressure,
                                  double sticking_coeff, double T_surf,
                                  double desorption_energy);

double sticking_probability_arrhenius(double T_surf, double E_act,
                                      double S0);

/*
 * L6: Reactive ion etching (RIE) — full process model
 */
void run_rie_process_model(EtchModel *model,
                           const char *material, const char *gas_mix,
                           double P_rf, double pressure, double flow,
                           double T_wall, double T_substrate,
                           double process_time_min);

/*
 * L7: Wafer-scale uniformity model
 */
double wafer_etch_uniformity(double *radial_rates, int n_points,
                             double wafer_radius);

double pecvd_thickness_uniformity(double *radial_thickness, int n_points,
                                   double wafer_radius);

/*
 * L8: Atomic Layer Etching (ALE) — self-limited, cycle-based process
 */
double ale_etch_per_cycle(double T_substrate, double ion_energy_eV,
                          double precursor_coverage, double E_sp);

int ale_synergy_parameter(double ER_ale, double ER_spontaneous,
                          double ER_ion_alone);

#endif /* SURFACE_H */
