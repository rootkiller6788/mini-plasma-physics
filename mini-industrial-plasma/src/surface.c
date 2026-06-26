/**
 * surface.c - Plasma-Surface Interactions for Industrial Processing
 *
 * Reference: Lieberman §9 (Etching), §14-16 (PECVD/Sputtering)
 *   Sigmund, Phys. Rev. 184, 383 (1969) - Sputtering theory
 *   Coburn & Winters, J. Appl. Phys. 50, 3189 (1979) - Ion-enhanced etching
 *   Steinbruchel, Appl. Phys. Lett. 55, 1960 (1989) - RIE modeling
 *   Yamamura & Tawara, At. Data Nucl. Data Tables 62, 149 (1996)
 *
 * Course: MIT 22.611, Berkeley EECS 245
 */

#include "surface.h"
#include "sheath.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

/* ============================================================
 * L1: Sputtering Yields
 * ============================================================ */

/*
 * Sigmund linear cascade sputtering theory (1969).
 *
 * For amorphous/polycrystalline targets with keV ion energies,
 * the sputtering yield (atoms/ion) in the linear cascade regime:
 *
 *   Y(E) = 0.042 * alpha * S_n(E) / U_s
 *
 * where:
 *   alpha = efficiency factor (depends on M_target/M_ion ratio)
 *   S_n(E) = nuclear stopping cross section
 *   U_s = surface binding energy [eV]
 *
 * The nuclear stopping power can be approximated:
 *   S_n(E) = 4*pi*a*Z1*Z2*e^2 * M1/(M1+M2) * s_n(epsilon)
 *
 * For simplicity, we use the empirical Yamamura formula
 * which is more practical for industrial applications.
 *
 * Complexity: O(1)
 */
double sigmund_sputtering_yield(double E_ion_eV, double M_ion,
                                double M_target, double U_surface_eV,
                                double alpha_efficiency)
{
    if (E_ion_eV <= 0.0 || U_surface_eV <= 0.0) return 0.0;

    /* Nuclear stopping power (simplified Lindhard-Scharff):
     * S_n ~ (Z1*Z2)^{2/3} * E^{1/2} / (M1/M2)^{1/2} */
    double mass_ratio = M_target / M_ion;
    if (mass_ratio <= 0.0) return 0.0;

    /* Simplified nuclear stopping using extrapolation */
    double S_n = 1.0e-15 * sqrt(E_ion_eV / mass_ratio); /* rough estimate */

    /* Sigmund formula: Y = 0.042 * alpha * S_n / U_s */
    double yield = 0.042 * alpha_efficiency * S_n / U_surface_eV;
    return yield > 0.0 ? yield : 0.0;
}

/*
 * Bohdansky empirical sputtering formula.
 *
 * Simple threshold-based fit widely used in fusion and
 * plasma processing:
 *
 *   Y(E) = Q * (1 - (E_th/E)^{2/3}) * (1 - E_th/E)^2  for E > E_th
 *
 * where Q is a material-dependent yield factor.
 * This captures the threshold behavior near E_th more
 * accurately than Sigmund theory at low energies.
 *
 * Reference: Bohdansky et al., Nucl. Instrum. Methods B2, 587 (1984)
 */
double bohdansky_sputtering_yield(double E_ion_eV, double E_threshold,
                                  double Q_yield, double alpha_empirical)
{
    if (E_ion_eV <= E_threshold || E_threshold <= 0.0) return 0.0;

    double r = E_threshold / E_ion_eV;
    double yield = Q_yield * (1.0 - pow(r, 2.0/3.0))
                 * (1.0 - r) * (1.0 - r);

    (void)alpha_empirical;
    return yield > 0.0 ? yield : 0.0;
}

/*
 * Yamamura semi-empirical sputtering formula.
 *
 * This is the most widely used formula for ion energies
 * from ~10 eV to ~10 keV. It combines the Lindhard-Scharff
 * electronic stopping with an empirical nuclear stopping
 * parametrization:
 *
 *   Y(E) = 0.042 * Q(Z2) * alpha*(M2/M1) * S_n(E) / (U_s*(1 + Gamma*k_e*epsilon^{0.3}))
 *
 * where epsilon is the reduced energy, S_n is the reduced
 * nuclear stopping cross section, and Q, alpha, Gamma, k_e
 * are empirical parameters.
 *
 * Reference: Yamamura & Tawara, ADNDT 62, 149 (1996)
 */
double yamamura_sputtering_yield(double E_ion_eV, double M_ion,
                                  double M_target, double U_s,
                                  double Q_param, double s_param)
{
    if (E_ion_eV <= 0.0 || U_s <= 0.0 || M_ion <= 0.0 || M_target <= 0.0)
        return 0.0;

    /* Reduced energy (Lindhard) */
    double a_L = 0.529e-10; /* Bohr radius — here just for scale */
    double Z1 = 18.0;      /* approximate Z for Ar+ */
    double Z2 = 14.0;      /* approximate Z for Si */
    double e2 = PLASMA_E_CHARGE * PLASMA_E_CHARGE;

    double epsilon = E_ion_eV * M_target
                   / ((M_ion + M_target) * Z1 * Z2 * e2 / a_L);

    /* Reduced nuclear stopping (Kr-C potential) */
    double sn;
    if (epsilon > 0.0 && epsilon < 30.0) {
        double log_eps = log(epsilon);
        sn = 0.5 * log(1.0 + 1.2288 * epsilon)
           / (epsilon + 0.1728 * sqrt(epsilon) + 0.008 * pow(epsilon, s_param));
        (void)log_eps;
    } else {
        sn = 0.0;
    }

    /* Mass-dependent efficiency */
    double mass_ratio = M_target / M_ion;
    double alpha_star;
    if (mass_ratio <= 0.5) {
        alpha_star = 0.2;
    } else {
        alpha_star = 0.3 * mass_ratio / (1.0 + 0.3 * mass_ratio);
    }

    /* Yamamura formula: Y = 0.042 * Q * alpha* * sn / (U_s * (1 + Gamma*k_e*eps^0.3)) */
    double Gamma = 0.35;
    double k_e = 0.3;
    double denom = U_s * (1.0 + Gamma * k_e * pow(epsilon, 0.3));

    if (denom <= 0.0) return 0.0;

    double yield = 0.042 * Q_param * alpha_star * sn / denom;
    return yield > 0.0 ? yield : 0.0;
}

/* ============================================================
 * L2: Ion-Enhanced Etching (Coburn-Winters Synergism)
 * ============================================================ */

/*
 * Compute etch rates including the ion-neutral synergy.
 *
 * The Coburn-Winters experiment (1979) demonstrated that
 * XeF2 gas (chemical etchant) + Ar+ bombardment (physical)
 * produce an etch rate ~10x greater than either process alone.
 *
 * Steinbruchel's model (1989):
 *   ER_total = ER_chem + ER_phys + beta * sqrt(J_ion) * theta
 *
 * where theta is the surface coverage of reactive species
 * and beta is the synergy coefficient.
 *
 * Reference: Coburn & Winters, J. Appl. Phys. 50, 3189 (1979)
 *            Steinbruchel, Appl. Phys. Lett. 55, 1960 (1989)
 */
void compute_etch_rates(EtchRates *rates, double T_e_eV, double n_e,
                        double V_sheath, double n_reactive,
                        double sticking_coeff, const char *material,
                        const char *etchant_gas)
{
    if (!rates) return;

    /* Ion flux from Bohm criterion: Gamma_i = n_e * c_s */
    double m_i = PLASMA_M_AR;
    double c_s = sqrt(PLASMA_E_CHARGE * T_e_eV / m_i);
    double ion_flux = n_e * c_s;

    /* Neutral flux: Gamma_n = (1/4) * n_reactive * v_thermal */
    double m_neutral = 88.0 * 1.66e-27; /* CF4 ~ 88 amu */
    double v_thermal = sqrt(8.0 * PLASMA_K_B * 300.0 / (M_PI * m_neutral));
    double neutral_flux = 0.25 * n_reactive * v_thermal;

    /* Ion energy = e * |V_sheath| (assuming collisionless sheath) */
    double E_ion = PLASMA_E_CHARGE * fabs(V_sheath);

    /* Surface coverage (Langmuir isotherm) */
    double theta = surface_coverage_langmuir(n_reactive * PLASMA_K_B * 300.0,
                                              sticking_coeff, 300.0, 0.5);

    /* Chemical etch rate (spontaneous, isotropic):
     * ER_chem = Gamma_n * sticking_coeff * yield_per_incident * theta */
    double ER_chemical = neutral_flux * sticking_coeff * 1.0 * theta;
    /* Convert to nm/min: 1e-9 m/min per atom removed */
    ER_chemical *= 1e9 * 60.0 * 3e-29; /* atomic volume ~Si */

    /* Physical sputtering rate (directional):
     * ER_phys = Gamma_i * Y(E_ion) * (1-theta) */
    double Y_sputter = sigmund_sputtering_yield(E_ion, m_i, 28.0 * 1.66e-27,
                                                 4.7, 0.3);
    double ER_physical = ion_flux * Y_sputter * (1.0 - theta);
    ER_physical *= 1e9 * 60.0 * 3e-29;

    /* Ion-enhanced rate (synergistic):
     * ER_ie = beta * sqrt(Gamma_i) * theta */
    double beta = 5e-14; /* synergy coefficient [m^3/s] */
    double ER_ion_enhanced = beta * sqrt(ion_flux) * theta;
    ER_ion_enhanced *= 1e9 * 60.0;

    rates->ER_chemical = ER_chemical;
    rates->ER_physical = ER_physical;
    rates->ER_ion_enhanced = ER_ion_enhanced;
    rates->ER_total = ER_chemical + ER_physical + ER_ion_enhanced;
    rates->ion_flux = ion_flux * 1e-4; /* convert to cm^-2 s^-1 */
    rates->neutral_flux = neutral_flux * 1e-4;

    (void)material;
    (void)etchant_gas;
}

/* ============================================================
 * L3: Feature Profile Evolution — ARDE & RIE Lag
 * ============================================================ */

/*
 * Aspect Ratio Dependent Etching (ARDE) / RIE lag.
 *
 * Narrower features etch slower due to:
 * 1. Neutral shadowing (depleted reactant flux at trench bottom)
 * 2. Ion shadowing (narrower angular distribution reaches bottom)
 * 3. Product re-deposition on sidewalls
 *
 * Empirical model:
 *   ER(AR) = ER(0) / (1 + AR / AR_char)
 *
 * where AR = depth/width is the aspect ratio and
 * AR_char is a characteristic ratio (~3-10 for SiO2 etching).
 */
double rie_lag_etch_rate(double open_area_rate, double aspect_ratio,
                         double sticking_coefficient)
{
    if (aspect_ratio <= 0.0) return open_area_rate;

    /* Effective transport through high-AR feature:
     * Knudsen transport regime when feature width << mean free path.
     * Flux reduction ~ (1 + 3*AR*S/8)^{-1} where S is sticking coeff. */
    double AR_char = 8.0 / (3.0 * sticking_coefficient);
    double reduction_factor = 1.0 / (1.0 + aspect_ratio / AR_char);

    return open_area_rate * reduction_factor;
}

/*
 * Neutral flux attenuation at trench bottom (Knudsen transport).
 *
 * In the molecular flow regime (Kn >> 1), the flux at the
 * trench bottom scales as:
 *
 *   F_bottom / F_top ~ (1 + 3*AR*s/4)^{-1}
 *
 * where s is the effective sticking coefficient
 * and AR = depth / width.
 */
double neutral_flux_at_trench_bottom(double aspect_ratio,
                                     double neutral_mfp, double trench_width)
{
    if (aspect_ratio <= 0.0) return 1.0;

    double Kn = neutral_mfp / trench_width; /* Knudsen number */
    if (Kn < 0.01) {
        /* Continuum regime: flux preserved */
        return 1.0;
    }

    /* Molecular flow regime */
    double transmission = 1.0 / (1.0 + 0.75 * aspect_ratio);
    (void)Kn;
    return (transmission > 0.0) ? transmission : 0.0;
}

/*
 * Ion flux at trench bottom considering angular distribution.
 *
 * Ions gain directed energy in the sheath (few eV to keV).
 * The angular spread depends on T_i/T_e and sheath collisionality:
 *
 *   sigma_theta ~ sqrt(T_i / (2*T_e))
 *
 * For collisionless sheaths: sigma_theta ~ 2-5 degrees.
 * For collisional sheaths: sigma_theta ~ 10-30 degrees.
 *
 * The fraction reaching the bottom through an AR feature:
 *   F_bottom ~ 1 - erf(AR * sigma_theta) for small angles.
 */
double ion_flux_at_trench_bottom(double aspect_ratio, double ion_energy_eV,
                                 double T_e_eV, double V_sheath)
{
    if (aspect_ratio <= 0.0) return 1.0;

    /* Angular spread of ion velocity distribution */
    double T_i_eV = 0.025; /* room temperature ions ~0.025 eV */
    double sigma_theta = sqrt(T_i_eV / (2.0 * T_e_eV));

    /* Collisional broadening */
    double E_ion = ion_energy_eV;
    double E_sheath = PLASMA_E_CHARGE * V_sheath;
    if (E_sheath > 0.0) {
        double collisionality = 1.0 - E_ion / E_sheath;
        sigma_theta += 0.1 * collisionality;
    }

    /* Fraction passing through AR feature:
     * Approximate as 1 - erf(AR * tan(sigma_theta)) */
    double arg = aspect_ratio * tan(sigma_theta);
    double transmission = 1.0 - erf(arg);
    return (transmission > 0.0) ? transmission : 0.0;
}

/* ============================================================
 * L4: PECVD Surface Kinetics
 * ============================================================ */

/*
 * PECVD deposition rate (Arrhenius model).
 *
 * In PECVD of a-Si:H from SiH4/H2 plasma:
 *
 *   R_dep = R0 * exp(-E_a / kT) * [SiH3]^n
 *
 * where E_a ~ 0.1-0.3 eV is the activation energy for
 * surface reaction, n ~ 0.5-1 is reaction order.
 *
 * Below T_sub ~ 300°C: surface reaction limited
 * Above T_sub ~ 350°C: mass transport limited
 *
 * Reference: Matsuda, J. Non-Cryst. Solids 59-60, 767 (1983)
 */
double pecvd_deposition_rate(double T_substrate, double n_SiH3,
                             double E_activation, double R0)
{
    if (T_substrate <= 0.0 || n_SiH3 <= 0.0) return 0.0;

    double kT = PLASMA_K_B * T_substrate;
    double arrhenius = exp(-E_activation / kT);

    /* Typical SiH3 reaction order ~0.5-1 for a-Si:H */
    double rate = R0 * arrhenius * sqrt(n_SiH3);
    return rate;
}

/*
 * Hydrogen content in a-Si:H films.
 *
 * The H content decreases with increasing substrate temperature
 * due to enhanced surface diffusion of H and preferential
 * etching of Si-H bonds by atomic H:
 *
 *   C_H[at%] ~ C_0 * exp(-T_sub / T_char)
 *
 * For device-quality a-Si:H (T_sub ~ 250°C): C_H ~ 8-15 at%.
 *
 * Reference: Street, "Hydrogenated Amorphous Silicon" (1991)
 */
double pecvd_hydrogen_content(double T_substrate, double deposition_rate)
{
    if (T_substrate <= 0.0) return 30.0;

    /* Empirical: H content decreases with T and increases with rate */
    double C_0 = 30.0; /* at% H at low T */
    double T_char = 250.0; /* characteristic temperature [K] */

    double thermal_factor = exp(-T_substrate / (T_char * 2.0));
    double rate_factor = (deposition_rate > 0.0) ?
                         1.0 + 0.5 * sqrt(deposition_rate) : 1.0;

    double C_H = C_0 * thermal_factor * rate_factor;
    return (C_H > 0.0) ? C_H : 1.0;
}

/*
 * PECVD film stress model.
 *
 * Total stress = thermal stress + intrinsic stress.
 *
 * Thermal stress:
 *   sigma_th = E_film * (alpha_sub - alpha_film) * (T_dep - T_meas) / (1-nu)
 *
 * Intrinsic stress: due to ion peening, hydrogen incorporation,
 * and microstructural evolution during growth.
 *   sigma_int ~ ion_energy * ion_flux / dep_rate
 *
 * For a-Si:H: intrinsic stress typically compressive (-100 to -500 MPa).
 */
double pecvd_film_stress(double T_substrate, double T_deposition,
                         double CTE_film, double CTE_substrate,
                         double E_modulus, double intrinsic_stress)
{
    double delta_alpha = CTE_substrate - CTE_film;
    double delta_T = T_deposition - T_substrate;

    /* Thermal stress (biaxial): sigma = E * dAlpha * dT / (1-nu)
     * Typical Poisson ratio nu ~ 0.25 */
    double nu = 0.25;
    double thermal_stress = E_modulus * delta_alpha * delta_T / (1.0 - nu);

    /* Total stress: positive = tensile, negative = compressive */
    return thermal_stress + intrinsic_stress;
}

/* ============================================================
 * L5: Surface Reaction Kinetics (Langmuir-Hinshelwood)
 * ============================================================ */

/*
 * Surface coverage (Langmuir isotherm).
 *
 * At steady state, the fractional coverage theta of adsorbed
 * species on the surface is:
 *
 *   theta = K * P / (1 + K * P)
 *
 * where K = s0 / (nu_des * sqrt(2*pi*m*kT)) * exp(E_des/kT)
 * is the equilibrium constant for adsorption/desorption.
 */
double surface_coverage_langmuir(double precursor_pressure,
                                  double sticking_coeff, double T_surf,
                                  double desorption_energy)
{
    if (precursor_pressure <= 0.0 || T_surf <= 0.0) return 0.0;

    double kT = PLASMA_K_B * T_surf / PLASMA_E_CHARGE; /* kT in eV */
    double K_eq = sticking_coeff * exp(desorption_energy / kT);
    double theta = K_eq * precursor_pressure / (1.0 + K_eq * precursor_pressure);

    return (theta >= 0.0 && theta <= 1.0) ? theta : (theta > 1.0 ? 1.0 : 0.0);
}

/*
 * Sticking probability (Arrhenius form).
 *
 * S = S0 * exp(-E_act / kT_surf)
 *
 * Most radicals (SiH3, CH3, CF3) have high sticking (~0.1-1)
 * on clean surfaces, decreasing with temperature as desorption
 * competes with chemisorption.
 */
double sticking_probability_arrhenius(double T_surf, double E_act, double S0)
{
    if (T_surf <= 0.0) return S0;
    double kT = PLASMA_K_B * T_surf;
    return S0 * exp(-E_act / kT);
}

/* ============================================================
 * L6: Reactive Ion Etching (RIE) - Full Process Model
 * ============================================================ */

/*
 * Run a complete RIE process model.
 *
 * This simulates a typical RIE recipe including:
 * - Gas phase chemistry (fragment formation)
 * - Ion flux and energy at the wafer
 * - Chemical and physical etch components
 * - ARDE effects for given feature geometry
 *
 * Output: EtchModel with predicted etch rates, selectivity,
 * anisotropy, and uniformity.
 */
void run_rie_process_model(EtchModel *model,
                           const char *material, const char *gas_mix,
                           double P_rf, double pressure, double flow,
                           double T_wall, double T_substrate,
                           double process_time_min)
{
    if (!model) return;

    /* Typical CCP RIE parameters:
     * n_e ~ 1e16 m^-3, T_e ~ 3 eV, V_sheath ~ -200 V */
    double n_e = 1e16;
    double T_e = 3.0;
    double V_sheath = -200.0;

    /* Scale with RF power: n_e ~ sqrt(P_rf) approximately */
    n_e *= sqrt(P_rf / 500.0);
    V_sheath *= sqrt(P_rf / 500.0);

    /* Ion flux at wafer */
    double m_i = PLASMA_M_AR;
    double c_s = sqrt(PLASMA_E_CHARGE * T_e / m_i);
    double ion_flux = n_e * c_s;

    /* Ion energy at wafer */
    double E_ion = PLASMA_E_CHARGE * fabs(V_sheath);

    /* Physical sputtering yield */
    double Y_phys;
    if (strstr(material, "Si") || strstr(material, "silicon")) {
        Y_phys = yamamura_sputtering_yield(E_ion, 40.0*1.66e-27,
                  28.0*1.66e-27, 4.7, 1.0, 2.5);
    } else if (strstr(material, "SiO") || strstr(material, "oxide")) {
        Y_phys = yamamura_sputtering_yield(E_ion, 40.0*1.66e-27,
                  60.0*1.66e-27, 5.0, 0.8, 2.5);
    } else {
        Y_phys = sigmund_sputtering_yield(E_ion, 40.0*1.66e-27,
                  28.0*1.66e-27, 5.0, 0.3);
    }

    /* Chemical etch component (depends on gas chemistry) */
    double ER_chem = 0.0;
    if (strstr(gas_mix, "SF6") || strstr(gas_mix, "CF4")) {
        /* Fluorine-based etching of Si: ~500 nm/min pure chemical */
        ER_chem = 500.0 * (pressure / 10.0); /* scales with etchant flux */
    } else if (strstr(gas_mix, "Cl2") || strstr(gas_mix, "BCl3")) {
        /* Chlorine-based: slower, more anisotropic */
        ER_chem = 100.0 * (pressure / 10.0);
    } else if (strstr(gas_mix, "O2")) {
        /* Oxygen plasma: photoresist ashing ~200 nm/min */
        ER_chem = 200.0 * (pressure / 10.0);
    }

    /* Physical sputter rate: Y_phys * Gamma_i * (1/atomic_density) */
    double n_atomic = 5e28; /* atoms/m^3 for Si */
    double ER_phys = Y_phys * ion_flux / n_atomic * 1e9 * 60.0; /* nm/min */

    /* Ion-enhanced synergy (Steinbruchel) */
    double beta = 1e-14;
    double ER_ion_enh = beta * sqrt(ion_flux) * 1e9 * 60.0;

    double ER_total = ER_chem + ER_phys + ER_ion_enh;

    /* Anisotropy: A = 1 - (lateral_rate / vertical_rate) */
    double lateral_rate = ER_chem * 0.3; /* chemical component is partially isotropic */
    double vertical_rate = ER_phys + ER_ion_enh + ER_chem * 0.7;
    double anisotropy;
    if (vertical_rate > 0.0) {
        anisotropy = 1.0 - lateral_rate / vertical_rate;
    } else {
        anisotropy = 0.0;
    }

    /* Selectivity to mask (photoresist):
     * Typically 2:1 to 5:1 for oxide etching with fluorocarbon */
    double sel_mask = 3.0;

    /* Uniformity: worse at higher pressure (more collisions -> less directional) */
    double uniformity = 3.0 + 2.0 * (pressure / 10.0);

    model->etch_rate_chemical = ER_chem;
    model->etch_rate_physical = ER_phys;
    model->etch_rate_total = ER_total;
    model->anisotropy = anisotropy;
    model->selectivity_mask = sel_mask;
    model->selectivity_etchstop = 10.0;
    model->uniformity_3sigma = uniformity;
    model->aspect_ratio = 5.0; /* typical feature AR */

    (void)flow;
    (void)T_wall;
    (void)T_substrate;
    (void)process_time_min;
}

/* ============================================================
 * L7: Wafer-Scale Uniformity
 * ============================================================ */

/*
 * Compute 3-sigma non-uniformity from radial measurements.
 */
double wafer_etch_uniformity(double *radial_rates, int n_points,
                             double wafer_radius)
{
    if (!radial_rates || n_points < 2) return 0.0;

    double mean = 0.0, std = 0.0;

    for (int i = 0; i < n_points; i++) {
        mean += radial_rates[i];
    }
    mean /= n_points;

    for (int i = 0; i < n_points; i++) {
        double diff = radial_rates[i] - mean;
        std += diff * diff;
    }
    std = sqrt(std / (n_points - 1));

    if (mean > 0.0) {
        return 300.0 * std / mean; /* 3-sigma percentage */
    }
    (void)wafer_radius;
    return 0.0;
}

double pecvd_thickness_uniformity(double *radial_thickness, int n_points,
                                   double wafer_radius)
{
    /* Same formula as etch uniformity, different name for semantic clarity */
    return wafer_etch_uniformity(radial_thickness, n_points, wafer_radius);
}

/* ============================================================
 * L8: Atomic Layer Etching (ALE)
 * ============================================================ */

/*
 * ALE etch depth per cycle.
 *
 * ALE is a cyclic process:
 * 1. Chemisorption of precursor (self-limited to ~1 monolayer)
 * 2. Ion bombardment to activate/remove the reacted layer
 *
 * EPC (etch per cycle) = theta_max * S * f(E_ion)
 *
 * where theta_max is the saturation coverage (~1 ML),
 * S is the number of atoms removed per incident ion,
 * and f(E_ion) captures the ion energy dependence.
 *
 * Typical ALE: ~0.1-1 nm/cycle, infinite selectivity to
 * underlying material (self-limited).
 *
 * Reference: George, Chem. Rev. 110, 111 (2010)
 *            Kanarik et al., JVST A 33, 020802 (2015)
 */
double ale_etch_per_cycle(double T_substrate, double ion_energy_eV,
                          double precursor_coverage, double E_sp)
{
    if (T_substrate <= 0.0 || ion_energy_eV <= 0.0) return 0.0;

    /* Saturation coverage ~ 1 monolayer */
    double theta_max = 1.0;

    /* Ion energy window: must be above sputtering threshold of
     * reacted layer but below threshold of underlying material */
    double E_reacted = E_sp * 0.5;  /* modified layer easier to remove */
    double E_underlying = E_sp;     /* bulk material threshold */

    double f_E;
    if (ion_energy_eV < E_reacted) {
        f_E = 0.0; /* no removal */
    } else if (ion_energy_eV < E_underlying) {
        f_E = (ion_energy_eV - E_reacted) / (E_underlying - E_reacted);
    } else {
        f_E = 1.0; /* saturated, but losing selectivity */
    }

    /* Synergy parameter:
     * S = 1 for ideal ALE (exactly 1 layer removed per cycle)
     * S < 1 for incomplete removal
     * S > 1 for partial sputtering of underlying material */
    double S = 0.3; /* typical: 0.1-0.5 ML/cycle */

    double EPC = theta_max * S * f_E * precursor_coverage;
    /* Convert to nm: 1 ML ~ 0.3 nm */
    double ML_to_nm = 0.3;
    return EPC * ML_to_nm;
}

/*
 * ALE synergy parameter.
 *
 * Synergy = ER_ALE / (ER_spontaneous + ER_ion_alone)
 *
 * For true ALE: Synergy >> 1 (the combined process is much
 * faster than either process alone).
 * For simple ion-enhanced etching: Synergy ~ 2-5.
 * For purely additive process: Synergy = 1.
 */
int ale_synergy_parameter(double ER_ale, double ER_spontaneous,
                          double ER_ion_alone)
{
    double denom = ER_spontaneous + ER_ion_alone;
    if (denom <= 0.0) return 0;

    double synergy = ER_ale / denom;
    /* Return:
     * 2 if clearly ALE regime (synergy >> 1)
     * 1 if moderate synergy
     * 0 if no synergy or additive */
    if (synergy > 10.0) return 2;
    if (synergy > 2.0) return 1;
    return 0;
}
