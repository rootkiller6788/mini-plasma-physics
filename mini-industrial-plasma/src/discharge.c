/**
 * discharge.c - Plasma Discharge Models for Industrial Processing
 *
 * Reference: Lieberman & Lichtenberg 2005, Ch.9-13
 *   Chabert & Braithwaite, "Physics of Radio-Frequency Plasmas" (2011)
 *   Paschen, Wied. Ann. 37, 69 (1889)
 *   Townsend, "Electricity in Gases" (1915)
 *
 * Course: MIT 22.611, Stanford EE 414, Berkeley EECS 245
 */

#include "discharge.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

/* ============================================================
 * L1: Townsend Discharge Theory
 * ============================================================ */

/*
 * Townsend first ionization coefficient alpha.
 * The number of ionizing collisions per unit length along the E-field:
 *   alpha = A * p * exp(-B*p / E)
 *
 * Where A [1/(m·Pa)] and B [V/(m·Pa)] are gas-dependent constants.
 * For Argon: A ≈ 12.2, B ≈ 180 (in units of 1/(m·Pa), V/(m·Pa)).
 * This empirical law was established by Townsend in 1902-1915.
 *
 * Reference: Townsend, Phil. Mag. 3, 557 (1902)
 */
double townsend_alpha(double pressure, double E_field, double A, double B)
{
    if (pressure <= 0.0 || E_field <= 0.0) return 0.0;
    double exponent = -B * pressure / E_field;
    if (exponent < -500.0) return 0.0;
    return A * pressure * exp(exponent);
}

/*
 * Effective secondary electron emission coefficient.
 * Accounts for electron multiplication in the gap:
 *   gamma_eff = gamma_se * (exp(alpha*d) - 1)
 *
 * where gamma_se is the secondary electron yield per ion impact.
 */
double townsend_gamma_effective(double gamma_se, double alpha, double d)
{
    if (gamma_se <= 0.0 || alpha <= 0.0 || d <= 0.0) return 0.0;
    double electron_mult = exp(alpha * d) - 1.0;
    if (electron_mult < 0.0) return 0.0;
    return gamma_se * electron_mult;
}

/*
 * Townsend breakdown condition:
 *   gamma_se * (exp(alpha * d) - 1) = 1
 *
 * Returns 1 if breakdown occurs, 0 otherwise.
 * The `breakdown_threshold` output gives the ratio
 * gamma*(exp(alpha*d)-1) which must be >= 1 for breakdown.
 */
int townsend_breakdown_condition(double gamma_se, double alpha, double d,
                                 double *breakdown_threshold)
{
    double electron_mult = exp(alpha * d);
    double threshold = gamma_se * (electron_mult - 1.0);

    if (breakdown_threshold) {
        *breakdown_threshold = threshold;
    }

    return (threshold >= 1.0) ? 1 : 0;
}

/* ============================================================
 * L2: Paschen Curve - Breakdown Voltage vs p·d
 * ============================================================ */

/*
 * Theorem (Paschen, 1889): The DC breakdown voltage V_br of a
 * parallel-plate discharge depends only on the product p·d
 * (pressure × gap distance) and the gas type:
 *
 *   V_br = B * p * d / ln(A * p * d / ln(1 + 1/gamma_se))
 *
 * This produces a characteristic U-shaped curve with a minimum
 * (the Paschen minimum). At very low p·d, few collisions occur;
 * at very high p·d, electrons lose energy too quickly to ionize.
 *
 * Derivation: Combine Townsend criterion gamma*(exp(alpha*d)-1) = 1
 * with alpha = A*p*exp(-B*p/E) and E = V/d.
 */
double paschen_breakdown_voltage(double pd, double A, double B,
                                 double gamma_se)
{
    if (pd <= 0.0 || gamma_se <= 0.0) return INFINITY;

    double denom = log(A * pd / log(1.0 + 1.0/gamma_se));
    if (denom <= 0.0) return INFINITY;

    return B * pd / denom;
}

/*
 * Paschen minimum voltage: V_min = (e*B/A) * ln(1 + 1/gamma_se)
 * where e = exp(1) = 2.71828...
 *
 * This occurs at pd_min = (e/A) * ln(1 + 1/gamma_se).
 *
 * For Ar (A≈13.6/Pa/m, B≈235 V/Pa/m, gamma≈0.1):
 *   V_min ≈ 137 V at pd_min ≈ 0.5 Pa·m.
 */
double paschen_minimum_voltage(double A, double B, double gamma_se)
{
    if (A <= 0.0 || gamma_se <= 0.0) return INFINITY;
    return (M_E * B / A) * log(1.0 + 1.0/gamma_se);
}

double paschen_pd_at_minimum(double A, double B, double gamma_se)
{
    if (A <= 0.0 || gamma_se <= 0.0) return INFINITY;
    (void)B;
    return (M_E / A) * log(1.0 + 1.0/gamma_se);
}

/*
 * Compute full Paschen curve over a range of p·d values.
 *
 * The curve reveals three discharge regimes:
 * - Left branch (low pd): electron mean free path >> d,
 *   few ionizing collisions. Used in vacuum electronics.
 * - Minimum: optimal pd for lowest breakdown voltage.
 *   Critical for designing gas discharge lamps and plasma sources.
 * - Right branch (high pd): frequent collisions impede
 *   electron acceleration. Used in high-pressure arcs.
 */
void compute_paschen_curve(PaschenCurve *curve, double A, double B,
                           double gamma_se, double pd_min_param, double pd_max,
                           int n_points)
{
    if (!curve || n_points <= 0) return;

    curve->n_points = n_points;
    curve->pd_values = (double*)malloc(n_points * sizeof(double));
    curve->V_breakdown = (double*)malloc(n_points * sizeof(double));

    if (!curve->pd_values || !curve->V_breakdown) {
        free(curve->pd_values);
        free(curve->V_breakdown);
        curve->n_points = 0;
        return;
    }

    double log_pd_min = log(pd_min_param);
    double log_pd_max = log(pd_max);
    double dlog = (log_pd_max - log_pd_min) / (n_points - 1);

    double V_min = INFINITY;
    double pd_at_V_min = 0.0;

    for (int i = 0; i < n_points; i++) {
        double pd = exp(log_pd_min + i * dlog);
        curve->pd_values[i] = pd;
        curve->V_breakdown[i] = paschen_breakdown_voltage(pd, A, B, gamma_se);

        if (curve->V_breakdown[i] < V_min) {
            V_min = curve->V_breakdown[i];
            pd_at_V_min = pd;
        }
    }

    curve->V_min = V_min;
    curve->pd_min = pd_at_V_min;
    curve->gamma_se = gamma_se;
    curve->A_coeff = A;
    curve->B_coeff = B;
}

/* ============================================================
 * L2: DC Discharge Regimes
 * ============================================================ */

/*
 * Normal glow current density scaling.
 * In the normal glow regime, the current density is approximately
 * constant (the discharge adjusts its area to carry the current):
 *
 *   j_n ≈ C_gas * p^2
 *
 * For argon: C_gas ≈ 2e-4 A/(m^2·Pa^2) (experimental).
 *
 * Reference: von Engel, "Ionized Gases" (1965)
 */
double normal_glow_current_density(double pressure, double gas_constant)
{
    if (pressure <= 0.0) return 0.0;
    return gas_constant * pressure * pressure;
}

/*
 * Abnormal glow voltage-current characteristic.
 *
 * In the abnormal glow regime (the entire cathode is covered),
 * voltage increases with current:
 *
 *   V = V_n * (j / j_n)^{1/n}
 *
 * where n ≈ 3-5 depending on gas/electrode material.
 *
 * This is due to the cathode fall adjusting to provide the
 * required current density via ion-induced secondary emission.
 */
double abnormal_glow_voltage(double current_density, double pressure,
                             double d, double gas_constant,
                             double exponent)
{
    double j_n = normal_glow_current_density(pressure, gas_constant);
    if (j_n <= 0.0 || exponent <= 0.0) return 0.0;
    (void)d;
    /* V ≈ V_normal * (j/j_n)^{1/exponent}, using scaling:
     * V_normal ~ several hundred volts. Use Paschen minimum as estimate. */
    double V_n = 200.0; /* approximate normal cathode fall for Ar */
    return V_n * pow(current_density / j_n, 1.0 / exponent);
}

/* ============================================================
 * L3: CCP - Capacitively Coupled Plasma
 * ============================================================ */

/*
 * DC self-bias in asymmetric CCP.
 *
 * In a geometrically asymmetric discharge (powered electrode
 * smaller than grounded), a negative DC self-bias develops:
 *
 *   V_dc ≈ -V_rf * (A_ground / A_powered)^k
 *
 * with k ≈ 1 experimentally, but theoretically between
 * 1 and 4 depending on sheath model.
 *
 * The self-bias arises because the smaller electrode has a
 * smaller sheath capacitance, requiring a larger voltage
 * drop to balance equal ion fluxes in steady state.
 *
 * This is the operating principle of RIE (Reactive Ion Etching):
 * the DC bias accelerates ions toward the wafer for anisotropic etching.
 */
double ccp_self_bias(double V_rf, double area_powered, double area_ground)
{
    if (area_powered <= 0.0 || area_ground <= 0.0) return 0.0;

    double area_ratio = area_ground / area_powered;
    /* Experimental scaling: k ~ 1 to 2.5.
     * Lieberman Eq. 11.2.40: V_dc/V_rf = sin(pi*theta/2)
     * with theta ~ A_powered/A_ground. Simplified: k=1.5 */
    double k = 1.5;
    double V_dc = -V_rf * pow(area_ratio, k);
    return V_dc;
}

/*
 * Sheath heating power per unit area.
 *
 * Lieberman Eq. 11.2.31:
 *   P_sh = 0.5 * n_s * e * c_s * V_sh * A
 *
 * Sheath heating (also called stochastic or Fermi acceleration)
 * is the primary electron heating mechanism in low-pressure CCP.
 * Electrons gain energy by reflecting from the moving sheath edge,
 * analogous to a ball bouncing off a moving wall.
 */
double ccp_sheath_heating_power(double n_sheath, double T_e_eV,
                                double V_sheath, double m_i, double area)
{
    double e = PLASMA_E_CHARGE;
    double c_s = sqrt(e * T_e_eV / m_i);
    double P = 0.5 * n_sheath * e * c_s * V_sheath * area;
    return P;
}

/*
 * Ohmic (collisional) heating power in bulk plasma.
 *
 * P_ohm = J^2 * d / sigma_dc
 * sigma_dc = e^2 * n_e / (m_e * nu_en)
 *
 * Dominates at high pressure where electron-neutral collisions
 * are frequent enough to randomize the directed energy gained
 * from the RF field.
 */
double ccp_ohmic_heating_power(double current_density, double gap,
                               double n_e, double nu_en)
{
    if (n_e <= 0.0 || nu_en <= 0.0) return 0.0;

    double sigma_dc = PLASMA_E_CHARGE * PLASMA_E_CHARGE * n_e
                    / (PLASMA_M_E * nu_en);
    if (sigma_dc <= 0.0) return 0.0;

    return current_density * current_density * gap / sigma_dc;
}

/*
 * CCP equilibrium plasma density from global power balance.
 *
 * The discharge reaches steady state when absorbed RF power
 * balances electron and ion losses to the walls. For a
 * cylindrical discharge (radius R, gap L):
 *
 *   P_abs = e * n_e * u_B * A_eff * (E_c + E_e + V_sheath)
 *
 * where u_B = sqrt(e*T_e/m_i) is the Bohm velocity,
 * A_eff is effective loss area,
 * E_c ~ 5-10*T_e is collisional energy loss per e-ion pair,
 * E_e ~ 2*T_e is kinetic energy carried to walls.
 *
 * Solving for n_e:
 *   n_e = P_abs / [e * u_B * A_eff * (E_c + E_e + V_sheath)]
 *
 * Reference: Lieberman §10.2 (Global model)
 */
double ccp_equilibrium_density(double P_abs, double V_rf,
                               double T_e_eV, double m_i,
                               double area_powered, double area_ground,
                               double pressure, double gap)
{
    double e = PLASMA_E_CHARGE;
    double u_B = sqrt(e * T_e_eV / m_i);

    /* Effective loss area: plasma contacts all surfaces.
     * Approximate for parallel-plate geometry */
    double total_area = area_powered + area_ground;
    double A_eff = 0.5 * total_area; /* sheath edge area ~ half total */

    /* Energy loss per e-ion pair created */
    double E_c = 5.0 * T_e_eV; /* collisional loss (approximation for Ar) */
    double E_e = 2.0 * T_e_eV; /* kinetic energy to walls */
    double E_loss = E_c + E_e + fabs(V_rf);

    double denom = e * u_B * A_eff * E_loss;
    if (denom <= 0.0) return 0.0;

    double n_e = P_abs / denom;

    (void)pressure;
    (void)gap;
    return n_e;
}

/* ============================================================
 * L3: ICP - Inductively Coupled Plasma
 * ============================================================ */

/*
 * Collisionless skin depth for RF field penetration.
 *
 * delta = c / omega_pe
 *
 * For typical ICP conditions (n_e ~ 1e17 m^-3):
 *   omega_pe ~ 1.8e10 rad/s
 *   delta ~ 1.7 cm
 *
 * When delta << chamber radius, the power is deposited
 * in a thin annular region near the chamber wall (skin effect).
 * When delta > chamber radius, the discharge is "weakly
 * penetrated" and volume-averaged heating occurs.
 */
double icp_skin_depth(double n_e)
{
    if (n_e <= 0.0) return INFINITY;
    double omega_pe = plasma_freq_electron(n_e);
    if (omega_pe <= 0.0) return INFINITY;
    return 299792458.0 / omega_pe; /* c/omega_pe */
}

/*
 * ICP plasma resistance (transformer model).
 *
 * The plasma acts as a single-turn secondary of a transformer
 * whose primary is the N-turn coil:
 *
 *   R_p = (2*pi*R)/ (sigma * delta * L) * N^2
 *
 * where sigma_dc = e^2*n_e/(m_e*nu_en) is the DC conductivity.
 * For collisionless regime, nu_en is replaced by an effective
 * collision frequency from anomalous skin effect.
 *
 * Reference: Piejak et al., Plasma Sources Sci. Technol. 1, 179 (1992)
 */
double icp_plasma_resistance(double n_e, double radius, double height,
                             int N_coils, double freq)
{
    if (n_e <= 0.0 || radius <= 0.0 || height <= 0.0) return 0.0;

    double e = PLASMA_E_CHARGE;
    double m_e = PLASMA_M_E;

    /* Estimate effective collision frequency */
    double nu_eff = 1e7; /* ~10 MHz for low-pressure Ar at 1 Pa */

    double sigma = e * e * n_e / (m_e * nu_eff);
    double delta = icp_skin_depth(n_e);
    if (delta <= 0.0) return 0.0;

    /* Cylindrical geometry */
    double N2 = (double)(N_coils * N_coils);
    double R_p = (2.0 * M_PI * radius) / (sigma * delta * height) * N2;

    (void)freq;
    return R_p;
}

/*
 * E-H mode transition threshold power.
 *
 * ICP discharges exhibit a hysteresis: at low power they operate
 * in E-mode (capacitive coupling, low density), transitioning
 * to H-mode (inductive coupling, high density) above a threshold:
 *
 *   P_th ∝ (pressure)^{1/2} * frequency^{-1}
 *
 * This scaling arises from the requirement that the induced
 * electric field exceeds a critical value for sustaining
 * ionization via inductive coupling.
 *
 * Reference: Kortshagen et al., Phys. Rev. Lett. 76, 1170 (1996)
 */
double icp_eh_mode_threshold(double pressure, double freq, double radius,
                             double height)
{
    if (pressure <= 0.0 || freq <= 0.0) return INFINITY;
    double volume = M_PI * radius * radius * height;
    double scaling = 1.0; /* [W · sqrt(Pa) · Hz] empirical prefactor */
    return scaling * volume * sqrt(pressure) / freq;
}

/* ============================================================
 * L4: Electron Heating Mechanisms
 * ============================================================ */

/*
 * Stochastic (Fermi) heating power density at oscillating sheath.
 *
 * Electrons gain energy by colliding with the moving sheath
 * boundary, analogous to a ball bouncing between fixed and
 * moving walls (Fermi acceleration mechanism):
 *
 *   S_stoc = 0.5 * m_e * n_e * v_e^3 * (v_sheath/v_e)^2
 *           = 0.25 * m_e * n_e * v_e * v_sheath^2
 *
 * where v_e = sqrt(8e*T_e/(pi*m_e)) is the mean electron speed
 * and v_sheath ~ omega*V_sheath is the sheath oscillation velocity.
 *
 * Reference: Lieberman, IEEE Trans. Plasma Sci. 16, 638 (1988)
 */
double stochastic_heating_power_density(double n_sheath, double T_e_eV,
                                        double V_sheath, double freq)
{
    if (n_sheath <= 0.0 || T_e_eV <= 0.0) return 0.0;

    double v_e = sqrt(8.0 * PLASMA_E_CHARGE * T_e_eV / (M_PI * PLASMA_M_E));
    double v_sh = 2.0 * M_PI * freq * V_sheath; /* approximate sheath velocity */

    /* P_stoc/A = 0.25 * m_e * n_e * v_e * v_sh^2 */
    double P_per_area = 0.25 * PLASMA_M_E * n_sheath * v_e * v_sh * v_sh;
    return P_per_area;
}

/*
 * Ion transit time heating in RF sheaths.
 *
 * Ions traversing the modulated sheath gain energy that depends
 * on the ratio of ion transit time to RF period. This determines
 * the ion energy distribution function (IEDF) width:
 *
 *   P_ion = n_i * m_i * v_B^3 / 2 (collisionless limit)
 *
 * Reference: Edelberg & Aydil, J. Appl. Phys. 86, 4799 (1999)
 */
double ion_transit_heating_power(double n_i, double m_i, double V_rf,
                                 double freq, double sheath_width)
{
    if (n_i <= 0.0 || m_i <= 0.0) return 0.0;

    double e = PLASMA_E_CHARGE;
    /* Ion transit time across sheath (collisionless):
     * tau_ion = 3*s/c_s (from Child law) */
    double c_s = sqrt(e * 3.0 / m_i); /* estimate T_e ~ 3 eV */
    double tau_ion = 3.0 * sheath_width / c_s;
    double tau_rf = 1.0 / freq;

    /* Energy modulation depends on tau_ion / tau_rf */
    double ratio = tau_ion / tau_rf;
    if (ratio > 1.0) ratio = 1.0; /* cap at full modulation */

    double E_mod = e * V_rf * ratio;
    double P = n_i * c_s * E_mod; /* flux * energy per ion */

    return P > 0.0 ? P : 0.0;
}

/* ============================================================
 * L5: Global (0-D) Discharge Model
 * ============================================================ */

/*
 * Global model electron temperature.
 *
 * In steady state, the rate of ionization must balance the rate
 * of particle loss to the walls:
 *
 *   K_iz(T_e) * n_gas * V = u_B(T_e) * A_eff / V
 *
 * or equivalently:
 *
 *   K_iz(T_e) / u_B(T_e) = 1 / (n_gas * d_eff)
 *
 * where d_eff = V/A_eff is the effective plasma size.
 * Since K_iz depends exponentially on T_e, this determines T_e
 * uniquely from gas pressure and chamber geometry.
 *
 * For argon, the ionization rate coefficient can be approximated:
 *   K_iz(T_e) ≈ 2.3e-14 * T_e^{0.59} * exp(-17.44/T_e) [m^3/s]
 *
 * Reference: Lieberman §10.2, Lee & Lieberman JVST A 13, 368 (1995)
 */
double global_model_electron_temp(double pressure, double R, double L,
                                  double n_gas, double u_B_energy)
{
    /* Effective size d_eff:
     * For a cylinder (R, L): A_eff = 2*pi*R*(R*h_L + L*h_R)
     * with edge-to-center density ratios h_L, h_R.
     * Simplified: d_eff ≈ R/2 for typical aspect ratios. */
    double V = M_PI * R * R * L;
    double A_eff = 2.0 * M_PI * R * R + 2.0 * M_PI * R * L;
    double d_eff = V / A_eff;

    /* Solve K_iz/u_B = 1/(n_gas * d_eff) numerically.
     * For Ar: approximate T_e using Lieberman Eq. 10.2.9:
     *   T_e ≈ epsilon_iz / ln(K_iz0 * n_gas * d_eff / u_B) */
    double E_iz = 15.76; /* Ar ionization energy [eV] */
    double K_iz0 = 2.3e-14; /* rate constant prefactor [m^3/s] */

    double u_B0 = sqrt(PLASMA_E_CHARGE * 1.0 / PLASMA_M_AR);
    double tau_loss = d_eff / u_B0;
    double product = K_iz0 * n_gas * tau_loss;

    if (product <= 1.0) return E_iz; /* no solution */

    double T_e = E_iz / log(product);
    /* Iterate to self-consistency (u_B depends on T_e) */
    for (int iter = 0; iter < 5; iter++) {
        double u_B = sqrt(PLASMA_E_CHARGE * T_e / PLASMA_M_AR);
        tau_loss = d_eff / u_B;
        product = K_iz0 * n_gas * tau_loss;
        if (product <= 1.0) break;
        T_e = E_iz / log(product);
    }

    (void)pressure;
    (void)u_B_energy;
    if (T_e < 0.5) T_e = 0.5;
    return T_e;
}

/*
 * Global model plasma density from power balance.
 *
 * Given T_e from particle balance, the power balance determines n_e:
 *
 *   P_abs = e * n_e * u_B * A_eff * E_T
 *
 * where E_T is total energy lost per e-ion pair:
 *   E_T ≈ E_iz + E_exc + E_elas + 2*T_e + V_sheath
 *
 * For Ar: E_T ≈ 5*T_e + E_iz + 2*T_e ~ 30-50 eV.
 *
 * Solving: n_e = P_abs / (e * u_B * A_eff * E_T)
 */
double global_model_plasma_density(double P_abs, double T_e,
                                   double R, double L, double u_B_energy)
{
    if (T_e <= 0.0) return 0.0;

    double V = M_PI * R * R * L;
    double A_eff = 2.0 * M_PI * R * R + 2.0 * M_PI * R * L;

    double u_B = sqrt(PLASMA_E_CHARGE * T_e / PLASMA_M_AR);
    double E_T = u_B_energy + 5.0 * T_e + 15.76 + 2.0 * T_e;

    double denom = PLASMA_E_CHARGE * u_B * A_eff * E_T;
    if (denom <= 0.0) return 0.0;
    (void)V;

    return P_abs / denom;
}
