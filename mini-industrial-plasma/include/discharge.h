#ifndef DISCHARGE_H
#define DISCHARGE_H

#include "plasma_types.h"
#include <stddef.h>

/*
 * discharge.h — Plasma Discharge Models for Industrial Processing
 *
 * Reference: Lieberman §9-13, Chabert & Braithwaite "Physics of
 *   Radio-Frequency Plasmas" (2011)
 * Course: MIT 22.611, Stanford EE 414, Berkeley EECS 245
 */

/*
 * L1: DC glow discharge — Townsend breakdown and normal glow regime
 */

/* Townsend first ionization coefficient: alpha = A*p*exp(-B*p/E) */
double townsend_alpha(double pressure, double E_field, double A, double B);

/* Townsend second coefficient: effective secondary emission */
double townsend_gamma_effective(double gamma_se, double alpha, double d);

/* Breakdown condition: gamma*(exp(alpha*d) - 1) = 1 */
int townsend_breakdown_condition(double gamma_se, double alpha, double d,
                                 double *breakdown_threshold);

/*
 * L2: Paschen curve — breakdown voltage vs p*d
 * V_br = B*p*d / ln(A*p*d / ln(1 + 1/gamma_se))
 */
double paschen_breakdown_voltage(double pd, double A, double B,
                                 double gamma_se);

double paschen_minimum_voltage(double A, double B, double gamma_se);

double paschen_pd_at_minimum(double A, double B, double gamma_se);

void compute_paschen_curve(PaschenCurve *curve, double A, double B,
                           double gamma_se, double pd_min, double pd_max,
                           int n_points);

/*
 * L2: DC discharge regimes — voltage-current characteristic
 */

typedef enum {
    DC_DARK_DISCHARGE = 0,
    DC_TOWNSEND       = 1,
    DC_CORONA         = 2,
    DC_NORMAL_GLOW    = 3,
    DC_ABNORMAL_GLOW  = 4,
    DC_ARC            = 5
} DCDischargeRegime;

/* Normal glow current density (experimental scaling):
 * j_n ~ p^2 * (Constant depending on gas/electrode) */
double normal_glow_current_density(double pressure, double gas_constant);

double abnormal_glow_voltage(double current_density, double pressure,
                             double d, double gas_constant, double exponent);

/*
 * L3: RF sheath models — Capacitively Coupled Plasma (CCP)
 */

/* CCP homogeneous discharge model (Lieberman §11.2):
 * Relates RF voltage, current, power to plasma density */
typedef struct {
    double V_rf;        /* applied RF voltage amplitude [V] */
    double I_rf;        /* RF current amplitude [A] */
    double V_dc_bias;   /* self-generated DC bias [V] */
    double P_abs;       /* absorbed power [W] */
    double n_e_bulk;    /* bulk plasma density [m^-3] */
    double s_m;         /* time-averaged sheath width [m] */
    double s_max;       /* maximum sheath width [m] */
    double T_e;         /* electron temperature [eV] */
} CCPState;

/* Compute DC self-bias in asymmetric CCP:
 * V_dc ≈ -V_rf * (A_g/A_p)^k  with k ≈ 1 (experimental) */
double ccp_self_bias(double V_rf, double area_powered, double area_ground);

/* Sheath heating power per unit area (Lieberman Eq. 11.2.31):
 * P_sh = 0.5 * n_s * e * c_s * V_sh * A */
double ccp_sheath_heating_power(double n_sheath, double T_e_eV,
                                double V_sheath, double m_i, double area);

/* Ohmic heating power in bulk plasma:
 * P_ohm = J^2 * d / sigma, sigma = e^2*n_e/(m_e*nu_en) */
double ccp_ohmic_heating_power(double current_density, double gap,
                               double n_e, double nu_en);

/* Global CCP power balance: solve for n_e given absorbed power */
double ccp_equilibrium_density(double P_abs, double V_rf,
                               double T_e_eV, double m_i,
                               double area_powered, double area_ground,
                               double pressure, double gap);

/*
 * L3: Inductively Coupled Plasma (ICP) — transformer model
 */

typedef struct {
    double P_abs;        /* absorbed power [W] */
    double freq;         /* driving frequency [Hz] */
    int    N_coils;      /* number of coil turns */
    double R_coil;       /* coil radius [m] */
    double R_chamber;    /* chamber radius [m] */
    double L_chamber;    /* chamber height [m] */
    double n_e;          /* plasma density [m^-3] */
    double T_e;          /* electron temperature [eV] */
    double skin_depth;   /* collisionless skin depth [m] */
    double R_plasma;     /* plasma resistance [Ohm] */
    double L_plasma;     /* plasma inductance [H] */
    double coupling_eff; /* power transfer efficiency [0-1] */
} ICPState;

/* Collisionless skin depth: delta = c/omega_pe */
double icp_skin_depth(double n_e);

/* Plasma resistance (collisionless, anomalous):
 * R_p ~ (2*pi*R_chamber*N^2) / (sigma_dc * delta * L_chamber) */
double icp_plasma_resistance(double n_e, double radius, double height,
                             int N_coils, double freq);

/* H-mode transition threshold power:
 * P_th ≈ scaling * pressure^exponent * freq^(-1) */
double icp_eh_mode_threshold(double pressure, double freq, double radius,
                             double height);

/*
 * L4: Electron heating mechanisms
 */

/* Stochastic (Fermi) heating in oscillating sheath:
 * P_stoc = 0.5*m_e*n_s*v_e^3 * (v_sh/v_e)^2  per unit area */
double stochastic_heating_power_density(double n_sheath, double T_e_eV,
                                        double V_sheath, double freq);

/* Ion transit time heating in RF sheaths */
double ion_transit_heating_power(double n_i, double m_i, double V_rf,
                                 double freq, double sheath_width);

/*
 * L5: Global (0-D) discharge model
 */
double global_model_electron_temp(double pressure, double R, double L,
                                  double n_gas, double u_B_energy);

double global_model_plasma_density(double P_abs, double T_e,
                                   double R, double L, double u_B_energy);

#endif /* DISCHARGE_H */
