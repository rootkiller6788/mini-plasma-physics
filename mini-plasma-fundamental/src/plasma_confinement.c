/**
 * plasma_confinement.c -- Magnetic Confinement Fusion Physics
 *
 * Tokamak, stellarator, and mirror machine physics.
 *
 * References:
 *   - Wesson, "Tokamaks" (2011)
 *   - Freidberg, "Plasma Physics and Fusion Energy" (2007)
 *
 * Knowledge Coverage:
 *   L4: MHD equilibrium, transport
 *   L7: ITER/SPARC parameters, energy confinement scaling
 *   L8: H-mode pedestal, mirror confinement
 */
#include "plasma_constants.h"
#include "plasma_params.h"
#include <math.h>

double aspect_ratio(double R0, double a) {
    if (a <= 0.0) return INFINITY;
    return R0 / a;
}

double toroidal_field(double B0, double R0, double R) {
    if (R <= 0.0) return INFINITY;
    return B0 * R0 / R;
}

double poloidal_field(double I_p, double a) {
    if (a <= 0.0) return 0.0;
    return MU_0 * I_p / (2.0 * M_PI * a);
}

double edge_safety_factor(double a, double R0, double B_phi, double I_p) {
    if (R0 <= 0.0 || I_p <= 0.0) return INFINITY;
    return (2.0 * M_PI * a * a * B_phi) / (MU_0 * R0 * I_p);
}

double plasma_inductance(double R0, double a, double li) {
    if (a <= 0.0) return INFINITY;
    return MU_0 * R0 * (log(8.0 * R0 / a) - 2.0 + 0.5 * li);
}

double volt_second_requirement(double I_p, double L_p, double R_p,
                                double t_ramp) {
    return L_p * I_p + R_p * I_p * t_ramp;
}

double greenwald_density_limit(double I_p_MA, double a) {
    if (a <= 0.0) return INFINITY;
    return I_p_MA / (M_PI * a * a) * 1.0e20;
}

double troyon_beta_limit(double beta_N_crit, double a, double B, double I_p) {
    if (I_p <= 0.0) return INFINITY;
    return beta_N_crit * I_p / (a * B);
}

int is_kink_stable(double q_edge) {
    return (q_edge > 1.0) ? 1 : 0;
}

typedef struct {
    double R0, a, B_t, I_p, n_e, T_e, T_i;
    double P_fusion, P_aux, Q, tau_E, beta_N, q_95;
} TokamakParams;

void iter_parameters(TokamakParams *p) {
    p->R0 = 6.2; p->a = 2.0; p->B_t = 5.3; p->I_p = 15.0;
    p->n_e = 1.0; p->T_e = 8.8; p->T_i = 8.0;
    p->P_fusion = 500.0; p->P_aux = 50.0; p->Q = 10.0;
    p->tau_E = 3.7; p->beta_N = 1.8; p->q_95 = 3.0;
}

void sparc_parameters(TokamakParams *p) {
    p->R0 = 1.85; p->a = 0.57; p->B_t = 12.2; p->I_p = 8.7;
    p->n_e = 3.0; p->T_e = 10.0; p->T_i = 10.0;
    p->P_fusion = 140.0; p->P_aux = 25.0; p->Q = 5.6;
    p->tau_E = 0.8; p->beta_N = 2.5; p->q_95 = 3.4;
}

double energy_confinement_time_lmode(double I_p_MA, double R0, double a,
                                     double kappa, double ne_20,
                                     double B_t, double A_mass,
                                     double P_loss_MW) {
    return 0.048 * pow(I_p_MA, 0.85) * pow(R0, 1.2)
           * pow(a, 0.3) * pow(kappa, 0.5) * pow(ne_20, 0.1)
           * pow(B_t, 0.2) * pow(A_mass, 0.5) * pow(P_loss_MW, -0.5);
}

double energy_confinement_time_hmode(double I_p_MA, double R0, double a,
                                     double kappa, double ne_20,
                                     double B_t, double A_mass,
                                     double P_loss_MW) {
    double epsilon = a / R0;
    return 0.0562 * pow(I_p_MA, 0.93) * pow(R0, 1.97)
           * pow(epsilon, 0.58) * pow(kappa, 0.78) * pow(ne_20, 0.41)
           * pow(B_t, 0.15) * pow(A_mass, 0.19) * pow(P_loss_MW, -0.69);
}

double lh_transition_power(double ne_20, double B_t, double S) {
    return 2.15 * pow(ne_20, 0.72) * pow(B_t, 0.80) * pow(S, 0.94);
}

double pedestal_pressure(double I_p_MA, double a) {
    if (a <= 0.0) return 0.0;
    double j = I_p_MA / (M_PI * a * a);
    return 0.1 * j * j * 1.0e3;
}

double pedestal_temperature(double p_ped, double n_ped) {
    if (n_ped <= 0.0) return 0.0;
    return p_ped / (2.0 * n_ped * 1.0e3 * E_CHARGE / 1.0e3);
}

/* Stellarator and mirror machine physics */

double rotational_transform(double iota_0, double iota_1,
                            double r_over_a) {
    return iota_0 + iota_1 * r_over_a * r_over_a;
}

typedef struct { double R0, a, B_t, n_e, T_e, P_ECRH; } StellaratorParams;

void w7x_parameters(StellaratorParams *p) {
    p->R0 = 5.5; p->a = 0.53; p->B_t = 3.0;
    p->n_e = 0.3; p->T_e = 3.0; p->P_ECRH = 10.0;
}

double mirror_ratio(double B_max, double B_min) {
    if (B_min <= 0.0) return INFINITY;
    return B_max / B_min;
}

double loss_cone_angle(double B_max, double B_min) {
    if (B_max <= B_min || B_min <= 0.0) return 0.0;
    return asin(sqrt(B_min / B_max));
}

double mirror_confinement_time(double tau_ii, double R_m) {
    if (R_m <= 1.0) return tau_ii;
    return tau_ii * R_m * exp(0.5 * R_m) * log10(R_m);
}
