/**
 * plasma_instabilities.h -- Instability Catalog
 * References: Mikhailovskii (1974), Weiland (2000), Kadomtsev (1965)
 */
#ifndef PLASMA_INSTABILITIES_H
#define PLASMA_INSTABILITIES_H
#include "waves_instabilities.h"
#ifdef __cplusplus
extern "C" {
#endif
double two_stream_growth_rate(double n_beam, double n_plasma, double omega_pe);
int two_stream_dispersion_roots(double k, double omega_pe, double omega_b, double v_b, ComplexOmega roots[4]);
double two_stream_growth_at_k(double k, double omega_pe, double omega_b, double v_b);
double buneman_growth_rate(double omega_pe, double v_d, double v_th_e);
double ion_two_stream_growth(double omega_pi, double v_b, double v_th_i);
double bump_on_tail_growth(double omega_pe, double n_b, double n_0, double v_b, double delta_v);
double weibel_growth_rate(double omega_pe, double v_th, double T_perp, double T_par, double beta);
double weibel_relativistic_growth(double omega_pe, double n_b, double n_0, double gamma_b);
int firehose_unstable(double p_parallel, double p_perp, double B);
double firehose_growth_rate_sq(double k_parallel, double p_parallel, double p_perp, double B, double rho);
int mirror_unstable(double p_parallel, double p_perp, double B);
double mirror_growth_rate(double k_parallel, double v_th, double beta_perp, double beta_par);
double current_driven_ia_growth(double k, double c_s, double v_d, double m_e, double m_i);
double rayleigh_taylor_growth(double g_eff, double k, double rho_heavy, double rho_light);
double magnetic_rt_growth(double g_eff, double k, double rho_heavy, double rho_light, double B_perp);
double kelvin_helmholtz_growth(double k, double v0, double rho1, double rho2);
double magnetic_kh_growth(double k, double v0, double rho1, double rho2, double B_parallel);
double interchange_growth_rate(double grad_ln_n, double cs, double R_curvature);
double ballooning_growth_rate(double cs, double grad_p_over_p, double R_curvature, double q);
double tearing_mode_growth_rate(double eta, double a, double v_A, double delta_prime);
int kink_unstable(double q, int m, int n);
int sausage_unstable(double pressure_gradient, double btheta_gradient);
double drift_wave_frequency(double k_y, double T_e, double grad_ln_n, double B);
double drift_wave_growth_rate(double omega_star, double k_parallel, double v_th_e);
double itg_mode_frequency(double omega_star, double eta_i, double eta_i_crit, double k_perp_rho_s);
double itg_growth_rate(double omega_star, double eta_i, double eta_i_crit, double k_perp_rho_s);
double tem_frequency(double omega_star, double eta_e, double epsilon);
double tem_growth_rate(double omega_star, double eta_e, double epsilon);
double etg_growth_rate(double omega_star_e, double eta_e, double eta_e_crit, double k_perp_rho_e);
double loss_cone_growth_rate(double omega_ce, double n_hot, double n_cold, double anisotropy);
double cyclotron_maser_growth(double omega_pe, double omega_ce, double v_perp_over_c, double n_hot_over_n);
#ifdef __cplusplus
}
#endif
#endif
