/**
 * nonlinear_waves.h -- Nonlinear Wave Phenomena
 *
 * Three-wave/four-wave interactions, parametric instabilities,
 * wave steepening, solitons, modulational instability, Zakharov equations.
 *
 * References: Sagdeev & Galeev (1969), Zakharov (1972), Davidson (1972)
 * Knowledge: L8 advanced topics
 */
#ifndef NONLINEAR_WAVES_H
#define NONLINEAR_WAVES_H
#ifdef __cplusplus
extern "C" {
#endif
int three_wave_resonance_check(double w1, double k1, double w2, double k2, double w3, double k3, int sum_freq, double tol);
double three_wave_coupling_coeff(double k1, double k2, double k3, double w1, double w2, double w3, double omega_pe);
double parametric_decay_growth(double V, double E_pump, double w1, double w2);
double srs_growth_rate(double k_L, double v_osc, double omega_pe, double omega_s);
double sbs_growth_rate(double omega_pi, double v_osc, double w0, double ws, double v_th_e);
double two_plasmon_decay_growth(double k_L, double v_osc, double w0, double omega_pe);
double modulational_instability_growth(double k_pert, double omega_pe, double k_De, double W, double n0, double T_e);
/* Three-wave dynamics */
typedef struct { double A1, A2, A3, phi, V, gamma1, gamma2, gamma3; } ThreeWaveState;
void three_wave_advance_rk4(ThreeWaveState *s, double dt, double w1, double w2, double w3);
void manley_rowe_invariants(const ThreeWaveState *s, double *I13, double *I23);
double nonlinear_frequency_shift(double A_sq, double Q);
/* Saturation mechanisms */
double quasilinear_saturation_time(double gamma, double omega_pe, double dv, double v_phase, double nb_over_n0);
double trapping_saturation_field(double gamma, double k, double m, double q);
double pump_depletion_fraction(double gamma_PDI, double gamma_damp);
double wave_steepening_distance(double c_s, double gamma_eff, double omega, double dn_over_n0);
int langmuir_collapse_threshold(double W, double n0, double T_e, double k, double lambda_De);
double nonlinear_landau_damping(double W_spec, double V_coupling_sq, double dw);
#ifdef __cplusplus
}
#endif
#endif
