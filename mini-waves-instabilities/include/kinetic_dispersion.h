/**
 * kinetic_dispersion.h -- Kinetic Dispersion Theory
 *
 * Vlasov-Maxwell kinetic treatment: plasma dispersion function Z(zeta),
 * Landau damping, cyclotron damping, Bernstein modes, hot plasma dielectric.
 *
 * References:
 *   Landau, J. Phys. USSR 10, 25 (1946)
 *   Fried & Conte, "The Plasma Dispersion Function" (1961)
 *   Stix, "Waves in Plasmas" (1992), Ch. 8-11
 *   Swanson, "Plasma Waves" (2003), Ch. 4-6
 *
 * Knowledge: L3 (Z function, Bessel sums), L4 (Vlasov-Maxwell),
 *             L5 (numerical Z), L6 (Landau, Bernstein)
 */
#ifndef KINETIC_DISPERSION_H
#define KINETIC_DISPERSION_H
#include <complex.h>
#ifdef __cplusplus
extern "C" {
#endif
/* Plasma dispersion function Z(zeta) -- Fried & Conte (1961) */
double complex plasma_dispersion_Z(double complex zeta);
double complex plasma_dispersion_Zp(double complex zeta);
double plasma_Z_real(double x);
double plasma_Z_imag(double x);
/* Kinetic dispersion relations */
ComplexOmega kinetic_langmuir_dispersion(double k, double omega_pe, double v_th_e);
ComplexOmega kinetic_ion_acoustic_dispersion(double k, double T_e, double T_i, double n0, double m_i);
ComplexOmega bernstein_mode_dispersion(double k_perp, double omega_pe, double omega_ce, double rho_e, int n_max);
ComplexOmega kinetic_alfven_full_dispersion(double k_parallel, double k_perp, double v_A, double rho_i, double rho_s);
double whistler_cyclotron_damping(double omega, double k_parallel, double omega_pe, double omega_ce, double v_th_e);
double complex hot_plasma_dielectric_es(double omega, double k_perp, double k_par, double omega_p, double omega_c, double v_th, int n_max);
/* Bessel function utilities */
double modified_bessel_I(int n, double x);
double gamma_n_bessel(int n, double x);
double bessel_J(int n, double x);
void gamma_0_and_derivative(double x, double *Gamma0, double *dGamma0);
#ifdef __cplusplus
}
#endif
#endif
