#ifndef SHEATH_H
#define SHEATH_H

#include "plasma_types.h"
#include <stddef.h>
#include <math.h>

/*
 * L2: Bohm Criterion — Ion entry velocity into plasma sheath
 * Reference: Bohm (1949), Lieberman §6.2
 *
 * c_s = sqrt(e * T_e / m_i) is the ion sound speed (Bohm velocity).
 * Ions must enter sheath with v_i >= c_s for stable sheath formation.
 */

double bohm_velocity(double T_e_eV, double m_i);

int bohm_criterion_check(double v_i, double T_e_eV, double m_i);

/*
 * Modified Bohm criterion for electronegative plasmas:
 * c_s_eff = c_s * sqrt((1+alpha_s)/(1+alpha_s*gamma))
 * alpha_s = n_-/n_e at sheath edge, gamma = T_e/T_-
 * Reference: Braithwaite & Allen, J. Phys. D 21, 1733 (1988)
 */
double bohm_velocity_electronegative(double T_e_eV, double m_i,
                                     double alpha_s, double gamma_ratio);

/*
 * L3: Debye Sheath — Thickness estimation
 * Matrix model: s = sqrt(2/3) * lam_D * (2*V_bias/T_e)^{3/4}
 * Lieberman Eq. 6.3.15
 */
double debye_sheath_thickness(double V_bias, double T_e_eV, double lam_D);

/*
 * Collisional sheath (mobility-limited regime):
 * s = ((9/8) * eps0 * mu_i * V_bias^2 / J_i)^{1/3}
 * Valid when mean free path << sheath width.
 */
double debye_sheath_thickness_collisional(double V_bias, double mu_i,
                                          double J_i);

/*
 * L4: Child-Langmuir Law — Space-charge-limited ion current
 * J = (4/9) * eps0 * sqrt(2*e/m_i) * V^{3/2} / d^2
 * References: Child (1911), Langmuir (1913, 1924)
 */
double child_langmuir_current_density(double V, double d, double m_i);

double child_langmuir_sheath_width(double J, double V, double m_i);

/*
 * Generalized Child-Langmuir with finite ion injection velocity v0.
 * J = (4*eps0/9) * sqrt(2*e/m_i) * (V^{3/2} + 3*V*v0^2*m_i/(4*e)) / d^2
 */
double child_langmuir_with_initial_velocity(double V, double d,
                                            double m_i, double v0);

/*
 * L4: Floating potential — electron-ion flux balance at insulating wall
 * V_f = -(T_e/2) * ln(m_i/(2*pi*m_e))
 * Reference: Lieberman Eq. 6.2.7
 */
double floating_potential(double T_e_eV, double m_i);

/*
 * Floating potential with RF modulation (I0 = modified Bessel function).
 * Uses series approximation for I0(x).
 */
double floating_potential_rf(double T_e_eV, double m_i, double V_rf);

/*
 * L4: Ion impact energy at electrode
 * E_ion = e*V_sheath (collisionless)
 * With collisions: E *= exp(-s/lambda_cx)
 */
double ion_impact_energy(double V_sheath, double sheath_width,
                         double lambda_cx);

/*
 * IEDF width from RF modulation of sheath voltage
 */
double iedf_energy_width(double V_rf, double freq, double tau_ion);

/*
 * L4: Sheath capacitance — nonlinear capacitance of plasma sheath
 * C/A = eps0 / s(V) with s ~ V^{3/4}
 */
double sheath_capacitance_per_area(double V_sheath, double T_e_eV,
                                   double lam_D);

void sheath_impedance(double V_sheath, double T_e_eV, double lam_D,
                      double freq, double n_e, double *R_s, double *C_s);

/*
 * L5: Numerical sheath solver — RK4 integration of Poisson equation
 */
int solve_planar_sheath(SheathSolution *result, double n0,
                        double T_e_eV, double m_i,
                        double V_wall, int n_points);

/*
 * L5: Matrix (uniform ion density) sheath model
 * phi(x) = V_wall * (x/s)^2, s = sqrt(2*eps0*|V_wall|/(e*n_i))
 */
double matrix_sheath_width(double V_wall, double n_i);

#endif /* SHEATH_H */
