/**
 * fusion_mhd.h — MHD Equilibrium and Stability Header
 *
 * Refs: Freidberg "Ideal MHD" (2014), Wesson (2011),
 *       Goedbloed & Poedts "Principles of MHD" (2004)
 *
 * L2: Flux surface geometry, Shafranov shift
 * L3: Grad-Shafranov equation, flux coordinates
 * L4: Ideal MHD force balance, virial theorem
 * L5: Equilibrium reconstruction (finite difference)
 * L8: Vertical stability, RWM, error fields
 */

#ifndef MINI_FUSION_MHD_H
#define MINI_FUSION_MHD_H

#include "fusion_plasma.h"

/* L3: Flux Surface Geometry */
double shafranov_shift(double r, double R0, double beta_p, double li);
double connection_length(double R, double q);
double rotational_transform(double q);
double toroidal_flux(double B_phi, double r);
double poloidal_flux_from_current(double r, double R, double B_phi, double q);

/* L4: MHD Equilibrium */
double force_balance_residual(Vec3 j, Vec3 B, Vec3 grad_p);
double virial_integrand(double p, double B);
double internal_inductance(double q_95);
double normalized_internal_inductance(double li);
double poloidal_beta_from_params(double p_avg, double B_p_a);
double edge_poloidal_field(double Ip, double a, double kappa);

/* L5: Equilibrium Reconstruction */
double cylindrical_equilibrium_beta(double p0, double B_phi0, double B_phi_a, double a);
double delta_star_psi(const double *psi_grid, const double *R_grid,
                       const double *Z_grid, int i, int j,
                       int nx, int ny, double dR, double dZ);
void magnetic_axis_find(const double *psi_grid, int nx, int ny,
                         double *R_axis, double *Z_axis, double *psi_axis,
                         const double *R_grid, const double *Z_grid);
double boundary_flux_find(const double *psi_grid, int nx, int ny,
                           double psi_axis, const double *R_grid,
                           const double *Z_grid, int axis_i, int axis_j);

/* L8: Advanced Stability */
double vertical_stability_growth_rate(double n_index, double tau_wall);
double resistive_wall_mode_growth_rate(double beta, double beta_no_wall,
                                        double beta_ideal_wall, double tau_wall);
double error_field_amplification(double n_index);
double plasma_inductance(double R, double a, double li);
double flux_consumption(double L_p, double Ip, double psi_bd);

#endif /* MINI_FUSION_MHD_H */