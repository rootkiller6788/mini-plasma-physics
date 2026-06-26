/**
 * dispersion_solvers.h -- Numerical Dispersion Relation Solvers
 *
 * Root finding (Brent), complex Newton, continuation method,
 * shooting for eigenvalue problems, Nyquist stability analysis,
 * CMA diagram topology, Penrose criterion.
 *
 * References: Numerical Recipes (2007), Stix (1992) Appendices
 * Knowledge: L5 computational methods
 */
#ifndef DISPERSION_SOLVERS_H
#define DISPERSION_SOLVERS_H
#include "waves_instabilities.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef double (*DispersionFunc)(double omega, const double params[]);
typedef void (*ComplexDispersionFunc)(double omega_r, double gamma, const double params[], double *D_real, double *D_imag);
int find_dispersion_roots_brent(DispersionFunc D, const double params[], double omega_min, double omega_max, int n_scan, double roots[], int max_roots, double eps);
int find_complex_root_newton(ComplexDispersionFunc D, const double params[], double wr0, double g0, double *wr, double *gamma, int max_iter, double tol);
typedef void (*ParamBuilder)(double k, double params[]);
int trace_dispersion_curve(DispersionFunc D, ParamBuilder build_params, double k_min, double k_max, int n_pts, double k_out[], double w_out[], double w_guess);
int mhd_stability_eigenvalue(const double rho[], const double F_diag[], const double F_off[], int n_grid, double *w_sq, double eigenvec[], int max_iter);
typedef double (*ODE_RHS)(double r, double psi, double w_sq, const double params[]);
int shooting_eigenvalue(ODE_RHS rhs, const double params[], int n_grid, const double r_grid[], double *w_sq, double psi[]);
int nyquist_unstable_count(ComplexDispersionFunc D, const double params[], double w_min, double w_max, double R, int n_pts);
int penrose_criterion(const double f[], const double v_grid[], int n_v, double omega_p);
void cma_diagram_grid(double a_min, double a_max, int n_a, double b_min, double b_max, int n_b, int regions[], int n_prop[]);
int special_frequencies(double omega_pe, double omega_ce, int n_sp, const double wp_ions[], const double wc_ions[], double special[], char labels[][32]);
#ifdef __cplusplus
}
#endif
#endif
