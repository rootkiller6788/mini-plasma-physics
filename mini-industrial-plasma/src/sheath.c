/**
 * sheath.c - Plasma Sheath Physics Implementation
 * Reference: Lieberman & Lichtenberg 2005, Ch.6,11
 *   Bohm (1949), Child (1911), Langmuir (1913,1924)
 * Course: MIT 22.611 / Berkeley EECS 245
 */

#include "sheath.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

/* ============================================================
 * L2: Bohm Criterion
 * Theorem (Bohm 1949): Ions must enter sheath at v_i >= c_s
 * c_s = sqrt(e*T_e/m_i) is the ion sound speed.
 * Derivation: Expand Poisson equation near sheath edge,
 * quasineutrality breaks down unless v_i >= c_s.
 * ============================================================ */

double bohm_velocity(double T_e_eV, double m_i)
{
    assert(T_e_eV > 0.0);
    assert(m_i > 0.0);
    return sqrt(PLASMA_E_CHARGE * T_e_eV / m_i);
}

int bohm_criterion_check(double v_i, double T_e_eV, double m_i)
{
    double c_s = bohm_velocity(T_e_eV, m_i);
    return (v_i >= c_s * 0.9999) ? 1 : 0;
}

/* Modified Bohm velocity for electronegative plasmas.
 * Braithwaite & Allen (1988): negative ions reduce space charge,
 * c_s_eff = c_s * sqrt((1+alpha_s)/(1+alpha_s*gamma)) */
double bohm_velocity_electronegative(double T_e_eV, double m_i,
                                     double alpha_s, double gamma_ratio)
{
    double c_s = bohm_velocity(T_e_eV, m_i);
    if (alpha_s <= 0.0) return c_s;
    return c_s * sqrt((1.0 + alpha_s) / (1.0 + alpha_s * gamma_ratio));
}

/* ============================================================
 * L3: Debye Sheath Thickness
 * Lieberman Eq. 6.3.15 (matrix model):
 * s = sqrt(2/3) * lam_D * (2*V/T_e)^{3/4}
 * ============================================================ */

double debye_sheath_thickness(double V_bias, double T_e_eV, double lam_D)
{
    double V = fabs(V_bias);
    if (V < 1e-9 || T_e_eV <= 0.0 || lam_D <= 0.0) return 0.0;
    double scaled = 2.0 * V / T_e_eV;
    if (scaled < 0.0) return 0.0;
    return sqrt(2.0/3.0) * lam_D * pow(scaled, 0.75);
}

/* Collisional sheath: s = ((9/8)*eps0*mu_i*V^2/J_i)^{1/3} */
double debye_sheath_thickness_collisional(double V_bias, double mu_i,
                                          double J_i)
{
    double V = fabs(V_bias);
    if (V < 1e-9 || mu_i <= 0.0 || J_i <= 0.0) return 0.0;
    double factor = (9.0/8.0) * PLASMA_EPS0 * mu_i * V * V / J_i;
    if (factor <= 0.0) return 0.0;
    return pow(factor, 1.0/3.0);
}

/* ============================================================
 * L4: Child-Langmuir Law
 * J = (4/9)*eps0*sqrt(2e/m_i)*V^{3/2}/d^2
 * References: Child, Phys. Rev. 32, 492 (1911)
 *             Langmuir, Phys. Rev. 2, 450 (1913)
 *             Langmuir & Blodgett, Phys. Rev. 24, 49 (1924)
 * ============================================================ */

double child_langmuir_current_density(double V, double d, double m_i)
{
    assert(V > 0.0);
    assert(d > 0.0);
    assert(m_i > 0.0);
    double prefactor = sqrt(2.0 * PLASMA_E_CHARGE / m_i);
    double J = (4.0/9.0) * PLASMA_EPS0 * prefactor * pow(V, 1.5) / (d * d);
    return (J > 0.0) ? J : 0.0;
}

double child_langmuir_sheath_width(double J, double V, double m_i)
{
    assert(J > 0.0);
    assert(V > 0.0);
    assert(m_i > 0.0);
    double prefactor = sqrt(2.0 * PLASMA_E_CHARGE / m_i);
    double d_squared = (4.0/9.0) * PLASMA_EPS0 * prefactor * pow(V, 1.5) / J;
    if (d_squared <= 0.0) return 0.0;
    return sqrt(d_squared);
}

/* Generalized Child-Langmuir with finite injection velocity v0.
 * Jaffe, Phys. Rev. 65, 91 (1944):
 * J = (4e0/9)*sqrt(2e/m_i)*[(V+V0)^{3/2} - V0^{3/2}]/d^2 */
double child_langmuir_with_initial_velocity(double V, double d,
                                            double m_i, double v0)
{
    assert(V > 0.0);
    assert(d > 0.0);
    assert(m_i > 0.0);
    double e = PLASMA_E_CHARGE;
    double V0 = 0.5 * m_i * v0 * v0 / e;
    if (V0 < 1e-12) {
        return child_langmuir_current_density(V, d, m_i);
    }
    double prefactor = sqrt(2.0 * e / m_i);
    double numerator = pow(V + V0, 1.5) - pow(V0, 1.5);
    double J = (4.0/9.0) * PLASMA_EPS0 * prefactor * numerator / (d * d);
    return (J > 0.0) ? J : 0.0;
}

/* ============================================================
 * L4: Floating Potential
 * V_f = -(T_e/2)*ln(m_i/(2*pi*m_e))
 * Langmuir probe theory, Lieberman Eq. 6.2.7
 * ============================================================ */

double floating_potential(double T_e_eV, double m_i)
{
    assert(T_e_eV > 0.0);
    assert(m_i > 0.0);
    double mass_ratio = m_i / (2.0 * M_PI * PLASMA_M_E);
    return -0.5 * T_e_eV * log(mass_ratio);
}

/* Modified Bessel function I0(x) series approximation.
 * For x<3.75: I0(x)=1 + x^2/4 + x^4/64 + x^6/2304 + ...
 * For x>=3.75: I0(x) ~ exp(x)/sqrt(2*pi*x)
 * Reference: Abramowitz & Stegun 9.8.1 */
static double bessel_I0_approx(double x)
{
    if (x < 0.0) x = -x;
    if (x > 15.0) {
        return exp(x) / sqrt(2.0 * M_PI * x);
    }
    double x2 = x * x;
    double x4 = x2 * x2;
    double x6 = x4 * x2;
    double x8 = x4 * x4;
    double x10 = x6 * x4;
    return 1.0 + x2/4.0 + x4/64.0 + x6/2304.0 + x8/147456.0
         + x10/14745600.0 + x10*x2/2123366400.0;
}

double floating_potential_rf(double T_e_eV, double m_i, double V_rf)
{
    assert(T_e_eV > 0.0);
    assert(m_i > 0.0);
    double V_float_dc = floating_potential(T_e_eV, m_i);
    double x = V_rf / T_e_eV;
    double I0_val = bessel_I0_approx(x);
    double delta_V = T_e_eV * log(I0_val);
    return V_float_dc - delta_V;
}

/* ============================================================
 * L4: Ion Impact Energy at Electrode
 * Davis & Vanderslice, Phys. Rev. 131, 219 (1963)
 * ============================================================ */

double ion_impact_energy(double V_sheath, double sheath_width,
                         double lambda_cx)
{
    double E_max = PLASMA_E_CHARGE * V_sheath;
    if (lambda_cx <= 0.0 || sheath_width <= 0.0) return E_max;
    return E_max * exp(-sheath_width / lambda_cx);
}

/* IEDF width from RF modulation
 * Kawamura et al., Plasma Sources Sci. Technol. 8, R45 (1999) */
double iedf_energy_width(double V_rf, double freq, double tau_ion)
{
    if (freq <= 0.0 || tau_ion <= 0.0 || V_rf <= 0.0) return 0.0;
    double ratio = tau_ion * freq;
    double delta_E = (2.0 / M_PI) * PLASMA_E_CHARGE * V_rf * ratio;
    return (delta_E > 0.0) ? delta_E : 0.0;
}

/* ============================================================
 * L4: Sheath Capacitance and Impedance
 * C/A = eps0 / s(V), s ~ V^{3/4} -> nonlinear capacitor
 * R_s models ion current limit + stochastic heating
 * ============================================================ */

double sheath_capacitance_per_area(double V_sheath, double T_e_eV,
                                   double lam_D)
{
    double V = fabs(V_sheath);
    if (V < 1e-9 || T_e_eV <= 0.0 || lam_D <= 0.0) return INFINITY;
    double s = debye_sheath_thickness(V, T_e_eV, lam_D);
    if (s <= 0.0) return INFINITY;
    return PLASMA_EPS0 / s;
}

void sheath_impedance(double V_sheath, double T_e_eV, double lam_D,
                      double freq, double n_e, double *R_s, double *C_s)
{
    double V = fabs(V_sheath);
    if (V < 1e-9 || T_e_eV <= 0.0 || n_e <= 0.0) {
        *R_s = INFINITY;
        *C_s = 0.0;
        return;
    }
    double s = debye_sheath_thickness(V, T_e_eV, lam_D);
    double cap_per_area = (s > 0.0) ? PLASMA_EPS0 / s : 0.0;

    /* Ion current limited resistance */
    double m_i = PLASMA_M_AR;
    double c_s = sqrt(PLASMA_E_CHARGE * T_e_eV / m_i);
    double J_sat = PLASMA_E_CHARGE * n_e * c_s;
    double res_per_area = (J_sat > 0.0) ? V / J_sat : INFINITY;

    /* Stochastic (Fermi) heating resistance */
    double v_e = sqrt(8.0 * PLASMA_E_CHARGE * T_e_eV / (M_PI * PLASMA_M_E));
    double R_stoc = PLASMA_M_E * v_e
                  / (PLASMA_E_CHARGE * PLASMA_E_CHARGE * n_e);

    (void)freq;
    *R_s = res_per_area + R_stoc;
    *C_s = cap_per_area;
}

/* ============================================================
 * L5: Planar Sheath RK4 ODE Solver
 * Solves: d2phi/dx2 = (e/eps0)*(n_e - n_i)
 *   n_e = n0*exp(e*phi/kT_e)
 *   n_i = n0/sqrt(1 - 2e*phi/(m_i*c_s^2))
 * BC: phi(0)=0, dphi/dx(0)=0 at sheath edge
 * Reference: Franklin, J. Phys. D 9, 1709 (1976)
 * ============================================================ */

typedef struct {
    double phi;
    double E;
} SheathODE;

static void sheath_ode_deriv(double x, const SheathODE *y, SheathODE *dydx,
                              double n0, double T_e_eV, double m_i)
{
    double e = PLASMA_E_CHARGE;
    double eps0 = PLASMA_EPS0;
    double kT = e * T_e_eV;
    double cs2 = kT / m_i;

    /* Boltzmann electrons */
    double n_e;
    if (y->phi / T_e_eV < -50.0) {
        n_e = 0.0;
    } else {
        n_e = n0 * exp(y->phi / T_e_eV);
    }

    /* Cold ions with flux conservation */
    double arg = 1.0 - 2.0 * e * y->phi / (m_i * cs2);
    double n_i;
    if (arg > 1e-6) {
        n_i = n0 / sqrt(arg);
    } else {
        n_i = 1e6 * n0;
    }

    dydx->phi = -y->E;
    dydx->E = (e / eps0) * (n_i - n_e);
    (void)x;
}

static void rk4_step_integrate(SheathODE *y, double *x, double *h,
                               double n0, double T_e_eV, double m_i)
{
    SheathODE k1, k2, k3, k4, yt;
    double hh = *h;
    double xx = *x;

    sheath_ode_deriv(xx, y, &k1, n0, T_e_eV, m_i);

    yt.phi = y->phi + 0.5*hh*k1.phi;
    yt.E   = y->E   + 0.5*hh*k1.E;
    sheath_ode_deriv(xx + 0.5*hh, &yt, &k2, n0, T_e_eV, m_i);

    yt.phi = y->phi + 0.5*hh*k2.phi;
    yt.E   = y->E   + 0.5*hh*k2.E;
    sheath_ode_deriv(xx + 0.5*hh, &yt, &k3, n0, T_e_eV, m_i);

    yt.phi = y->phi + hh*k3.phi;
    yt.E   = y->E   + hh*k3.E;
    sheath_ode_deriv(xx + hh, &yt, &k4, n0, T_e_eV, m_i);

    y->phi += (hh/6.0) * (k1.phi + 2.0*k2.phi + 2.0*k3.phi + k4.phi);
    y->E   += (hh/6.0) * (k1.E   + 2.0*k2.E   + 2.0*k3.E   + k4.E);
    *x += hh;
}

int solve_planar_sheath(SheathSolution *result, double n0,
                        double T_e_eV, double m_i,
                        double V_wall, int n_points)
{
    if (!result || n0 <= 0.0 || T_e_eV <= 0.0 || m_i <= 0.0 || n_points < 3) {
        return -1;
    }

    double V_target = -fabs(V_wall);
    int max_steps = n_points * 120;

    double *x_traj = (double*)malloc(max_steps * sizeof(double));
    double *phi_traj = (double*)malloc(max_steps * sizeof(double));
    double *ni_traj = (double*)malloc(max_steps * sizeof(double));
    double *ne_traj = (double*)malloc(max_steps * sizeof(double));
    double *E_traj = (double*)malloc(max_steps * sizeof(double));

    if (!x_traj || !phi_traj || !ni_traj || !ne_traj || !E_traj) {
        free(x_traj); free(phi_traj); free(ni_traj);
        free(ne_traj); free(E_traj);
        return -2;
    }

    SheathODE yobj;
    yobj.phi = 0.0;
    yobj.E   = 0.0;

    double lam_D = debye_length_electron(n0, T_e_eV);
    double h = lam_D * 0.005;
    double x = 0.0;

    int n_traj = 0;
    x_traj[n_traj] = x;
    phi_traj[n_traj] = yobj.phi;
    ni_traj[n_traj] = n0;
    ne_traj[n_traj] = n0;
    E_traj[n_traj] = yobj.E;
    n_traj++;

    int converged = 0;
    for (int i = 0; i < max_steps - 1; i++) {
        rk4_step_integrate(&yobj, &x, &h, n0, T_e_eV, m_i);

        x_traj[n_traj] = x;
        phi_traj[n_traj] = yobj.phi;
        E_traj[n_traj] = yobj.E;

        double ne_val;
        if (yobj.phi / T_e_eV < -50.0) {
            ne_val = 0.0;
        } else {
            ne_val = n0 * exp(yobj.phi / T_e_eV);
        }
        ne_traj[n_traj] = ne_val;

        double cs2 = PLASMA_E_CHARGE * T_e_eV / m_i;
        double arg = 1.0 - 2.0 * PLASMA_E_CHARGE * yobj.phi / (m_i * cs2);
        ni_traj[n_traj] = (arg > 1e-9) ? n0 / sqrt(arg) : 1e10 * n0;
        n_traj++;

        if (yobj.phi <= V_target) {
            converged = 1;
            break;
        }

        if (x > 15.0 * lam_D) h = lam_D * 0.2;
        if (yobj.phi < 0.95 * V_target) h = lam_D * 0.002;
        if (yobj.E > 1e7) h = lam_D * 0.001;
    }

    if (converged && n_traj > n_points) {
        int skip = n_traj / n_points;
        if (skip < 1) skip = 1;

        result->x = (double*)malloc(n_points * sizeof(double));
        result->potential = (double*)malloc(n_points * sizeof(double));
        result->ion_density = (double*)malloc(n_points * sizeof(double));
        result->electron_density = (double*)malloc(n_points * sizeof(double));
        result->electric_field = (double*)malloc(n_points * sizeof(double));

        if (!result->x || !result->potential || !result->ion_density
            || !result->electron_density || !result->electric_field) {
            free(result->x); free(result->potential);
            free(result->ion_density); free(result->electron_density);
            free(result->electric_field);
            converged = -3;
        } else {
            int out = 0;
            for (int i = 0; i < n_traj && out < n_points; i += skip) {
                result->x[out] = x_traj[i];
                result->potential[out] = phi_traj[i];
                result->ion_density[out] = ni_traj[i];
                result->electron_density[out] = ne_traj[i];
                result->electric_field[out] = E_traj[i];
                out++;
            }
            if (out < n_points) {
                int last = n_points - 1;
                result->x[last] = x_traj[n_traj-1];
                result->potential[last] = phi_traj[n_traj-1];
                result->ion_density[last] = ni_traj[n_traj-1];
                result->electron_density[last] = ne_traj[n_traj-1];
                result->electric_field[last] = E_traj[n_traj-1];
                out = n_points;
            }
            result->n_points = out;
        }
    } else {
        result->n_points = 0;
    }

    if (converged == 1) {
        result->sheath_width = x;
        result->V_wall = yobj.phi;
        double c_s = bohm_velocity(T_e_eV, m_i);
        result->ion_flux = ni_traj[n_traj-1] * c_s;
        result->ion_energy = PLASMA_E_CHARGE * fabs(yobj.phi);
    } else if (converged >= 0) {
        result->sheath_width = 0.0;
        result->V_wall = 0.0;
        result->ion_flux = 0.0;
        result->ion_energy = 0.0;
        result->n_points = 0;
    }

    free(x_traj); free(phi_traj); free(ni_traj);
    free(ne_traj); free(E_traj);
    return (converged == 1) ? 0 : -4;
}

/* ============================================================
 * L5: Matrix (Uniform Ion Density) Sheath Model
 * phi(x) = V_wall * (x/s)^2, s = sqrt(2*eps0*|V_wall|/(e*n_i))
 * Simplest analytic model widely used in feature-scale simulation.
 * ============================================================ */

double matrix_sheath_width(double V_wall, double n_i)
{
    double V = fabs(V_wall);
    if (V < 1e-9 || n_i <= 0.0) return 0.0;
    double term = 2.0 * PLASMA_EPS0 * V / (PLASMA_E_CHARGE * n_i);
    if (term <= 0.0) return 0.0;
    return sqrt(term);
}
