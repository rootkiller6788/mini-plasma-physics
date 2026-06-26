/**
 * electrostatic_waves.c -- Electrostatic Waves in Plasma
 *
 * L1: Core plasma parameter computations (omega_p, omega_c, v_th, lambda_D, etc.)
 * L4: Langmuir wave (Bohm-Gross), ion acoustic wave dispersions
 * L5: Plasma dispersion function Z(zeta), Bessel functions for kinetic theory
 * L6: Canonical Langmuir, ion acoustic, Bernstein mode systems
 *
 * References:
 *   Bohm & Gross, Phys. Rev. 75, 1851 (1949)
 *   Landau, J. Phys. USSR 10, 25 (1946)
 *   Bernstein, Phys. Rev. 109, 10 (1958)
 *   Stix, "Waves in Plasmas" (1992)
 */

#include "waves_instabilities.h"
#include "kinetic_dispersion.h"
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <stdio.h>

/* ===============================================================
 * L1: Core Plasma Parameter Computations
 * =============================================================== */

void compute_plasma_params(const PlasmaSpecies *sp, PlasmaParams *p)
{
    if (!sp || !p) return;
    double abs_q = fabs(sp->q);
    p->omega_p = (sp->n0 > 0 && sp->m > 0)
        ? sqrt(sp->n0 * abs_q * abs_q / (EPSILON0 * sp->m)) : 0.0;
    p->omega_c = (sp->B0 > 0 && abs_q > 0)
        ? (abs_q * sp->B0 / sp->m) : 0.0;
    p->v_th = (sp->T > 0 && sp->m > 0)
        ? sqrt(2.0 * K_BOLTZMANN * sp->T / sp->m) : 0.0;
    p->lambda_D = (sp->n0 > 0 && sp->T > 0)
        ? sqrt(EPSILON0 * K_BOLTZMANN * sp->T
               / (sp->n0 * E_CHARGE * E_CHARGE)) : 0.0;
    p->r_L = (p->omega_c > 0)
        ? (sp->m * p->v_th / (abs_q * sp->B0)) : 0.0;
    p->beta = (sp->B0 > 0)
        ? (2.0 * MU0 * sp->n0 * K_BOLTZMANN * sp->T
           / (sp->B0 * sp->B0)) : 0.0;
    p->c_s = (sp->m > 0 && sp->T > 0)
        ? sqrt(K_BOLTZMANN * sp->T / sp->m) : 0.0;
    p->v_A = (sp->B0 > 0 && sp->n0 > 0 && sp->m > 0)
        ? (sp->B0 / sqrt(MU0 * sp->n0 * sp->m)) : 0.0;
}

double plasma_frequency(double n, double q, double m)
{
    if (n <= 0.0 || m <= 0.0) return 0.0;
    double abs_q = fabs(q);
    return sqrt(n * abs_q * abs_q / (EPSILON0 * m));
}

double electron_plasma_frequency(double n_e)
{ return plasma_frequency(n_e, E_CHARGE, M_ELECTRON); }

double ion_plasma_frequency(double n_i, double Z_i, double m_i)
{
    if (Z_i <= 0.0) return 0.0;
    return plasma_frequency(n_i, Z_i * E_CHARGE, m_i);
}

double cyclotron_frequency(double q, double B, double m)
{
    if (m <= 0.0 || B <= 0.0) return 0.0;
    return fabs(q) * B / m;
}

double debye_length(double n, double T)
{
    if (n <= 0.0 || T <= 0.0) return 0.0;
    return sqrt(EPSILON0 * K_BOLTZMANN * T
                / (n * E_CHARGE * E_CHARGE));
}

double thermal_velocity(double T, double m)
{
    if (T <= 0.0 || m <= 0.0) return 0.0;
    return sqrt(2.0 * K_BOLTZMANN * T / m);
}

double ion_sound_speed(double T_e, double m_i)
{
    if (T_e <= 0.0 || m_i <= 0.0) return 0.0;
    return sqrt(K_BOLTZMANN * T_e / m_i);
}

double alfven_speed(double B, double n, double m_i)
{
    if (B <= 0.0 || n <= 0.0 || m_i <= 0.0) return 0.0;
    return B / sqrt(MU0 * n * m_i);
}

double plasma_beta(double n, double T, double B)
{
    if (B <= 0.0 || n <= 0.0) return 0.0;
    return 2.0 * MU0 * n * K_BOLTZMANN * T / (B * B);
}

double larmor_radius(double m, double v_perp, double q, double B)
{
    if (q == 0.0 || B <= 0.0) return 0.0;
    return m * v_perp / (fabs(q) * B);
}

double coulomb_logarithm(double n, double lambda_D)
{
    if (n <= 0.0 || lambda_D <= 0.0) return 1.0;
    double arg = 12.0 * PLASMA_PI * n
                 * lambda_D * lambda_D * lambda_D;
    if (arg < 1.0) return 1.0;
    return log(arg);
}

double ei_collision_frequency(double n, double T_e, double ln_Lambda)
{
    if (n <= 0.0 || T_e <= 0.0) return 0.0;
    double e4 = E_CHARGE * E_CHARGE * E_CHARGE * E_CHARGE;
    double num = sqrt(2.0) * n * e4 * ln_Lambda;
    double den = 12.0 * pow(PLASMA_PI, 1.5) * EPSILON0 * EPSILON0
                 * sqrt(M_ELECTRON) * pow(T_e, 1.5);
    return num / den;
}

/* ===============================================================
 * L4: Langmuir Wave Dispersion (Bohm-Gross)
 * =============================================================== */

double langmuir_wave_dispersion(double k, double omega_pe, double v_th_e)
{
    if (omega_pe <= 0.0) return 0.0;
    return sqrt(omega_pe * omega_pe + 3.0 * v_th_e * v_th_e * k * k);
}

ComplexOmega langmuir_wave_damped(double k, double omega_pe, double v_th_e)
{
    ComplexOmega r = {0.0, 0.0};
    if (omega_pe <= 0.0 || k <= 0.0 || v_th_e <= 0.0) return r;
    double lambda_De = v_th_e / (sqrt(2.0) * omega_pe);
    double kl = k * lambda_De;
    r.omega_r = omega_pe * sqrt(1.0 + 3.0 * kl * kl);
    if (kl < 0.4 && kl > 1e-10) {
        double arg = -1.0 / (2.0 * kl * kl) - 1.5;
        if (arg > -700.0)
            r.gamma = -omega_pe * sqrt(PLASMA_PI / 8.0)
                      / (kl * kl * kl) * exp(arg);
    }
    return r;
}

double ion_acoustic_wave_dispersion(double k, double c_s, double lambda_De)
{
    if (c_s <= 0.0 || k <= 0.0) return 0.0;
    return k * c_s / sqrt(1.0 + k * k * lambda_De * lambda_De);
}

ComplexOmega ion_acoustic_wave_damped(double k, double c_s, double lambda_De,
                                       double T_e, double T_i, double m_i)
{
    ComplexOmega r = {0.0, 0.0};
    if (k <= 0.0 || c_s <= 0.0) return r;
    double kl = k * lambda_De;
    r.omega_r = k * c_s / sqrt(1.0 + kl * kl);
    if (r.omega_r <= 0.0) return r;
    double ion_d = sqrt(M_ELECTRON / m_i), elec_d = 0.0;
    if (T_i > 0.0 && T_e > 0.0) {
        double rt = T_e / T_i;
        if (rt < 50.0) {
            double arg = -rt / 2.0 - 1.5;
            if (arg > -700.0) elec_d = pow(rt, 1.5) * exp(arg);
        }
    }
    r.gamma = -fabs(r.omega_r) * sqrt(PLASMA_PI / 8.0)
              * (ion_d + elec_d);
    return r;
}

/* ===============================================================
 * L6: Bernstein Mode Dispersion
 * =============================================================== */

ComplexOmega bernstein_mode_dispersion(double k_perp, double omega_pe,
                                        double omega_ce, double rho_e,
                                        int n_max)
{
    ComplexOmega r = {0.0, 0.0};
    if (omega_pe <= 0.0 || omega_ce <= 0.0 || k_perp <= 0.0
        || n_max < 2 || rho_e <= 0.0) return r;
    if (n_max > 30) n_max = 30;
    double x = k_perp * k_perp * rho_e * rho_e;
    double Gam[31];
    for (int n = 0; n <= n_max; n++)
        Gam[n] = gamma_n_bessel(n, x);
    double wlo = 1.001 * omega_ce, whi = 1.999 * omega_ce;
    double fac = omega_pe * omega_pe / (omega_ce * omega_ce);
    for (int it = 0; it < 60; it++) {
        double wm = 0.5 * (wlo + whi), sum = 0.0;
        double o2 = wm * wm / (omega_ce * omega_ce);
        for (int n = 1; n <= n_max; n++) {
            double dn = o2 - (double)(n * n);
            if (fabs(dn) < 1e-12) dn = (dn > 0) ? 1e-12 : -1e-12;
            sum += 2.0 * n * n / dn * Gam[n];
        }
        double fm = 1.0 + fac * sum;
        if (fabs(fm) < 1e-10) { r.omega_r = wm; break; }
        double suml = 0.0, o2l = wlo * wlo / (omega_ce * omega_ce);
        for (int n = 1; n <= n_max; n++) {
            double dn = o2l - (double)(n * n);
            if (fabs(dn) < 1e-12) dn = (dn > 0) ? 1e-12 : -1e-12;
            suml += 2.0 * n * n / dn * Gam[n];
        }
        double fl = 1.0 + fac * suml;
        if (fl * fm < 0.0) whi = wm; else wlo = wm;
        if (fabs(whi - wlo) < 1e-10 * omega_ce) {
            r.omega_r = wm; break;
        }
    }
    return r;
}

/* ===============================================================
 * L5: Plasma Dispersion Function Z(zeta)
 * =============================================================== */

double complex plasma_dispersion_Z(double complex zeta)
{
    double complex i_unit = _Complex_I;
    double az = cabs(zeta);
    if (az < 0.5) {
        double complex z2 = zeta * zeta;
        double complex s = 1.0 - 2.0*z2/3.0 + 4.0*z2*z2/15.0
                           - 8.0*z2*z2*z2/105.0;
        return i_unit * sqrt(PLASMA_PI) * cexp(-z2) - 2.0 * zeta * s;
    }
    if (az > 6.0) {
        double complex z2 = zeta * zeta;
        double complex a = 1.0 + 0.5/z2 + 0.75/(z2*z2);
        double sigma = (cimag(zeta) > 0.0) ? 0.0
                     : (cimag(zeta) < 0.0) ? 2.0 : 1.0;
        return -1.0/zeta * a
               + i_unit * sigma * sqrt(PLASMA_PI) * cexp(-z2);
    }
    double x = creal(zeta), y = cimag(zeta);
    double ex = exp(-x*x);
    double im = sqrt(PLASMA_PI) * ex
                * (1.0 + y*y*(x*x - 0.5)/(1.0 + fabs(y)));
    double re = -2.0*x*(1.0 - 2.0*x*x/3.0 + 4.0*x*x*x*x/15.0);
    if (fabs(x) > 1.5)
        re = -1.0/x*(1.0 + 0.5/(x*x) + 0.75/(x*x*x*x));
    return re + i_unit * im;
}

double complex plasma_dispersion_Zp(double complex zeta)
{ return -2.0 * (1.0 + zeta * plasma_dispersion_Z(zeta)); }

double plasma_Z_real(double x)
{
    if (fabs(x) < 0.5)
        return -2.0*x*(1.0 - 2.0*x*x/3.0 + 4.0*x*x*x*x/15.0);
    if (fabs(x) > 6.0)
        return -1.0/x*(1.0 + 0.5/(x*x) + 0.75/(x*x*x*x));
    double x2 = x*x;
    return -2.0*x*(1.0 - 0.6667*x2 + 0.2667*x2*x2)/(1.0 + 0.1*x2);
}

double plasma_Z_imag(double x)
{ return sqrt(PLASMA_PI) * exp(-x*x); }

/* ===============================================================
 * L5: Kinetic Dispersion Solvers
 * =============================================================== */

ComplexOmega kinetic_langmuir_dispersion(double k, double omega_pe,
                                          double v_th_e)
{
    ComplexOmega r = {0.0, 0.0};
    if (k <= 0.0 || omega_pe <= 0.0 || v_th_e <= 0.0) return r;
    double kl = k * v_th_e / (sqrt(2.0) * omega_pe);
    if (kl < 0.3) {
        r.omega_r = omega_pe * (1.0 + 1.5 * kl * kl);
        if (kl > 1e-10) {
            double arg = -1.0/(2.0*kl*kl) - 1.5;
            if (arg > -700.0)
                r.gamma = -omega_pe*sqrt(PLASMA_PI/8.0)
                          /(kl*kl*kl)*exp(arg);
        }
    } else {
        r.omega_r = omega_pe*sqrt(1.0 + 3.0*kl*kl);
        r.gamma = -omega_pe * fmax(0.0, (kl-0.3)*0.4);
    }
    return r;
}

ComplexOmega kinetic_ion_acoustic_dispersion(double k, double T_e, double T_i,
                                              double n0, double m_i)
{
    ComplexOmega r = {0.0, 0.0};
    if (k <= 0.0 || T_e <= 0.0 || n0 <= 0.0 || m_i <= 0.0) return r;
    double lDe = sqrt(EPSILON0*K_BOLTZMANN*T_e/(n0*E_CHARGE*E_CHARGE));
    double cs = sqrt(K_BOLTZMANN*T_e/m_i), kl = k*lDe;
    r.omega_r = k*cs/sqrt(1.0 + kl*kl);
    double rt = (T_i > 0.0) ? T_e/T_i : 100.0;
    double it = sqrt(M_ELECTRON/m_i), et = 0.0;
    if (rt < 50.0) {
        double arg = -rt/2.0 - 1.5;
        if (arg > -700.0) et = pow(rt, 1.5)*exp(arg);
    }
    r.gamma = -fabs(r.omega_r)*sqrt(PLASMA_PI/8.0)*(it+et);
    return r;
}

ComplexOmega kinetic_alfven_full_dispersion(double k_par, double k_perp,
                                             double v_A, double rho_i,
                                             double rho_s)
{
    ComplexOmega r = {0.0, 0.0};
    if (v_A <= 0.0 || k_par <= 0.0) return r;
    double kp2 = k_perp * k_perp;
    double flr = 1.0 + kp2*(rho_s*rho_s + 0.75*rho_i*rho_i);
    r.omega_r = fabs(k_par)*v_A*sqrt(flr);
    if (kp2 > 0 && rho_s > 0)
        r.gamma = -fabs(r.omega_r)*0.1*kp2*rho_s*rho_s
                  /(1.0 + kp2*rho_s*rho_s);
    return r;
}

double whistler_cyclotron_damping(double omega, double k_par,
                                   double omega_pe, double omega_ce,
                                   double v_th_e)
{
    if (omega <= 0 || omega_ce <= 0 || v_th_e <= 0) return 0.0;
    if (fabs(k_par) < 1e-10) return 0.0;
    double vr = (omega - omega_ce)/k_par, vrn = vr/v_th_e;
    if (fabs(vrn) > 6.0) return 0.0;
    return -PLASMA_PI*omega_pe*omega_pe*omega_ce
            /(k_par*k_par*C_LIGHT*C_LIGHT)*exp(-vrn*vrn);
}

/* ===============================================================
 * L3: Bessel Functions for Kinetic Theory
 * =============================================================== */

double modified_bessel_I(int n, double x)
{
    if (n < 0) n = -n;
    if (x < 0.0) x = -x;
    if (x == 0.0) return (n == 0) ? 1.0 : 0.0;
    if (x > 20.0 + 2.0*n) {
        double n2 = (double)(4*n*n);
        return exp(x)/sqrt(2.0*PLASMA_PI*x)
               *(1.0 - (n2-1.0)/(8.0*x)
                     + (n2-1.0)*(n2-9.0)/(128.0*x*x));
    }
    double sum = 0.0, term = 1.0;
    for (int k = 2; k <= n; k++) term /= k;
    sum = term; double xh2 = (x/2.0)*(x/2.0);
    for (int k = 1; k < 200; k++) {
        term *= xh2/(k*(n+k)); sum += term;
        if (fabs(term) < 1e-16*fabs(sum)) break;
    }
    double pf = 1.0;
    for (int k = 0; k < n; k++) pf *= x/2.0;
    return pf*sum;
}

double gamma_n_bessel(int n, double x)
{
    if (n < 0) n = -n;
    if (x < 0.0) x = -x;
    if (x == 0.0) return (n == 0) ? 1.0 : 0.0;
    if (x < 0.01) {
        if (n == 0) return 1.0 - x + 0.5*x*x;
        if (n == 1) return 0.5*x*(1.0 - x);
        return pow(x/2.0, n)/tgamma(n+1.0)*(1.0 - x);
    }
    return modified_bessel_I(n, x)*exp(-x);
}

double bessel_J(int n, double x)
{
    if (n < 0) {
        double jn = bessel_J(-n, x);
        return (n%2==0) ? jn : -jn;
    }
    if (fabs(x) < 1e-10) return (n==0) ? 1.0 : 0.0;
    double sum = 0.0, term = 1.0;
    for (int k = 2; k <= n; k++) term /= k;
    double xh2 = (x/2.0)*(x/2.0); int sgn = 1;
    for (int k = 0; k < 200; k++) {
        sum += sgn*term; sgn = -sgn;
        term *= xh2/((k+1)*(n+k+1));
        if (fabs(term) < 1e-16*fabs(sum) + 1e-30) break;
    }
    double pf = 1.0;
    for (int k = 0; k < n; k++) pf *= x/2.0;
    return pf*sum;
}

void gamma_0_and_derivative(double x, double *G0, double *dG0)
{
    if (!G0 || !dG0) return;
    if (x < 0.0) x = -x;
    if (x < 1e-10) { *G0 = 1.0; *dG0 = -1.0; return; }
    double I0 = modified_bessel_I(0, x), I1 = modified_bessel_I(1, x);
    double ex = exp(-x);
    *G0 = I0*ex; *dG0 = (I1-I0)*ex;
}

double complex hot_plasma_dielectric_es(double omega,
                                         double k_perp, double k_par,
                                         double omega_p, double omega_c,
                                         double v_th, int n_max)
{
    double k2 = k_perp*k_perp + k_par*k_par;
    if (k2 < 1e-30 || omega_p <= 0.0 || v_th <= 0.0)
        return 1.0 - omega_p*omega_p/(omega*omega + 1e-30);
    if (omega_c < omega_p*1e-3 || n_max <= 0) {
        double z = omega/(sqrt(2.0)*sqrt(k2)*v_th);
        double complex Zp = plasma_dispersion_Zp(z);
        return 1.0 + omega_p*omega_p/(k2*v_th*v_th)*Zp;
    }
    double x = k_perp*k_perp*v_th*v_th/(omega_c*omega_c);
    double complex sum = 0.0;
    for (int n = -n_max; n <= n_max; n++) {
        double zn = (omega - n*omega_c)
                    /(sqrt(2.0)*fabs(k_par)*v_th + 1e-30);
        sum += gamma_n_bessel(n, x)*plasma_dispersion_Z(zn);
    }
    double z0 = omega/(sqrt(2.0)*fabs(k_par)*v_th + 1e-30);
    double pf = omega_p*omega_p/(k2*v_th*v_th);
    return 1.0 + pf*(1.0 + z0*sum);
}
