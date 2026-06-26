/*
 * mhd_equilibrium.c -- MHD Equilibrium
 * Reference: Freidberg (2014) Ch.8, Wesson (2011) Ch.3
 * Knowledge: L6 (pinch equilibria), L8 (Grad-Shafranov, Taylor)
 */

#include "mhd_equilibrium.h"
#include "mhd_defs.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void mhd_zpinch_equilibrium(double (*J_z)(double), double r_max, int nr,
                             double p_edge,
                             double *r, double *B_phi, double *p) {
    if (!J_z || !r || !B_phi || !p || nr < 2) return;
    double dr = r_max / (double)(nr - 1);
    double integral = 0.0;
    for (int i = 0; i < nr; i++) {
        r[i] = i * dr;
        if (i == 0) { B_phi[i] = 0.0; }
        else {
            integral += 0.5 * (r[i-1]*J_z(r[i-1]) + r[i]*J_z(r[i])) * dr;
            B_phi[i] = MHD_MU0 * integral / r[i];
        }
    }
    p[nr-1] = p_edge;
    for (int i = nr-2; i >= 0; i--) {
        double Jz_m = 0.5*(J_z(r[i]) + J_z(r[i+1]));
        double Bm = 0.5*(B_phi[i] + B_phi[i+1]);
        p[i] = p[i+1] + Jz_m * Bm * dr;
    }
}

void mhd_thetapinch_equilibrium(double (*B_z)(double), double r_max, int nr,
                                 double p_edge,
                                 double *r, double *Bz, double *p) {
    if (!B_z || !r || !Bz || !p || nr < 2) return;
    double dr = r_max / (double)(nr - 1);
    for (int i = 0; i < nr; i++) { r[i] = i * dr; Bz[i] = B_z(r[i]); }
    double ptot = p_edge + Bz[nr-1]*Bz[nr-1]/(2.0*MHD_MU0);
    for (int i = 0; i < nr; i++) {
        p[i] = ptot - Bz[i]*Bz[i]/(2.0*MHD_MU0);
        if (p[i] < 0.0) p[i] = 0.0;
    }
}
void mhd_screwpinch_equilibrium(double (*q)(double), double R0,
                                 double r_max, int nr,
                                 double p0, double B0,
                                 double *r, double *B_phi, double *B_z, double *p) {
    if (!q || !r || !B_phi || !B_z || !p || nr < 2) return;
    double dr = r_max / (double)(nr - 1);
    for (int i = 0; i < nr; i++) {
        r[i] = i * dr;
        double ra = r[i]/r_max;
        B_z[i] = B0 * (1.0 - ra*ra);
        if (B_z[i] < 0.0) B_z[i] = 0.0;
        double qv = q(r[i]);
        if (fabs(qv) < 1e-10) qv = 1e-10;
        B_phi[i] = (r[i] < 1e-10) ? 0.0 : r[i]*B_z[i]/(R0*qv);
    }
    p[nr-1] = 0.0;
    double im = 1.0/MHD_MU0;
    for (int i = nr-2; i >= 0; i--) {
        double dBz = (B_z[i+1]-B_z[i])/dr;
        double Jp = -dBz * im;
        double rB1 = r[i]*B_phi[i], rB2 = r[i+1]*B_phi[i+1];
        double drB = (rB2 - rB1)/dr;
        double rm = 0.5*(r[i]+r[i+1]);
        double Jz = drB * im / rm;
        double Bzm = 0.5*(B_z[i]+B_z[i+1]);
        double Bpm = 0.5*(B_phi[i]+B_phi[i+1]);
        p[i] = p[i+1] - (Jp*Bzm - Jz*Bpm) * dr;
        if (p[i] < 0.0) p[i] = 0.0;
    }
    if (p[0] > 0.0) { double s = p0/p[0]; for (int i=0; i<nr; i++) p[i] *= s; }
}

double mhd_bennett_temperature(double I, double N) {
    if (N < 1e-40) return INFINITY;
    return MHD_MU0 * I * I / (16.0 * M_PI * MHD_KB * N);
}

double mhd_bennett_radius(double I, double N, double T) {
    if (I < 1e-40) return INFINITY;
    return sqrt(8.0 * M_PI * N * MHD_KB * T / (MHD_MU0 * I * I));
}

double mhd_grad_shafranov_operator(const double **psi,
                                    int iR, int iZ,
                                    double dR, double dZ, double R) {
    if (!psi) return 0.0;
    double dpp = (psi[iR+1][iZ] - psi[iR][iZ]) / dR;
    double dpm = (psi[iR][iZ] - psi[iR-1][iZ]) / dR;
    double Rp = R + 0.5*dR, Rm = R - 0.5*dR;
    double tR = R * (dpp/Rp - dpm/Rm) / dR;
    double tZ = (psi[iR][iZ+1] - 2.0*psi[iR][iZ] + psi[iR][iZ-1])/(dZ*dZ);
    return tR + tZ;
}

double mhd_soloviev_psi(double R, double Z, double R0, double kappa, double psi0) {
    double t1 = (R*R - R0*R0)*(R*R - R0*R0)/(4.0*R0*R0);
    double t2 = kappa*kappa*Z*Z/(R0*R0);
    return psi0 * (t1 + t2);
}

void mhd_soloviev_field(double R, double Z, double R0, double k,
                         double psi0, double F0,
                         double *BR, double *Bphi, double *BZ) {
    double dpZ = psi0*2.0*k*k*Z/(R0*R0);
    *BR = -dpZ/R;
    double dpR = psi0*(R*R-R0*R0)*R/(R0*R0);
    *BZ = dpR/R;
    *Bphi = F0/R;
}

double mhd_safety_factor_q(double r, double R0, double Bp, double Bt) {
    if (R0 < 1e-40 || fabs(Bt) < 1e-40) return INFINITY;
    return r * Bp / (R0 * Bt);
}

double mhd_safety_factor_from_flux(double psi, double F, double dV, double R0) {
    (void)psi;
    if (R0 < 1e-40 || dV < 1e-40) return INFINITY;
    return (F/R0) * dV / (2.0*M_PI);
}
/* ---- Bessel function approximations (Abramowitz & Stegun) ---- */

static double bessel_J0(double x) {
    double ax = fabs(x);
    if (ax < 1e-6) return 1.0 - 0.25*x*x;
    if (ax < 3.0) {
        double y = (x/3.0)*(x/3.0);
        return 1.0 - 2.2499997*y + 1.2656208*y*y - 0.3163866*y*y*y
             + 0.0444479*y*y*y*y - 0.0039444*y*y*y*y*y + 0.0002100*y*y*y*y*y*y;
    } else {
        double y = 3.0/ax;
        double th = ax - 0.78539816 + y*(-0.04166397 + y*(-0.00003954
                 + y*(0.00262573 + y*(-0.00054125))));
        return sqrt(0.63661977/ax) * cos(th);
    }
}

static double bessel_J1(double x) {
    double ax = fabs(x);
    double s = (x >= 0) ? 1.0 : -1.0;
    if (ax < 1e-6) return 0.5*x;
    if (ax < 3.0) {
        double y = (x/3.0)*(x/3.0);
        return s*x*(0.5 - 0.56249985*y + 0.21093573*y*y
             - 0.03954289*y*y*y + 0.00443319*y*y*y*y
             - 0.00031761*y*y*y*y*y + 0.00001109*y*y*y*y*y*y);
    } else {
        double y = 3.0/ax;
        double th = ax - 2.35619449 + y*(0.12499612 + y*(0.00005650
                 + y*(-0.00637879 + y*(0.00074348))));
        return s * sqrt(0.63661977/ax) * cos(th);
    }
}

static double taylor_helicity(double lam, double a, double B0, double L, int n) {
    double dr = a / (double)n, integ = 0.0;
    for (int i = 0; i < n; i++) {
        double r = (i + 0.5)*dr;
        double J0 = bessel_J0(lam*r), J1 = bessel_J1(lam*r);
        integ += (J0*J0 + J1*J1) * r * dr;
    }
    return 2.0*M_PI*L*B0*B0*integ/lam;
}

int mhd_taylor_relaxation(double a, double H0, double B0,
                           int nr, double *r, double *lam_out,
                           double *Bz, double *Bphi) {
    if (!r||!lam_out||!Bz||!Bphi||nr<2) return -1;
    double lmin=0.5/a, lmax=2.4048/a;
    double Hl=taylor_helicity(lmin,a,B0,1.0,200);
    double Hh=taylor_helicity(lmax,a,B0,1.0,200);
    double lam;
    if (H0<Hl) lam=lmin;
    else if (H0>Hh) lam=lmax;
    else {
        double lo=lmin, hi=lmax;
        for (int it=0; it<60; it++) {
            double m=0.5*(lo+hi);
            double Hm=taylor_helicity(m,a,B0,1.0,200);
            if (Hm<H0) lo=m; else hi=m;
            if (hi-lo<1e-8) break;
        }
        lam=0.5*(lo+hi);
    }
    *lam_out=lam;
    double dr=a/(double)(nr-1);
    for (int i=0; i<nr; i++) {
        r[i]=i*dr;
        double lr=lam*r[i];
        Bz[i]=B0*bessel_J0(lr);
        Bphi[i]=B0*bessel_J1(lr);
    }
    return 0;
}

double mhd_beltrami_check(double alpha, double Bx, double By, double Bz,
                           double Jx, double Jy, double Jz) {
    double dx=Jx-alpha*Bx, dy=Jy-alpha*By, dz=Jz-alpha*Bz;
    double dn=sqrt(dx*dx+dy*dy+dz*dz);
    double Bn=sqrt(Bx*Bx+By*By+Bz*Bz);
    if (Bn<1e-40) return 0.0;
    return dn/Bn;
}
