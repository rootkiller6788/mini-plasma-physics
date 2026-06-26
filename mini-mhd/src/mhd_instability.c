/* mhd_instability.c -- MHD Instability Analysis */
/* Reference: Freidberg (2014) Ch.9-12, Bateman (1978) */
#include "mhd_instability.h"
#include "mhd_defs.h"
#include <math.h>
#include <stdlib.h>

double mhd_energy_principle_fluid(const MHDState *eq, double xix, double xiy, double xiz, double g) {
    if(!eq)return 0.0;
    double B2=eq->Bx*eq->Bx+eq->By*eq->By+eq->Bz*eq->Bz;
    double Bmag=sqrt(B2);
    double inv=1.0/MHD_MU0;
    double d2=xix*xix+xiy*xiy+xiz*xiz;
    double Qx=Bmag*xix, Qy=Bmag*xiy, Qz=Bmag*xiz;
    double Q2=Qx*Qx+Qy*Qy+Qz*Qz;
    return 0.5*(Q2*inv+g*eq->p*d2);
}

double mhd_suydam_criterion(double r, double Bz, double q, double dq, double dp) {
    if(fabs(q)<1e-40)return -1e308;
    return r*Bz*Bz/(8.0*MHD_MU0)*(dq*dq)/(q*q)+dp;
}

int mhd_kruskal_shafranov_condition(double a,double R0,double Bp,double Bt){
    if(R0<1e-40||fabs(Bt)<1e-40)return 0;
    double qa=a*Bp/(R0*Bt);
    return (qa>1.0)?1:((fabs(qa-1.0)<1e-10)?0:-1);
}

double mhd_kruskal_shafranov_margin(double qa){return qa-1.0;}

double mhd_sausage_growth_rate(double r,double Bp,double dp,double rho){
    if(rho<1e-40||r<1e-40||dp>=0.0)return 0.0;
    double g2=-(2.0/MHD_MU0)*(Bp*Bp/r)*dp/(rho*fabs(Bp));
    return (g2<0.0)?0.0:sqrt(g2);
}

double mhd_kink_growth_rate(double va,double R0,double q){
    if(R0<1e-40)return 0.0;
    return (va/R0)*fabs(1.0-1.0/q);
}

double mhd_interchange_criterion(double kx,double ky,double kz,double gx,double gy,double gz){
    return kx*gx+ky*gy+kz*gz;
}

double mhd_ballooning_parameter(double q,double R,double dp,double B){
    if(fabs(B)<1e-40)return INFINITY;
    return -q*q*R*dp*(2.0*MHD_MU0/(B*B));
}

double mhd_mercier_criterion(double q,double dq,double dp,double R,double r,double B0){
    if(fabs(B0)<1e-40||R<1e-40)return INFINITY;
    double DS=r*r/(q*q)*dq*dq/(R*R);
    double DI=2.0*MHD_MU0*r*dp/(B0*B0);
    return DS+DI;
}

double mhd_tearing_growth_rate(double a,double va,double eta,double Dp){
    if(a<1e-40)return 0.0;
    double g=pow(eta/(a*a),0.6)*pow(va/a,0.4);
    if(Dp>0)g*=pow(Dp,0.8);
    return g;
}

double mhd_rayleigh_taylor_growth(double g,double k,double At,double kB,double rho){
    if(rho<1e-40)return 0.0;
    double g2=g*k*At-(kB*kB)/(MHD_MU0*rho);
    return (g2<=0.0)?0.0:sqrt(g2);
}

double mhd_kelvin_helmholtz_growth(double k,double dv,double kB,double rho){
    if(rho<1e-40)return 0.0;
    double g2=(k*dv/2.0)*(k*dv/2.0)-(kB*kB)/(MHD_MU0*rho);
    return (g2<=0.0)?0.0:sqrt(g2);
}

int mhd_linear_stability_scan(const MHDState *eq,double q,double R0,double a,double eta,double gamma){
    (void)R0; (void)a; (void)gamma;
    if(!eq)return 0;
    int flags=0;
    double Bm=sqrt(eq->Bx*eq->Bx+eq->By*eq->By+eq->Bz*eq->Bz);
    if(q<1.0)flags|=1;
    if(q<0.5)flags|=2;
    if(q>3.0)flags|=4;
    double beta=mhd_plasma_beta(eq->p,Bm);
    if(beta>0.1)flags|=8;
    if(eta>1e-40)flags|=16;
    return flags;
}
