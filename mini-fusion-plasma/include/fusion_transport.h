/**
 * fusion_transport.h — Plasma Transport Theory Header
 *
 * Refs: Braginskii "Transport Processes in a Plasma" (1965),
 *       Helander & Sigmar "Collisional Transport in Magnetized Plasmas" (2002),
 *       ITER Physics Basis (1999)
 *
 * L2: Collision frequencies, energy equilibration
 * L3: Coulomb logarithm, Spitzer resistivity, 
 *     neoclassical regimes (banana/plateau/Pfirsch-Schluter)
 * L5: Transport coefficient models (GLF23, MMM, TGLF)
 * L8: Turbulent transport, zonal flows, geodesic acoustic modes
 */

#ifndef MINI_FUSION_TRANSPORT_H
#define MINI_FUSION_TRANSPORT_H

#include "fusion_plasma.h"

/* ================================================================
 * L2: Collisionality and Regime Classification
 * ================================================================ */

/** CollisionalityRegime — neoclassical regime classification */
typedef enum {
    COLLISIONALITY_BANANA = 0,
    COLLISIONALITY_PLATEAU = 1,
    COLLISIONALITY_PFIRSCH_SCHLUTER = 2,
    COLLISIONALITY_TRANSITION = 3
} CollisionalityRegime;

/**
 * classify_collisionality — Classify neoclassical regime
 *
 * nu* < 0.1: Banana regime
 * 0.1 <= nu* <= 10: Plateau regime
 * nu* > 10: Pfirsch-Schluter regime
 *
 * @param nu_star  normalized collisionality
 * @return         regime classification
 */
CollisionalityRegime classify_collisionality(double nu_star);

/**
 * neoclassical_diffusivity_unified — Unified neoclassical diffusivity
 *
 * Smooth transition between banana, plateau, and Pfirsch-Schluter
 * regimes using interpolation formula.
 *
 * @param epsilon  inverse aspect ratio a/R
 * @param q        safety factor
 * @param rho_i    ion gyroradius [m]
 * @param nu_ii    ion-ion collision frequency [1/s]
 * @param v_thi    ion thermal velocity [m/s]
 * @param R        major radius [m]
 * @return         neoclassical diffusivity [m^2/s]
 */
double neoclassical_diffusivity_unified(double epsilon, double q, double rho_i,
                                         double nu_ii, double v_thi, double R);

/**
 * electron_neoclassical_diffusivity — Electron neoclassical diffusivity
 *
 * Similar to ion but with electron parameters and an additional
 * factor from the electron temperature gradient.
 */
double electron_neoclassical_diffusivity(double epsilon, double q, double rho_e,
                                          double nu_ei, double v_the, double R);

/* ================================================================
 * L5: Turbulent Transport Models
 * ================================================================ */

/**
 * itg_turbulent_diffusivity — ITG turbulent diffusivity [m^2/s]
 *
 * Ion Temperature Gradient (ITG) driven turbulence.
 * Mixing length estimate: chi_ITG ~ (c_s * rho_i) * (R/L_Ti)^(3/2)
 * where L_Ti = |Ti / grad(Ti)| is the temperature gradient scale length.
 *
 * @param c_s         ion sound speed [m/s]
 * @param rho_i       ion gyroradius [m]
 * @param R_over_LTi  normalized temperature gradient R/L_Ti
 * @param q           safety factor
 * @param shear       magnetic shear s
 * @return            ITG heat diffusivity [m^2/s]
 */
double itg_turbulent_diffusivity(double c_s, double rho_i, double R_over_LTi,
                                  double q, double shear);

/**
 * etg_turbulent_diffusivity — ETG turbulent diffusivity [m^2/s]
 *
 * Electron Temperature Gradient (ETG) driven turbulence.
 * chi_ETG ~ (v_the * rho_e) * (R/L_Te)^(3/2) * f(q, shear)
 * ETG produces fine-scale (rho_e scale) turbulence.
 *
 * @param v_the       electron thermal velocity [m/s]
 * @param rho_e       electron gyroradius [m]
 * @param R_over_LTe  normalized electron temperature gradient
 * @param q           safety factor
 * @return            ETG heat diffusivity [m^2/s]
 */
double etg_turbulent_diffusivity(double v_the, double rho_e, double R_over_LTe,
                                  double q);

/**
 * tem_turbulent_diffusivity — TEM turbulent diffusivity [m^2/s]
 *
 * Trapped Electron Mode (TEM) turbulence.
 * chi_TEM ~ (c_s * rho_i^2 / L_n) * f(epsilon, collisionality)
 * TEM is driven by density gradient and trapped electron dynamics.
 */
double tem_turbulent_diffusivity(double c_s, double rho_i, double L_n,
                                  double epsilon, double nu_ei, double v_the, double R);

/**
 * critical_gradient_threshold — Critical gradient for ITG/TEM onset
 *
 * Definition: (R/L_T)_crit ~ 4/(1 + tau) * f(q, shear, epsilon)
 *          tau = Te/Ti
 * Physics: instability threshold. Below critical gradient,
 *          turbulence is suppressed and internal transport
 *          barriers (ITBs) can form.
 *
 * @param tau_ratio  Te/Ti
 * @param q          safety factor
 * @param shear      magnetic shear s
 * @param epsilon    inverse aspect ratio
 * @return           critical R/L_T
 */
double critical_gradient_threshold(double tau_ratio, double q, double shear, double epsilon);

/* ================================================================
 * L8: Zonal Flows and GAMs
 * ================================================================ */

/**
 * geodesic_acoustic_mode_frequency — GAM frequency [rad/s]
 *
 * Definition: omega_GAM ~ sqrt(2) * c_s / R
 * Physics: finite-frequency zonal flow oscillation unique to
 *          toroidal geometry. GAMs regulate turbulence.
 *          omega_GAM >> turbulence correlation rate suppresses
 *          transport.
 */
double geodesic_acoustic_mode_frequency(double c_s, double R);

/**
 * zonal_flow_shearing_rate — Zonal flow shearing rate [1/s]
 *
 * Definition: omega_E = (r/q) * d(Er/B)/dr
 * Physics: sheared ExB flows tear apart turbulent eddies
 *          when omega_E > gamma_linear (turbulence decorrelation).
 *          Key mechanism for L-H transition and ITB formation.
 */
double zonal_flow_shearing_rate(double Er, double B, double r, double q,
                                 double dq_dr);

/* ================================================================
 * L5: Transport Matrix
 * ================================================================ */

/** TransportMatrix — coupled particle and heat fluxes
 *
 * (Gamma/Q) = [[D, DT], [DT, chi]] * (-grad(n), -grad(T))
 * Off-diagonal terms couple particle and heat transport.
 */
typedef struct {
    double D;       /* particle diffusivity [m^2/s] */
    double chi_e;   /* electron heat diffusivity [m^2/s] */
    double chi_i;   /* ion heat diffusivity [m^2/s] */
    double v_pinch; /* particle pinch velocity [m/s] */
    double DT;      /* off-diagonal thermo-diffusion [m^2/s] */
} TransportMatrix;

/**
 * transport_matrix_empirical — Empirical transport matrix
 *
 * Fills a TransportMatrix with typical H-mode transport coefficients
 * based on experimental database scaling.
 *
 * @param tm    output: transport matrix
 * @param rho   normalized minor radius
 * @param p     plasma parameters
 */
void transport_matrix_empirical(TransportMatrix *tm, double rho,
                                 const PlasmaParameters *p);

#endif /* MINI_FUSION_TRANSPORT_H */