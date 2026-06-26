/*
 * wakefield.c -- Laser wakefield acceleration (LWFA)
 *
 * Implements the physics of electron plasma wave excitation
 * by ultra-short laser pulses, including operating regime
 * calculations for state-of-the-art LWFA.
 *
 * References:
 *   - Tajima & Dawson (1979), PRL 43, 267
 *   - Esarey, Schroeder & Leemans (2009), Rev. Mod. Phys. 81, 1229
 *   - Lu et al. (2006), Phys. Rev. ST Accel. Beams 10, 061301
 *   - Pukhov & Meyer-ter-Vehn (2002), Appl. Phys. B 74, 355
 *
 * Knowledge Layers: L1, L2, L4, L5
 */

#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "plasma_constants.h"
#include "plasma_params.h"
#include "wakefield.h"

static const double c  = PLASMA_C;
static const double e  = PLASMA_E;
static const double me = PLASMA_ME;

/* ============================================================
 *  L1: plasma_wavelength
 *
 *  lambda_p = 2 pi c / omega_p
 *
 * The spatial period of a relativistic electron plasma wave.
 * For ne = 10^18 cm^-3: lambda_p ~ 33 um
 * For ne = 10^19 cm^-3: lambda_p ~ 10.5 um
 *
 * Resonant wakefield excitation requires pulse length
 * L_pulse ~ lambda_p / 2 (resonance condition).
 *
 * Source: Esarey et al. Sec II.A
 *
 * Complexity: O(1)
 * ============================================================ */
double plasma_wavelength(double wp)
{
    if (wp <= 0.0) return DBL_MAX;
    return 2.0 * M_PI * c / wp;
}

/* ============================================================
 *  L1: cold_wavebreaking_field -- Tajima-Dawson field
 *
 *  E_wb = me c omega_p / e
 *       ~ 96 * sqrt(ne[cm^-3])  V/cm
 *
 * This is the maximum electric field that a non-relativistic,
 * cold plasma wave can support before the wave "breaks":
 * thermal electrons become trapped in the wave potential
 * and the ordered oscillation transitions to turbulence.
 *
 * For ne = 10^18 cm^-3: E_wb ~ 96 GV/m
 * (cf. conventional RF accelerators: ~100 MV/m => 1000x gain)
 *
 * Source: Tajima & Dawson 1979; Esarey et al. Eq. 9
 * Theorem: 1D cold fluid theory: maximum density perturbation
 *          delta n/n = 1 gives E_max = (me c wp)/e
 *
 * Complexity: O(1)
 * ============================================================ */
double cold_wavebreaking_field(double wp)
{
    if (wp <= 0.0) return 0.0;
    return me * c * wp / e;
}

/* ============================================================
 *  L1: relativistic_wavebreaking_field
 *
 *  E_wb^rel = sqrt(2 (gamma_phi - 1)) * me c omega_p / e
 *
 * For relativistic phase velocities (gamma_phi >> 1), the
 * wave-breaking threshold increases because the wave potential
 * is deeper in the moving frame.
 *
 * For typical LWFA (gamma_phi ~ 10-100), E_wb^rel is 4-14x higher
 * than the cold limit.
 *
 * Source: Esarey et al. Eq. 12
 *
 * Complexity: O(1)
 * ============================================================ */
double relativistic_wavebreaking_field(double wp, double gamma_phi)
{
    if (wp <= 0.0 || gamma_phi < 1.0) return 0.0;
    double factor = sqrt(2.0 * (gamma_phi - 1.0));
    return factor * cold_wavebreaking_field(wp);
}

/* ============================================================
 *  L1: wakefield_gradient_1D_linear -- Linear regime
 *
 *  E_acc = (a0^2 / sqrt(1 + a0^2/2)) E_wb
 *
 * For a0 << 1: the wakefield amplitude is proportional to a0^2
 * (the ponderomotive force is proportional to intensity).
 *
 * This formula assumes a Gaussian pulse profile and includes
 * the first relativistic correction factor.
 *
 * Source: Esarey et al. Eq. 33
 *
 * Complexity: O(1)
 * ============================================================ */
double wakefield_gradient_1D_linear(double a0, double E_wb)
{
    if (a0 < 0.0 || E_wb < 0.0) return 0.0;
    double denom = sqrt(1.0 + a0 * a0 / 2.0);
    return (a0 * a0 / denom) * E_wb;
}

/* ============================================================
 *  L1: wakefield_gradient_1D_nonlinear -- Nonlinear regime
 *
 *  E_max / E_wb ~ a0^2 / (2 sqrt(1 + a0^2))  for a0 > 1
 *
 * In the nonlinear regime, the plasma wave profile becomes
 * spiked (sawtooth waveform) and the peak field can exceed
 * the linear prediction.  This formula captures the 1D
 * nonlinear wave solution.
 *
 * For a0 = 2: E_max ~ 0.45 E_wb
 * For a0 = 4: E_max ~ 0.97 E_wb
 *
 * Source: Esarey et al. Eq. 53 (Akhiezer-Polovin solution)
 *
 * Complexity: O(1)
 * ============================================================ */
double wakefield_gradient_1D_nonlinear(double a0, double E_wb)
{
    if (a0 < 0.0 || E_wb < 0.0) return 0.0;
    return (a0 * a0) / (2.0 * sqrt(1.0 + a0 * a0)) * E_wb;
}

/* ============================================================
 *  L2: bubble_accelerating_field -- Blowout regime
 *
 *  E_acc ~ 0.5 * (me c omega_p / e)
 *
 * In the bubble (blowout) regime (a0 >= 2), the laser pulse
 * completely expels electrons from its path, creating a
 * pure ion cavity.  The accelerating field inside this
 * bubble is roughly uniform and approximately half the
 * wave-breaking field.
 *
 * This is the regime used in all GeV-class LWFA experiments.
 *
 * Source: Lu et al. (2006); Pukhov & Meyer-ter-Vehn (2002)
 *
 * Complexity: O(1)
 * ============================================================ */
double bubble_accelerating_field(double wp)
{
    if (wp <= 0.0) return 0.0;
    return 0.5 * cold_wavebreaking_field(wp);
}

/* ============================================================
 *  L2: dephasing_length_1D -- Relativistic dephasing length
 *
 *  L_deph = (lambda_p / 2) * (nc / ne)
 *         = (lambda_p / 2) * (omega^2 / omega_p^2)
 *
 * This is the distance over which a relativistic electron
 * (v ~ c) outruns the accelerating phase of a plasma wave
 * (v_phi ~ v_g < c) by half a plasma wavelength, moving
 * from the accelerating to the decelerating phase.
 *
 * For ne = 10^18 cm^-3, lambda=0.8 um:
 *   L_deph ~ 8.3 mm  -> theoretical max energy ~ 0.8 GeV
 *
 * Source: Esarey et al. Sec IV.A, Eq. 88
 * Theorem: Phase slip rate = c - v_g ~ c wp^2/(2 omega^2)
 *          => L_deph = (lambda_p/2) / (1 - v_g/c) = (lambda_p/2)*(omega^2/wp^2)
 *
 * Complexity: O(1)
 * ============================================================ */
double dephasing_length_1D(double wp, double omega, double nc, double ne)
{
    if (wp <= 0.0 || omega <= 0.0 || nc <= 0.0 || ne <= 0.0)
        return 0.0;

    double lambda_p = plasma_wavelength(wp);
    return (lambda_p / 2.0) * (nc / ne);
}

/* ============================================================
 *  L2: pump_depletion_length
 *
 *  L_pump = (omega^2 / omega_p^2) * c * tau_fwhm / a0^2
 *
 * Distance over which the laser pulse transfers an order-unity
 * fraction of its energy to the plasma wave.  The scaling
 * omega^2/wp^2 reflects the fact that most laser energy is EM
 * field energy, while only the ponderomotive component couples
 * to the plasma.
 *
 * For matched LWFA operation, L_deph ~ L_pump, which sets the
 * optimal a0 for a given density:
 *   a0_opt ~ (omega_p tau_fwhm)^{1/2}
 *
 * Source: Esarey et al. Sec IV.B, Eq. 97
 *
 * Complexity: O(1)
 * ============================================================ */
double pump_depletion_length(double omega, double wp, double a0,
                              double tau_fwhm)
{
    if (omega <= 0.0 || wp <= 0.0 || a0 <= 0.0 || tau_fwhm <= 0.0)
        return DBL_MAX;

    double omega_ratio_sq = (omega * omega) / (wp * wp);
    return omega_ratio_sq * c * tau_fwhm / (a0 * a0);
}

/* ============================================================
 *  L2: acceleration_length
 *
 *  L_acc = min(L_deph, L_pump, L_plasma)
 *
 * The achievable acceleration length is limited by the most
 * restrictive of three constraints:
 *
 *   1. Dephasing: electron outruns the accelerating wave
 *   2. Pump depletion: laser runs out of energy
 *   3. Plasma extent: physical end of the plasma medium
 *
 * In an optimized LWFA, L_deph ~ L_pump => operation at
 * the "matched" condition.
 *
 * Complexity: O(1)
 * ============================================================ */
double acceleration_length(double L_deph, double L_pump, double L_plasma)
{
    double L_min = L_deph;
    if (L_pump < L_min) L_min = L_pump;
    if (L_plasma < L_min) L_min = L_plasma;
    return L_min;
}

/* ============================================================
 *  L2: maximum_energy_gain
 *
 *  W_max = e E_acc L_acc  ->  W_max[eV] = E_acc * L_acc
 *
 * For dephasing-limited operation:
 *   W_max ~ me c^2 * (omega^2 / omega_p^2) * (E_acc / E_wb)
 *        ~ 0.511 MeV * (nc/ne) * (E_acc / E_wb)
 *
 * For ne = 10^18 cm^-3, lambda = 0.8 um:
 *   nc/ne ~ 1.7e3, W_max ~ 0.8 GeV (assuming E_acc = E_wb)
 *
 * This is why LWFA can achieve GeV energies in cm-scale plasmas.
 *
 * Source: Esarey et al. Sec IV.D
 *
 * Complexity: O(1)
 * ============================================================ */
double maximum_energy_gain(double E_acc, double L_acc)
{
    return E_acc * L_acc;  /* returns energy in eV if E_acc in V/m */
}

/* ============================================================
 *  L2: optimal_pulse_duration -- Resonant condition
 *
 *  tau_opt = lambda_p / (2 c) = pi / omega_p
 *
 * A laser pulse with FWHM duration matching half the plasma
 * period most efficiently excites a wakefield.
 *
 * For ne = 10^18 cm^-3: omega_p ~ 1.8e13 rad/s
 *   => tau_opt ~ 170 fs
 *
 * This is why LWFA requires femtosecond laser pulses.
 *
 * Complexity: O(1)
 * ============================================================ */
double optimal_pulse_duration(double wp)
{
    if (wp <= 0.0) return DBL_MAX;
    return M_PI / wp;
}

/* ============================================================
 *  L2: optimal_spot_size -- Matched laser spot
 *
 *  w0_match = (2 sqrt(a0) c) / omega_p
 *
 * For a given a0, this spot size produces a bubble radius
 * that matches the laser spot, enabling efficient guiding
 * and energy transfer.
 *
 * Source: Lu et al. (2006), Eq. 4
 *
 * Complexity: O(1)
 * ============================================================ */
double optimal_spot_size(double a0, double wp)
{
    if (a0 < 0.0 || wp <= 0.0) return 0.0;
    return 2.0 * sqrt(a0) * c / wp;
}

/* ============================================================
 *  L2: beam_loading_limit
 *
 *  n_bunch_max ~ E_acc / (e lambda_p) * (rb / lambda_p)^2
 *
 * Maximum charge density that can be accelerated without
 * significantly reducing the wakefield amplitude.  Beam
 * loading flattens the accelerating field, which can be
 * beneficial for producing monoenergetic electron beams.
 *
 * Source: Esarey et al. Sec VI
 *
 * Complexity: O(1)
 * ============================================================ */
double beam_loading_limit(double E_acc, double wp, double rb)
{
    if (E_acc <= 0.0 || wp <= 0.0 || rb <= 0.0) return 0.0;

    double lambda_p = plasma_wavelength(wp);
    double aspect_ratio = rb / lambda_p;
    return (E_acc / (e * lambda_p)) * aspect_ratio * aspect_ratio;
}

/* ============================================================
 *  L4: compute_wakefield_params -- Complete LWFA analysis
 *
 * Fills the WakefieldParams structure with all relevant
 * acceleration parameters for the given plasma and laser
 * conditions.
 *
 * Returns 0 on success, -1 for unphysical input.
 *
 * Complexity: O(1)
 * ============================================================ */
int compute_wakefield_params(double ne, double nc, double omega,
                              double wp, double a0, double tau_fwhm,
                              double L_plasma, WakefieldParams *w)
{
    if (!w) return -1;
    if (ne <= 0.0 || nc <= 0.0 || omega <= 0.0 || wp <= 0.0
        || a0 < 0.0 || tau_fwhm <= 0.0)
        return -1;

    w->lambda_p = plasma_wavelength(wp);
    w->E_wb     = cold_wavebreaking_field(wp);

    /* Choose accelerating field based on regime */
    if (a0 < 0.5) {
        w->E_acc = wakefield_gradient_1D_linear(a0, w->E_wb);
    } else if (a0 < 2.0) {
        w->E_acc = wakefield_gradient_1D_nonlinear(a0, w->E_wb);
    } else {
        w->E_acc = bubble_accelerating_field(wp);
    }

    w->L_deph   = dephasing_length_1D(wp, omega, nc, ne);
    w->L_pump   = pump_depletion_length(omega, wp, a0, tau_fwhm);
    w->L_acc    = acceleration_length(w->L_deph, w->L_pump, L_plasma);
    w->W_max    = maximum_energy_gain(w->E_acc, w->L_acc);
    w->tau_opt  = optimal_pulse_duration(wp);
    w->a0_opt   = sqrt(wp * tau_fwhm / (2.0 * M_PI));
    w->n_bunch  = beam_loading_limit(w->E_acc, wp,
                                      optimal_spot_size(a0, wp));

    return 0;
}
