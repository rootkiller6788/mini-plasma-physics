#ifndef WAKEFIELD_H
#define WAKEFIELD_H

/*
 * wakefield.h -- Laser wakefield acceleration (LWFA)
 *
 * Models for electron plasma wave excitation by an ultra-short laser
 * pulse and subsequent electron acceleration in the wakefield.
 *
 * References:
 *   - Tajima & Dawson (1979), "Laser Electron Accelerator", PRL 43, 267
 *   - Esarey, Schroeder & Leemans (2009), Rev. Mod. Phys. 81, 1229
 *   - Gibbon (2005), Ch. 7
 *
 * Knowledge Layers:
 *   L1: Wakefield gradient, plasma wavelength, dephasing length
 *   L2: Bubble regime, beam loading, self-injection
 *   L4: Wakefield excitation equations, 1D nonlinear theory
 *   L5: PIC-based wakefield simulation
 *
 * Courses: MIT 22.611, Stanford PHYSICS 450, Berkeley PHYS 242
 */

#include "plasma_constants.h"
#include "plasma_params.h"

/* ============================================================
 *  L1: Wakefield Acceleration Parameters
 * ============================================================ */

/**
 * WakefieldParams -- LWFA operating regime parameters
 */
typedef struct {
    double lambda_p;    /* plasma wavelength 2pi c/wp [m]       */
    double E_wb;        /* cold wave-breaking field [V/m]        */
    double E_acc;       /* accelerating gradient (operating)     */
    double L_deph;      /* dephasing length [m]                  */
    double L_pump;      /* pump depletion length [m]             */
    double L_acc;       /* acceleration length [m]               */
    double W_max;       /* maximum energy gain [eV]              */
    double tau_opt;     /* optimal pulse duration (resonant) [s]  */
    double a0_opt;      /* optimal a0 for blowout                 */
    double n_bunch;     /* beam loading density [m^-3]            */
} WakefieldParams;

/**
 * BubbleState -- Descriptor for the bubble (blowout) regime
 *
 * rb      : bubble radius [m]
 * E_peak  : peak accelerating field inside bubble [V/m]
 * Ez_slope: longitudinal field gradient (focusing) [V/m^2]
 */
typedef struct {
    double rb;
    double E_peak;
    double Ez_slope;
} BubbleState;

/* ============================================================
 *  L1: Core Wakefield Quantities
 * ============================================================ */

/**
 * plasma_wavelength -- Plasma wavelength
 *
 *   lambda_p = 2pi c / omega_p ~ 3.34e6 / sqrt(ne)  [m, ne in cm^-3]
 *
 * The spatial period of a relativistic plasma wave. Resonant
 * wakefield excitation requires pulse length L_pulse ~ lambda_p/2.
 * [Esarey et al. Sec II.A]
 */
double plasma_wavelength(double wp);

/**
 * cold_wavebreaking_field -- 1D cold non-relativistic limit
 *
 *   E_wb = me c omega_p / e ~ 96 sqrt(ne[cm^-3])  [V/cm]
 *
 * Maximum electric field a non-relativistic, cold plasma wave
 * can sustain before wave-breaking (trapping). This is the
 * "Tajima-Dawson field".
 * [Tajima & Dawson 1979]
 */
double cold_wavebreaking_field(double wp);

/**
 * relativistic_wavebreaking_field -- Relativistic correction
 *
 *   E_wb^rel = sqrt(2 (gamma_phi - 1)) me c omega_p / e
 *
 * For mildly relativistic phase velocities, the wave-breaking
 * limit increases.
 * [Esarey et al. Sec II.B]
 */
double relativistic_wavebreaking_field(double wp, double gamma_phi);

/**
 * wakefield_gradient_1D_linear -- Linear regime wakefield
 *
 *   E_acc = (a0^2 / sqrt(1 + a0^2)) E_wb
 *
 * For a0 << 1, the wakefield amplitude scales as a0^2.
 * For a0 ~ 1, nonlinear corrections become important.
 * [Esarey et al. Sec III.A]
 */
double wakefield_gradient_1D_linear(double a0, double E_wb);

/**
 * wakefield_gradient_1D_nonlinear -- 1D nonlinear wakefield
 *
 *   E_max = (a0^2 / 2) * (sqrt(1 + a0^2) - 1) * E_wb / a0
 *
 * In the nonlinear regime (a0 > 1), the wakefield shape becomes
 * spiked and the peak field can exceed the linear prediction.
 * [Esarey et al. Sec III.B]
 */
double wakefield_gradient_1D_nonlinear(double a0, double E_wb);

/**
 * bubble_accelerating_field -- Blowout/bubble regime field
 *
 *   E_acc ~ sqrt(ne[cm^-3]) * (me c omega_p / e) / 2
 *
 * In the bubble regime (a0 >= 2), the accelerating field is
 * roughly half the wave-breaking field and uniform in radius.
 * [Lu et al. 2006]
 */
double bubble_accelerating_field(double wp);

/* ============================================================
 *  L2: Acceleration Length Scales
 * ============================================================ */

/**
 * dephasing_length_1D -- Electron dephasing length
 *
 *   L_deph = (lambda_p / 2) * nc / ne
 *          ~ (lambda_p / 2) * (omega^2 / omega_p^2)
 *
 * Distance over which an electron outruns the accelerating phase
 * of the wakefield by half a plasma wavelength (relativistic).
 * [Esarey et al. Sec IV.A]
 */
double dephasing_length_1D(double wp, double omega, double nc, double ne);

/**
 * pump_depletion_length -- Laser pump depletion length
 *
 *   L_pump = (omega^2 / omega_p^2) * c * tau_fwhm / a0^2
 *
 * Distance over which the laser pulse loses a significant fraction
 * of its energy to the wakefield.
 * [Esarey et al. Sec IV.B]
 */
double pump_depletion_length(double omega, double wp, double a0,
                              double tau_fwhm);

/**
 * acceleration_length -- Achievable acceleration distance
 *
 *   L_acc = min(L_deph, L_pump, L_plasma)
 *
 * Net acceleration is limited by the shortest of: dephasing,
 * pump depletion, and the physical plasma length.
 */
double acceleration_length(double L_deph, double L_pump, double L_plasma);

/**
 * maximum_energy_gain -- Maximum electron energy gain in LWFA
 *
 *   W_max = e E_acc L_acc
 *
 * Product of accelerating gradient and available length.
 * For L_deph-limited: W_max ~ (nc/ne) me c^2 ~ (omega^2/wp^2) me c^2
 * [Esarey et al. Sec IV.D]
 */
double maximum_energy_gain(double E_acc, double L_acc);

/* ============================================================
 *  L2/L5: Optimal Pulse Parameters
 * ============================================================ */

/**
 * optimal_pulse_duration -- Resonant pulse duration
 *
 *   tau_opt = lambda_p / (2 c)
 *
 * A laser pulse of this duration most efficiently drives a
 * wakefield (resonant condition).
 */
double optimal_pulse_duration(double wp);

/**
 * optimal_spot_size -- Matched laser spot size
 *
 *   w0_opt = (2 sqrt(a0) c) / omega_p
 *
 * For efficient guiding in the bubble regime, the laser spot
 * size should be matched to the bubble radius.
 * [Lu et al. 2006]
 */
double optimal_spot_size(double a0, double wp);

/**
 * beam_loading_limit -- Maximum charge density for beam loading
 *
 *   n_bunch_max ~ (E_acc / (e lambda_p)) * (rb^2 / lambda_p^2)
 *
 * When the trapped bunch charge approaches this value, its
 * self-field significantly reduces the accelerating gradient.
 * [Esarey et al. Sec VI]
 */
double beam_loading_limit(double E_acc, double wp, double rb);

/**
 * compute_wakefield_params -- Compute all LWFA parameters
 *
 * Fills WakefieldParams from given plasma density, laser wavelength,
 * laser a0, and pulse duration.
 *
 * Returns 0 on success, -1 on invalid input.
 */
int compute_wakefield_params(double ne, double nc, double omega,
                              double wp, double a0, double tau_fwhm,
                              double L_plasma, WakefieldParams *wp_out);

#endif /* WAKEFIELD_H */
