/**
 * example_magnetosphere.c -- Earth's Magnetosphere Demo
 *
 * Demonstrates:
 *   1. Chapman-Ferraro magnetopause standoff distance
 *   2. Bow shock location and MHD shock jumps
 *   3. Plasmapause location (Carpenter-Anderson model)
 *   4. Ring current energy and Dst relation
 *   5. Particle drifts in dipole field
 *
 * Reference: Kivelson & Russell S7-S10
 *            MIT 22.611 Lectures 9-14
 */
#include "../include/space_plasma.h"
#include "../include/plasma_parameters.h"
#include "../include/magnetosphere.h"
#include <stdio.h>
#include <math.h>

int main(void) {
    printf("=== Earth Magnetosphere Model ===\n\n");

    /* Solar wind conditions */
    solar_wind_t sw_q, sw_fast, sw_cme;
    solar_wind_quiet_sun(&sw_q);
    solar_wind_fast(&sw_fast);
    solar_wind_cme(&sw_cme);

    /* Chapman-Ferraro magnetopause */
    printf("Magnetopause Standoff Distance:\n");
    printf("  %-20s %-10s %-10s %-10s\n",
           "Condition", "n[cm^-3]", "v[km/s]", "R_mp[R_E]");

    solar_wind_t *sw_list[] = {&sw_q, &sw_fast, &sw_cme};
    const char *names[] = {"Quiet Sun", "Fast Wind", "CME"};

    for (int i = 0; i < 3; i++) {
        double R_mp = chapman_ferraro_standoff(sw_list[i]->n_sw,
                                                sw_list[i]->v_sw, 2.44);
        printf("  %-20s %-10.1f %-10.0f %-10.1f\n",
               names[i], sw_list[i]->n_sw/1e6,
               sw_list[i]->v_sw/1e3, R_mp);
    }
    printf("\n");

    /* Bow shock */
    printf("Bow Shock Location:\n");
    for (int i = 0; i < 3; i++) {
        double R_mp = chapman_ferraro_standoff(sw_list[i]->n_sw,
                                                sw_list[i]->v_sw, 2.44);
        double R_bs = bow_shock_standoff(R_mp, sw_list[i]->M_ms, 5.0/3.0);
        printf("  %-20s R_bs = %.2f R_E  (standoff = %.2f R_E)\n",
               names[i], R_bs, R_bs - R_mp);
    }
    printf("\n");

    /* MHD Shock Jump Conditions */
    printf("MHD Shock Jump (perpendicular):\n");
    double M_A = sw_q.Ma;
    printf("  Upstream Mach M_A = %.2f\n", M_A);
    double rho_r, p_r, B_r;
    mhd_shock_jump(M_A, 5.0/3.0, &rho_r, &p_r, &B_r);
    printf("  Density jump:  rho2/rho1 = %.2f\n", rho_r);
    printf("  Pressure jump: p2/p1 = %.2f\n", p_r);
    printf("  B-field jump:  B2/B1 = %.2f\n\n", B_r);

    /* Plasmasphere */
    printf("Plasmasphere Profile (Carpenter-Anderson):\n");
    double Kp_vals[] = {1.0, 3.0, 6.0};
    const char *kp_desc[] = {"Quiet", "Unsettled", "Storm"};

    for (int k = 0; k < 3; k++) {
        double L_pp = plasmapause_L(Kp_vals[k]);
        printf("  Kp=%.0f (%s): L_pp = %.2f R_E\n", Kp_vals[k], kp_desc[k], L_pp);
        printf("    Density: L=2: %.0f cm^-3, L=%d: %.1f cm^-3, L=6: %.1f cm^-3\n",
               plasmasphere_density(2.0, L_pp, 1000.0, 1.0, 3.5),
               (int)L_pp,
               plasmasphere_density((double)(int)L_pp, L_pp, 1000.0, 1.0, 3.5),
               plasmasphere_density(6.0, L_pp, 1000.0, 1.0, 3.5));
    }
    printf("\n");

    /* Ring Current and Dst */
    printf("Ring Current Dst Relation:\n");
    double U_R = 1.0e22;  /* ring current energy [J] */
    double Dst = dessler_parker_sckopke(U_R);
    printf("  U_R = %.1e J -> Dst = %.0f nT\n", U_R, Dst);
    double U_density = ring_current_energy_density(1.0e6, 3.0e4, 3.0e4);
    printf("  Energy density at L=4 (isotropic, 30 keV): %.2e J/m^3\n\n", U_density);

    /* Charged particle drifts */
    printf("Particle Drifts in Dipole Field:\n");
    printf("  %-15s %-12s %-12s %-12s %-15s\n",
           "Particle", "E_kin", "L-shell", "Drift Dir", "Bounce Period");

    struct { const char *name; double m; double q; double E; double L; } parts[] = {
        {"1 MeV e-",     SP_ME, -SP_EC, 1e6*SP_EC, 4.0},
        {"100 keV H+",   SP_MP,  SP_EC, 1e5*SP_EC, 6.0},
        {"10 keV O+",    16.0*SP_MP, SP_EC, 1e4*SP_EC, 3.0},
        {"50 MeV proton", SP_MP, SP_EC, 5e7*SP_EC, 2.0},
    };

    for (int i = 0; i < 4; i++) {
        double v_d[3];
        dipole_drift_speed(parts[i].E, parts[i].L, parts[i].q, v_d);
        double tau_b = dipole_bounce_period(parts[i].E, parts[i].m,
                                             parts[i].L, M_PI/2.0);
        const char *dir = (v_d[1] > 0) ? "Eastward" : "Westward";
        printf("  %-15s %-12.1e %-12.1f %-12s %-12.2f s\n",
               parts[i].name, parts[i].E/SP_EC, parts[i].L, dir, tau_b);
    }
    printf("\n");

    /* E x B convection */
    printf("Magnetospheric Convection:\n");
    double E_conv[3] = {0.0, 0.5e-3, 0.0};  /* 0.5 mV/m dawn-dusk */
    double B_dip[3];
    double x_eq[3] = {6.0, 0.0, 0.0};  /* equatorial at L=6 */

    earth_dipole_field(x_eq, SP_MU_EARTH, B_dip);
    double v_conv[3];
    magnetosphere_exb_drift(E_conv, B_dip, v_conv);

    printf("  At L=6 equator, E_dusk=0.5 mV/m:\n");
    printf("  B = (%.1f, %.1f, %.1f) nT\n",
           B_dip[0]*1e9, B_dip[1]*1e9, B_dip[2]*1e9);
    printf("  E x B drift v = (%.1f, %.1f, %.1f) km/s\n",
           v_conv[0]/1e3, v_conv[1]/1e3, v_conv[2]/1e3);

    /* Cross-polar cap potential */
    double cpcp = cross_polar_cap_potential(400.0, 5.0);
    printf("  Cross-polar cap potential (v=400km/s, Bz=-5nT): %.0f kV\n\n", cpcp);

    printf("=== Done ===\n");
    return 0;
}
