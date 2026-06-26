#include "sheath.h"
#include "discharge.h"
#include "plasma_types.h"
#include <stdio.h>
#include <math.h>

extern void industrial_plasma_defaults(const char*,PlasmaState*,DischargeParams*,DischargeGeometry*);

int main(void) {
    printf("=== Example 1: CCP RIE Etch Process ===\n");
    PlasmaState ps; DischargeParams dp; DischargeGeometry g;
    industrial_plasma_defaults("RIE", &ps, &dp, &g);
    printf("Process: Reactive Ion Etching (CCP)\n");
    printf("RF Power: %.0f W @ %.2f MHz\n", dp.power, dp.frequency/1e6);
    printf("DC Self-Bias: %.0f V\n", dp.dc_bias);
    double Vdc_calc = ccp_self_bias(dp.voltage_amplitude, g.electrode_area, g.ground_area);
    printf("Calculated V_dc = %.0f V (area ratio = %.1f)\n", Vdc_calc, g.ground_area/g.electrode_area);
    printf("Plasma: n_e=%.1e m^-3, T_e=%.1f eV\n", ps.n_e, ps.T_e);
    printf("Debye length = %.4f mm\n", ps.lam_D*1000);
    double cs = bohm_velocity(ps.T_e, PLASMA_M_AR);
    printf("Bohm velocity = %.0f m/s\n", cs);
    double s = debye_sheath_thickness(fabs(dp.dc_bias), ps.T_e, ps.lam_D);
    printf("Sheath width = %.3f mm\n", s*1000);
    double J = child_langmuir_current_density(fabs(dp.dc_bias), s, PLASMA_M_AR);
    printf("Ion current density = %.1f A/m^2\n", J);
    double Eion = ion_impact_energy(fabs(dp.dc_bias), s, s*10);
    printf("Mean ion energy = %.0f eV\n", Eion/PLASMA_E_CHARGE);
    printf("Wafer: %.0f mm diameter\n", g.wafer_diameter*1000);
    printf("\nTypical RIE: 500W, 5 Pa, SF6 50 sccm, 13.56 MHz\n");
    printf("Expected Si etch rate: ~500 nm/min, anisotropy >0.9\n");
    return 0;
}
