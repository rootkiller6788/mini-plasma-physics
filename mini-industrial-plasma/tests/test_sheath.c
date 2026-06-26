#include "sheath.h"
#include <stdio.h>
#include <assert.h>
#include <math.h>

int main(void) {
    double cs = bohm_velocity(3.0, PLASMA_M_AR);
    assert(cs > 0);
    printf("PASS: Bohm velocity = %.1f m/s\n", cs);
    assert(bohm_criterion_check(cs*1.01, 3.0, PLASMA_M_AR) == 1);
    assert(bohm_criterion_check(cs*0.5, 3.0, PLASMA_M_AR) == 0);
    double lamD = debye_length_electron(1e16, 3.0);
    double s = debye_sheath_thickness(200.0, 3.0, lamD);
    assert(s > 0);
    double J = child_langmuir_current_density(200.0, s, PLASMA_M_AR);
    assert(J > 0);
    double Vf = floating_potential(3.0, PLASMA_M_AR);
    assert(Vf < 0);
    double C = sheath_capacitance_per_area(200.0, 3.0, lamD);
    assert(C > 0 && C < INFINITY);
    double sm = matrix_sheath_width(200.0, 1e16);
    assert(sm > 0);
    printf("All sheath tests passed.\n");
    return 0;
}
