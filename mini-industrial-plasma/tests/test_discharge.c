#include "discharge.h"
#include <stdio.h>
#include <assert.h>
#include <math.h>

int main(void) {
    double a = townsend_alpha(10.0, 10000.0, 13.6, 235.0);
    assert(a > 0);
    double Vbr = paschen_breakdown_voltage(1.0, 13.6, 235.0, 0.1);
    assert(Vbr > 0 && Vbr < INFINITY);
    double Vmin = paschen_minimum_voltage(13.6, 235.0, 0.1);
    assert(Vmin > 0);
    double Vdc = ccp_self_bias(500.0, 0.05, 0.2);
    assert(Vdc < 0);
    double d = icp_skin_depth(5e17);
    assert(d > 0);
    double Te = global_model_electron_temp(10.0, 0.15, 0.1, 2.4e20, 5.0);
    assert(Te > 0.5 && Te < 20.0);
    printf("All discharge tests passed.\n");
    return 0;
}
