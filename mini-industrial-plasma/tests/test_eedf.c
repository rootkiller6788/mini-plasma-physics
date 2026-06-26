#include "eedf.h"
#include <stdio.h>
#include <assert.h>
#include <math.h>

int main(void) {
    double fM = maxwellian_eedf(3.0, 3.0); assert(fM > 0);
    double fD = druyvesteyn_eedf(3.0, 3.0); assert(fD > 0);
    double fBM = bimaxwellian_eedf(5.0, 2.0, 8.0, 0.1); assert(fBM > 0);
    CrossSectionModel cs = {2.5e-20, 25.0, 1.8, 1e-20, 10.0};
    double k = rate_coefficient_maxwellian(&cs, 3.0); assert(k >= 0);
    CrossSectionModel cs2[5]; int nc;
    argon_default_cross_sections(cs2, &nc); assert(nc == 3);
    printf("All EEDF tests passed.\n");
    return 0;
}
