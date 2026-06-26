#include "surface.h"
#include <stdio.h>
#include <assert.h>

int main(void) {
    double Y = sigmund_sputtering_yield(500.0, 6.63e-26, 4.65e-26, 4.7, 0.3); assert(Y>=0);
    double Yb = bohdansky_sputtering_yield(500.0, 25.0, 2.0, 0.5); assert(Yb>0);
    EtchRates er; compute_etch_rates(&er, 3.0, 1e16, -200.0, 1e20, 0.1, "Si", "CF4"); assert(er.ER_total>0);
    double ER = rie_lag_etch_rate(500.0, 5.0, 0.1); assert(ER>0 && ER<500);
    double R = pecvd_deposition_rate(523.0, 1e20, 0.2*1.602e-19, 1e-10); assert(R>=0);
    double th = surface_coverage_langmuir(1.0, 0.5, 300.0, 0.3); assert(th>0 && th<=1);
    EtchModel m; run_rie_process_model(&m, "Si", "SF6", 500.0, 5.0, 100.0, 300.0, 300.0, 1.0); assert(m.etch_rate_total>0);
    double epc = ale_etch_per_cycle(300.0, 50.0, 0.8, 5.0); assert(epc>=0);
    printf("All surface tests passed.\n");
    return 0;
}
