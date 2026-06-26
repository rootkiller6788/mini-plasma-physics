/**
 * fusion_equilibrium.h — Operating Points and Power Balance Header
 *
 * Refs: ITER Technical Basis (2002),
 *       Sorbom et al. (2015) — SPARC,
 *       Maisonnier et al. (2005) — DEMO
 *
 * L6: ITER, DEMO, SPARC, JET operating points
 * L7: Fusion power plant design, tritium fuel cycle
 */

#ifndef MINI_FUSION_EQUILIBRIUM_H
#define MINI_FUSION_EQUILIBRIUM_H

#include "fusion_plasma.h"

/* L6: Operating Points */
void iter_operating_point(PlasmaParameters *p);
void demo_operating_point(PlasmaParameters *p);
void sparc_operating_point(PlasmaParameters *p);
void jet_dt_record(PlasmaParameters *p);
void iter_power_balance(const PlasmaParameters *p, EnergyBalance *eb);

/* L7: Reactor Scaling */
double fusion_power_scaling(const PlasmaParameters *p);
double optimal_fusion_temperature_dt(void);
double required_triple_product_for_Q(double Q_target);
double capital_cost_scaling(const PlasmaParameters *p);

/* L7: Tritium Fuel Cycle */
double lithium_enrichment_for_tbr(double TBR_target, double natural_capability);
double tritium_inventory_required(double P_fusion, double t_burn, double f_burn);
double tritium_doubling_time(double inventory, double production, double consumption);

#endif /* MINI_FUSION_EQUILIBRIUM_H */