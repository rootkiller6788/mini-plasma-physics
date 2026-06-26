#include "mhd_defs.h"
#include "mhd_eqns.h"
#include "mhd_waves.h"
#include "mhd_equilibrium.h"
#include "mhd_instability.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#define GAMMA 1.6666666666666667
#define EPS 1e-6
static int tp=0, tf=0;
#define T(n) printf("  %-45s ", n)
#define C(c) do { if (c) { printf("PASS\n"); tp++; } else { printf("FAIL\n"); tf++; } } while(0)
#define CF(a,b,e) C(fabs((a)-(b)) < (e))

void t_constants(void){
  T("Constants non-zero"); C(MHD_MU0 > 0.0);
  T("Speed of light"); CF(MHD_C, 299792458.0, 1.0);
}
void t_state(void){
  T("State roundtrip");
  MHDState s; mhd_state_init_uniform(&s, 1,0,0,0, 0,0,1, 1);
  MHDConserved c; mhd_primitive_to_conserved(&s, GAMMA, &c);
  MHDState s2; mhd_conserved_to_primitive(&c, GAMMA, &s2);
  CF(s2.rho, s.rho, EPS); CF(s2.p, s.p, EPS);
}
void t_core(void){
  T("Sound speed"); CF(mhd_sound_speed(1,1,GAMMA), sqrt(GAMMA), EPS);
  T("Alfven speed"); CF(mhd_alfven_speed(1,1), 1.0/sqrt(MHD_MU0), EPS);
  T("Mag pressure"); CF(mhd_magnetic_pressure(1), 0.5/MHD_MU0, EPS);
  T("Plasma beta"); CF(mhd_plasma_beta(1,1), 2.0*MHD_MU0, EPS);
  T("Beta=INF at B=0"); C(isinf(mhd_plasma_beta(1,0)));
  T("Mag Reynolds"); CF(mhd_magnetic_reynolds(1e3,1,0.1), 1e4, 1e-4*1e4);
  T("Lundquist > 0"); C(mhd_lundquist(1,1,1,0.1) > 0);
  T("Mag Prandtl"); CF(mhd_magnetic_prandtl(0.01,0.1), 0.1, EPS);
}
void t_vector(void){
  T("|v|"); CF(mhd_vector_magnitude(3,4,0), 5, EPS);
  T("Dot"); CF(mhd_dot_product(1,2,3,4,5,6), 32, EPS);
  T("Cross"); double cx,cy,cz;
  mhd_cross_product(1,0,0, 0,1,0, &cx,&cy,&cz);
  CF(cx,0,EPS); CF(cy,0,EPS); CF(cz,1,EPS);
}
void t_eqns(void){
  T("div(B)=0 uniform");
  MHDState s[3][3][3]; int i,j,k;
  for(i=0;i<3;i++)for(j=0;j<3;j++)for(k=0;k<3;k++)
    mhd_state_init_uniform(&s[i][j][k], 1, 0,0,0, 0,0,1, 1);
  CF(mhd_divB(s,0.1,0.1,0.1), 0, EPS);
  T("Continuity uniform"); CF(mhd_continuity_residual(s,0.1,0.1,0.1), 0, EPS);
  T("Cross helicity"); CF(mhd_cross_helicity(1,2,3,4,5,6), 32, EPS);
  T("Mag helicity"); CF(mhd_magnetic_helicity(1,2,3,4,5,6), 32, EPS);
}
void t_waves(void){
  T("Wave speeds > 0");
  MHDState s; mhd_state_init_uniform(&s, 1, 0,0,0, 1,0,0, 1);
  MHDWaveSpeeds ws; mhd_wavespeeds_compute(&s, 1,0,0, GAMMA, &ws);
  C(ws.c_fast > 0); C(ws.c_slow <= ws.c_fast);
  T("Alfven parallel"); C(mhd_alfven_wave_speed_projected(&s,1,0,0) > 0);
  T("Alfven perp=0"); CF(mhd_alfven_wave_speed_projected(&s,0,1,0), 0, EPS);
}
void t_equil(void){
  T("Bennett T > 0"); C(mhd_bennett_temperature(1e6,1e20) > 0);
  T("Safety factor q"); CF(mhd_safety_factor_q(0.5,2,5,1), 1.25, EPS);
  T("Soloviev psi=0"); CF(mhd_soloviev_psi(2,0,2,1.5,1), 0, EPS);
}
void t_instab(void){
  T("KS q>1 stable"); C(mhd_kruskal_shafranov_condition(0.5,2,5,1) == 1);
  T("KS q<1 unstable"); C(mhd_kruskal_shafranov_condition(0.5,2,1,1) == -1);
  T("Sausage stable"); C(mhd_sausage_growth_rate(0.1,1,1,1) == 0);
  T("Kink zero q=1"); CF(mhd_kink_growth_rate(1000,1,1), 0, EPS);
}
void t_adv(void){
  T("Beltrami exact"); CF(mhd_beltrami_check(2,1,0,0,2,0,0), 0, EPS);
}
int main(void){
  printf("=== mini-mhd Tests ===\n");
  t_constants(); t_state(); t_core(); t_vector();
  t_eqns(); t_waves(); t_equil(); t_instab(); t_adv();
  printf("=== %d pass, %d fail ===\n", tp, tf);
  return tf > 0 ? 1 : 0;
}
