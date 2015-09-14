#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "test.h"
#include "collision.h"
#include "waterwave2d.h"
#include "implicit.h"

const  real a[6]={1,1,1,1,1,1};
const  real b[6]={1,1,1,1,1,1};
const  real c[6]={1,1,1,1,1,1};
  
void TestSteady_Wave_ImposedData(const real *x, const real t, real *w);
void TestSteady_Wave_InitData(real *x, real *w);
void TestSteady_Wave_Source(const real *xy, const real t, const real *w, 
			    real *S);
void Wave_Upwind_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
			      real *flux);


int main() 
{
  // unit tests
    
  int resu = Test_Wave_Steady();
	 
  if (resu) printf("wave steady  test OK !\n");
  else printf("wave steady test failed !\n");

  return !resu;
} 

int Test_Wave_Steady() 
{
  bool test = true;


  field f;
  init_empty_field(&f);

  ReadMacroMesh(&f.macromesh,"../test/testcube.msh");
  Detect2DMacroMesh(&f.macromesh);
  
  real A[3][3] = {{_LENGTH_DOMAIN, 0, 0}, {0, _LENGTH_DOMAIN, 0}, {0, 0,1}};
  real x0[3] = {0, 0, 0};
  AffineMapMacroMesh(&f.macromesh,A,x0);
  BuildConnectivity(&f.macromesh);

  f.model.m=3; 
  f.model.NumFlux=Wave_Upwind_NumFlux;
  f.model.InitData = TestSteady_Wave_InitData; 
  f.model.ImposedData = TestSteady_Wave_ImposedData; 
  f.model.BoundaryFlux = Wave_Upwind_BoundaryFlux; 
  f.model.Source = TestSteady_Wave_Source; 

  f.deg[0] = 3;
  f.deg[1] = 3;
  f.deg[2] = 0;

  f.raf[0] = 2;
  f.raf[1] = 2;
  f.raf[2] = 1;
  
  CheckMacroMesh(&f.macromesh, f.deg, f.raf);

  // FIXME: init f

  real tmax = 0.01;
  real dt = 1e-5; // FIXME
  f.model.cfl=0.2;
  f.vmax=_SPEED_WAVE;
  RK4(&f,tmax, dt);
 
  real dd = 0;
  dd = L2error(&f);

  printf("erreur explicit L2=%.12e\n", dd);

  /* PlotFields(0,false, &simu, "p", "dgvisu_exp.msh"); */
  /* PlotFields(1,false, &simu, "u", "dgvisu_exu.msh"); */
  /* PlotFields(2,false, &simu, "v", "dgvisu_exv.msh"); */

  real tolerance = _SMALL;

  test = test && (dd < tolerance);

  // FIXME: init f

  ThetaTimeScheme(&f, tmax, dt);
  
  dd = L2error(&f);

  printf("erreur implicit L2=%.12e\n", dd);

  /* PlotFields(0,false, &simu2, "p", "dgvisu_imp.msh"); */
  /* PlotFields(1,false, &simu2, "u", "dgvisu_imu.msh"); */
  /* PlotFields(2,false, &simu2, "v", "dgvisu_imv.msh"); */

  test = test && (dd < tolerance);
  
  return test;
}



void TestSteady_Wave_ImposedData(const real *xy, const real t, real *w) {

  real x=xy[0];
  real y=xy[1];


  w[0] = x*(1-x)*y*(1-y)+1;
  w[1] = 2*x*(1-x)*y*(1-y)+2;
  w[2] = 3*x*(1-x)*y*(1-y)+3;


}


void TestSteady_Wave_Source(const real *xy, const real t, const real *w, real *S){
  
  real x=xy[0];
  real y=xy[1];

  S[0] = 2*(1-2*x)*(y*(1-y))+3*(1-2*y)*(x*(1-x));
  S[1] = (1-2*x)*(y*(1-y));
  S[2] = (1-2*y)*(x*(1-x));

  S[0] *= _SPEED_WAVE;
  S[1] *= _SPEED_WAVE;
  S[2] *= _SPEED_WAVE;

}

void TestSteady_Wave_InitData(real *x, real *w) {
  real t = 0;
  TestSteady_Wave_ImposedData(x, t, w);
}




void Wave_Upwind_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
				       real *flux) {
  real wR[3];
  TestSteady_Wave_ImposedData(x , t, wR);
  Wave_Upwind_NumFlux(wL, wR, vnorm, flux);
}
