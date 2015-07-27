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
void TestSteady_Wave_Source(const real *xy, const real t, const real *w, real *S);


void Wave_Upwind_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
			      real *flux);




int main(void) {
  
  // unit tests
    
  int resu = Test_Wave_Steady();
	 
  if (resu) printf("wave steady  test OK !\n");
  else printf("wave steady test failed !\n");

  return !resu;
} 

int Test_Wave_Steady(void) {

  bool test = true;

  MacroMesh mesh;
  ReadMacroMesh(&mesh,"test/testcube.msh");
  Detect2DMacroMesh(&mesh);
  
  real A[3][3] = {{_LENGTH_DOMAIN, 0, 0}, {0, _LENGTH_DOMAIN, 0}, {0, 0,1}};
  real x0[3] = {0, 0, 0};
  AffineMapMacroMesh(&mesh,A,x0);

  BuildConnectivity(&mesh);

  Model model;

  model.m=3;
  model.NumFlux=Wave_Upwind_NumFlux;
 
  model.InitData = TestSteady_Wave_InitData;
  model.ImposedData = TestSteady_Wave_ImposedData;
  model.BoundaryFlux = Wave_Upwind_BoundaryFlux;
  model.Source = TestSteady_Wave_Source;

  /* model.m=1; */
  /* model.NumFlux = TransNumFlux2d; */
  /* model.BoundaryFlux = TransBoundaryFlux2d; */
  /* model.InitData = TransInitData2d; */
  /* model.ImposedData = TransImposedData2d; */
  /* model.Source = NULL; */

  int deg[]={3, 3, 0};
  int raf[]={2, 2, 1};


  assert(mesh.is2d);

  CheckMacroMesh(&mesh, deg, raf);
  Simulation simu;

  InitSimulation(&simu, &mesh, deg, raf, &model);
 
  real tmax = .1;
  simu.cfl=0.2;
  simu.vmax=_SPEED_WAVE;
  RK4(&simu,tmax);

  real dd = 0;
  dd = L2error(&simu);

  printf("erreur L2=%f\n", dd);

  real tolerance = 0.001;

  test = dd < tolerance;

  // Save the results and the error
  PlotFields(0,false, &simu, "p", "dgvisup.msh");
  PlotFields(1,false, &simu, "u1", "dgvisuu1.msh");
  PlotFields(2,false, &simu, "u2", "dgvisuu2.msh");

  LinearSolver solver;
  
  InitImplicitLinearSolver(&simu, &solver);
  

  return test;
}



void TestSteady_Wave_ImposedData(const real *xy, const real t, real *w) {

  real x=xy[0];
  real y=xy[1];


  w[0] = a[3] * x * x + a[5] * x * y + a[4] * y * y + a[1] * x + a[2] * y + a[0];
  w[1] = b[3] * x * x + b[5] * x * y + b[4] * y * y + b[1] * x + b[2] * y + b[0];
  w[2] = c[3] * x * x + c[5] * x * y + c[4] * y * y + c[1] * x + c[2] * y + c[0];


}

void TestSteady_Wave_Source(const real *xy, const real t, const real *w, real *S){
  
  real x=xy[0];
  real y=xy[1];

  S[0] = 2 * x * b[3] + x * c[5] + y * b[5] + 2 * y * c[4] + b[1] + c[2];
  S[1] = 2 * x * a[3] + y * a[5] + a[1];
  S[2] = x * a[5] + 2 * y * a[4] + a[2];

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


