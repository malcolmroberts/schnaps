#include "implicit.h"
#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "test.h"
#include "collision.h"
#include "waterwave2d.h"

void TestSteady_Wave_ImposedData(const real *x, const real t, real *w);
void TestSteady_Wave_InitData(real *x, real *w);
void TestSteady_Wave_Source(const real *xy, const real t, const real *w, real *S);
void Wave_Steady_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
			      real *flux);
void Wave_Periodic_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
				real *flux);


int main(void) {
  
  // unit tests
    
  int resu1 = 0;
  int resu2 = 0;
  resu1=Test_Wave_Periodic();
  resu2=Test_Wave_Steady();
	 
  if (resu1 + resu2 >1) printf("wave periodic  test OK !\n");
  else printf("wave periodic test failed !\n");

  return !resu1;
} 

int Test_Wave_Periodic(void) {

  bool test = true;

#ifdef PARALUTION 
  paralution_begin();
#endif 

  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testcube.msh");
  Detect2DMacroMesh(&mesh);
  
  real A[3][3] = {{_LENGTH_DOMAIN, 0, 0}, {0, _LENGTH_DOMAIN, 0}, {0, 0,1}};
  real x0[3] = {0, 0, 0};
  AffineMapMacroMesh(&mesh,A,x0);

  BuildConnectivity(&mesh);

  Model model;

  model.m = 3; 
  model.NumFlux=Wave_Upwind_NumFlux;
  model.InitData = TestPeriodic_Wave_InitData;
  model.ImposedData = TestPeriodic_Wave_ImposedData;
  model.BoundaryFlux = Wave_Periodic_BoundaryFlux;
  model.Source = NULL;

  int deg[]={4, 4, 0};
  int raf[]={4, 4, 1};


  assert(mesh.is2d);

  CheckMacroMesh(&mesh, deg, raf);
  Simulation simu;

  InitSimulation(&simu, &mesh, deg, raf, &model);
 
  real tmax = 0.025;
  simu.cfl=0.2;
  simu.vmax=_SPEED_WAVE;
  RK4(&simu,tmax);

  real dd = 0;
  dd = L2error(&simu);

  printf("erreur L2=%.12e\n", dd);

  real tolerance = 0.002;

  test = test && (dd < tolerance);

  PlotFields(0,false, &simu, "p", "dgvisu_exp.msh");
  PlotFields(1,false, &simu, "u", "dgvisu_exu.msh");
  PlotFields(2,false, &simu, "v", "dgvisu_exv.msh");

  Simulation simu2;

  InitSimulation(&simu2, &mesh, deg, raf, &model);
  ThetaTimeScheme(&simu2, tmax, simu.dt);

  
  dd = L2error(&simu2);

  printf("erreur implicit L2=%.12e\n", dd);

  test = test && (dd < tolerance);

  PlotFields(0,false, &simu2, "p", "dgvisu_imp.msh");
  PlotFields(1,false, &simu2, "u", "dgvisu_imu.msh");
  PlotFields(2,false, &simu2, "v", "dgvisu_imv.msh");
 
#ifdef PARALUTION 
  paralution_end();
#endif

  FreeMacroMesh(&mesh);

  return test;
}

int Test_Wave_Steady(void) {

  bool test = true;

#ifdef PARALUTION 
  paralution_begin();
#endif 

  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testcube.msh");
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
  model.BoundaryFlux = Wave_Steady_BoundaryFlux; 
  model.Source = TestSteady_Wave_Source; 

  int deg[]={3, 3, 0};
  int raf[]={2, 2, 1};
  
  CheckMacroMesh(&mesh, deg, raf);
  Simulation simu, simu2;

  InitSimulation(&simu, &mesh, deg, raf, &model);

  real tmax = 0.01;
  simu.cfl=0.2;
  simu.vmax=_SPEED_WAVE;
  RK4(&simu,tmax);
 
  real dd = 0;
  dd = L2error(&simu);

  printf("erreur explicit L2=%.12e\n", dd);

  PlotFields(0,false, &simu, "p", "dgvisu_exp.msh");
  PlotFields(1,false, &simu, "u", "dgvisu_exu.msh");
  PlotFields(2,false, &simu, "v", "dgvisu_exv.msh");

  real tolerance = 0.0000001;

  test = test && (dd < tolerance);

  InitSimulation(&simu2, &mesh, deg, raf, &model);
  ThetaTimeScheme(&simu2, tmax, simu.dt);
  
  dd = L2error(&simu2);

  printf("erreur implicit L2=%.12e\n", dd);

  PlotFields(0,false, &simu2, "p", "dgvisu_imp.msh");
  PlotFields(1,false, &simu2, "u", "dgvisu_imu.msh");
  PlotFields(2,false, &simu2, "v", "dgvisu_imv.msh");

  test = test && (dd < tolerance);
   
#ifdef PARALUTION 
  paralution_end();
#endif 

  FreeMacroMesh(&mesh);
  
  return test;
}

void TestPeriodic_Wave_ImposedData(const real *x, const real t, real *w) {
  real pi=4.0*atan(1.0);
  real L=_LENGTH_DOMAIN;
  real Coef=(2.0*pi)/L;
  real a=_SPEED_WAVE;

  w[0] = -a*Coef*sqrt(2.0)*sin(a*Coef*sqrt(2.0)*t)*cos(Coef*x[0])*cos(Coef*x[1]);
  w[1] = a*Coef*cos(Coef*a*sqrt(2.0)*t)*sin(Coef*x[0])*cos(Coef*x[1]);
  w[2] = a*Coef*cos(Coef*a*sqrt(2.0)*t)*cos(Coef*x[0])*sin(Coef*x[1]);
  

}

void TestPeriodic_Wave_InitData(real *x, real *w) {
  real t = 0;
  TestPeriodic_Wave_ImposedData(x, t, w);
}


void Wave_Periodic_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
				       real *flux) {
  real wR[3];
  TestPeriodic_Wave_ImposedData(x , t, wR);
  Wave_Upwind_NumFlux(wL, wR, vnorm, flux);
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


void Wave_Steady_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
				       real *flux) {
  real wR[3];
  TestSteady_Wave_ImposedData(x , t, wR);
  Wave_Upwind_NumFlux(wL, wR, vnorm, flux);
}
