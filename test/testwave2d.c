#include "implicit.h"
#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "test.h"
#include "collision.h"
#include "waterwave2d.h"



int main(void) {
  
  // unit tests
    
  int resu = Test_Wave_Periodic();
	 
  if (resu) printf("wave periodic  test OK !\n");
  else printf("wave periodic test failed !\n");

  return !resu;
} 

int Test_Wave_Periodic(void) {

  bool test = true;

  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testcube.msh");
  Detect2DMacroMesh(&mesh);
  
  real A[3][3] = {{_LENGTH_DOMAIN, 0, 0}, {0, _LENGTH_DOMAIN, 0}, {0, 0,1}};
  real x0[3] = {0, 0, 0};
  AffineMapMacroMesh(&mesh,A,x0);

  BuildConnectivity(&mesh);

  Model model;

  //model.m = 7;
  //model.NumFlux = Maxwell2DNumFlux_upwind;
  //model.BoundaryFlux = Maxwell2DBoundaryFlux_upwind;
  //model.InitData = Maxwell2DInitData;
  //model.ImposedData = Maxwell2DImposedData;
  //model.Source = Maxwell2DSource;
  //model.Source = NULL;

  model.m = 3; 
  model.NumFlux=Wave_Upwind_NumFlux;
  model.InitData = TestPeriodic_Wave_InitData;
  model.ImposedData = TestPeriodic_Wave_ImposedData;
  model.BoundaryFlux = Wave_Upwind_BoundaryFlux;
  model.Source = NULL;

  //model.m = 2;
  //model.NumFlux = VecTransNumFlux2d;
  //model.BoundaryFlux = VecTransBoundaryFlux2d;
  //model.InitData = VecTransInitData2d;
  //model.ImposedData = VecTransImposedData2d;
  //model.Source = NULL;

  int deg[]={3, 3, 0};
  int raf[]={16, 16, 1};


  assert(mesh.is2d);

  CheckMacroMesh(&mesh, deg, raf);
  Simulation simu;

  InitSimulation(&simu, &mesh, deg, raf, &model);
 
  real tmax = 0.5;
  simu.cfl=0.2;
  simu.vmax=_SPEED_WAVE;
  RK4(&simu,tmax);

  real dd = 0;
  dd = L2error(&simu);

  printf("erreur L2=%.12e\n", dd);

  real tolerance = 0.001;

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


void Wave_Upwind_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
				       real *flux) {
  real wR[3];
  TestPeriodic_Wave_ImposedData(x , t, wR);
  Wave_Upwind_NumFlux(wL, wR, vnorm, flux);
}


