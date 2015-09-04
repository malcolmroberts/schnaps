#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "test.h"
#include "collision.h"
#include "quantities_vp.h"
#include "diagnostics_vp.h"
#include "solverpoisson.h"

void TestPeriodic_ImposedData(const real *x, const real t, real *w);
void TestPeriodic_InitData(real *x, real *w);
void TestPeriodic_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
			       real *flux);

int main(void) {
  // unit tests
    
  int resu = TestPeriodic();
	 
  if (resu) printf("periodic test OK !\n");
  else printf("periodic test failed !\n");

  return !resu;
} 

int TestPeriodic(void) {
  bool test = true;

  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testcube.msh");
  Detect1DMacroMesh(&mesh);
  assert(mesh.is1d);
  // periodic mesh
  mesh.period[0] = 1;
  BuildConnectivity(&mesh);
  int deg[]={2, 0, 0};
  int raf[]={16, 1, 1};

  CheckMacroMesh(&mesh,deg,raf);
  PrintMacroMesh(&mesh);


  Model model;

  model.m = 1; // only one conservative variable
  model.NumFlux = TransNumFlux2d;
  model.BoundaryFlux = TransBoundaryFlux2d;
  model.InitData = TransInitData2d;
  model.ImposedData = TransImposedData2d;

  model.m = _INDEX_MAX; // num of conservative variables
  model.NumFlux = VlasovP_Lagrangian_NumFlux;
  model.Source = NULL;
  
  model.BoundaryFlux = TestPeriodic_BoundaryFlux;
  model.InitData = TestPeriodic_InitData;
  model.ImposedData = TestPeriodic_ImposedData;



  Simulation simu;

  InitSimulation(&simu, &mesh, deg, raf, &model);
  printf("cfl param =%f\n",simu.hmin);

  simu.vmax = _VMAX; // maximal wave speed 
  simu.cfl = 0.05;
  real tmax = 0.5;
 
  RK2(&simu,tmax);
 
  // save the results and the error
  PlotFields(0, false, &simu, "sol","dgvisu.msh");
  PlotFields(0, true, &simu, "error","dgerror.msh");

  real dd = L2error(&simu);
  //real dd_Kinetic = L2_Kinetic_error(&simu);
  //printf("erreur kinetic L2: %lf\n", dd_Kinetic);

  printf("erreur L2: %lf\n", dd);
  test = test && (dd<2e-1);

  //SolvePoisson(&f);

  return test;
}

void TestPeriodic_ImposedData(const real x[3], const real t,real w[])
{
  real pi = 4 * atan(1.0);
  for(int i = 0; i < _INDEX_MAX_KIN + 1 ; ++i) {
    int j = i % _DEG_V; // local connectivity put in function
    int nel = i / _DEG_V; // element num (TODO : function)

    real vi = (-_VMAX + nel * _DV + _DV * glop(_DEG_V, j));

    w[i] = cos(2 * pi * ( x[0] - vi * t) );
  }
  // exact value of the potential and electric field
  w[_INDEX_PHI] = 0;
  w[_INDEX_EX] = 0;
  w[_INDEX_RHO] = 2.0; //rho init
  w[_INDEX_VELOCITY] = 0; // u init
  w[_INDEX_PRESSURE] = 0; // p init
  w[_INDEX_TEMP] = 0; // e ou T init
}

void TestPeriodic_InitData(real x[3], real w[])
{
  real t = 0;
  TestPeriodic_ImposedData(x, t, w);
}

void TestPeriodic_BoundaryFlux(real x[3], real t, real wL[], real *vnorm,
			       real* flux)
{
  real wR[_MV + 6];
  TestPeriodic_ImposedData(x, t, wR);
  VlasovP_Lagrangian_NumFlux(wL, wR, vnorm, flux);
  assert(false);
}


