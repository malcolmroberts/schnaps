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


  field f;
  init_empty_field(&f);  

  f.model.m = _INDEX_MAX; // num of conservative variables
  f.vmax = _VMAX; // maximal wave speed 
  f.model.NumFlux = VlasovP_Lagrangian_NumFlux;
  f.model.Source = NULL;
  
  f.model.BoundaryFlux = TestPeriodic_BoundaryFlux;
  f.model.InitData = TestPeriodic_InitData;
  f.model.ImposedData = TestPeriodic_ImposedData;
 
  f.varindex=GenericVarindex;
  f.pre_dtfield = NULL;
  f.post_dtfield = NULL;
  f.update_after_rk = NULL; 
  f.model.cfl = 0.05;

  f.deg[0] = 2;  // x direction degree
  f.deg[1] = 2;  // y direction degree
  f.deg[2] = 0;  // z direction degree
  f.raf[0] = 4;  // x direction refinement
  f.raf[1] = 4;  // y direction refinement
  f.raf[2] = 1;  // z direction refinement

  // read the gmsh file
  ReadMacroMesh(&f.macromesh, "test/testcube.msh");
  // try to detect a 2d mesh
  Detect1DMacroMesh(&f.macromesh);
  assert(f.macromesh.is1d);

  // mesh preparation
  f.macromesh.period[0] = 1;

  BuildConnectivity(&f.macromesh);

  PrintMacroMesh(&f.macromesh);
 
  // prepare the initial fields
  Initfield(&f);
  f.nb_diags = 0;

  CheckMacroMesh(&f.macromesh, f.deg, f.raf);

  printf("minimum grid spacing: %f\n", f.hmin);

  f.vmax = _VMAX;
  real dt = set_dt(&f);
  real tmax = 0.5;
  RK2(&f, tmax, dt);
 
  // save the results and the error
  Plotfield(0, false, &f, "sol","dgvisu.msh");
  Plotfield(0, true, &f, "error","dgerror.msh");

  real dd = L2error(&f);
  real dd_Kinetic = L2_Kinetic_error(&f);
  
  printf("kinetic L2 error: %lf\n", dd_Kinetic);
  printf("L2 error: %lf\n", dd);
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


