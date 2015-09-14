#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <assert.h>
#include <math.h>

int TestfieldRK2_2D(void) {
  bool test = true;
  field f;
  init_empty_field(&f);

  f.model.cfl = 0.05;
  f.model.m = 1; // only one conservative variable
  f.model.NumFlux = TransNumFlux2d;
  f.model.BoundaryFlux = TransBoundaryFlux2d;
  f.model.InitData = TransInitData2d;
  f.model.ImposedData = TransImposedData2d;
  f.varindex = GenericVarindex;


  f.deg[0] = 2;  // x direction degree
  f.deg[1] = 2;  // y direction degree
  f.deg[2] = 0;  // z direction degree
  f.raf[0] = 1;  // x direction refinement
  f.raf[1] = 1;  // y direction refinement
  f.raf[2] = 1;  // z direction refinement

  ReadMacroMesh(&f.macromesh, "../test/testdisque2d.msh");
  Detect2DMacroMesh(&f.macromesh);
  // require a 2d computation
  assert(f.macromesh.is2d);
  BuildConnectivity(&f.macromesh);

  Initfield(&f);

  CheckMacroMesh(&f.macromesh, f.deg, f.raf);

  printf("cfl param =%f\n",f.hmin);

  real tmax = 0.2;
  f.vmax=1;
  real dt = 0;
  RK2(&f, tmax, dt);
 
  Plotfield(0, false, &f, NULL, "dgvisu.msh");
  Plotfield(0, true, &f, "error", "dgerror.msh");

  real dd = L2error(&f);

  printf("erreur L2=%f\n", dd);

  test = test && (dd < 0.01);

  return test;
}

int main(void) {
  int resu = TestfieldRK2_2D();
  if (resu) 
    printf("field RK2 2D test OK !\n");
  else 
    printf("field RK2 2D test failed !\n");
  return !resu;
} 
