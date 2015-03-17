#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "test.h"
#include "collision.h"

int TestCollision(void) {
  bool test = true;
  field f;
  int vec = 1;
  
  f.model.m = _MV; // num of conservative variables
  f.model.cfl = 0.05;
  f.model.NumFlux = CollisionNumFlux;
  f.model.BoundaryFlux = CollisionBoundaryFlux;
  f.model.InitData = CollisionInitData;
  f.model.ImposedData = CollisionImposedData;
  f.varindex = GenericVarindex;
  
  f.interp.interp_param[0] = _MV;  // _M
  f.interp.interp_param[1] = 2;  // x direction degree
  f.interp.interp_param[2] = 1;  // y direction degree
  f.interp.interp_param[3] = 0;  // z direction degree
  f.interp.interp_param[4] = 8;  // x direction refinement
  f.interp.interp_param[5] = 1;  // y direction refinement
  f.interp.interp_param[6]=  1;  // z direction refinement
  // read the gmsh file
  ReadMacroMesh(&(f.macromesh), "test/testcube.msh");
  // try to detect a 2d mesh
  Detect2DMacroMesh(&(f.macromesh));
  assert(f.macromesh.is2d);

  // mesh preparation
  BuildConnectivity(&(f.macromesh));

  //AffineMapMacroMesh(&(f.macromesh));
 
  // prepare the initial fields
  Initfield(&f);
  f.is2d = true;

  // prudence...
  CheckMacroMesh(&(f.macromesh), f.interp.interp_param + 1);

  printf("cfl param =%f\n", f.hmin);

  // time derivative
  //dtfield(&f);
  //Displayfield(&f);
  //assert(1==2);
  // apply the DG scheme
  // time integration by RK2 scheme 
  // up to final time = 1.
  double tmax = 1.0;
  RK2(&f, tmax);
 
  // save the results and the error
  Plotfield(0, false, &f, "dgvisu.msh", "dgvisu.msh");
  Plotfield(0, true, &f, "dgerror.msh", "dgerror.msh");

  double dd = L2error(&f);

  printf("L2 error: %lf\n", dd);
  test = test && (dd < 2e-4);
  return test;
}

int main(void) {
  int resu = TestCollision();
  
  if (resu) printf("collision test OK !\n");
  else printf("collision test failed !\n");

  return !resu;
} 


