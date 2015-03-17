#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "test.h"
#include "collision.h"

int TestCollision(void) {
  bool test = true;
  field f;
  int vec = 1;
  
  f.model.m=_MV+1; // num of conservative variables
  f.model.vmax = _VMAX; // maximal wave speed 

  f.model.cfl = 0.5; // maximal wave speed 
  f.model.NumFlux=Collision_Lagrangian_NumFlux;
  f.model.BoundaryFlux=Collision_Lagrangian_BoundaryFlux;
  f.model.InitData=CollisionInitData;
  f.model.ImposedData=CollisionImposedData;
  f.varindex=GenericVarindex;
    
  f.vmax = 1;
    
  f.interp.interp_param[0]=_MV+1;  // _M
  f.interp.interp_param[1]=3;  // x direction degree
  f.interp.interp_param[2]=0;  // y direction degree
  f.interp.interp_param[3]=0;  // z direction degree
  f.interp.interp_param[4]=10;  // x direction refinement
  f.interp.interp_param[5]=1;  // y direction refinement
  f.interp.interp_param[6]=1;  // z direction refinement

  // read the gmsh file
  ReadMacroMesh(&(f.macromesh), "test/testcube.msh");

  Detect1DMacroMesh(&(f.macromesh));
  assert(f.macromesh.is1d);

  // mesh preparation
  BuildConnectivity(&(f.macromesh));

  //AffineMapMacroMesh(&(f.macromesh));
 

  Initfield(&f);
  f.macromesh.is2d = true;
  f.macromesh.is1d = true;

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

  double dd=L2error(&f);
  double dd_Kinetic=L2_Kinetic_error(&f);
  
  printf("erreur kinetic L2=%lf\n",dd_Kinetic);
  printf("erreur L2=%lf\n",dd);
  test= test && (dd<2e-4);


  SolvePoisson(&f);

  return test;
}

int main(void) {
  int resu = TestCollision();
  
  if (resu) printf("collision test OK !\n");
  else printf("collision test failed !\n");

  return !resu;
} 





