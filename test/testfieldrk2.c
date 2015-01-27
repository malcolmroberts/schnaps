#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <assert.h>
#include <math.h>

int main(void) {
  int resu = TestFieldRK2();
  if(resu) 
    printf("Field RK2 test OK !\n");
  else 
    printf("Field RK2 test failed !\n");
  return !resu;
} 

int TestFieldRK2(void){
  int test = true;

  Field f;
  
  // 2D meshes:
  // test/disque2d.msh
  // test/testdisque2d.msh
  // test/testmacromesh.msh
  // test/unit-cube.msh

  // 3D meshes"
  // test/testdisque.msh

  char *mshname =  "test/disque2d.msh";
  
  ReadMacroMesh(&(f.macromesh), mshname);
  Detect2DMacroMesh(&(f.macromesh));
  BuildConnectivity(&(f.macromesh));

#if 1
  // 2D version
  f.model.cfl = 0.05;
  f.model.m = 1;

  f.model.NumFlux = TransNumFlux2d;
  f.model.BoundaryFlux = TransBoundaryFlux2d;
  f.model.InitData = TransInitData2d;
  f.model.ImposedData = TransImposedData2d;
  f.varindex = GenericVarindex;

  f.interp.interp_param[0] = f.model.m;
  f.interp.interp_param[1] = 2; // x direction degree
  f.interp.interp_param[2] = 2; // y direction degree
  f.interp.interp_param[3] = 0; // z direction degree
  f.interp.interp_param[4] = 4; // x direction refinement
  f.interp.interp_param[5] = 4; // y direction refinement
  f.interp.interp_param[6] = 1; // z direction refinement

  assert(f.macromesh.is2d);
#else
  // 3D version
  f.model.cfl = 0.05;
  f.model.m = 1;
  f.model.NumFlux = TransNumFlux;
  f.model.BoundaryFlux = TestTransBoundaryFlux;
  f.model.InitData = TestTransInitData;
  f.model.ImposedData = TestTransImposedData;
  f.varindex = GenericVarindex;

  f.interp.interp_param[0] = f.model.m;
  f.interp.interp_param[1] = 2; // x direction degree
  f.interp.interp_param[2] = 2; // y direction degree
  f.interp.interp_param[3] = 2; // z direction degree
  f.interp.interp_param[4] = 3; // x direction refinement
  f.interp.interp_param[5] = 3; // y direction refinement
  f.interp.interp_param[6] = 3; // z direction refinement
#endif

  // 2015-01-19: the below parameters fail with testmacrocellinterface
  // but pass the test here (perhaps because the error is hidden by
  // the RK error?)
  /*
  f.model.cfl = 0.05;
  f.model.m = 1;
  f.model.NumFlux = TransNumFlux;
  f.model.BoundaryFlux = TestTransBoundaryFlux;
  f.model.InitData = TestTransInitData;
  f.model.ImposedData = TestTransImposedData;
  f.varindex = GenericVarindex;

  f.interp.interp_param[0] = f.model.m;
  f.interp.interp_param[1] = 3; // x direction degree
  f.interp.interp_param[2] = 3; // y direction degree
  f.interp.interp_param[3] = 3; // z direction degree
  f.interp.interp_param[4] = 1; // x direction refinement
  f.interp.interp_param[5] = 1; // y direction refinement
  f.interp.interp_param[6] = 1; // z direction refinement
  */

  //AffineMapMacroMesh(&(f.macromesh));
  InitField(&f);

  CheckMacroMesh(&(f.macromesh), f.interp.interp_param + 1);
 
  double tmax = 0.1;
  double dt = 0.0;
  RK2(&f, tmax, dt);
 
  PlotField(0, false, &f, NULL, "dgvisu.msh");
  PlotField(0, true , &f, "error", "dgerror.msh");

  double dd = L2error(&f);

  printf("erreur L2=%f\n", dd);

  double tolerance = 0.001;

  test = dd < tolerance;
  
  return test;
};
