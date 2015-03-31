#include "../test.h"
#include "schnaps.h"
#include<stdio.h>
#include <assert.h>
#include <math.h>


int TestFieldSource(void){

  int test = true;

  Field f;
  //f.model.cfl = 0.05;
  f.model.m = 1; // only one conservative variable
  f.model.NumFlux = TransportNumFlux;
  f.model.BoundaryFlux = TestSourceBoundaryFlux;
  f.model.InitData = TestSourceInitData;
  f.model.ImposedData = TestSourceImposedData;
  f.model.Source = TestSourceSource;
  f.varindex = GenericVarindex;

  f.interp.interp_param[0] = 1; // _M
  f.interp.interp_param[1] = 3; // x direction degree
  f.interp.interp_param[2] = 3; // y direction degree
  f.interp.interp_param[3] = 3; // z direction degree
  f.interp.interp_param[4] = 2; // x direction refinement
  f.interp.interp_param[5] = 2; // y direction refinement
  f.interp.interp_param[6] = 2; // z direction refinement

  ReadMacroMesh(&(f.macromesh), "test/testcube2.msh");
  //ReadMacroMesh(&(f.macromesh),"test/testmacromesh.msh");
  BuildConnectivity(&(f.macromesh));

  PrintMacroMesh(&(f.macromesh));
  //AffineMapMacroMesh(&(f.macromesh));
  PrintMacroMesh(&(f.macromesh));
  
  InitField(&f);
  CheckMacroMesh(&(f.macromesh), f.interp.interp_param + 1);

  dtField(&f);
  
  DisplayField(&f);

  int yes_compare = 1;
  int no_compare = 0;

  PlotField(0,no_compare,&f,"visu.msh");
  PlotField(0,yes_compare,&f,"error.msh");

  // Test the time derivative with the exact solution
  for(int i = 0; 
      i < f.model.m * f.macromesh.nbelems * NPG(f.interp.interp_param+1); 
      i++){
    test = test && fabs(f.dtwn[i]) < 1e-10;
    assert(test);
  }
  
  return test;
};

int main(void) {
  // Unit tests
  int resu = TestFieldSource();
  if (resu) printf("field source test OK !\n");
  else printf("field source test failed !\n");
  return !resu;
} 
