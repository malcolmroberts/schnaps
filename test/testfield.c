#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <assert.h>
#include <math.h>

int Testfield(void){
  int test = true;

  field f;
  init_empty_field(&f);
  
  f.model.cfl = 0.05;
  f.model.m = 1; // only one conservative variable
  f.model.NumFlux = TransNumFlux;
  f.model.BoundaryFlux = TestTransBoundaryFlux;
  f.model.InitData = TestTransInitData;
  f.model.ImposedData = TestTransImposedData;
  f.varindex = GenericVarindex;

  f.deg[0] = 2; // x direction degree
  f.deg[1] = 2; // y direction degree
  f.deg[2] = 2; // z direction degree
  f.raf[0] = 2; // x direction refinement
  f.raf[1] = 2; // y direction refinement
  f.raf[2] = 2; // z direction refinement

  ReadMacroMesh(&f.macromesh,"../test/testmacromesh.msh");
  BuildConnectivity(&f.macromesh);

  Initfield(&f);
  CheckMacroMesh(&f.macromesh, f.deg, f.raf);

  Plotfield(0, false, &f, NULL, "testvisufield.msh");
  
  return test;
}

int main(void) {
  // Unit tests
  int resu = Testfield();
  if (resu)
    printf("field test OK !\n");
  else 
    printf("field test failed !\n");
  return !resu;
} 
