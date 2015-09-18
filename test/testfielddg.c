#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <assert.h>
#include <math.h>

int TestfieldDG()
{
  int test = true;

  field f;
  init_empty_field(&f);
  
  f.model.cfl = 0.05;
  f.model.m = 1; // only one conservative variable
  f.model.NumFlux = TransNumFlux;
  f.model.BoundaryFlux = TestTransBoundaryFlux;
  f.model.InitData = TestTransInitData;
  f.model.ImposedData = TestTransImposedData;
  f.model.Source = NULL;
  f.varindex = GenericVarindex;

  f.deg[0] = 2; // x direction degree
  f.deg[1] = 2; // y direction degree
  f.deg[2] = 2; // z direction degree
  f.raf[0] = 1; // x direction refinement
  f.raf[1] = 1; // y direction refinement
  f.raf[2] = 1; // z direction refinement
  // FIXME: set back to 2.

  ReadMacroMesh(&f.macromesh, "../test/testcube2.msh");
  BuildConnectivity(&f.macromesh);

  PrintMacroMesh(&f.macromesh);
  PrintMacroMesh(&f.macromesh);
  
  Initfield(&f);
  CheckMacroMesh(&f.macromesh, f.deg, f.raf);

  real tnow = 0.0;
  dtfield(&f, f.wn, f.dtwn, tnow);
  
  Displayfield(&f);

  /* Plotfield(0, false, &f, NULL, "visu.msh"); */
  /* Plotfield(0, true, &f, "error", "error.msh"); */

  real maxerr = 0.0;

  // Test the time derivative with the exact solution
  const int wsize = f.model.m * f.macromesh.nbelems * NPG(f.deg, f.raf); 
  for(int i = 0; i < wsize; i++) {
    real err = fabs(4 * f.wn[i] - pow(f.dtwn[i], 2));
    if(err > maxerr)
      maxerr = err;
    printf("i: %d\terr: %f\tdtw: %f\tw: %f\n", i, err, f.dtwn[i], f.wn[i]);
  }

  test = maxerr < 1e-2;
  printf("maxerr: %f\n", maxerr);
  
  return test;
}

int main(void) {
  // Unit tests
  int resu = TestfieldDG();
  if (resu) printf("field DG test OK !\n");
  else printf("field DG test failed !\n");
  return !resu;
} 
