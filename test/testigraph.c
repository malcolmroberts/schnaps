#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "test.h"
#include "schnaps.h"

int main(void) {
  // Unit tests
  int resu=TestIGraph();
  if (resu) printf("IGraph test OK !\n");
  else printf("IGraph test failed !\n");
  return !resu;
} 

// some unit tests of the macromesh code
int TestIGraph(void)
{
  MacroMesh m;

  int deg[]={2, 2, 2};
  int raf[]={2, 2, 2};
  
  // test gmsh file reading
  //ReadMacroMesh(&m, "../test/testmacromesh.msh");
  ReadMacroMesh(&m, "cubegros.msh");
  BuildConnectivity(&m);
  CheckMacroMesh(&m, deg, raf);
  //PrintMacroMesh(&m);
  real vit[3] = {1, 1, 0};
  BuildMacroMeshGraph(&m, vit, deg, raf);
  
  int test = true;


  return test;
}
