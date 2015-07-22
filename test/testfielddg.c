#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <assert.h>
#include <math.h>


int TestfieldDG(void){

  int test = true;

  Model model;
  
  model.cfl = 0.05;
  model.m = 1; // only one conservative variable
  model.NumFlux = TransNumFlux;
  model.BoundaryFlux = TestTransBoundaryFlux;
  model.InitData = TestTransInitData;
  model.ImposedData = TestTransImposedData;

  int deg[]={1, 1, 1};
  int raf[]={1, 1, 1};
  
  MacroMesh mesh;
  //ReadMacroMesh(&mesh,"test/testmacromesh.msh");
  ReadMacroMesh(&mesh,"test/testcube2.msh");
  BuildConnectivity(&mesh);

  /* real A[3][3] = {{10,2 , 0}, {0, 1, -0.1}, {0, 0.1,1}}; */
  /* real x0[3] = {1, 2, 3}; */
  /* AffineMapMacroMesh(&mesh,A,x0); */

  CheckMacroMesh(&mesh, deg, raf);

  //PrintMacroMesh(&mesh);

  Simulation simu;

  InitSimulation(&simu, &mesh, deg, raf, &model);

  real *w = simu.fd[0].wn;
  real *dtw = simu.fd[0].dtwn;

  DtFields(&simu, w, dtw);
  
  DisplaySimulation(&simu);

  PlotFields(0, false, &simu, NULL, "visu.msh");
  //PlotFields(0, true, &simu, "error", "error.msh");

  // Test the time derivative with the exact solution
  for(int i = 0; 
      i < model.m * mesh.nbelems * NPG(deg,raf); 
      i++){
    test = test && fabs(4 * w[i] - pow(dtw[i], 2)) < 1e-2;
    printf("i=%d err=%f \n",i,4 * w[i] - pow(dtw[i], 2));
    //assert(test);
  }
  
  return test;
};

int main(void) {
  // Unit tests
  int resu = TestfieldDG();
  if (resu) printf("field DG test OK !\n");
  else printf("field DG test failed !\n");
  return !resu;
} 
