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
  model.Source = NULL;

  int deg[]={4, 4, 4};
  int raf[]={4, 4, 4};
  
  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testmacromesh.msh");
  //ReadMacroMesh(&mesh,"../test/testcube2.msh");
  BuildConnectivity(&mesh);

  /* real A[3][3] = {{10,2 , 0}, {0, 1, -0.1}, {0, 0.1,1}}; */
  /* real x0[3] = {1, 2, 3}; */
  /* AffineMapMacroMesh(&mesh,A,x0); */

  CheckMacroMesh(&mesh, deg, raf);

  //PrintMacroMesh(&mesh);

  Simulation simu;

  InitSimulation(&simu, &mesh, deg, raf, &model);


  DtFields(&simu, simu.w, simu.dtw);
  
  DisplaySimulation(&simu);

  //PlotFields(0, false, &simu, NULL, "visu.msh");
  //PlotFields(0, true, &simu, "error", "error.msh");

  // Test the time derivative with the exact solution
  real test2 = 0;
  for(int i = 0; 
      i < model.m * mesh.nbelems * NPG(deg,raf); 
      i++){
    real errloc = fabs(4 * simu.w[i] - pow(simu.dtw[i], 2));
    test2 += errloc * errloc;
    test = test && errloc < 1e-2;
    //printf("i=%d err=%f \n",i,4 * w[i] - pow(dtw[i], 2));
    //assert(test);
  }

  printf("error=%f\n",sqrt(test2/ (mesh.nbelems * NPG(deg,raf)) ));
  return test;
};

int main(void) {
  // Unit tests
  int resu = TestfieldDG();
  if (resu) printf("field DG test OK !\n");
  else printf("field DG test failed !\n");
  return !resu;
} 
