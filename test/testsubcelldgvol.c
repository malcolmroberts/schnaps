#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <assert.h>
#include <math.h>

int TestfieldSubCellDGVol()
{

  int test = true;

  Model model;
  model.cfl = 0.05;
  model.m = 1; // only one conservative variable
   
  field f;
  init_empty_field(&f);

  model.NumFlux = TransNumFlux;
  model.BoundaryFlux = TestTransBoundaryFlux;
  model.InitData = TestTransInitData;
  model.ImposedData = TestTransImposedData;
  model.Source = NULL;
  f.varindex = GenericVarindex;

  int deg[] = {2, 2, 2};
  int raf[] = {2, 2, 1};
  
  MacroMesh mesh;
  ReadMacroMesh(&mesh, "../test/testcube.msh");
  BuildConnectivity(&mesh);
  CheckMacroMesh(&mesh, deg, raf);
  
  Simulation simu;
  InitSimulation(&simu, &mesh, deg, raf, &model);

  // number of points in a macrocell.
  int mcell_size =  simu.wsize / simu.macromesh.nbelems;
  
  for(int ifa = 0; ifa < simu.macromesh.nbfaces; ifa++){
    //DGMacroCellInterface(ifa, &f, f.wn, f.dtwn);
  }
  for(int ie = 0; ie < simu.macromesh.nbelems; ie++) {
    /* DGSubCellInterface(ie, &f, f.wn, f.dtwn); */
    /* DGVolume(ie, &f, f.wn, f.dtwn); */
    /* DGMass(ie, &f, f.dtwn); */
    /* DGSource(ie, &f, f.wn, f.dtwn); */
  }
  
  /* Displayfield(&f); */

  /* Plotfield(0, false, &f, NULL, "visu.msh"); */
  /* Plotfield(0, true, &f, "error", "error.msh"); */
  
  // test the time derivative with the exact solution
  real maxerr = 0.0;
  int stop = model.m * mesh.nbelems * NPG(deg,raf);
  for(int i = 0; i < stop; i++) {
    real errloc = fabs(4 * simu.w[i] - pow(simu.dtw[i], 2));
    test = test && errloc < 1e-2;
    //printf("i=%d err=%f \n",i,4 * w[i] - pow(dtw[i], 2));
    //assert(test);
  }
  
  /* int stop = model.m * simu.macromesh.nbelems * NPG(deg, raf); */
  /* for(int i=0; i < stop; i++) { */
  /*   real err = fabs(4 * f.wn[i] - pow(f.dtwn[i] , 2)); */
    
  /*   printf("i: %d\terr: %f\tw: %f\tdtw: %f \n",i, err, ); */
  /*   if(err > maxerr) */
  /*     maxerr = err; */
  /* } */
  test = maxerr < 1e-2;

  printf("\nmaxerr: %f\n", maxerr);
  return test;
}

int main()
{
  int resu = TestfieldSubCellDGVol();
  if(resu) 
    printf("field DG Subcell Vol test OK !\n");
  else 
    printf("field DG Subcell Vol test failed !\n");
  return !resu;
} 
