#include "implicit.h"
#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "test.h"
#include "collision.h"
#include "solvercontinuous.h"
#include "waterwave2d.h"
#include "continuouspc.h"



void TestExtraction_ImposedData(const real *x, const real t, real *w);
void TestExtraction_InitData(real *x, real *w);

int main(void) {
  
  // unit tests
    
  int resu = Test_Extraction();
	 
  if (resu) printf("wave periodic  test OK !\n");
  else printf("wave periodic test failed !\n");

  return !resu;
} 

int Test_Extraction(void) {

  bool test = true;
  real dd;

  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testcube.msh");
  Detect2DMacroMesh(&mesh);
  
  real A[3][3] = {{_LENGTH_DOMAIN, 0, 0}, {0, _LENGTH_DOMAIN, 0}, {0, 0,1}};
  real x0[3] = {0, 0, 0};
  AffineMapMacroMesh(&mesh,A,x0);

  BuildConnectivity(&mesh);

  Model model;

  model.m = 3; 
  model.NumFlux=Wave_Upwind_NumFlux;
  model.InitData = TestExtraction_InitData;
  model.ImposedData = TestExtraction_ImposedData;
  model.BoundaryFlux = Wave_Upwind_BoundaryFlux;
  model.Source = NULL;

  int deg[]={2, 2, 0};
  int raf[]={4, 4, 1};


  assert(mesh.is2d);

  CheckMacroMesh(&mesh, deg, raf);
  Simulation simu;

  InitSimulation(&simu, &mesh, deg, raf, &model);
  
  // Apply division by the mass matrix
  field *f = &simu.fd[0];
  int m = f->model.m;
  
  for(int ipg = 0; ipg < NPG(f->deg, f->raf); ipg++) {
    real dtau[3][3], codtau[3][3], xpgref[3], xphy[3], wpg;
    ref_pg_vol(f->deg, f->raf, ipg, xpgref, &wpg, NULL);
    Ref2Phy(f->physnode, // phys. nodes
      xpgref, // xref
      NULL, -1, // dpsiref, ifa
      xphy, dtau, // xphy, dtau
      codtau, NULL, NULL); // codtau, dpsi, vnds
    real det = dot_product(dtau[0], codtau[0]);
    for(int iv = 0; iv < f->model.m; iv++) {
      int imem = f->varindex(f->deg, f->raf, f->model.m, ipg, iv);
      simu.fd[0].wn[imem] *= (wpg * det);
    }
  }
  
  physicPC_wave(&simu,simu.fd[0].dtwn,simu.fd[0].wn);
 
  dd = L2error(&simu);

  //printf("erreur L2=%.12e\n", dd);

 // PlotFields(0,false, &simu, "p", "dgvisu_exp.msh");
 // PlotFields(1,false, &simu, "u", "dgvisu_exu.msh");
 // PlotFields(2,false, &simu, "v", "dgvisu_exv.msh");

  return test;
}



void TestExtraction_ImposedData(const real *x, const real t, real *w) {

  w[0] = 5;
  w[1] = 2; 
  w[2] = 3;
  

}

void TestExtraction_InitData(real *x, real *w) {
  real t = 0;
  TestExtraction_ImposedData(x, t, w);
}


void Wave_Upwind_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
				       real *flux) {
  real wR[3];
  TestExtraction_ImposedData(x , t, wR);
  Wave_Upwind_NumFlux(wL, wR, vnorm, flux);
}


