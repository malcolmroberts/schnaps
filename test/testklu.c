#include "implicit.h"
#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "test.h"
#include "collision.h"
#include "waterwave2d.h"
#include "suitesparse/klu.h"

schnaps_real vit[3] = {-1.4, -0.7, 0};

void KLU_Upwind_NumFlux(schnaps_real wL[],schnaps_real wR[],schnaps_real* vnorm,schnaps_real* flux);
  void TestSteady_KLU_ImposedData(const schnaps_real *x, const schnaps_real t, schnaps_real *w);
void TestSteady_KLU_InitData(schnaps_real *x, schnaps_real *w);
void TestSteady_KLU_Source(const schnaps_real *xy, const schnaps_real t, const schnaps_real *w, schnaps_real *S);
void KLU_Steady_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,
			      schnaps_real *flux);

int Test_KLU_Steady(void);

int main(void) {

  int n = 5 ;
  int Ap [ ] = {0, 2, 5, 9, 10, 12} ;
  int Ai [ ] = { 0, 1, 0, 2, 4, 1, 2, 3, 4, 2, 1, 4} ;
  double Ax [ ] = {2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.} ;
  double b [ ] = {8., 45., -3., 3., 19.} ;
  klu_symbolic *Symbolic ;
  klu_numeric *Numeric ;
  klu_common Common ;
  int i ;
  klu_defaults (&Common) ;
  Symbolic = klu_analyze (n, Ap, Ai, &Common) ;
  Numeric = klu_factor (Ap, Ai, Ax, Symbolic, &Common) ;
  klu_solve (Symbolic, Numeric, 5, 1, b, &Common) ;
  klu_free_symbolic (&Symbolic, &Common) ;
  klu_free_numeric (&Numeric, &Common) ;
  for (i = 0 ; i < n ; i++) printf ("x [%d] = %g\n", i, b [i]) ;
  assert(1==2);
  // unit tests
    
  int resu1 = 0;
  resu1=Test_KLU_Steady();
	 
  if (resu1 >=1) printf("wave periodic  test OK !\n");
  else printf("wave periodic test failed !\n");

  return !resu1;
} 


int Test_KLU_Steady(void) {

  bool test = true;

  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testcube.msh");
  Detect2DMacroMesh(&mesh);
  
  schnaps_real A[3][3] = {{_LENGTH_DOMAIN, 0, 0}, {0, _LENGTH_DOMAIN, 0}, {0, 0,1}};
  schnaps_real x0[3] = {0, 0, 0};
  AffineMapMacroMesh(&mesh,A,x0);

  BuildConnectivity(&mesh);

  Model model;

  model.m=1; 
  model.NumFlux=KLU_Upwind_NumFlux;
  model.InitData = TestSteady_KLU_InitData; 
  model.ImposedData = TestSteady_KLU_ImposedData; 
  model.BoundaryFlux = KLU_Steady_BoundaryFlux; 
  model.Source = TestSteady_KLU_Source; 

  int deg[]={1, 1, 0};
  int raf[]={3, 3, 1};
  
  CheckMacroMesh(&mesh, deg, raf);

  schnaps_real tmax = 0.01;
  
  Simulation simu2;
  EmptySimulation(&simu2);
  InitSimulation(&simu2, &mesh, deg, raf, &model);
  simu2.dt = tmax / 10;
  ThetaTimeScheme(&simu2, tmax, simu2.dt);
  
  schnaps_real dd = L2error(&simu2);

  printf("erreur implicit L2=%.12e\n", dd);

  PlotFields(0,false, &simu2, "f0", "dgvisu_f0.msh");
  PlotFields(1,false, &simu2, "f1", "dgvisu_f1.msh");
  schnaps_real tolerance = _SMALL;
  test = test && (dd < tolerance);

  FreeMacroMesh(&mesh);
  
  return test;
}



void TestSteady_KLU_ImposedData(const schnaps_real *xy, const schnaps_real t, schnaps_real *w) {

  schnaps_real x=xy[0];
  schnaps_real y=xy[1];


  w[0] = x + y;
  //w[1] = x - y;


}


void TestSteady_KLU_Source(const schnaps_real *xy, const schnaps_real t, const schnaps_real *w, schnaps_real *S){
  
  schnaps_real x=xy[0];
  schnaps_real y=xy[1];
  schnaps_real z=xy[2];

  S[0] = vit[0] * 1 + vit[1] * 1 + vit[2] * 0;
  //S[1] = 0;

}

void TestSteady_KLU_InitData(schnaps_real *x, schnaps_real *w) {
  schnaps_real t = 0;
  TestSteady_KLU_ImposedData(x, t, w);
}


void KLU_Steady_BoundaryFlux(schnaps_real *x, schnaps_real t, schnaps_real *wL, schnaps_real *vnorm,
				       schnaps_real *flux) {
  schnaps_real wR[3];
  TestSteady_KLU_ImposedData(x , t, wR);
  KLU_Upwind_NumFlux(wL, wR, vnorm, flux);
}

void KLU_Upwind_NumFlux(schnaps_real wL[],schnaps_real wR[],schnaps_real* vnorm,schnaps_real* flux){
  schnaps_real flux_temp=0;

  schnaps_real vn = vit[0] * vnorm[0] + vit[1] * vnorm[1] + vit[2] * vnorm[2];
  schnaps_real vnp = vn > 0.0 ? vn : 0.0;
  schnaps_real vnm = vn - vnp;
  flux[0] = vnp * wL[0] + vnm * wR[0];
  //flux[1] = 0;
  
};
