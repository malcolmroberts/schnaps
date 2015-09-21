#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "test.h"
#include "collision.h"
#include "waterwave2d.h"
#include "implicit.h"

const  real a[6]={1,1,1,1,1,1};
const  real b[6]={1,1,1,1,1,1};
const  real c[6]={1,1,1,1,1,1};
  
void TestSteady_Transport_NumFlux(real *wL, real *wR, real *vnorm, real *flux);
void TestSteady_Transport_ImposedData(const real *x, const real t, real *w);
void TestSteady_Transport_InitData(real *x, real *w);
void TestSteady_Transport_Source(const real *xy, const real t, const real *w, real *S);
void Transport_Upwind_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
                              real *flux);


void AssemblyImplicitLinearSolver(Simulation *simu, LinearSolver *solver,real theta,real dt);




int main(void) {
  
  // unit tests
    
  int resu = Test_Transport_Steady();
	 
  if (resu) printf("transport steady  test OK !\n");
  else printf("transport steady test failed !\n");

  return !resu;
} 

int Test_Transport_Steady(void) {

  bool test = true;

  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testcube.msh");
  //ReadMacroMesh(&mesh,"../test/testdisque2d.msh");
  //ReadMacroMesh(&mesh,"../test/testmacromesh.msh");
  Detect2DMacroMesh(&mesh);
  
  real A[3][3] = {{_LENGTH_DOMAIN, 0, 0}, {0, _LENGTH_DOMAIN, 0}, {0, 0,1}};
  real x0[3] = {0, 0, 0};
  AffineMapMacroMesh(&mesh,A,x0);

  BuildConnectivity(&mesh);

  Model model;

  model.m=1;
  model.NumFlux = TestSteady_Transport_NumFlux;
  model.BoundaryFlux = Transport_Upwind_BoundaryFlux;
  model.InitData = TestSteady_Transport_InitData;
  model.ImposedData = TestSteady_Transport_ImposedData;
  model.Source = TestSteady_Transport_Source;

  //int deg[]={4, 4, 0};
  //int raf[]={2, 2, 1};
  int deg[]={1, 1, 0};
  int raf[]={1, 1, 1};
  
  CheckMacroMesh(&mesh, deg, raf);
  Simulation simu;
  EmptySimulation(&simu);

  InitSimulation(&simu, &mesh, deg, raf, &model);

  real tmax = 1000.0;
  simu.cfl=0.2;
  simu.vmax=_SPEED_WAVE;
 
  real dd = 0;

  real tolerance = 3e-5;

  test = test && (dd < tolerance);

  ThetaTimeScheme(&simu, tmax, 1);

  
  
  dd = L2error(&simu);

  printf("erreur implicit L2=%.12e\n", dd);

  test = test && (dd < tolerance);
   
  return test;
}

void TestSteady_Transport_ImposedData(const real *xy, const real t, real *w) {

  real x=xy[0];
  real y=xy[1];

  w[0] = x * (1 - x) * y * (1-y) + 1;
}

void TestSteady_Transport_Source(const real *xy, const real t, const real *w, real *S){

  real x=xy[0];
  real y=xy[1];

  const real v2[] = {sqrt(0.5), sqrt(0.5), 0};

  S[0] = v2[0] * (1 - 2 * x) * y * (1 - y) +
    v2[1] * (1 - 2 * y) * x * (1 - x);

}

void TestSteady_Transport_InitData(real *x, real *w) {
  real t = 0;
  TestSteady_Transport_ImposedData(x, t, w);
}


void Transport_Upwind_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
                                       real *flux) {
  real wR[3];
  TestSteady_Transport_ImposedData(x , t, wR);
  TestSteady_Transport_NumFlux(wL, wR, vnorm, flux);
}
 
void TestSteady_Transport_NumFlux(real *wL, real *wR, real *vnorm, real *flux)
{
  //const real transport_v2d[] = {sqrt(0.5), sqrt(0.5), 0};
  const real transport_v2d[] = {1,0, 0};
  real vn 
    = transport_v2d[0] * vnorm[0]
    + transport_v2d[1] * vnorm[1]
    + transport_v2d[2] * vnorm[2];
  real vnp = vn > 0 ? vn : 0;
  real vnm = vn - vnp;
  flux[0] = vnp * wL[0] + vnm * wR[0];
  /* if (fabs(vnorm[2])>1e-6) { */
  /*   printf("vnds %lf %lf %lf \n", vnorm[0], vnorm[1], vnorm[2]); */
  /* } */
  // verify that 2d computations are actually
  // activated
  //assert(fabs(vnorm[2]) < 1e-8);
}
