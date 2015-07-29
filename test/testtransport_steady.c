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
  Detect2DMacroMesh(&mesh);
  
  real A[3][3] = {{_LENGTH_DOMAIN, 0, 0}, {0, _LENGTH_DOMAIN, 0}, {0, 0,1}};
  real x0[3] = {0, 0, 0};
  AffineMapMacroMesh(&mesh,A,x0);

  BuildConnectivity(&mesh);

  Model model;

   model.m=1;
  model.NumFlux = TransNumFlux2d;
  model.BoundaryFlux = Transport_Upwind_BoundaryFlux;
  model.InitData = TestSteady_Transport_InitData;
  model.ImposedData = TestSteady_Transport_ImposedData;
  model.Source = TestSteady_Transport_Source;

  int deg[]={4, 4, 0};
  int raf[]={2, 2, 1};
  
  CheckMacroMesh(&mesh, deg, raf);
  Simulation simu;

  InitSimulation(&simu, &mesh, deg, raf, &model);

  real tmax = .1;
  simu.cfl=0.2;
  simu.vmax=_SPEED_WAVE;
  RK4(&simu,tmax);
 
  real dd = 0;
  dd = L2error(&simu);

  printf("erreur explicit L2=%f\n", dd);

  real tolerance = 0.0000001;

  test = test && (dd < tolerance);

  LinearSolver solver_implicit;
  LinearSolver solver_explicit;  

  real theta=0.5;
  //real dt=simu.dt;

  real dt = 1;
  InitImplicitLinearSolver(&simu, &solver_implicit);
  InitImplicitLinearSolver(&simu, &solver_explicit);
  AssemblyImplicitLinearSolver(&simu, &solver_implicit,theta,dt);
  AssemblyImplicitLinearSolver(&simu, &solver_explicit,-(1.0-theta),dt);

  real *res = calloc(simu.wsize, sizeof(real));

  for(int tstep=0;tstep<10;tstep++){

    MatVect(&solver_explicit, simu.w, res);

    for(int i=0;i<solver_implicit.neq;i++){
      solver_implicit.rhs[i]=solver_explicit.rhs[i]+res[i];
    }
  
    SolveLinearSolver(&solver_implicit);

    for(int i=0;i<solver_implicit.neq;i++){
      simu.w[i]=solver_implicit.sol[i];
    }
   

 }
  
  dd = L2error(&simu);

  printf("erreur implicit L2=%f\n", dd);

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
  TransNumFlux2d(wL, wR, vnorm, flux);
}
 
