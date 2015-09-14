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
void TestSteady_Transport_Source(const real *xy, const real t, const real *w, 
				 real *S);
void Transport_Upwind_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
				   real *flux);

int main() 
{
  // unit tests
    
  int resu = Test_Transport_Steady();
  
  if (resu) printf("transport steady  test OK !\n");
  else printf("transport steady test failed !\n");

  return !resu;
} 

int Test_Transport_Steady()
{

  bool test = true;

  field f;
  init_empty_field(&f);

  ReadMacroMesh(&f.macromesh, "../test/testdisque2d.msh");
  //ReadMacroMesh(&mesh,"../test/testmacromesh.msh");
  Detect2DMacroMesh(&f.macromesh);
  
  real A[3][3] = {{_LENGTH_DOMAIN, 0, 0}, {0, _LENGTH_DOMAIN, 0}, {0, 0,1}};
  real x0[3] = {0, 0, 0};
  AffineMapMacroMesh(&f.macromesh,A,x0);
  BuildConnectivity(&f.macromesh);

  f.model.m=1;
  f.model.NumFlux = TransNumFlux2d;
  f.model.BoundaryFlux = Transport_Upwind_BoundaryFlux;
  f.model.InitData = TestSteady_Transport_InitData;
  f.model.ImposedData = TestSteady_Transport_ImposedData;
  f.model.Source = TestSteady_Transport_Source;

  f.deg[0] = 4;
  f.deg[1] = 4;
  f.deg[2] = 0;

  f.raf[0] = 2;
  f.raf[1] = 2;
  f.raf[2] = 1;

  f.model.cfl = 0.2;
  f.vmax=_SPEED_WAVE;

  CheckMacroMesh(&f.macromesh, f.deg, f.raf);
  
  // FIXME: init

  real tmax = 1000.0;
  real dd = 0;

  real tolerance = 3e-5;

  test = test && (dd < tolerance);

  ThetaTimeScheme(&f, tmax, 10);
  
  dd = L2error(&f);

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
  TransNumFlux2d(wL, wR, vnorm, flux);
}
 
