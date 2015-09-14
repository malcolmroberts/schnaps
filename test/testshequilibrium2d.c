#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "test.h"
#include "collision.h"
#include "waterwave2d.h"

int main(void) {
  // unit tests
    
  int resu = Test_SH_equilibrium();
	 
  if (resu ==1) printf("shallow water equilibrium  test OK !\n");
  else printf("shallow water equilibrium test failed !\n");

  return !resu;
} 

int Test_SH_equilibrium(void) 
{
  int test1 = 0;
  int test2 = 0;
  int test3=0;
  int test = 0;
  real tolerance=0.0; 
  real dd=0.0;

  field f;
  init_empty_field(&f);

  ReadMacroMesh(&f.macromesh, "../test/testcube.msh");
  Detect2DMacroMesh(&f.macromesh);

  real A[3][3] = {{_LENGTH_DOMAIN, 0, 0}, {0, _LENGTH_DOMAIN, 0}, {0, 0,1}};
  real x0[3] = {0, 0, 0};
  AffineMapMacroMesh(&f.macromesh, A, x0);

  BuildConnectivity(&f.macromesh);

  Model model;

  int deg[]={2, 2, 0};
  int raf[]={8, 8, 1};

  assert(f.macromesh.is2d);

  CheckMacroMesh(&f.macromesh, deg, raf);

  f.deg[0] = deg[0];
  f.deg[1] = deg[1];
  f.deg[2] = deg[2];

  f.raf[0] = raf[0];
  f.raf[1] = raf[1];
  f.raf[2] = raf[2];

  f.model.cfl = 0.2;

  real tmax = 0.01;
  real dt = 1e-5; // FIXME: make variable
  f.vmax = _SPEED_WAVE;

  /******* Test for Rusanov ******/
  f.model.m=6; 
  f.model.NumFlux=ShallowWater_Rusanov_NumFlux;
  f.model.InitData = TestSH_equilibrium_InitData;
  f.model.ImposedData = TestSH_equilibrium_ImposedData;
  f.model.BoundaryFlux = ShallowWater_Rusanov_BoundaryFlux;
  f.model.Source = ShallowWater_classical_SourceTerm;

  RK2(&f, tmax, dt);

  dd = L2error(&f);
  tolerance = 1e-3;
  
  if(dd < tolerance){
    test1=1;     
  }
  printf("L2 Rusanov error %.8e\n", dd);
  
  /******* Test for hLL ******/
  f.model.m=6; 
  f.model.NumFlux=ShallowWater_HLL_NumFlux;
  f.model.InitData = TestSH_equilibrium_InitData;
  f.model.ImposedData = TestSH_equilibrium_ImposedData;
  f.model.BoundaryFlux = ShallowWater_HLL_BoundaryFlux;
  f.model.Source = ShallowWater_classical_SourceTerm;

  // FIXME: init!
  RK2(&f, tmax, dt);

  dd = L2error(&f);
  tolerance = 1e-3;
  
  if(dd < tolerance){
    test2=1;     
  }
  printf("L2 HLL error %.8e\n", dd);

  /******* Test for Roe ******/
  f.model.m=6; 
  f.model.NumFlux=ShallowWater_Roe_NumFlux;
  f.model.InitData = TestSH_equilibrium_InitData;
  f.model.ImposedData = TestSH_equilibrium_ImposedData;
  f.model.BoundaryFlux = ShallowWater_Roe_BoundaryFlux;
  f.model.Source = ShallowWater_classical_SourceTerm;

  // FIXME: init field
  RK2(&f, tmax, dt);

  dd = L2error(&f);
  tolerance = 1e-3;
  
  if(dd < tolerance){
    test3=1;     
  }
  printf("L2 HLL error %.8e\n", dd);

  /******* Test WB HLL ******/
  /* f.model.NumFlux=ShallowWater_HLLWB_NumFlux;
     f.model.BoundaryFlux = ShallowWater_HLLWB_BoundaryFlux;
     Initfield(&f);
     // maximal wave speed
     f.nb_diags = 0;
     f.pre_dtfield = NULL;
     f.update_after_rk = NULL;
     f.model.Source = ShallowWater_HLLWB_SourceTerm;
  
     tmax = 0.01;
     dt = set_dt(&f);
     RK2(&f, tmax, dt);

     dd = L2error(&f);
     tolerance = 1e-10;
  
     if(dd < tolerance){
     test4=1;     
     }
     printf("L2 HLL WB error %.8e\n", dd); */

  // Save the results and the error
  /* PlotFields(0,false, &simu3, "h", "dgvisu_h.msh"); */
  /* PlotFields(1,false, &simu3, "u1", "dgvisu_u1.msh"); */
  /* PlotFields(2,false, &simu3, "u2", "dgvisu_u2.msh"); */

  if(test1 +test2+test3 > 2){
    test=1;     
  }

  return test;
}


void TestSH_equilibrium_ImposedData(const real *x, const real t, real *w) 
{

  w[0] = 1.0-0.8*exp(-50*(pow(x[0]-0.5*_LENGTH_DOMAIN,2.0)+pow(x[1]-0.5*_LENGTH_DOMAIN,2.0)));
  w[1] = 0.0;
  w[2] = 0.0;
  w[3] = 0.8*exp(-50*(pow(x[0]-0.5*_LENGTH_DOMAIN,2.0)+pow(x[1]-0.5*_LENGTH_DOMAIN,2.0)));
  w[4] = -0.8*50*(2.0*x[0]-_LENGTH_DOMAIN)*exp(-50*(pow(x[0]-0.5*_LENGTH_DOMAIN,2.0)+pow(x[1]-0.5*_LENGTH_DOMAIN,2.0)));
  w[5] = -0.8*50*(2.0*x[1]-_LENGTH_DOMAIN)*exp(-50*(pow(x[0]-0.5*_LENGTH_DOMAIN,2.0)+pow(x[1]-0.5*_LENGTH_DOMAIN,2.0)));
  

}

void TestSH_equilibrium_InitData(real *x, real *w) 
{
  real t = 0;
  TestSH_equilibrium_ImposedData(x, t, w);
}


void ShallowWater_Rusanov_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
				       real *flux) 
{
  real wR[6];
  TestSH_equilibrium_ImposedData(x , t, wR);
  ShallowWater_Rusanov_NumFlux(wL, wR, vnorm, flux);
}

void ShallowWater_HLL_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
				   real *flux) 
{
  real wR[6];
  TestSH_equilibrium_ImposedData(x , t, wR);
  ShallowWater_HLL_NumFlux(wL, wR, vnorm, flux);
}


void ShallowWater_Roe_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
				   real *flux) 
{
  real wR[6];
  TestSH_equilibrium_ImposedData(x , t, wR);
  ShallowWater_Roe_NumFlux(wL, wR, vnorm, flux);
}


void ShallowWater_HLLWB_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
				     real *flux) 
{
  real wR[6];
  TestSH_equilibrium_ImposedData(x , t, wR);
  ShallowWater_HLLWB_NumFlux(wL, wR, vnorm, flux);
}
