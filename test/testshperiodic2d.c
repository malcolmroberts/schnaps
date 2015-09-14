#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "test.h"
#include "collision.h"
#include "waterwave2d.h"

int main(void) 
{
  // unit tests
    
  int resu1 = Test_SH_periodic();
	 
  if (resu1==1) printf("shallow water periodic  test OK !\n");
  else printf("shallow water periodic test failed !\n");

  return !resu1;
} 


int Test_SH_periodic() 
{
  int test1 = 0,test2 = 0,test3=0,test = 0;
  real tmax=0.01; 
  real dt=1e-4; // FIXME: make variable 
  real tolerance=0.0;
  real dd=0.0;
  real dd2=0;
  real dd3=0;

  field f;
  init_empty_field(&f);
  
  f.vmax = _SPEED_WAVE;
  f.model.cfl = 0.2;

  ReadMacroMesh(&f.macromesh, "../test/testcube.msh");
  Detect2DMacroMesh(&f.macromesh);
  assert(f.macromesh.is2d);

  real A[3][3] = {{_LENGTH_DOMAIN, 0, 0}, {0, _LENGTH_DOMAIN, 0}, {0, 0,1}};
  real x0[3] = {0, 0, 0};
  AffineMapMacroMesh(&f.macromesh,A,x0);
  BuildConnectivity(&f.macromesh);

  CheckMacroMesh(&f.macromesh, f.deg, f.raf);

  int deg[]={2, 2, 0};
  int raf[]={8, 8, 1};

  f.deg[0] = deg[0];
  f.deg[1] = deg[1];
  f.deg[2] = deg[2];

  f.raf[0] = raf[0];
  f.raf[1] = raf[1];
  f.raf[2] = raf[2];

  f.model.cfl = 0.2;
  f.vmax = _SPEED_WAVE;

  /******* Test for Rusanov ******/
  f.model.m=6; 
  f.model.NumFlux=ShallowWater_Rusanov_NumFlux;
  f.model.InitData = TestSH_periodic_InitData;
  f.model.ImposedData = TestSH_periodic_ImposedData;
  f.model.BoundaryFlux = ShallowWater_Rusanov_BoundaryFlux;
  f.model.Source = ShallowWater_periodic_SourceTerm;

  // FIXME: init

  RK2(&f, tmax, dt);

  dd = L2error_onefield(&f,0);
  dd2 = L2error_onefield(&f,1);
  dd3 = L2error_onefield(&f,2);
  tolerance = 5e-3;
  
  if(dd+dd2+dd3 < tolerance){
    test1=1;     
  }
  printf("L2 Rusanov error %.8e\n", dd+dd2+dd3);
  
  /******* Test for hLL ******/
  f.model.m=6; 
  f.model.NumFlux=ShallowWater_HLL_NumFlux;
  f.model.InitData = TestSH_periodic_InitData;
  f.model.ImposedData = TestSH_periodic_ImposedData;
  f.model.BoundaryFlux = ShallowWater_HLL_BoundaryFlux;
  f.model.Source = ShallowWater_periodic_SourceTerm;

  // FIXME: init
  
  RK2(&f, tmax, dt);

  dd = L2error_onefield(&f,0);
  dd2 = L2error_onefield(&f,1);
  dd3 = L2error_onefield(&f,2);
  tolerance = 5e-3;
  
  if(dd+dd2+dd3 < tolerance){
    test2=1;     
  }
  printf("L2 HLL error %.8e\n", dd+dd2+dd3);

  /******* Test for Roe ******/
  f.model.m=6; 
  f.model.NumFlux=ShallowWater_Roe_NumFlux;
  f.model.InitData = TestSH_periodic_InitData;
  f.model.ImposedData = TestSH_periodic_ImposedData;
  f.model.BoundaryFlux = ShallowWater_Roe_BoundaryFlux;
  f.model.Source = ShallowWater_periodic_SourceTerm;

  // FIXME: init
  RK2(&f, tmax, dt);

  dd = L2error_onefield(&f,0);
  dd2 = L2error_onefield(&f,1);
  dd3 = L2error_onefield(&f,2);
  tolerance = 5e-3;
  
  if(dd+dd2+dd3 < tolerance){
    test3=1;     
  }
  printf("L2 HLL error %.8e\n", dd+dd2+dd3);

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
  /* PlotFields(0,false, &f, "h", "dgvisu_h.msh"); */
  /* PlotFields(1,false, &f, "u1", "dgvisu_u1.msh"); */
  /* PlotFields(2,false, &f, "u2", "dgvisu_u2.msh"); */

  if(test1 +test2+test3 > 2){
    test=1;     
  }

  return test;
}

void TestSH_periodic_ImposedData(const real *x, const real t, real *w) {
  real u0=2;
  real v0=2;
  real pi=4.0*atan(1.0);
  
  w[0] = 1.0+0.2*sin(((2.*pi)/_LENGTH_DOMAIN)*(x[0]+x[1]-u0*t-v0*t));
  w[1] = u0*w[0];
  w[2] = v0*w[0];
  w[3] = 1.0-(1.0+0.2*sin(((2.*pi)/_LENGTH_DOMAIN)*(x[0]+x[1]-u0*t-v0*t)));
  w[4] = -0.2*((2.*pi)/_LENGTH_DOMAIN)*cos(((2.*pi)/_LENGTH_DOMAIN)*(x[0]+x[1]-u0*t-v0*t));
  w[5] = -0.2*((2.*pi)/_LENGTH_DOMAIN)*cos(((2.*pi)/_LENGTH_DOMAIN)*(x[0]+x[1]-u0*t-v0*t));
}

void TestSH_periodic_InitData(real *x, real *w) 
{
  real t = 0;
  TestSH_periodic_ImposedData(x, t, w);
}

void ShallowWater_periodic_SourceTerm(const real *x, const real t, 
				      const real *w, real *source)
{
  real g=_GRAVITY;
  real hL=0, hR=0, uL=0, uR=0, vL=0, vR=0;
  real S=0,b=0,bx=0,by=0;
  real wexact[6];
  
  hL = w[0];
  uL=w[1]/w[0];
  vL=w[2]/w[0];

  TestSH_periodic_ImposedData(x , t, wexact);

  source[0]= 0.0;
  source[1]= -g*hL*wexact[4];
  source[2]= -g*hL*wexact[5];
  source[3]= 0.0;
  source[4]= 0.0;
  source[5]= 0.0;
}

void ShallowWater_Rusanov_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
				       real *flux) 
{
  real wR[6];
  TestSH_periodic_ImposedData(x , t, wR);
  ShallowWater_Rusanov_NumFlux(wL, wR, vnorm, flux);
}

void ShallowWater_HLL_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
				   real *flux) 
{
  real wR[6];
  TestSH_periodic_ImposedData(x , t, wR);
  ShallowWater_HLL_NumFlux(wL, wR, vnorm, flux);
}


void ShallowWater_Roe_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
				   real *flux) 
{
  real wR[6];
  TestSH_periodic_ImposedData(x , t, wR);
  ShallowWater_Roe_NumFlux(wL, wR, vnorm, flux);
}


void ShallowWater_HLLWB_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
				     real *flux) 
{
  real wR[6];
  TestSH_periodic_ImposedData(x , t, wR);
  ShallowWater_HLLWB_NumFlux(wL, wR, vnorm, flux);
}
