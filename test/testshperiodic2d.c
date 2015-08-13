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
    
  int resu1 = Test_SH_periodic();
	 
   if (resu1==1) printf("shallow water periodic  test OK !\n");
   else printf("shallow water periodic test failed !\n");

  return !resu1;
} 


int Test_SH_periodic(void) {

  int test1 = 0,test2 = 0,test3=0,test = 0;
  real tmax=0.0,dt=0.0,tolerance=0.0,dd=0.0,dd2=0,dd3=0;

  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testcube.msh");
  Detect2DMacroMesh(&mesh);
  
  real A[3][3] = {{_LENGTH_DOMAIN, 0, 0}, {0, _LENGTH_DOMAIN, 0}, {0, 0,1}};
  real x0[3] = {0, 0, 0};
  AffineMapMacroMesh(&mesh,A,x0);

  BuildConnectivity(&mesh);

  Model model;

  int deg[]={2, 2, 0};
  int raf[]={8, 8, 1};


  assert(mesh.is2d);

  CheckMacroMesh(&mesh, deg, raf);
  

  /******* Test for Rusanov ******/
  model.m=6; 
  model.NumFlux=ShallowWater_Rusanov_NumFlux;
  model.InitData = TestSH_periodic_InitData;
  model.ImposedData = TestSH_periodic_ImposedData;
  model.BoundaryFlux = ShallowWater_Rusanov_BoundaryFlux;
  model.Source = ShallowWater_periodic_SourceTerm;

  Simulation simu;

  InitSimulation(&simu, &mesh, deg, raf, &model);
  
  tmax = 0.01;
  simu.vmax = 1;//_SPEED_WAVE;
  simu.cfl = 0.025;
  RK2(&simu, tmax);

  dd = L2error_onefield(&simu,0);
  dd2 = L2error_onefield(&simu,1);
  dd3 = L2error_onefield(&simu,2);
  tolerance = 5e-3;
  
  if(dd+dd2+dd3 < tolerance){
    test1=1;     
  }
  printf("L2 Rusanov error %.8e\n", dd+dd2+dd3);
  
  /******* Test for hLL ******/
  model.m=6; 
  model.NumFlux=ShallowWater_HLL_NumFlux;
  model.InitData = TestSH_periodic_InitData;
  model.ImposedData = TestSH_periodic_ImposedData;
  model.BoundaryFlux = ShallowWater_HLL_BoundaryFlux;
  model.Source = ShallowWater_periodic_SourceTerm;

  Simulation simu2;

  InitSimulation(&simu2, &mesh, deg, raf, &model);
  
  tmax = 0.01;
  simu2.vmax = 1;//_SPEED_WAVE;
  simu2.cfl = 0.025;
  RK2(&simu2, tmax);

  dd = L2error_onefield(&simu2,0);
  dd2 = L2error_onefield(&simu2,1);
  dd3 = L2error_onefield(&simu2,2);
  tolerance = 5e-3;
  
  if(dd+dd2+dd3 < tolerance){
    test2=1;     
  }
  printf("L2 HLL error %.8e\n", dd+dd2+dd3);

  /******* Test for Roe ******/
   model.m=6; 
  model.NumFlux=ShallowWater_Roe_NumFlux;
  model.InitData = TestSH_periodic_InitData;
  model.ImposedData = TestSH_periodic_ImposedData;
  model.BoundaryFlux = ShallowWater_Roe_BoundaryFlux;
  model.Source = ShallowWater_periodic_SourceTerm;

  Simulation simu3;

  InitSimulation(&simu3, &mesh, deg, raf, &model);
  
  tmax = 0.01;
  simu3.vmax = 1;//_SPEED_WAVE;
  simu3.cfl = 0.025;
  RK2(&simu3, tmax);

  dd = L2error_onefield(&simu3,0);
  dd2 = L2error_onefield(&simu3,1);
  dd3 = L2error_onefield(&simu3,2);
  tolerance = 5e-3;

  real dd4 = L2error_onefield(&simu2,3);
  real dd5 = L2error_onefield(&simu2,4);
  real dd6 = L2error_onefield(&simu2,5);

  printf("erreur h L2=%.12e\n", dd);
  printf("erreur u1 L2=%.12e\n", dd2);
  printf("erreur u2 L2=%.12e\n", dd3);
  printf("erreur b L2=%.12e\n", dd4);
  printf("erreur bx L2=%.12e\n", dd5);
  printf("erreur by L2=%.12e\n", dd6);
  
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
  PlotFields(0,false, &simu3, "h", "dgvisu_h.msh");
  PlotFields(1,false, &simu3, "u1", "dgvisu_u1.msh");
  PlotFields(2,false, &simu3, "u2", "dgvisu_u2.msh");

  if(test1 +test2+test3 > 2){
    test=1;     
  }

    FreeMacroMesh(&mesh);

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

void TestSH_periodic_InitData(real *x, real *w) {
  real t = 0;
  TestSH_periodic_ImposedData(x, t, w);
}



void ShallowWater_periodic_SourceTerm(const real *x, const real t, const real *w, real *source){
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
 

};

void ShallowWater_Rusanov_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
				       real *flux) {
  real wR[6];
  TestSH_periodic_ImposedData(x , t, wR);
  ShallowWater_Rusanov_NumFlux(wL, wR, vnorm, flux);
}

void ShallowWater_HLL_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
				       real *flux) {
  real wR[6];
  TestSH_periodic_ImposedData(x , t, wR);
  ShallowWater_HLL_NumFlux(wL, wR, vnorm, flux);
}


void ShallowWater_Roe_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
				       real *flux) {
  real wR[6];
  TestSH_periodic_ImposedData(x , t, wR);
  ShallowWater_Roe_NumFlux(wL, wR, vnorm, flux);
}


void ShallowWater_HLLWB_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
				       real *flux) {
  real wR[6];
  TestSH_periodic_ImposedData(x , t, wR);
  ShallowWater_HLLWB_NumFlux(wL, wR, vnorm, flux);
}
