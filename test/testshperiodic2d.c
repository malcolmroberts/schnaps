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

  int test1 = 0,test2=0,test3=0,test=0;
  real tmax=0.0,dt=0.0,tolerance=0.0,dd=0.0;

  field f;
  init_empty_field(&f);

  int vec=1;
  f.model.m=6; 
  f.model.BoundaryFlux = NULL;

 
  f.model.InitData = TestSH_periodic_InitData;
  f.model.ImposedData = TestSH_periodic_ImposedData;
 

  f.varindex = GenericVarindex;
    
  f.interp.interp_param[0] = f.model.m;  // _M
  f.interp.interp_param[1] = 2;  // x direction degree
  f.interp.interp_param[2] = 2;  // y direction degree
  f.interp.interp_param[3] = 0;  // z direction degree
  f.interp.interp_param[4] = 8;  // x direction refinement
  f.interp.interp_param[5] = 8;  // y direction refinement
  f.interp.interp_param[6] = 1;  // z direction refinement
 // read the gmsh file

  ReadMacroMesh(&(f.macromesh), "test/testcube.msh");

  Detect2DMacroMesh(&(f.macromesh));
  assert(f.macromesh.is2d);
  real A[3][3] = {{_LENGTH_DOMAIN, 0, 0}, {0, _LENGTH_DOMAIN, 0}, {0, 0,1}};
  real x0[3] = {0, 0, 0};
  AffineMapMacroMesh(&(f.macromesh),A,x0);

  f.macromesh.period[0]=_LENGTH_DOMAIN;
  f.macromesh.period[1]=_LENGTH_DOMAIN;
  BuildConnectivity(&(f.macromesh));


  CheckMacroMesh(&(f.macromesh), f.interp.interp_param + 1);

  printf("cfl param =%f\n", f.hmin);

  f.vmax = _SPEED_WAVE;
  f.model.cfl = 0.2;

  /******* Test for Rusanov ******/
  Initfield(&f);
   // maximal wave speed
  f.nb_diags = 0;
  f.pre_dtfield = NULL;
  f.update_after_rk = NULL;
  f.model.Source = ShallowWater_periodic_SourceTerm;
  f.model.NumFlux=ShallowWater_Rusanov_NumFlux;
  
  tmax = 0.01;
  dt = set_dt(&f);
  RK2(&f, tmax, dt);

  double ddp = L2error_onefield(&f,0);
  double ddu1 = L2error_onefield(&f,1);
  double ddu2 = L2error_onefield(&f,2);
  printf("L2 Rusanov h, u ,v  error %.8e %.8e %.8e \n", ddp,ddu1,ddu2);
  tolerance = 5e-3;
  
  if(ddp+ddu1+ddu2 < tolerance){
    test1=1;     
  }

   /******* Test for Roe ******/
  Initfield(&f);
   // maximal wave speed
  f.nb_diags = 0;
  f.pre_dtfield = NULL;
  f.update_after_rk = NULL;
  f.model.Source = ShallowWater_periodic_SourceTerm;
   f.model.NumFlux=ShallowWater_Roe_NumFlux;
  
  tmax = 0.01;
  dt = set_dt(&f);
  RK2(&f, tmax, dt);
  
  ddp = L2error_onefield(&f,0);
  ddu1 = L2error_onefield(&f,1);
  ddu2 = L2error_onefield(&f,2);
  printf("L2 Roe h, u ,v  error %.8e %.8e %.8e \n", ddp,ddu1,ddu2);
  tolerance = 5e-3;
  
  if(ddp+ddu1+ddu2 < tolerance){
    test2=1;     
  }

   /******* Test for hLL ******/
  Initfield(&f);
   // maximal wave speed
  f.nb_diags = 0;
  f.pre_dtfield = NULL;
  f.update_after_rk = NULL;
  f.model.Source = ShallowWater_periodic_SourceTerm;
  f.model.NumFlux=ShallowWater_HLL_NumFlux;
  
  tmax = 0.01;
  dt = set_dt(&f);
  RK2(&f, tmax, dt);

  ddp = L2error_onefield(&f,0);
  ddu1 = L2error_onefield(&f,1);
  ddu2 = L2error_onefield(&f,2);
  printf("L2 HLL h, u ,v  error %.8e %.8e %.8e \n", ddp,ddu1,ddu2);
  tolerance = 5e-3;
  
  if(ddp+ddu1+ddu2 < tolerance){
    test3=1;     
  }

  // Save the results and the error
  Plotfield(0,false, &f, "h", "dgvisu_h.msh");
  Plotfield(1,false, &f, "u", "dgvisu_u.msh");
  Plotfield(2,false, &f, "v", "dgvisu_v.msh");

  if(test1+test2+test3== 3){
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


