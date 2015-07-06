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

int Test_SH_equilibrium(void) {

  int test1 = 1,test2 = 1,test3=1,test4=1,test = 1;
  real tmax=0.0,dt=0.0,tolerance=0.0,dd=0.0;

  field f;
  init_empty_field(&f);

  int vec=1;
  f.model.m=6; 
  f.model.NumFlux=ShallowWater_Rusanov_NumFlux;
 
  //f.model.Source = NULL;
 
  f.model.InitData = TestSH_equilibrium_InitData;
  f.model.ImposedData = TestSH_equilibrium_ImposedData;
  f.model.BoundaryFlux = ShallowWater_Rusanov_BoundaryFlux;

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

  // prepare the initial fields
  f.vmax = _SPEED_WAVE;
  f.model.cfl = 0.1;


  /******* Test for Rusanov ******/
  Initfield(&f);
   // maximal wave speed
  f.nb_diags = 0;
  f.pre_dtfield = NULL;
  f.update_after_rk = NULL;
  f.model.Source = ShallowWater_classical_SourceTerm;
  
  tmax = 0.01;
  dt = set_dt(&f);
  RK2(&f, tmax, dt);

  dd = L2error(&f);
  tolerance = 1e-3;
  
  if(dd < tolerance){
    test1=1;     
  }
  printf("L2 Rusanov error %.8e\n", dd);
  
  /******* Test for hLL ******/
  f.model.NumFlux=ShallowWater_HLL_NumFlux;
  f.model.BoundaryFlux = ShallowWater_HLL_BoundaryFlux;
  Initfield(&f);
   // maximal wave speed
  f.nb_diags = 0;
  f.pre_dtfield = NULL;
  f.update_after_rk = NULL;
  f.model.Source = ShallowWater_classical_SourceTerm;
  
  tmax = 0.01;
  dt = set_dt(&f);
  RK2(&f, tmax, dt);

  dd = L2error(&f);
  tolerance = 1e-3;
  
  if(dd < tolerance){
    test2=1;     
  }
  printf("L2 HLL error %.8e\n", dd);

  /******* Test for Roe ******/
   f.model.NumFlux=ShallowWater_Roe_NumFlux;
  f.model.BoundaryFlux = ShallowWater_Roe_BoundaryFlux;
  Initfield(&f);
   // maximal wave speed
  f.nb_diags = 0;
  f.pre_dtfield = NULL;
  f.update_after_rk = NULL;
  f.model.Source = ShallowWater_classical_SourceTerm;
  
  tmax = 0.01;
  dt = set_dt(&f);
  RK2(&f, tmax, dt);

  dd = L2error(&f);
  tolerance = 1e-3;
  
  if(dd < tolerance){
    test3=1;     
  }
  printf("L2 HLL error %.8e\n", dd);

  // Save the results and the error
  Plotfield(0,false, &f, "h", "dgvisu_h.msh");
  Plotfield(1,false, &f, "u1", "dgvisu_u1.msh");
  Plotfield(2,false, &f, "u2", "dgvisu_u2.msh");

  /******* Test WB HLL ******/
   f.model.NumFlux=ShallowWater_HLLWB_NumFlux;
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
  printf("L2 HLL WB error %.8e\n", dd);

  // Save the results and the error
  Plotfield(0,false, &f, "h", "dgvisu_h.msh");
  Plotfield(1,false, &f, "u1", "dgvisu_u1.msh");
  Plotfield(2,false, &f, "u2", "dgvisu_u2.msh");


  if(test1 +test2+test3+test4 > 3){
    test=1;     
  }

  return test;
}



void TestSH_equilibrium_ImposedData(const real *x, const real t, real *w) {

  w[0] = 1.0-0.8*exp(-50*(pow(x[0]-0.5*_LENGTH_DOMAIN,2.0)+pow(x[1]-0.5*_LENGTH_DOMAIN,2.0)));
  w[1] = 0.0;
  w[2] = 0.0;
  w[3] = 0.8*exp(-50*(pow(x[0]-0.5*_LENGTH_DOMAIN,2.0)+pow(x[1]-0.5*_LENGTH_DOMAIN,2.0)));
  w[4] = -0.8*50*(2.0*x[0]-_LENGTH_DOMAIN)*exp(-50*(pow(x[0]-0.5*_LENGTH_DOMAIN,2.0)+pow(x[1]-0.5*_LENGTH_DOMAIN,2.0)));
  w[5] = -0.8*50*(2.0*x[1]-_LENGTH_DOMAIN)*exp(-50*(pow(x[0]-0.5*_LENGTH_DOMAIN,2.0)+pow(x[1]-0.5*_LENGTH_DOMAIN,2.0)));
  

}

void TestSH_equilibrium_InitData(real *x, real *w) {
  real t = 0;
  TestSH_equilibrium_ImposedData(x, t, w);
}


void ShallowWater_Rusanov_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
				       real *flux) {
  real wR[3];
  TestSH_equilibrium_ImposedData(x , t, wR);
  ShallowWater_Rusanov_NumFlux(wL, wR, vnorm, flux);
}

void ShallowWater_HLL_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
				       real *flux) {
  real wR[3];
  TestSH_equilibrium_ImposedData(x , t, wR);
  ShallowWater_HLL_NumFlux(wL, wR, vnorm, flux);
}


void ShallowWater_Roe_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
				       real *flux) {
  real wR[3];
  TestSH_equilibrium_ImposedData(x , t, wR);
  ShallowWater_Roe_NumFlux(wL, wR, vnorm, flux);
}


void ShallowWater_HLLWB_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
				       real *flux) {
  real wR[3];
  TestSH_equilibrium_ImposedData(x , t, wR);
  ShallowWater_HLLWB_NumFlux(wL, wR, vnorm, flux);
}
