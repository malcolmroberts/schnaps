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
    
   int resu1 = Test_SH_periodic_HLL();
  int resu2 = Test_SH_periodic_Rusanov();
  int resu3 = Test_SH_periodic_Roe();
  //int resu4 = Test_SH_equilibrium_HLLWB();
	 
  if (resu1==1 && resu2==1 && resu3==1) printf("shallow water equilibrium  test OK !\n");
  else printf("shallow water equilibrium test failed !\n");

  return !resu1;
} 

int Test_SH_periodic_Rusanov(void) {

  int test = 1;
  real tmax=0.0,dt=0.0,tolerance=0.0,dd=0.0;

  field f;
  init_empty_field(&f);

  int vec=1;
  f.model.m=6; 
  f.model.NumFlux=ShallowWater_Rusanov_NumFlux;
  f.model.BoundaryFlux = ShallowWater_periodic_BoundaryFlux;

 
  f.model.InitData = TestSH_periodic_InitData;
  f.model.ImposedData = TestSH_periodic_ImposedData;
 

  f.varindex = GenericVarindex;
    
  f.interp.interp_param[0] = f.model.m;  // _M
  f.interp.interp_param[1] = 2;  // x direction degree
  f.interp.interp_param[2] = 2;  // y direction degree
  f.interp.interp_param[3] = 0;  // z direction degree
  f.interp.interp_param[4] = 16;  // x direction refinement
  f.interp.interp_param[5] = 16;  // y direction refinement
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
    test=1;     
  }
  printf("L2 Rusanov error %.8e\n", dd);

  // Save the results and the error
  Plotfield(0,false, &f, "h", "dgvisu_h.msh");
  Plotfield(1,false, &f, "u1", "dgvisu_u1.msh");
  Plotfield(2,false, &f, "u2", "dgvisu_u2.msh");
  Plotfield(3,false, &f, "b", "dgvisu_b.msh");

  return test;
}


int Test_SH_periodic_HLL(void) {

  int test1 = 1,test2 = 1,test3=1,test4=1,test = 1;
  real tmax=0.0,dt=0.0,tolerance=0.0,dd=0.0;

  field f;
  init_empty_field(&f);

  int vec=1;
  f.model.m=6;
  f.model.NumFlux=ShallowWater_HLL_NumFlux;
  f.model.BoundaryFlux = ShallowWater_periodic_BoundaryFlux;
 
  //f.model.Source = NULL;
 
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
  f.model.cfl = 0.1;
  
  /******* Test for hLL ******/
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
    test=1;     
  }
  printf("L2 HLL error %.8e\n", dd);

  // Save the results and the error
  Plotfield(0,false, &f, "h", "dgvisu_h.msh");
  Plotfield(1,false, &f, "u1", "dgvisu_u1.msh");
  Plotfield(2,false, &f, "u2", "dgvisu_u2.msh");
  Plotfield(3,false, &f, "b", "dgvisu_b.msh");

  return test;
}


int Test_SH_periodic_Roe(void) {
  int test = 1;
  real tmax=0.0,dt=0.0,tolerance=0.0,dd=0.0;

  field f;
  init_empty_field(&f);

  int vec=1;
  f.model.m=6; 
  f.model.NumFlux=ShallowWater_Roe_NumFlux;
  f.model.BoundaryFlux = ShallowWater_periodic_BoundaryFlux;
  //f.model.Source = NULL;
 
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
  f.model.cfl = 0.1;

  /******* Test for Roe ******/
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
    test=1;     
  }
  printf("L2 Roe error %.8e\n", dd);

  // Save the results and the error
  Plotfield(0,false, &f, "h", "dgvisu_h.msh");
  Plotfield(1,false, &f, "u1", "dgvisu_u1.msh");
  Plotfield(2,false, &f, "u2", "dgvisu_u2.msh");
  Plotfield(3,false, &f, "b", "dgvisu_b.msh");

  return test;
}



void TestSH_periodic_ImposedData(const real *x, const real t, real *w) {
  real u0=2;
  real v0=2;
  real pi=4.0*atan(1.0);
  
  w[0] = 1.0+0.2*sin(pi*(x[0]+x[1]-u0*t-v0*t));
  w[1] = u0;
  w[2] = v0;
  w[3] = -(1.0+0.2*sin(pi*(x[0]+x[1]-u0*t-v0*t)));
  w[4] = -0.2*pi*cos(pi*(x[0]+x[1]-u0*t-v0*t));
  w[5] = -0.2*pi*cos(pi*(x[0]+x[1]-u0*t-v0*t));
   

}

void TestSH_periodic_InitData(real *x, real *w) {
  real t = 0;
  TestSH_periodic_ImposedData(x, t, w);
}


void ShallowWater_periodic_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
				       real *flux) {
  real wR[3];

  
}



