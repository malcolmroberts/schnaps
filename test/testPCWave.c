#include "implicit.h"
#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "test.h"
#include "solvercontinuous.h"
#include "waterwave2d.h"
#include "physBased_PC.h"


void SteadyStateOne_ImposedData(const real *x, const real t, real *w);
void SteadyStateOne_InitData(real *x, real *w);
void SteadyStateOne_Source(const real *xy, const real t, const real *w, real *S);
void SteadyStateOne_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
				                                                        real *flux);
void SteadyStateTwo_ImposedData(const real *x, const real t, real *w);
void SteadyStateTwo_InitData(real *x, real *w);
void SteadyStateTwo_Source(const real *xy, const real t, const real *w, real *S);
void SteadyStateTwo_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
				                                                        real *flux);

int main(void) {
  
  // unit tests
    
  int resu = TestPCWave();
	 
  if (resu) printf("wave periodic  test OK !\n");
  else printf("wave periodic test failed !\n");

  return !resu;
} 

int TestPCWave(void) {

  bool test = true;
  real dd;


#ifdef PARALUTION 
  paralution_begin();
#endif 

  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testcube.msh");
  Detect2DMacroMesh(&mesh);
  
  real A[3][3] = {{_LENGTH_DOMAIN, 0, 0}, {0, _LENGTH_DOMAIN, 0}, {0, 0,1}};
  real x0[3] = {0, 0, 0};
  AffineMapMacroMesh(&mesh,A,x0);

  BuildConnectivity(&mesh);

  /////////////////////////// Test One: No splitting error
  Model model;

  model.m = 3; 
  model.NumFlux=Wave_Upwind_NumFlux;
  model.InitData = SteadyStateOne_InitData;
  model.ImposedData = SteadyStateOne_ImposedData;
  model.BoundaryFlux = SteadyStateOne_BoundaryFlux;
  model.Source = SteadyStateOne_Source;

  int deg[]={3, 3, 0};
  int raf[]={2, 2, 1};


  assert(mesh.is2d);

  CheckMacroMesh(&mesh, deg, raf);
  Simulation simu;

  InitSimulation(&simu, &mesh, deg, raf, &model);


  LinearSolver solver_implicit;
  LinearSolver solver_explicit;  

  real theta=0.5;
  simu.theta=theta;
  simu.dt=0.1;
  simu.vmax=_SPEED_WAVE;
  real tmax=0.5;
  
  int itermax=tmax/simu.dt;
  simu.itermax_rk=itermax;
  InitImplicitLinearSolver(&simu, &solver_implicit);
  InitImplicitLinearSolver(&simu, &solver_explicit);
  real *res = calloc(simu.wsize, sizeof(real));

  PB_PC pb_pc;
  int mat2assemble[6] = {1, 1, 1, 1, 1, 1};
  InitPhy_Wave(&simu, &pb_pc, mat2assemble);

  simu.tnow=0;
  for(int ie=0; ie < simu.macromesh.nbelems; ++ie){
    simu.fd[ie].tnow=simu.tnow;
  } 

  for(int tstep=0;tstep<simu.itermax_rk;tstep++){
  

    if(tstep==0){ 
      solver_implicit.mat_is_assembly=false;
      solver_explicit.mat_is_assembly=false;
    } 
    else 
      { 
      solver_implicit.mat_is_assembly=true;
      solver_explicit.mat_is_assembly=true;
    } 

    solver_implicit.rhs_is_assembly=false;
    solver_explicit.rhs_is_assembly=false;

    
    AssemblyImplicitLinearSolver(&simu, &solver_explicit,-(1.0-theta),simu.dt);
    simu.tnow=simu.tnow+simu.dt;
    for(int ie=0; ie < simu.macromesh.nbelems; ++ie){
      simu.fd[ie].tnow=simu.tnow;
    } 
    AssemblyImplicitLinearSolver(&simu, &solver_implicit,theta,simu.dt);
  

    MatVect(&solver_explicit, simu.w, res);

    for(int i=0;i<solver_implicit.neq;i++){
      solver_implicit.rhs[i]=solver_implicit.rhs[i]+res[i]-solver_explicit.rhs[i];
    }
    solvePhy(&pb_pc,&simu,simu.w,solver_implicit.rhs);

    
    int freq = (1 >= simu.itermax_rk / 10)? 1 : simu.itermax_rk / 10;
    if (tstep % freq == 0)
      printf("t=%f iter=%d/%d dt=%f\n", simu.tnow, tstep, simu.itermax_rk, simu.dt);
  }
  dd = L2error(&simu);

  printf("erreur L2=%.12e\n", dd);

  test = test && (dd<2.e-2);
  freePB_PC(&pb_pc);
  freeSimulation(&simu);

  ///////////////////////////////// Test TWO: Splitting error
  printf("//////////////////////////////////////\n");
  printf("TEST TWO!\n");
  printf("//////////////////////////////////////\n");

  Model model2;

  model2.m = 3; 
  model2.NumFlux=Wave_Upwind_NumFlux;
  model2.InitData = SteadyStateTwo_InitData;
  model2.ImposedData = SteadyStateTwo_ImposedData;
  model2.BoundaryFlux = SteadyStateTwo_BoundaryFlux;
  model2.Source = SteadyStateTwo_Source;

  int deg2[]={4, 4, 0};
  int raf2[]={2, 2, 1};


  CheckMacroMesh(&mesh, deg2, raf2);
  Simulation simu2;

  InitSimulation(&simu2, &mesh, deg2, raf2, &model2);


  LinearSolver solver_implicit2;
  LinearSolver solver_explicit2;  

  real theta2=0.5;
  simu2.theta=theta2;
  simu2.dt=0.02;
  simu2.vmax=_SPEED_WAVE;
  real tmax2 = 0.08;
  
  real itermax2=tmax2/simu2.dt;
  simu2.itermax_rk=itermax2;
  InitImplicitLinearSolver(&simu2, &solver_implicit2);
  InitImplicitLinearSolver(&simu2, &solver_explicit2);
  real *res2 = calloc(simu2.wsize, sizeof(real));

  PB_PC pb_pc2;
  int mat2assemble2[6] = {1, 1, 1, 1, 1, 1};
  InitPhy_Wave(&simu2, &pb_pc2, mat2assemble2);

  simu2.tnow=0;
  for(int ie=0; ie < simu2.macromesh.nbelems; ++ie){
    simu2.fd[ie].tnow=simu2.tnow;
  } 

  for(int tstep=0;tstep<simu2.itermax_rk;tstep++){
  

    if(tstep==0){ 
      solver_implicit2.mat_is_assembly=false;
      solver_explicit2.mat_is_assembly=false;
    } 
    else 
      { 
      solver_implicit2.mat_is_assembly=true;
      solver_explicit2.mat_is_assembly=true;
    } 

    solver_implicit2.rhs_is_assembly=false;
    solver_explicit2.rhs_is_assembly=false;

    
    AssemblyImplicitLinearSolver(&simu2, &solver_explicit2,-(1.0-theta2),simu2.dt);
    simu2.tnow=simu2.tnow+simu2.dt;
    for(int ie=0; ie < simu2.macromesh.nbelems; ++ie){
      simu2.fd[ie].tnow=simu2.tnow;
    } 
    AssemblyImplicitLinearSolver(&simu2, &solver_implicit2,theta2,simu2.dt);
  

    MatVect(&solver_explicit2, simu2.w, res2);

    for(int i=0;i<solver_implicit2.neq;i++){
      solver_implicit2.rhs[i]=solver_implicit2.rhs[i]+res2[i]-solver_explicit2.rhs[i];
    }
    solvePhy(&pb_pc2,&simu2,simu2.fd[0].wn,solver_implicit2.rhs);

    
    int freq = (1 >= simu2.itermax_rk / 10)? 1 : simu2.itermax_rk / 10;
    if (tstep % freq == 0)
      printf("t=%f iter=%d/%d dt=%f\n", simu2.tnow, tstep, simu2.itermax_rk, simu2.dt);
  }
  dd = L2error(&simu2);

  printf("erreur L2=%.12e\n", dd);

  test = test && (dd<1.e-2);
  freePB_PC(&pb_pc2);
  freeSimulation(&simu2);

  ///////////////////////////////// Test THREE: Time-dependent problem
  printf("//////////////////////////////////////\n");
  printf("TEST THREE!\n");
  printf("//////////////////////////////////////\n");

  Model model3;

  model3.m = 3; 
  model3.NumFlux=Wave_Upwind_NumFlux;
  model3.InitData = TestPeriodic_Wave_InitData;
  model3.ImposedData = TestPeriodic_Wave_ImposedData;
  model3.BoundaryFlux = Wave_Upwind_BoundaryFlux;
  model3.Source = NULL;

  int deg3[]={4, 4, 0};
  int raf3[]={2, 2, 1};


  CheckMacroMesh(&mesh, deg3, raf3);
  Simulation simu3;

  InitSimulation(&simu3, &mesh, deg3, raf3, &model3);


  LinearSolver solver_implicit3;
  LinearSolver solver_explicit3;  

  real theta3=0.5;
  simu3.theta=theta3;
  simu3.dt=0.001667;//0.005;
  simu3.vmax=_SPEED_WAVE;
  real tmax3 = 0.1;
  
  real itermax3=tmax3/simu3.dt;
  simu3.itermax_rk=itermax3;
  InitImplicitLinearSolver(&simu3, &solver_implicit3);
  InitImplicitLinearSolver(&simu3, &solver_explicit3);
  real *res3 = calloc(simu3.wsize, sizeof(real));

  PB_PC pb_pc3;
  int mat2assemble3[6] = {1, 1, 1, 1, 1, 1};
  InitPhy_Wave(&simu3, &pb_pc3, mat2assemble3);

  simu3.tnow=0;
  for(int ie=0; ie < simu3.macromesh.nbelems; ++ie){
    simu3.fd[ie].tnow=simu3.tnow;
  } 

  for(int tstep=0;tstep<simu3.itermax_rk;tstep++){
  

    if(tstep==0){ 
      solver_implicit3.mat_is_assembly=false;
      solver_explicit3.mat_is_assembly=false;
    } 
    else 
    { 
      solver_implicit3.mat_is_assembly=true;
      solver_explicit3.mat_is_assembly=true;
    }

    solver_implicit3.rhs_is_assembly=false;
    solver_explicit3.rhs_is_assembly=false;

    
    AssemblyImplicitLinearSolver(&simu3, &solver_explicit3,-(1.0-theta3),simu3.dt);
    simu3.tnow=simu3.tnow+simu3.dt;
    for(int ie=0; ie < simu3.macromesh.nbelems; ++ie){
      simu3.fd[ie].tnow=simu3.tnow;
    } 
    AssemblyImplicitLinearSolver(&simu3, &solver_implicit3,theta3,simu3.dt);
  

    MatVect(&solver_explicit3, simu3.w, res3);

    for(int i=0;i<solver_implicit3.neq;i++){
      solver_implicit3.rhs[i]=solver_implicit3.rhs[i]+res3[i]-solver_explicit3.rhs[i];
    }
    
    solvePhy(&pb_pc3,&simu3,simu3.fd[0].wn,solver_implicit3.rhs);

    //
    int freq = (1 >= simu3.itermax_rk / 10)? 1 : simu3.itermax_rk / 10;
    if (tstep % freq == 0)
      printf("t=%f iter=%d/%d dt=%f\n", simu3.tnow, tstep, simu3.itermax_rk, simu3.dt);
  }
  dd = L2error(&simu3);

  printf("erreur L2=%.13e\n", dd);

  test = test && (dd<1.e-2);

  freePB_PC(&pb_pc3);
  freeSimulation(&simu3);

#ifdef PARALUTION 
  paralution_end();
#endif 

  return test;
}



void SteadyStateOne_ImposedData(const real *xy, const real t, real *w) {

  real x=xy[0];
  real y=xy[1];

  w[0] = 10+x*x+y*y*y;
  w[1] = x*y-y+2;
  w[2] = x-y*y*0.5+5.6;
}

void SteadyStateTwo_ImposedData(const real *xy, const real t, real *w) {

  real x=xy[0];
  real y=xy[1];


  w[0] = x*(1-x)*y*(1-y)+1;
  w[1] = 2*x*(1-x)*y*(1-y);
  w[2] = 3*x*(1-x)*y*(1-y);

}

void SteadyStateOne_Source(const real *xy, const real t, const real *w, real *S){
  
  real x=xy[0];
  real y=xy[1];

  S[0] = 0;
  S[1] = 2*x;
  S[2] = 3*y*y;

  S[0] *= _SPEED_WAVE;
  S[1] *= _SPEED_WAVE;
  S[2] *= _SPEED_WAVE;

}

void SteadyStateTwo_Source(const real *xy, const real t, const real *w, real *S){
  
  real x=xy[0];
  real y=xy[1];

  S[0] = 2*(1-2*x)*(y*(1-y))+3*(1-2*y)*(x*(1-x));
  S[1] = (1-2*x)*(y*(1-y));
  S[2] = (1-2*y)*(x*(1-x));

  S[0] *= _SPEED_WAVE;
  S[1] *= _SPEED_WAVE;
  S[2] *= _SPEED_WAVE;

}

void SteadyStateOne_InitData(real *x, real *w) {
  real t = 0;
  SteadyStateOne_ImposedData(x, t, w);
}

void SteadyStateTwo_InitData(real *x, real *w) {
  real t = 0;
  SteadyStateTwo_ImposedData(x, t, w);
}

void SteadyStateOne_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
				       real *flux) {
  real wR[3];
  SteadyStateOne_ImposedData(x , t, wR);
  Wave_Upwind_NumFlux(wL, wR, vnorm, flux);
}

void SteadyStateTwo_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
				       real *flux) {
  real wR[3];
  SteadyStateTwo_ImposedData(x , t, wR);
  Wave_Upwind_NumFlux(wL, wR, vnorm, flux);
}

void TestPeriodic_Wave_ImposedData(const real *x, const real t, real *w) {
  real pi=4.0*atan(1.0);
  real L=_LENGTH_DOMAIN;
  real Coef=(2.0*pi)/L;
  real a=_SPEED_WAVE;

  w[0] = -a*Coef*sqrt(2.0)*sin(a*Coef*sqrt(2.0)*t)*cos(Coef*x[0])*cos(Coef*x[1]);
  w[1] = a*Coef*cos(Coef*a*sqrt(2.0)*t)*sin(Coef*x[0])*cos(Coef*x[1]);
  w[2] = a*Coef*cos(Coef*a*sqrt(2.0)*t)*cos(Coef*x[0])*sin(Coef*x[1]);
  

}

void TestPeriodic_Wave_InitData(real *x, real *w) {
  real t = 0;
  TestPeriodic_Wave_ImposedData(x, t, w);
}


void Wave_Upwind_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
				       real *flux) {
  real wR[3];
  TestPeriodic_Wave_ImposedData(x , t, wR);
  Wave_Upwind_NumFlux(wL, wR, vnorm, flux);
}
