#include "implicit.h"
#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "../test/test.h"
#include "solverwave.h"
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
void Generic(ContinuousSolver *cs);

int main(void) {

  // unit tests

  int resu = Testrealpc();

  if (resu) printf("wave periodic  test OK !\n");
  else printf("wave periodic test failed !\n");

  return !resu;
} 

int Testrealpc(void) {

  bool test = true;
  real dd;
  int test1_ok=1,test2_ok=1;

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
  
  printf("//////////////////////////////////////\n");
  printf("TEST ONE : Stationnary CG !\n");
  printf("//////////////////////////////////////\n");
  if(test1_ok==1){
    Model model;

    model.m = 3; 
    model.NumFlux=Wave_Upwind_NumFlux;
    model.InitData = SteadyStateOne_InitData;
    model.ImposedData = SteadyStateOne_ImposedData;
    model.BoundaryFlux = SteadyStateOne_BoundaryFlux;
    model.Source = SteadyStateOne_Source;

    int deg[]={4, 4, 0};
    int raf[]={8, 8, 1};

    assert(mesh.is2d);
    CheckMacroMesh(&mesh, deg, raf);
    Simulation simu;
    EmptySimulation(&simu);
    InitSimulation(&simu, &mesh, deg, raf, &model);

    ContinuousSolver cs;
    ContinuousSolver csSolve;
    int nbvar=3;
    int *listvar=calloc(nbvar,sizeof(int));
    listvar[0]=0;
    listvar[1]=1;
    listvar[2]=2;
    InitContinuousSolver(&cs,&simu,1,nbvar,listvar);
    InitContinuousSolver(&csSolve,&simu,1,nbvar,listvar);

    real theta=0.5;
    simu.theta=theta;
    simu.dt=10;
    simu.vmax=_SPEED_WAVE;
    real tmax=5*simu.dt;
    int itermax=tmax/simu.dt;
    simu.itermax_rk=itermax;
    int size = cs.nb_fe_dof;
    real *resCG = calloc(size, sizeof(real));
    real *wCG = calloc(size, sizeof(real));
  
    csSolve.lsol.solver_type=LU;//GMRES;
    csSolve.lsol.tol=1.e-9;
    csSolve.lsol.pc_type=NONE;//PHY_BASED_U2;//PHY_BASED_U1//JACOBI;
    csSolve.lsol.iter_max=2000;
    csSolve.lsol.restart_gmres=30;
    csSolve.lsol.is_CG=true;

    cs.bc_flux=Wave_BC_pressure_imposed;
    cs.bc_assembly=BoundaryConditionFriedrichsAssembly;
    csSolve.bc_flux=Wave_BC_pressure_imposed;
    csSolve.bc_assembly=BoundaryConditionFriedrichsAssembly;
    //csSolve.rhs_assembly=Source_Assembly;
    csSolve.type_bc=2;
    
    Wave_test(&cs,-(1.0-simu.theta),simu.dt);
    GenericOperator_Continuous(&cs);
    Wave_test(&csSolve,simu.theta,simu.dt);
    GenericOperator_Continuous(&csSolve);
  
    PiDgToCg(&cs,simu.w,wCG);

    simu.tnow=0;
    for(int ie=0; ie < simu.macromesh.nbelems; ++ie){
      simu.fd[ie].tnow=simu.tnow;
    } 

    for(int tstep=0;tstep<simu.itermax_rk;tstep++){

      for (int i=0; i<size; i++){
	cs.lsol.rhs[i]=0;
      }
      cs.bc_assembly(&cs);
      MatVect(&cs.lsol,wCG,resCG);
  
      simu.tnow=simu.tnow+simu.dt;
      for(int ie=0; ie < simu.macromesh.nbelems; ++ie){
        simu.fd[ie].tnow=simu.tnow;
      }

      for (int i=0; i<size; i++){
	csSolve.lsol.rhs[i]=0;
      }
      csSolve.bc_assembly(&csSolve);
      //csSolve.rhs_assembly(&csSolve);
      
      for (int i=0; i<size; i++){
	csSolve.lsol.rhs[i]=csSolve.lsol.rhs[i]+resCG[i]-cs.lsol.rhs[i];
      }
      
       Advanced_SolveLinearSolver(&csSolve.lsol,&simu);
      
      for (int i=0; i<size; i++){
        wCG[i] = csSolve.lsol.sol[i];
      }
       
      int freq = (1 >= simu.itermax_rk / 10)? 1 : simu.itermax_rk / 10;
      if (tstep % freq == 0)
        printf("t=%f iter=%d/%d dt=%f\n", simu.tnow, tstep, simu.itermax_rk, simu.dt);
    }
    PiInvertCgToDg(&cs,wCG,simu.w);
    dd = L2error(&simu);

    printf("erreur L2=%.12e\n", dd);

    test = test && (dd<2.e-2);
    freeSimulation(&simu);
  }
  
  printf("//////////////////////////////////////\n");
  printf("TEST TWO : Time dependant CG !\n");
  printf("//////////////////////////////////////\n");
  if(test2_ok==1){
    Model model2;

    model2.m = 3; 
    model2.NumFlux=Wave_Upwind_NumFlux;
    model2.InitData = TestPeriodic_Wave_InitData;
    model2.ImposedData = TestPeriodic_Wave_ImposedData;
    model2.BoundaryFlux = Wave_Upwind_BoundaryFlux;
    model2.Source = NULL;

    int deg2[]={4, 4, 0};
    int raf2[]={16, 16 , 1};
    assert(mesh.is2d);
    CheckMacroMesh(&mesh, deg2, raf2);
    
    Simulation simu2;
    EmptySimulation(&simu2);
    InitSimulation(&simu2, &mesh, deg2, raf2, &model2); 
    ContinuousSolver cs;
    ContinuousSolver csSolve;
    int nbvar=3;
    int *listvar=calloc(nbvar,sizeof(int));
    listvar[0]=0;
    listvar[1]=1;
    listvar[2]=2;
    InitContinuousSolver(&cs,&simu2,1,nbvar,listvar);
    InitContinuousSolver(&csSolve,&simu2,1,nbvar,listvar);

    simu2.theta=0.5;
    simu2.dt=10;
    simu2.vmax=_SPEED_WAVE;
    real tmax2=1*simu2.dt;//;0.5;
    int itermax2=tmax2/simu2.dt;
    simu2.itermax_rk=itermax2;
    int size = cs.nb_fe_dof;
    real *resCG = calloc(size, sizeof(real));
    real *wCG = calloc(size, sizeof(real));

    csSolve.lsol.solver_type=GMRES;
    csSolve.lsol.tol=1.0e-9;
    csSolve.lsol.pc_type=PHY_BASED_P1;
    csSolve.lsol.iter_max=100;
    csSolve.lsol.restart_gmres=30;
    csSolve.lsol.is_CG=true;

    cs.bc_flux=Wave_BC_normalvelocity_null;
    cs.bc_assembly=BoundaryConditionFriedrichsAssembly;
    csSolve.bc_flux=Wave_BC_normalvelocity_null;
    csSolve.bc_assembly=BoundaryConditionFriedrichsAssembly;
    csSolve.type_bc=1;

    Wave_test(&cs,-(1.0-simu2.theta),simu2.dt);
    GenericOperator_Continuous(&cs);
    Wave_test(&csSolve,simu2.theta,simu2.dt);
    GenericOperator_Continuous(&csSolve);
  
    PiDgToCg(&cs,simu2.w,wCG);

    simu2.tnow=0;
    for(int ie=0; ie < simu2.macromesh.nbelems; ++ie){
      simu2.fd[ie].tnow=simu2.tnow;
    }

    for(int i=0;i<size;i++){
      csSolve.lsol.sol[i]=wCG[i];
    }

    for(int tstep=0;tstep<simu2.itermax_rk;tstep++){
      for (int i=0; i<size; i++){
	cs.lsol.rhs[i]=0;
      }
      cs.bc_assembly(&cs);
      MatVect(&cs.lsol,wCG,resCG);

      
      simu2.tnow=simu2.tnow+simu2.dt;
      for(int ie=0; ie < simu2.macromesh.nbelems; ++ie){
        simu2.fd[ie].tnow=simu2.tnow;
      }
      
      for (int i=0; i<size; i++){
	csSolve.lsol.rhs[i]=0;
      }      
      csSolve.bc_assembly(&csSolve);
      for (int i=0; i<size; i++){
	csSolve.lsol.rhs[i]=csSolve.lsol.rhs[i]+resCG[i]-cs.lsol.rhs[i];
      }
   
      Advanced_SolveLinearSolver(&csSolve.lsol,&simu2);
      
      for (int i=0; i<size; i++){
        wCG[i] = csSolve.lsol.sol[i];
      }
       
      int freq = (1 >= simu2.itermax_rk / 10)? 1 : simu2.itermax_rk / 10;
      if (tstep % freq == 0)
        printf("t=%f iter=%d/%d dt=%f\n", simu2.tnow, tstep, simu2.itermax_rk, simu2.dt);
    }
    PiInvertCgToDg(&cs,wCG,simu2.w);
    dd = L2error(&simu2);

    printf("erreur L2= %.12e\n", dd);

    test = test && (dd<5.e-2);
    freeSimulation(&simu2);
  }

#ifdef PARALUTION 
  paralution_end();
#endif 

  return test;
}



void SteadyStateOne_ImposedData(const real *xy, const real t, real *w) {

  real x=xy[0];
  real y=xy[1];

  w[0] = 10.0;//10
  w[1] = x*y;
  w[2] = -y*y*0.5;

}


void SteadyStateOne_Source(const real *xy, const real t, const real *w, real *S){

  real x=xy[0];
  real y=xy[1];

  S[0] = 0;
  S[1] = 0;//exp(x);//2*x;
  S[2] = 0;//2*exp(2*y);//3*y*y;

  S[0] *= _SPEED_WAVE;
  S[1] *= _SPEED_WAVE;
  S[2] *= _SPEED_WAVE;

}

void SteadyStateOne_InitData(real *x, real *w) {
  real t = 0;
  SteadyStateOne_ImposedData(x, t, w);
}

void SteadyStateOne_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
    real *flux) {
  real wR[3];
  SteadyStateOne_ImposedData(x , t, wR);
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





