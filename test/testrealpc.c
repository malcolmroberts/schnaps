#include "implicit.h"
#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "test.h"
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
void Interface(ContinuousSolver *cs,LinearSolver *solver,real theta, real dt);

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
  int test1_ok=0,test2_ok=0,test3_ok=1;


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

  ///////////////////////////////// Test THREE: Time-dependent problem
  printf("//////////////////////////////////////\n");
  printf("TEST ONE : Stationnary DG !\n");
  printf("//////////////////////////////////////\n");
  if(test1_ok==1){
    Model model3;

    model3.m = 3; 
    model3.NumFlux=Wave_Upwind_NumFlux;
    model3.InitData = SteadyStateOne_InitData;
    model3.ImposedData = SteadyStateOne_ImposedData;
    model3.BoundaryFlux = SteadyStateOne_BoundaryFlux;
    model3.Source = SteadyStateOne_Source;

    int deg3[]={4, 4, 0};
    int raf3[]={4, 4, 1};
    CheckMacroMesh(&mesh, deg3, raf3);
    Simulation simu3;

    InitSimulation(&simu3, &mesh, deg3, raf3, &model3);

    LinearSolver solver_implicit3;
    LinearSolver solver_explicit3;  

    real theta3=0.5;
    simu3.theta=theta3;
    simu3.dt=0.0001;
    simu3.vmax=_SPEED_WAVE;
    real tmax3 = 1*simu3.dt;

    real itermax3=tmax3/simu3.dt;
    simu3.itermax_rk=itermax3;
    InitImplicitLinearSolver(&simu3, &solver_implicit3);
    InitImplicitLinearSolver(&simu3, &solver_explicit3);
    
    real *res3 = calloc(simu3.wsize, sizeof(real));
    solver_implicit3.solver_type=GMRES;//LU;//GMRES;
    solver_implicit3.tol=1.e-8;
    solver_implicit3.pc_type=NONE;//PHY_BASED_P2;
    solver_implicit3.iter_max=2000;
    solver_implicit3.restart_gmres=15;
    solver_implicit3.is_CG=false;

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

      for(int i=0;i<solver_implicit3.neq;i++){
        solver_implicit3.rhs[i]=0;
      }
      
      AssemblyImplicitLinearSolver(&simu3, &solver_implicit3,theta3,simu3.dt);

      MatVect(&solver_explicit3, simu3.w, res3);

      for(int i=0;i<solver_implicit3.neq;i++){
        solver_implicit3.rhs[i]=solver_implicit3.rhs[i]+res3[i]-solver_explicit3.rhs[i];
      }

      Advanced_SolveLinearSolver(&solver_implicit3,&simu3);

      for(int i=0;i<simu3.wsize;i++){ 
        simu3.fd[0].wn[i]=solver_implicit3.sol[i];
      }
      for(int i=0;i<simu3.wsize;i++){ 
        solver_implicit3.sol[i]=0;
      } 

      int freq = (1 >= simu3.itermax_rk / 10)? 1 : simu3.itermax_rk / 10;
      if (tstep % freq == 0)
        printf("t=%f iter=%d/%d dt=%f\n", simu3.tnow, tstep, simu3.itermax_rk, simu3.dt);
    }
    dd = L2error(&simu3);

    printf("erreur L2=%.13e\n", dd);

    test = test && (dd<1.e-2);

    freeSimulation(&simu3);
  }
  
  printf("//////////////////////////////////////\n");
  printf("TEST TWO : Stationnary CG !\n");
  printf("//////////////////////////////////////\n");
  if(test2_ok==1){
    Model model4;

    model4.m = 3; 
    model4.NumFlux=Wave_Upwind_NumFlux;
    model4.InitData = SteadyStateOne_InitData;
    model4.ImposedData = SteadyStateOne_ImposedData;
    model4.BoundaryFlux = SteadyStateOne_BoundaryFlux;
    model4.Source = SteadyStateOne_Source;

    int deg4[]={4, 4, 0};
    int raf4[]={8, 8, 1};

    assert(mesh.is2d);
    CheckMacroMesh(&mesh, deg4, raf4);
    Simulation simu4;

    InitSimulation(&simu4, &mesh, deg4, raf4, &model4);

    ContinuousSolver cs;
    ContinuousSolver csSolve;
    int nbvar=3;
    int *listvar=calloc(nbvar,sizeof(int));
    listvar[0]=0;
    listvar[1]=1;
    listvar[2]=2;
    InitContinuousSolver(&cs,&simu4,1,nbvar,listvar);
    InitContinuousSolver(&csSolve,&simu4,1,nbvar,listvar);

    real theta4=0.5;
    simu4.theta=theta4;
    simu4.dt=100;
    simu4.vmax=_SPEED_WAVE;
    real tmax4=5*simu4.dt;
    int itermax4=tmax4/simu4.dt;
    simu4.itermax_rk=itermax4;
    int size = cs.nb_fe_dof;
    real *resCG = calloc(size, sizeof(real));
    real *wCG = calloc(size, sizeof(real));
  
    csSolve.lsol.solver_type=GMRES;
    csSolve.lsol.tol=1.e-8;
    csSolve.lsol.pc_type=PHY_BASED_U2;//PHY_BASED_U1//JACOBI;
    csSolve.lsol.iter_max=2000;
    csSolve.lsol.restart_gmres=15;
    csSolve.lsol.is_CG=true;

    cs.bc_flux=Wave_BC_pressure_imposed;
    cs.bc_assembly=BoundaryConditionFriedrichsAssembly;
    csSolve.bc_flux=Wave_BC_pressure_imposed;
    csSolve.bc_assembly=BoundaryConditionFriedrichsAssembly;
    csSolve.type_bc=2;
    
    Wave_test(&cs,-(1.0-simu4.theta),simu4.dt);
    GenericOperator_Continuous(&cs,&cs.lsol);
    Wave_test(&csSolve,simu4.theta,simu4.dt);
    GenericOperator_Continuous(&csSolve,&csSolve.lsol);
  
    PiDgToCg(&cs,simu4.w,wCG);

    simu4.tnow=0;
    for(int ie=0; ie < simu4.macromesh.nbelems; ++ie){
      simu4.fd[ie].tnow=simu4.tnow;
    } 

    for(int tstep=0;tstep<simu4.itermax_rk;tstep++){

      for (int i=0; i<size; i++){
	cs.lsol.rhs[i]=0;
      }
      cs.bc_assembly(&cs, &cs.lsol);
      MatVect(&cs.lsol,wCG,resCG);
  
      simu4.tnow=simu4.tnow+simu4.dt;
      for(int ie=0; ie < simu4.macromesh.nbelems; ++ie){
        simu4.fd[ie].tnow=simu4.tnow;
      }

      for (int i=0; i<size; i++){
	csSolve.lsol.rhs[i]=0;
      }
      csSolve.bc_assembly(&csSolve, &csSolve.lsol);
      
      for (int i=0; i<size; i++){
	csSolve.lsol.rhs[i]=csSolve.lsol.rhs[i]+resCG[i]-cs.lsol.rhs[i];
      }
      
       Advanced_SolveLinearSolver(&csSolve.lsol,&simu4);
      
      for (int i=0; i<size; i++){
        wCG[i] = csSolve.lsol.sol[i];
      }
       
      int freq = (1 >= simu4.itermax_rk / 10)? 1 : simu4.itermax_rk / 10;
      if (tstep % freq == 0)
        printf("t=%f iter=%d/%d dt=%f\n", simu4.tnow, tstep, simu4.itermax_rk, simu4.dt);
    }
    PiInvertCgToDg(&cs,wCG,simu4.w);
    dd = L2error(&simu4);

    printf("erreur L2=%.12e\n", dd);

    test = test && (dd<2.e-2);
    freeSimulation(&simu4);
  }
  
  printf("//////////////////////////////////////\n");
  printf("TEST THREE : Time dependant CG !\n");
  printf("//////////////////////////////////////\n");
  if(test3_ok==1){
    Model model5;

    model5.m = 3; 
    model5.NumFlux=Wave_Upwind_NumFlux;
    model5.InitData = TestPeriodic_Wave_InitData;
    model5.ImposedData = TestPeriodic_Wave_ImposedData;
    model5.BoundaryFlux = Wave_Upwind_BoundaryFlux;
    model5.Source = NULL;

    int deg5[]={3, 3, 0};
    int raf5[]={1, 1, 1};

    assert(mesh.is2d);
    CheckMacroMesh(&mesh, deg5, raf5);
    Simulation simu5;

    InitSimulation(&simu5, &mesh, deg5, raf5, &model5); 
    ContinuousSolver cs;
    ContinuousSolver csSolve;
    int nbvar=3;
    int *listvar=calloc(nbvar,sizeof(int));
    listvar[0]=0;
    listvar[1]=1;
    listvar[2]=2;
    InitContinuousSolver(&cs,&simu5,1,nbvar,listvar);
    InitContinuousSolver(&csSolve,&simu5,1,nbvar,listvar);

    real theta5=1.0;
    simu5.theta=theta5;
    simu5.dt=0.1;//0.001667;
    simu5.vmax=_SPEED_WAVE;
    real tmax5=5*simu5.dt;//;0.5;

    int itermax5=tmax5/simu5.dt;
    simu5.itermax_rk=itermax5;
    int size = cs.nb_fe_dof;
    real *resCG = calloc(size, sizeof(real));
    real *wCG = calloc(size, sizeof(real));

    csSolve.lsol.solver_type=GMRES;
    csSolve.lsol.tol=1.e-8;
    csSolve.lsol.pc_type=PHY_BASED_U1;
    csSolve.lsol.iter_max=5000;
    csSolve.lsol.restart_gmres=30;
    csSolve.lsol.is_CG=true;

    cs.bc_flux=Wave_BC_normalvelocity_null;
    cs.bc_assembly=BoundaryConditionFriedrichsAssembly;
    csSolve.bc_flux=Wave_BC_normalvelocity_null;
    csSolve.bc_assembly=BoundaryConditionFriedrichsAssembly;
    csSolve.type_bc=1;

    Wave_test(&cs,-(1.0-simu5.theta),simu5.dt);
    GenericOperator_Continuous(&cs,&cs.lsol);
    Wave_test(&csSolve,simu5.theta,simu5.dt);
    GenericOperator_Continuous(&csSolve,&csSolve.lsol);
  
    PiDgToCg(&cs,simu5.w,wCG);

    simu5.tnow=0;
    for(int ie=0; ie < simu5.macromesh.nbelems; ++ie){
      simu5.fd[ie].tnow=simu5.tnow;
    }

    for(int i=0;i<size;i++){
      csSolve.lsol.sol[i]=wCG[i];
    }

    for(int tstep=0;tstep<simu5.itermax_rk;tstep++){

      for (int i=0; i<size; i++){
	cs.lsol.rhs[i]=0;
      }
      cs.bc_assembly(&cs, &cs.lsol);
      MatVect(&cs.lsol,wCG,resCG);

      
      simu5.tnow=simu5.tnow+simu5.dt;
      for(int ie=0; ie < simu5.macromesh.nbelems; ++ie){
        simu5.fd[ie].tnow=simu5.tnow;
      }
      
      for (int i=0; i<size; i++){
	csSolve.lsol.rhs[i]=0;
      }      
      csSolve.bc_assembly(&csSolve, &csSolve.lsol);
      for(int i=0;i<size/3;i++){
        printf("CACA resCGP[%d] = %8e, resCGU[%d] = %8e, resCGV[%d] = %8e\n", i, resCG[3*i], i, resCG[3*i+1], i, resCG[3*i+2]);
        printf("CACA csERHS[%d] = %8e, csERHS[%d] = %8e,  csERHS[%d] = %8e\n", i, cs.lsol.rhs[3*i], i, cs.lsol.rhs[3*i+1], i, cs.lsol.rhs[3*i+2]);
        printf("CACA csIRHS[%d] = %8e, csIRHS[%d] = %8e,  csIRHS[%d] = %8e\n", i, csSolve.lsol.rhs[3*i], i, csSolve.lsol.rhs[3*i+1], i, csSolve.lsol.rhs[3*i+2]);
      }
      for (int i=0; i<size; i++){
	csSolve.lsol.rhs[i]=csSolve.lsol.rhs[i]+resCG[i]-cs.lsol.rhs[i];
      }

      Advanced_SolveLinearSolver(&csSolve.lsol,&simu5);
      
      for (int i=0; i<size; i++){
        wCG[i] = csSolve.lsol.sol[i];
      }
       
      int freq = (1 >= simu5.itermax_rk / 10)? 1 : simu5.itermax_rk / 10;
      if (tstep % freq == 0)
        printf("t=%f iter=%d/%d dt=%f\n", simu5.tnow, tstep, simu5.itermax_rk, simu5.dt);
    }
    PiInvertCgToDg(&cs,wCG,simu5.w);
    dd = L2error(&simu5);

    printf("erreur L2= %.12e\n", dd);

    test = test && (dd<5.e-2);
    freeSimulation(&simu5);
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





