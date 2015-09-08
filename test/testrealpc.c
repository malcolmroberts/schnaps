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
  int test1_ok=0,test2_ok=0,test3_ok=0,test4_ok=0,test5_ok=1;


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
  if(test1_ok==1){
    Model model;

    model.m = 3; 
    model.NumFlux=Wave_Upwind_NumFlux;
    model.InitData = SteadyStateOne_InitData;
    model.ImposedData = SteadyStateOne_ImposedData;
    model.BoundaryFlux = SteadyStateOne_BoundaryFlux;
    model.Source = SteadyStateOne_Source;

    real k=1;

    int deg[]={4, 4, 0};
    int raf[]={4*k, 4*k, 1};


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
    real tmax=5*simu.dt;

    int itermax=tmax/simu.dt;
    simu.itermax_rk=itermax;
    InitImplicitLinearSolver(&simu, &solver_implicit);
    InitImplicitLinearSolver(&simu, &solver_explicit);
    real *res = calloc(simu.wsize, sizeof(real));

    solver_implicit.solver_type=GMRES;
    solver_implicit.pc_type=PHY_BASED;
    solver_implicit.tol=1.e-10;
    solver_implicit.iter_max=200;
    solver_implicit.restart_gmres=30;


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

      Advanced_SolveLinearSolver(&solver_implicit,&simu);

      for(int i=0;i<simu.wsize;i++){ 
        simu.w[i]=solver_implicit.sol[i];
      }
      for(int i=0;i<simu.wsize;i++){ 
        solver_implicit.sol[i]=0; // sans cela ca marche pas avec pc for 4*4, avec ca marche pour 4*4 mais pas 8*8
      } 

      int freq = (1 >= simu.itermax_rk / 10)? 1 : simu.itermax_rk / 10;
      if (tstep % freq == 0)
        printf("t=%f iter=%d/%d dt=%f\n", simu.tnow, tstep, simu.itermax_rk, simu.dt);
    }
    dd = L2error(&simu);

    printf("erreur L2=%.12e\n", dd);

    test = test && (dd<2.e-2);
    freeSimulation(&simu);
  }
  ///////////////////////////////// Test TWO: Splitting error
  printf("//////////////////////////////////////\n");
  printf("TEST TWO!\n");
  printf("//////////////////////////////////////\n");
  if(test2_ok==1){
    Model model2;

    model2.m = 3; 
    model2.NumFlux=Wave_Upwind_NumFlux;
    model2.InitData = SteadyStateTwo_InitData;
    model2.ImposedData = SteadyStateTwo_ImposedData;
    model2.BoundaryFlux = SteadyStateTwo_BoundaryFlux;
    model2.Source = SteadyStateTwo_Source;

    int deg2[]={4, 4, 0};
    int raf2[]={10, 10, 1};


    CheckMacroMesh(&mesh, deg2, raf2);
    Simulation simu2;

    InitSimulation(&simu2, &mesh, deg2, raf2, &model2);


    LinearSolver solver_implicit2;
    LinearSolver solver_explicit2;  

    real theta2=0.5;
    simu2.theta=theta2;
    simu2.dt=0.01;
    simu2.vmax=_SPEED_WAVE;
    real tmax2 = 5*simu2.dt;//0.08;

    real itermax2=tmax2/simu2.dt;
    simu2.itermax_rk=itermax2;
    InitImplicitLinearSolver(&simu2, &solver_implicit2);
    InitImplicitLinearSolver(&simu2, &solver_explicit2);
    real *res2 = calloc(simu2.wsize, sizeof(real));

    solver_implicit2.solver_type=GMRES;//LU;//GMRES;
    solver_implicit2.tol=1.e-8;
    solver_implicit2.pc_type=PHY_BASED;
    solver_implicit2.iter_max=500;

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

      Advanced_SolveLinearSolver(&solver_implicit2,&simu2);

      for(int i=0;i<simu2.wsize;i++){ 
        simu2.fd[0].wn[i]=solver_implicit2.sol[i];
      }
      for(int i=0;i<simu2.wsize;i++){ 
        solver_implicit2.sol[i]=0; // sans cela ca marche pas le PC pk ?????????
      } 

      int freq = (1 >= simu2.itermax_rk / 10)? 1 : simu2.itermax_rk / 10;
      if (tstep % freq == 0)
        printf("t=%f iter=%d/%d dt=%f\n", simu2.tnow, tstep, simu2.itermax_rk, simu2.dt);
    }
    dd = L2error(&simu2);

    printf("erreur L2=%.12e\n", dd);

    test = test && (dd<1.e-2);
    freeSimulation(&simu2);
  } 

  ///////////////////////////////// Test THREE: Time-dependent problem
  printf("//////////////////////////////////////\n");
  printf("TEST THREE!\n");
  printf("//////////////////////////////////////\n");
  if(test3_ok==1){
    Model model3;

    model3.m = 3; 
    model3.NumFlux=Wave_Upwind_NumFlux;
    model3.InitData = TestPeriodic_Wave_InitData;
    model3.ImposedData = TestPeriodic_Wave_ImposedData;
    model3.BoundaryFlux = Wave_Upwind_BoundaryFlux;
    model3.Source = NULL;

    int deg3[]={4, 4, 0};
    int raf3[]={4, 4, 1};


    CheckMacroMesh(&mesh, deg3, raf3);
    Simulation simu3;

    InitSimulation(&simu3, &mesh, deg3, raf3, &model3);


    LinearSolver solver_implicit3;
    LinearSolver solver_explicit3;  

    real theta3=0.5;
    simu3.theta=theta3;
    simu3.dt=0.05;
    simu3.vmax=_SPEED_WAVE;
    real tmax3 = 10*simu3.dt;

    real itermax3=tmax3/simu3.dt;
    simu3.itermax_rk=itermax3;
    InitImplicitLinearSolver(&simu3, &solver_implicit3);
    InitImplicitLinearSolver(&simu3, &solver_explicit3);
    real *res3 = calloc(simu3.wsize, sizeof(real));

    solver_implicit3.solver_type=GMRES;//LU;//GMRES;
    solver_implicit3.tol=1.e-8;
    solver_implicit3.pc_type=PHY_BASED;
    solver_implicit3.iter_max=2000;
    solver_implicit3.restart_gmres=30;

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


      Advanced_SolveLinearSolver(&solver_implicit3,&simu3);

      for(int i=0;i<simu3.wsize;i++){ 
        simu3.fd[0].wn[i]=solver_implicit3.sol[i];
      }
      for(int i=0;i<simu3.wsize;i++){ 
        solver_implicit3.sol[i]=0;
      } 

      //
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
  printf("TEST FOUR!\n");
  printf("//////////////////////////////////////\n");
  if(test4_ok==1){
    Model model4;

    model4.m = 3; 
    model4.NumFlux=Wave_Upwind_NumFlux;
    model4.InitData = SteadyStateOne_InitData;
    model4.ImposedData = SteadyStateOne_ImposedData;
    model4.BoundaryFlux = SteadyStateOne_BoundaryFlux;
    model4.Source = SteadyStateOne_Source;

    int deg4[]={4, 4, 0};
    int raf4[]={16, 16, 1};


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
    simu4.dt=10;//0.001667;
    simu4.vmax=_SPEED_WAVE;
    real tmax4=10*simu4.dt;//10*simu.dt;//;0.5;
    int itermax4=tmax4/simu4.dt;
    simu4.itermax_rk=itermax4;
    int size = cs.nb_fe_dof;
    real *resCG = calloc(size, sizeof(real));
    real *wCG = calloc(size, sizeof(real));
  

    csSolve.lsol.solver_type=GMRES;//LU;
    csSolve.lsol.tol=1.e-10;
    csSolve.lsol.pc_type=PHY_BASED;//PHY_BASED;//;NONE;//EXACT;//PHY_BASED;
    csSolve.lsol.iter_max=10000;
    csSolve.lsol.restart_gmres=5;
    csSolve.lsol.is_CG=true;
    csSolve.bc_assembly=ExactDirichletContinuousMatrix;

    //////////////////////////////////
    PB_PC pb_pc;
     int mat2assemble[6] = {1, 1, 1, 1, 1, 1};
     Init_PhyBasedPC_SchurPressure_Wave(&simu4, &pb_pc, mat2assemble);
     Init_Parameters_PhyBasedPC(&pb_pc);
       real *solpc = calloc(size, sizeof(real));
     ////////////////////////////

    Wave_test(&cs,-(1.0-simu4.theta),simu4.dt);
    Generic(&cs);
    Wave_test(&csSolve,simu4.theta,simu4.dt);
    Generic(&csSolve);

    PiDgToCg(&cs,simu4.w,wCG);

    simu4.tnow=0;
    for(int ie=0; ie < simu4.macromesh.nbelems; ++ie){
      simu4.fd[ie].tnow=simu4.tnow;
    } 

    for(int tstep=0;tstep<simu4.itermax_rk;tstep++){

      MatVect(&cs.lsol,wCG,resCG);
      simu4.tnow=simu4.tnow+simu4.dt;
      for(int ie=0; ie < simu4.macromesh.nbelems; ++ie){
        simu4.fd[ie].tnow=simu4.tnow;
      } 
      for(int i=0;i<size;i++){
        csSolve.lsol.rhs[i]=resCG[i];
      }
      csSolve.bc_assembly(&csSolve, &csSolve.lsol);     
      Advanced_SolveLinearSolver(&csSolve.lsol,&simu4);
      
      ///////////////////////////////////////
      /*PhyBased_PC_Full(&pb_pc,&simu4,solpc,csSolve.lsol.rhs);
       real error=0;
       for (int i=0; i<size; i++){
	 error=error+fabs((solpc[i]-csSolve.lsol.sol[i])*(solpc[i]-csSolve.lsol.sol[i]));
	 //printf("pppp %d %.12e %.12e %.12e\n",i,solpc[i],csSolve.lsol.sol[i],wCG[i]-csSolve.lsol.sol[i]);
       }
       printf("pppp %.12e\n",sqrt(error));*/
	 /////////////////////////////////// */
      for (int i=0; i<size; i++){
        wCG[i] = csSolve.lsol.sol[i];
      }
      
      /*real error2=0;
      for (int i=0; i<size/3; i++){
	 error2=error2+(10.0-csSolve.lsol.sol[3*i])*(10.0-csSolve.lsol.sol[3*i]);
	 error2=error2+(2.0-csSolve.lsol.sol[3*i+1])*(2.0-csSolve.lsol.sol[3*i+1]);
	 error2=error2+(3.0-csSolve.lsol.sol[3*i+2])*(3.0-csSolve.lsol.sol[3*i+2]);
       }
       printf("pppp  sol %.12e\n",sqrt(error2));*/
      
      
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
  printf("TEST FIVE!\n");
  printf("//////////////////////////////////////\n");
  if(test5_ok==1){
    Model model5;

    model5.m = 3; 
    model5.NumFlux=Wave_Upwind_NumFlux;
    model5.InitData = TestPeriodic_Wave_InitData;
    model5.ImposedData = TestPeriodic_Wave_ImposedData;
    model5.BoundaryFlux = Wave_Upwind_BoundaryFlux;
    model5.Source = NULL;

    int deg5[]={4, 4, 0};
    int raf5[]={16, 16, 1};

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

    real theta5=0.5;
    simu5.theta=theta5;
    simu5.dt=10;//0.001667;
    simu5.vmax=_SPEED_WAVE;
    real tmax5=1*simu5.dt;//;0.5;

    int itermax5=tmax5/simu5.dt;
    simu5.itermax_rk=itermax5;
    int size = cs.nb_fe_dof;
    real *resCG = calloc(size, sizeof(real));
    real *wCG = calloc(size, sizeof(real));

    csSolve.lsol.solver_type=GMRES;
    csSolve.lsol.tol=1.e-10;
    csSolve.lsol.pc_type=PHY_BASED;//PHY_BASED;
    csSolve.lsol.iter_max=500;
    csSolve.lsol.restart_gmres=30;
    csSolve.lsol.is_CG=true;
    csSolve.bc_assembly=ExactDirichletContinuousMatrix;

     //////////////////////////////////
    PB_PC pb_pc;
     int mat2assemble[6] = {1, 1, 1, 1, 1, 1};
     Init_PhyBasedPC_SchurPressure_Wave(&simu5, &pb_pc, mat2assemble);
     //Init_PhyBasedPC_SchurFull_Wave(&simu5, &pb_pc, mat2assemble);
     Init_Parameters_PhyBasedPC(&pb_pc);
     real *solpc = calloc(size, sizeof(real));
     ////////////////////////////

    Wave_test(&cs,-(1.0-simu5.theta),simu5.dt);
    Generic(&cs);
    Wave_test(&csSolve,simu5.theta,simu5.dt);
    Generic(&csSolve);

    PiDgToCg(&cs,simu5.w,wCG);

    simu5.tnow=0;
    for(int ie=0; ie < simu5.macromesh.nbelems; ++ie){
      simu5.fd[ie].tnow=simu5.tnow;
    }

    for(int i=0;i<size;i++){
      csSolve.lsol.sol[i]=wCG[i];
    }

    for(int tstep=0;tstep<simu5.itermax_rk;tstep++){

      MatVect(&cs.lsol,wCG,resCG);
      simu5.tnow=simu5.tnow+simu5.dt;
      for(int ie=0; ie < simu5.macromesh.nbelems; ++ie){
        simu5.fd[ie].tnow=simu5.tnow;
      } 
      for(int i=0;i<size;i++){
        csSolve.lsol.rhs[i]=resCG[i];
      }
      csSolve.bc_assembly(&csSolve, &csSolve.lsol);
      Advanced_SolveLinearSolver(&csSolve.lsol,&simu5);
      
      for (int i=0; i<size; i++){
        wCG[i] = csSolve.lsol.sol[i];
      }

          ///////////////////////////////////////
      //PhyBased_PC_Full(&pb_pc,&simu5,solpc,csSolve.lsol.rhs);
      PhyBased_PC_InvertSchur_CG(&pb_pc,&simu5,solpc,csSolve.lsol.rhs);
       real error=0;
       for (int i=0; i<size; i++){
	 error=error+fabs((solpc[i]-csSolve.lsol.sol[i])*(solpc[i]-csSolve.lsol.sol[i]));
       }
       printf("pppp %.12e\n",sqrt(error));
	 /////////////////////////////////// */
       /*for (int i=0; i<size/3; i++){
	 printf("iter=%d p=%.12e %.12e %.12e\n",i,solpc[3*i],wCG[3*i],solpc[3*i]-wCG[3*i]);
       }

       for (int i=0; i<size/3; i++){
	 printf("iter=%d u=%.12e %.12e %.12e\n",i,solpc[3*i+1],wCG[3*i+1],solpc[3*i+1]-wCG[3*i+1]);
       }

       for (int i=0; i<size/3; i++){
	 printf("iter=%d v=%.12e %.12e %.12e\n",i,solpc[3*i+2],wCG[3*i+2],solpc[3*i+2]-wCG[3*i+2]);
	 }*/
       
      int freq = (1 >= simu5.itermax_rk / 10)? 1 : simu5.itermax_rk / 10;
      if (tstep % freq == 0)
        printf("t=%f iter=%d/%d dt=%f\n", simu5.tnow, tstep, simu5.itermax_rk, simu5.dt);
    }
    PiInvertCgToDg(&cs,wCG,simu5.w);
    dd = L2error(&simu5);

    printf("erreur L2=%.12e\n", dd);

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

  w[0] = 10.0;//+exp(x)+exp(2*y); // 10+x*x+y*y*y
  w[1] = 0.2*x*x*x*x*x*y-exp(y)+2;
  w[2] = exp(x)-x*x*x*x*y*y*0.5+5.6;

}

void SteadyStateTwo_ImposedData(const real *xy, const real t, real *w) {

  real x=xy[0];
  real y=xy[1];


  /*w[0] = x*(1-x)*y*(1-y);
    w[1] = 2*x*(1-x)*y*(1-y);
    w[2] = 3*x*(1-x)*y*(1-y);*/

  real x5=x*x*x;
  real y5=y*y*y;

  w[0] = x5*(1-x5)*y5*(1-y5);
  w[1] = 2.0*x5*(1-x5)*y5*(1-y5);
  w[2] = 2.0*x5*(1-x5)*y5*(1-y5);

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

void SteadyStateTwo_Source(const real *xy, const real t, const real *w, real *S){

  real x=xy[0];
  real y=xy[1];

  /*S[0] = 2*(1-2*x)*(y*(1-y))+3*(1-2*y)*(x*(1-x));
    S[1] = (1-2*x)*(y*(1-y));
    S[2] = (1-2*y)*(x*(1-x));*/


  real x5=x*x*x;
  real x4=3.0*x*x;
  real y5=y*y*y;
  real y4=3.0*y*y;

  S[0] = 2.0*(x5*(-x4)+x4*(1.0-x5))*y5*(1.0-y5)+2.0*(y5*(-y4)+y4*(1.0-y5))*x5*(1.0-x5);
  S[1] = (x5*(-x4)+x4*(1.0-x5))*y5*(1.0-y5);
  S[2] = (y5*(-y4)+y4*(1.0-y5))*x5*(1.0-x5);

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



void Generic(ContinuousSolver *cs){

  int nnodes = cs->nnodes;
  field* f0 = &cs->simu->fd[0];
  int nbel = cs->nbel;
  int dg_to_fe_index[cs->nb_dg_nodes];
  for (int i=0; i<cs->nb_dg_nodes; i++){
    dg_to_fe_index[i]=cs->dg_to_fe_index[i];
  }

 
  for(int ie = 0; ie < nbel; ie++){

     real aloc[cs->nnodes*cs->nb_phy_vars][cs->nnodes*cs->nb_phy_vars];
     for(int iloc = 0; iloc < cs->nnodes*cs->nb_phy_vars; iloc++){
       for(int jloc = 0; jloc < cs->nnodes*cs->nb_phy_vars; jloc++){
	 aloc[iloc][jloc] = 0;
       }
     }


    int iemacro = ie / (f0->raf[0] * f0->raf[1] * f0->raf[2]);
    int isubcell = ie % (f0->raf[0] * f0->raf[1] * f0->raf[2]);

    for(int ipg = 0;ipg < nnodes; ipg++){
      real wpg;
      real xref[3];
      int ipgmacro = ipg + isubcell * nnodes;

      ref_pg_vol(f0->deg,f0->raf,ipgmacro,xref,&wpg,NULL);

      for(int iloc = 0; iloc < nnodes; iloc++){
        real dtau[3][3],codtau[3][3];
        real dphiref_i[3],dphiref_j[3];
        real dphi_i[3],dphi_j[3];
        real basisPhi_i[4], basisPhi_j[4];
        int ilocmacro = iloc + isubcell * nnodes;
        int ino_dg = iloc + ie * nnodes;
        int ino_fe = dg_to_fe_index[ino_dg];
        grad_psi_pg(f0->deg,f0->raf,ilocmacro,ipgmacro,dphiref_i);
        Ref2Phy(cs->simu->fd[iemacro].physnode,
            xref,dphiref_i,0,NULL,
            dtau,codtau,dphi_i,NULL);

        real det = dot_product(dtau[0], codtau[0]);
        if (ilocmacro==ipgmacro){
          basisPhi_i[0]=1.;
        }
        else
        {
          basisPhi_i[0]=0.;
        }
        basisPhi_i[1]=dphi_i[0]/det;
        basisPhi_i[2]=dphi_i[1]/det;
        basisPhi_i[3]=dphi_i[2]/det;
        for(int jloc = 0; jloc < nnodes; jloc++){
          int jlocmacro = jloc + isubcell * nnodes;
          int jno_dg = jloc + ie * nnodes;
          int jno_fe = dg_to_fe_index[jno_dg];
          grad_psi_pg(f0->deg,f0->raf,jlocmacro,ipgmacro,dphiref_j);
          Ref2Phy(cs->simu->fd[iemacro].physnode,
              xref,dphiref_j,0,NULL,
              dtau,codtau,dphi_j,NULL);
          if (jlocmacro==ipgmacro){
            basisPhi_j[0]=1.;
          }
          else
          {
            basisPhi_j[0]=0.;
          }
          basisPhi_j[1]=dphi_j[0]/det;
          basisPhi_j[2]=dphi_j[1]/det;
          basisPhi_j[3]=dphi_j[2]/det;
          real val;
          real res[4];

          // Building Schur Matrix
          for (int iv1=0; iv1<cs->nb_phy_vars; iv1++){
            for (int iv2=0; iv2<cs->nb_phy_vars; iv2++){
              for (int i=0; i<4; i++){
                res[i]=0;
                for (int j=0; j<4; j++){
                  res[i]+=basisPhi_j[j]*cs->diff_op3vec[cs->nb_phy_vars*iv1+iv2][i][j];
                }
              }
              val = dot_product(basisPhi_i, res) * wpg * det;
	      aloc[iloc*cs->nb_phy_vars+iv1][jloc*cs->nb_phy_vars+iv2]+=val;
            } // end for iv1
          } // end for iv2

        } // end for jloc
      } // end for iloc
    } // end for ipg

    for(int iloc = 0; iloc < cs->nnodes; iloc++){
      for(int jloc = 0; jloc < cs->nnodes; jloc++){
          int ino_dg = iloc + ie * cs->nnodes;
          int jno_dg = jloc + ie * cs->nnodes;
          int ino_fe = cs->dg_to_fe_index[ino_dg];
          int jno_fe = cs->dg_to_fe_index[jno_dg];
          for (int iv1=0;iv1<cs->nb_phy_vars;iv1++){
            for (int iv2=0;iv2<cs->nb_phy_vars;iv2++){
              real val = aloc[iloc*cs->nb_phy_vars+iv1][jloc*cs->nb_phy_vars+iv2];
              AddLinearSolver(&cs->lsol,ino_fe*cs->nb_phy_vars+iv1,jno_fe*cs->nb_phy_vars+iv2,val);
            }
          }
      }
    }  
    } // end for ie 

  /* for(int ie = 0; ie < nbel; ie++){

     real aloc[cs->nnodes*cs->nb_phy_vars][cs->nnodes*cs->nb_phy_vars];
     for(int iloc = 0; iloc < cs->nnodes*cs->nb_phy_vars; iloc++){
       for(int jloc = 0; jloc < cs->nnodes*cs->nb_phy_vars; jloc++){
	 aloc[iloc][jloc] = 0;
       }
     }


    int iemacro = ie / (f0->raf[0] * f0->raf[1] * f0->raf[2]);
    int isubcell = ie % (f0->raf[0] * f0->raf[1] * f0->raf[2]);

    for(int ipg = 0;ipg < nnodes; ipg++){
      real wpg;
      real xref[3];
      int ipgmacro = ipg + isubcell * nnodes;

      ref_pg_vol(f0->deg,f0->raf,ipgmacro,xref,&wpg,NULL);

      for(int iloc = 0; iloc < nnodes; iloc++){
        real dtau[3][3],codtau[3][3];
        real dphiref_i[3],dphiref_j[3];
        real dphi_i[3],dphi_j[3];
        real basisPhi_i[4], basisPhi_j[4];
        int ilocmacro = iloc + isubcell * nnodes;
        int ino_dg = iloc + ie * nnodes;
        int ino_fe = dg_to_fe_index[ino_dg];
        grad_psi_pg(f0->deg,f0->raf,ilocmacro,ipgmacro,dphiref_i);
        Ref2Phy(cs->simu->fd[iemacro].physnode,
            xref,dphiref_i,0,NULL,
            dtau,codtau,dphi_i,NULL);

        real det = dot_product(dtau[0], codtau[0]);
        if (ilocmacro==ipgmacro){
          basisPhi_i[0]=1.;
        }
        else
        {
          basisPhi_i[0]=0.;
        }
        basisPhi_i[1]=dphi_i[0]/det;
        basisPhi_i[2]=dphi_i[1]/det;
        basisPhi_i[3]=dphi_i[2]/det;
        for(int jloc = 0; jloc < nnodes; jloc++){
          int jlocmacro = jloc + isubcell * nnodes;
          int jno_dg = jloc + ie * nnodes;
          int jno_fe = dg_to_fe_index[jno_dg];
          grad_psi_pg(f0->deg,f0->raf,jlocmacro,ipgmacro,dphiref_j);
          Ref2Phy(cs->simu->fd[iemacro].physnode,
              xref,dphiref_j,0,NULL,
              dtau,codtau,dphi_j,NULL);
          if (jlocmacro==ipgmacro){
            basisPhi_j[0]=1.;
          }
          else
          {
            basisPhi_j[0]=0.;
          }
          basisPhi_j[1]=dphi_j[0]/det;
          basisPhi_j[2]=dphi_j[1]/det;
          basisPhi_j[3]=dphi_j[2]/det;
          real val;
          real res[4];

          // Building Schur Matrix
          for (int iv1=0; iv1<cs->nb_phy_vars; iv1++){
            for (int iv2=0; iv2<cs->nb_phy_vars; iv2++){
              for (int i=0; i<4; i++){
                res[i]=0;
                for (int j=0; j<4; j++){
                  res[i]+=basisPhi_j[j]*cs->diff_op3vec[cs->nb_phy_vars*iv1+iv2][i][j];
                }
              }
              val = dot_product(basisPhi_i, res) * wpg * det;
	      AddLinearSolver(&cs->lsol,ino_fe*cs->nb_phy_vars+iv1,jno_fe*cs->nb_phy_vars+iv2,val);
            } // end for iv1
          } // end for iv2

        } // end for jloc
      } // end for iloc
    } // end for ipg

    
  }*/
  
}


