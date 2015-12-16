#include "implicit.h"
#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "../test/test.h"
#include "solverwave.h"
#include "waterwave2d.h"
#include "physBased_PC.h"
#include <math.h>


//! \brief boundary rusanov flux based for a Steady State with velocity imposed
//! \param[in] t : current time
//! \param[in] x : current position
//! \param[in] wL : states
//! \param[in] vnorm : normal vector
//! \param[out] flux : the flux
void SteadyState_U_Rusanov_BoundaryFlux(real *x, real t, real *wL, real *vnorm,real *flux);

//! \brief boundary HLL flux based for a Steady State with velocity imposed
//! \param[in] t : current time
//! \param[in] x : current position
//! \param[in] wL : states
//! \param[in] vnorm : normal vector
//! \param[out] flux : the flux
void SteadyState_U_HLL_BoundaryFlux(real *x, real t, real *wL, real *vnorm,real *flux);

//! \brief boundary Roe flux based for a Steady State with velocity imposed
//! \param[in] t : current time
//! \param[in] x : current position
//! \param[in] wL : states
//! \param[in] vnorm : normal vector
//! \param[out] flux : the flux
void SteadyState_U_Roe_BoundaryFlux(real *x, real t, real *wL, real *vnorm,real *flux);

//! \brief boundary rusanov flux based for a Steady State with pressure imposed
//! \param[in] t : current time
//! \param[in] x : current position
//! \param[in] wL : states
//! \param[in] vnorm : normal vector
//! \param[out] flux : the flux
void SteadyState_P_Rusanov_BoundaryFlux(real *x, real t, real *wL, real *vnorm,real *flux);

//! \brief boundary HLL flux based for a Steady State with pressure imposed
//! \param[in] t : current time
//! \param[in] x : current position
//! \param[in] wL : states
//! \param[in] vnorm : normal vector
//! \param[out] flux : the flux
void SteadyState_P_HLL_BoundaryFlux(real *x, real t, real *wL, real *vnorm,real *flux);

//! \brief boundary Roe flux based for a Steady State with pressure imposed
//! \param[in] t : current time
//! \param[in] x : current position
//! \param[in] wL : states
//! \param[in] vnorm : normal vector
//! \param[out] flux : the flux
void SteadyState_P_Roe_BoundaryFlux(real *x, real t, real *wL, real *vnorm,real *flux);

int main(void) {

  // unit tests

  int resu = TestpcSW();

  if (resu) printf("Steady SW test OK !\n");
  else printf("Steady SW test failed !\n");

  return !resu;
} 

int TestpcSW(void) {

  bool test = true;
  real dd;
  int test1_ok=1,test2_ok=0;

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
  printf("TEST ONE : Steady U !\n");
  printf("//////////////////////////////////////\n");
  if(test1_ok==1){
    Model model;

    model.m=6; 
    model.NumFlux=ShallowWater_Rusanov_NumFlux;
    model.InitData = TestSH_SteadyState_U_InitData;
    model.ImposedData = TestSH_SteadyState_U_ImposedData;
    model.BoundaryFlux = SteadyState_U_Rusanov_BoundaryFlux;
    model.Source = ShallowWater_SteadyState_U_SourceTerm;

    int deg[]={4, 4, 0};
    int raf[]={12, 12, 1};

    assert(mesh.is2d);
    CheckMacroMesh(&mesh, deg, raf);
    Simulation simu;
    EmptySimulation(&simu);
    InitSimulation(&simu, &mesh, deg, raf, &model);

    ContinuousSolver csSolve;
    int nbvar=3;
    int *listvar=calloc(nbvar,sizeof(int));
    listvar[0]=0;
    listvar[1]=1;
    listvar[2]=2;
    InitContinuousSolver(&csSolve,&simu,1,nbvar,listvar);

    real theta=0.5;
    simu.theta=theta;
    //simu.dt=0.002/8;
    simu.vmax=1;//_SPEED_WAVE;
    simu.cfl = 0.025;
    simu.dt = 0.0000001/2;//Get_Dt_RK(&simu)/16;
    real tmax=0.0000001;//1*simu.dt;//0.002;
    int itermax=tmax/simu.dt;
    simu.itermax_rk=itermax;

    csSolve.lsol.solver_type=LU;
    csSolve.lsol.tol=1.e-14;
    csSolve.lsol.pc_type=JACOBI;//PHY_BASED_U2;
    csSolve.lsol.iter_max=2000;
    csSolve.lsol.restart_gmres=20;
    csSolve.lsol.is_CG=true;

    csSolve.bc_flux=Wave_BC_normalvelocity_null;
    csSolve.bc_assembly=BoundaryConditionFriedrichsAssembly;
    csSolve.rhs_assembly=Source_Assembly;
    csSolve.type_bc=2;
    
    int size1varDG = csSolve.nb_dg_nodes;
    int size1varCG = csSolve.nb_fe_nodes;
    int sizeDG = csSolve.nb_dg_dof;
    int sizeCG = csSolve.nb_fe_dof;

    real *wDG = calloc(sizeDG, sizeof(real));
    real *wCG = calloc(sizeCG, sizeof(real));
    real *resCG = calloc(sizeCG, sizeof(real));
    //for(int i=0;i<size1varDG;i++){
    //  printf("i=%d, h=%.8e, u=%.8e, v=%.8e, a=%.8e, b=%8.e, c=%8.e\n",i,simu.w[i*6],simu.w[i*6+1],simu.w[i*6+2],simu.w[i*6+3],simu.w[i*6+4],simu.w[i*6+5]);
    //}
    PiDgToCg(&csSolve,model.m,simu.w,wCG);

    simu.tnow=0;
    for(int ie=0; ie < simu.macromesh.nbelems; ++ie){
      simu.fd[ie].tnow=simu.tnow;
    } 

    for(int tstep=0;tstep<simu.itermax_rk;tstep++){

      csSolve.reset_dt(&csSolve.lsol);

      simu.tnow=simu.tnow+simu.dt;
      for(int ie=0; ie < simu.macromesh.nbelems; ++ie){
        simu.fd[ie].tnow=simu.tnow;
      }

      SW_test(&csSolve,wCG,simu.theta,simu.dt);
      csSolve.rhs_assembly(&csSolve);
      csSolve.bc_assembly(&csSolve);
      Advanced_SolveLinearSolver(&csSolve.lsol,&simu);

      for (int i=0; i<sizeCG; i++){
        wCG[i] += csSolve.lsol.sol[i];
      }

      int freq = (1 >= simu.itermax_rk / 10)? 1 : simu.itermax_rk / 10;
      if (tstep % freq == 0)
        printf("t=%f iter=%d/%d dt=%f\n", simu.tnow, tstep, simu.itermax_rk, simu.dt);
    }

    PiInvertCgToDg(&csSolve,model.m,wCG,simu.w);
    dd = L2error(&simu);
    //for(int i=0;i<size1varDG;i++){
    //  printf("i=%d, h=%.8e, u=%.8e, v=%.8e, a=%.8e, b=%8.e, c=%8.e\n",i,simu.w[i*6],simu.w[i*6+1],simu.w[i*6+2],simu.w[i*6+3],simu.w[i*6+4],simu.w[i*6+5]);
    //}

    printf("erreur L2=%.12e\n", dd);
    dd = L2error_onefield(&simu,0);
    real dd2 = L2error_onefield(&simu,1);
    real dd3 = L2error_onefield(&simu,2);
    printf("erreur L2 sur h=%.12e\n", dd);
    printf("erreur L2 sur u=%.12e\n", dd2);
    printf("erreur L2 sur v=%.12e\n", dd3);

    test = test && (dd<2.e-2);
    freeSimulation(&simu);
  }
  
  printf("//////////////////////////////////////\n");
  printf("TEST TWO : Steady P !\n");
  printf("//////////////////////////////////////\n");
  if(test2_ok==1){
    Model model2;

    model2.m=6; 
    model2.NumFlux=ShallowWater_Rusanov_NumFlux;
    model2.InitData = TestSH_SteadyState_P_InitData;
    model2.ImposedData = TestSH_SteadyState_P_ImposedData;
    model2.BoundaryFlux = SteadyState_P_Rusanov_BoundaryFlux;
    model2.Source = ShallowWater_SteadyState_P_SourceTerm;

    int deg2[]={4, 4, 0};
    int raf2[]={12, 12, 1};

    assert(mesh.is2d);
    CheckMacroMesh(&mesh, deg2, raf2);
    Simulation simu2;
    EmptySimulation(&simu2);
    InitSimulation(&simu2, &mesh, deg2, raf2, &model2);

    ContinuousSolver csSolve2;
    int nbvar=3;
    int *listvar=calloc(nbvar,sizeof(int));
    listvar[0]=0;
    listvar[1]=1;
    listvar[2]=2;
    InitContinuousSolver(&csSolve2,&simu2,1,nbvar,listvar);

    real theta=0.5;
    simu2.theta=theta;
    //simu2.dt=0.002/8;
    simu2.vmax=1;//_SPEED_WAVE;
    simu2.cfl = 0.025;
    simu2.dt = 0.0000001/2;//Get_Dt_RK(&simu2)/16;
    real tmax=0.0000001;//1*simu2.dt;//0.002;
    int itermax=tmax/simu2.dt;
    simu2.itermax_rk=itermax;

    csSolve2.lsol.solver_type=LU;
    csSolve2.lsol.tol=1.e-14;
    csSolve2.lsol.pc_type=JACOBI;//PHY_BASED_U2;
    csSolve2.lsol.iter_max=2000;
    csSolve2.lsol.restart_gmres=20;
    csSolve2.lsol.is_CG=true;

    csSolve2.bc_flux=Wave_BC_pressure_imposed;
    csSolve2.bc_assembly=BoundaryConditionFriedrichsAssembly;
    csSolve2.rhs_assembly=Source_Assembly;
    csSolve2.type_bc=2;
    
    int size1varDG = csSolve2.nb_dg_nodes;
    int size1varCG = csSolve2.nb_fe_nodes;
    int sizeDG = csSolve2.nb_dg_dof;
    int sizeCG = csSolve2.nb_fe_dof;

    real *wDG = calloc(sizeDG, sizeof(real));
    real *wCG = calloc(sizeCG, sizeof(real));
    real *resCG = calloc(sizeCG, sizeof(real));
    //for(int i=0;i<size1varDG;i++){
    //  printf("i=%d, h=%.8e, u=%.8e, v=%.8e, a=%.8e, b=%8.e, c=%8.e\n",i,simu2.w[i*6],simu2.w[i*6+1],simu2.w[i*6+2],simu2.w[i*6+3],simu2.w[i*6+4],simu2.w[i*6+5]);
    //}
    PiDgToCg(&csSolve2,model2.m,simu2.w,wCG);

    simu2.tnow=0;
    for(int ie=0; ie < simu2.macromesh.nbelems; ++ie){
      simu2.fd[ie].tnow=simu2.tnow;
    } 

    for(int tstep=0;tstep<simu2.itermax_rk;tstep++){

      csSolve2.reset_dt(&csSolve2.lsol);

      simu2.tnow=simu2.tnow+simu2.dt;
      for(int ie=0; ie < simu2.macromesh.nbelems; ++ie){
        simu2.fd[ie].tnow=simu2.tnow;
      }

      SW_test(&csSolve2,wCG,simu2.theta,simu2.dt);
      csSolve2.rhs_assembly(&csSolve2);
      csSolve2.bc_assembly(&csSolve2);
      Advanced_SolveLinearSolver(&csSolve2.lsol,&simu2);
      for (int i=0; i<sizeCG; i++){
        wCG[i] += csSolve2.lsol.sol[i];
      }

      int freq = (1 >= simu2.itermax_rk / 10)? 1 : simu2.itermax_rk / 10;
      if (tstep % freq == 0)
        printf("t=%f iter=%d/%d dt=%f\n", simu2.tnow, tstep, simu2.itermax_rk, simu2.dt);
    }

    PiInvertCgToDg(&csSolve2,model2.m,wCG,simu2.w);
    dd = L2error(&simu2);
    //for(int i=0;i<size1varDG;i++){
    //  printf("i=%d, h=%.8e, u=%.8e, v=%.8e, a=%.8e, b=%8.e, c=%8.e\n",i,simu2.w[i*6],simu2.w[i*6+1],simu2.w[i*6+2],simu2.w[i*6+3],simu2.w[i*6+4],simu2.w[i*6+5]);
    //}

    printf("erreur L2=%.12e\n", dd);
    dd = L2error_onefield(&simu2,0);
    real dd2 = L2error_onefield(&simu2,1);
    real dd3 = L2error_onefield(&simu2,2);
    printf("erreur L2 sur h=%.12e\n", dd);
    printf("erreur L2 sur u=%.12e\n", dd2);
    printf("erreur L2 sur v=%.12e\n", dd3);

    test = test && (dd<2.e-2);
    freeSimulation(&simu2);
  }

#ifdef PARALUTION 
  paralution_end();
#endif 

  return test;
}


void TestSH_SteadyState_U_ImposedData(const real *x, const real t, real *w) {
  real alpha=1.0; 
  real p = 2;
  w[0] = 1.0+alpha*(pow(x[0],p)-pow(x[0],2*p))*(pow(x[1],p)-pow(x[1],2*p));
  w[1] =  (x[0]-x[0]*x[0])*(1.0-2.0*x[1]);
  w[2] = -(x[1]-x[1]*x[1])*(1.0-2.0*x[0]);
  w[3] = 0.0;
  w[4] = 0.0;
  w[5] = 0.0;
}

void TestSH_SteadyState_U_InitData(real *x, real *w) {
  real t = 0;
  TestSH_SteadyState_U_ImposedData(x, t, w);
}

void ShallowWater_SteadyState_U_SourceTerm(const real *x, const real t, const real *w, real *source){
  real g=_GRAVITY;
  real alpha=1.0;
  real p = 2;
  real S_11, S_12,S_13, S_21, S_22, S_23, S_factor;
  real wexact[6];

  real h  = 1.0+alpha*(x[0]-x[0]*x[0])*(x[1]-x[1]*x[1]);
  real hx = alpha * p * ( pow(x[1],p)-pow(x[1],2*p) ) * ( pow(x[0],p-1)-2.0*pow(x[0],2*p-1) );
  real hy = alpha * p * ( pow(x[0],p)-pow(x[0],2*p) ) * ( pow(x[1],p-1)-2.0*pow(x[1],2*p-1) );
  real u  = (x[0]-x[0]*x[0])*(1.0-2.0*x[1]);
  real ux = (1.0-2.0*x[0])*(1.0-2.0*x[1]);
  real uy = -2.0*(x[0]-x[0]*x[0]);
  real v  = -(x[1]-x[1]*x[1])*(1.0-2.0*x[0]);
  real vx = 2.0*(x[1]-x[1]*x[1]);
  real vy = -(1.0-2.0*x[0])*(1.0-2.0*x[1]);

  source[0]= h*(ux+vy) + (hx*u+hy*v);
  source[1]= h*(u*ux+v*uy) + g*h*hx;
  source[2]= h*(u*vx+v*vy) + g*h*hy;
  source[3]= 0.0;
  source[4]= 0.0;
  source[5]= 0.0;
};

void SteadyState_U_Rusanov_BoundaryFlux(real *x, real t, real *wL, real *vnorm,real *flux) {
  real wR[6];
  TestSH_SteadyState_U_ImposedData(x , t, wR);
  ShallowWater_Rusanov_NumFlux(wL, wR, vnorm, flux);
}

void SteadyState_U_HLL_BoundaryFlux(real *x, real t, real *wL, real *vnorm,real *flux) {
  real wR[6];
  TestSH_SteadyState_U_ImposedData(x , t, wR);
  ShallowWater_HLL_NumFlux(wL, wR, vnorm, flux);
}

void SteadyState_U_Roe_BoundaryFlux(real *x, real t, real *wL, real *vnorm,real *flux) {
  real wR[6];
  TestSH_SteadyState_U_ImposedData(x , t, wR);
  ShallowWater_Roe_NumFlux(wL, wR, vnorm, flux);
}



void TestSH_SteadyState_P_ImposedData(const real *x, const real t, real *w) {
  real g=_GRAVITY;
  real u01=2.0,u02=2.0; 
  
  w[0] = sqrt(2.0/g*(1.0+(x[0]-x[0]*x[0])*(x[1]-x[1]*x[1])));
  w[1] = u01;//*w[0];
  w[2] = u02;//*w[0];
  w[3] = 0.0;
  w[4] = 0.0;
  w[5] = 0.0;
}

void TestSH_SteadyState_P_InitData(real *x, real *w) {
  real t = 0;
  TestSH_SteadyState_P_ImposedData(x, t, w);
}


void ShallowWater_SteadyState_P_SourceTerm(const real *x, const real t, const real *w, real *source){
  real g=_GRAVITY;
  real S_h1, S_h2;
  real wexact[6];
  real u01=2.0,u02=2.0; 

  S_h1=1.0/sqrt(2.0*g*(1.0+(x[0]-x[0]*x[0])*(x[1]-x[1]*x[1])));
  S_h2=u01*(1-2.0*x[0])*(x[1]-x[1]*x[1])+u02*(1-2.0*x[1])*(x[0]-x[0]*x[0]);
  
  source[0]= S_h1*S_h2;
  source[1]= (1-2.0*x[0])*(x[1]-x[1]*x[1])+u01*source[0];
  source[2]= (1-2.0*x[1])*(x[0]-x[0]*x[0])+u02*source[0];
  source[3]= 0.0;
  source[4]= 0.0;
  source[5]= 0.0;
};

void SteadyState_P_Rusanov_BoundaryFlux(real *x, real t, real *wL, real *vnorm,real *flux) {
  real wR[6];
  TestSH_SteadyState_P_ImposedData(x , t, wR);
  ShallowWater_Rusanov_NumFlux(wL, wR, vnorm, flux);
}

void SteadyState_P_HLL_BoundaryFlux(real *x, real t, real *wL, real *vnorm,real *flux) {
  real wR[6];
  TestSH_SteadyState_P_ImposedData(x , t, wR);
  ShallowWater_HLL_NumFlux(wL, wR, vnorm, flux);
}


void SteadyState_P_Roe_BoundaryFlux(real *x, real t, real *wL, real *vnorm, real *flux) {
  real wR[6];
  TestSH_SteadyState_P_ImposedData(x , t, wR);
  ShallowWater_Roe_NumFlux(wL, wR, vnorm, flux);
}





