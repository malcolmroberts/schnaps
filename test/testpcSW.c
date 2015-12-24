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
#include "pcSW.h"

void SteadyStateOne_ImposedData(const real *x, const real t, real *w);
void SteadyStateOne_InitData(real *x, real *w);
void SteadyStateOne_Source(const real *xy, const real t, const real *w, real *S);
void SteadyStateOne_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
    real *flux);

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
  int test1_ok=1,test2_ok=1,test3_ok=1,test4_ok=1;

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
    int raf[]={4, 4, 1};

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
    simu.dt = 0.001/4;//Get_Dt_RK(&simu)/16;
    real tmax=0.001;
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
    csSolve.type_bc=1;
    
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
      //for (int i=0; i<size1varCG; i++){
      //  printf("i=%d, h=%.8e, u=%.8e, v=%.8e\n",i,csSolve.lsol.rhs[3*i],csSolve.lsol.rhs[3*i+1],csSolve.lsol.rhs[3*i+2]);
      //}
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

    test = test && (dd<1.e-5);
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
    int raf2[]={4, 4, 1};

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
    simu2.dt = 0.001/4;//Get_Dt_RK(&simu2)/16;
    real tmax= 0.001;//4*simu2.dt;//0.002;
    int itermax=tmax/simu2.dt;
    simu2.itermax_rk=itermax;

    csSolve2.lsol.solver_type=LU;//GMRES;
    csSolve2.lsol.tol=1.e-14;
    csSolve2.lsol.pc_type=JACOBI;//PHY_BASED_U2;
    csSolve2.lsol.iter_max=2000;
    csSolve2.lsol.restart_gmres=20;
    csSolve2.lsol.is_CG=true;

    csSolve2.bc_flux=Wave_BC_pressure_null;
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

    test = test && (dd<1.e-5);
    freeSimulation(&simu2);
  }

  printf("//////////////////////////////////////\n");
  printf("TEST THREE: SW U ASF !\n");
  printf("//////////////////////////////////////\n");
  if(test3_ok==1){
    Model model3;

    model3.m = 6; 
    model3.NumFlux=ShallowWater_Rusanov_NumFlux;
    model3.InitData = TestSH_SteadyState_U_InitData;
    model3.ImposedData = TestSH_SteadyState_U_ImposedData;
    model3.BoundaryFlux = SteadyState_U_Rusanov_BoundaryFlux;
    model3.Source = ShallowWater_SteadyState_U_SourceTerm;

    int deg3[]={4, 4, 0};
    int raf3[]={4, 4, 1};

    assert(mesh.is2d);
    CheckMacroMesh(&mesh, deg3, raf3);
    Simulation simu3;
    EmptySimulation(&simu3);
    InitSimulation(&simu3, &mesh, deg3, raf3, &model3);

    ContinuousSolver csSolve3;
    int nbvar=3;
    int *listvar=calloc(nbvar,sizeof(int));
    listvar[0]=0;
    listvar[1]=1;
    listvar[2]=2;
    InitContinuousSolver(&csSolve3,&simu3,1,nbvar,listvar);

    real theta=0.5;
    simu3.theta=theta;
    simu3.dt=0.001/4;
    simu3.vmax=_SPEED_WAVE;
    real tmax=0.001;//1*simu3.dt;
    int itermax=tmax/simu3.dt;
    simu3.itermax_rk=itermax;
    
    csSolve3.lsol.solver_type=LU;//GMRES;
    csSolve3.lsol.tol=1.e-8;
    csSolve3.lsol.pc_type=NONE;//PHY_BASED_U2;
    csSolve3.lsol.iter_max=2000;
    csSolve3.lsol.restart_gmres=10;
    csSolve3.lsol.is_CG=true;

    csSolve3.bc_flux=Wave_BC_normalvelocity_null;
    csSolve3.bc_assembly=BoundaryConditionFriedrichsAssembly;
    csSolve3.rhs_assembly=Source_Assembly;
    csSolve3.type_bc=2;

    int size1varDG = csSolve3.nb_dg_nodes;
    int size1varCG = csSolve3.nb_fe_nodes;
    int sizeDG = csSolve3.nb_dg_dof;
    int sizeCG = csSolve3.nb_fe_dof;

    real *wDG = calloc(sizeDG, sizeof(real));
    real *wCG = calloc(sizeCG, sizeof(real));
    real *resCG = calloc(sizeCG, sizeof(real));

    PB_PC pb_pc;    
    int mat2assemble[6] = {1, 1, 1, 1, 1, 1};
    Init_PBPC_SW_SchurVelocity_BCVelocity(&simu3, &pb_pc, mat2assemble);
    Init_Parameters_PhyBasedPC(&pb_pc);
  
    real h=simu3.vmax*simu3.dt*simu3.theta;
    csSolve3.FluxMatrix = calloc(csSolve3.nb_phy_vars,sizeof(real));
    for (int i=0; i<csSolve3.nb_phy_vars; i++){
      csSolve3.FluxMatrix[i] = calloc(csSolve3.nb_phy_vars,sizeof(real));
    }
    for (int i=0; i<csSolve3.nb_phy_vars; i++){
      for (int j=0; j<csSolve3.nb_phy_vars; j++){
        csSolve3.FluxMatrix[i][j] = 0.0;
      }
    }

    csSolve3.FluxMatrix[0][1]=h;
    csSolve3.FluxMatrix[0][2]=h;
    csSolve3.FluxMatrix[1][0]=h;
    csSolve3.FluxMatrix[2][0]=h;
    csSolve3.lsol.mat_is_assembly=true;

    PiDgToCg(&csSolve3,model3.m,simu3.w,wCG);

    //for(int i=0;i<size1varDG;i++){
    //  printf("i=%d, h=%.8e, u=%.8e, v=%.8e, a=%.8e, b=%8.e, c=%8.e\n",i,simu3.w[i*6],simu3.w[i*6+1],simu3.w[i*6+2],simu3.w[i*6+3],simu3.w[i*6+4],simu3.w[i*6+5]);
    //}
    simu3.tnow=0;
    for(int ie=0; ie < simu3.macromesh.nbelems; ++ie){
      simu3.fd[ie].tnow=simu3.tnow;
    } 

    for(int tstep=0;tstep<simu3.itermax_rk;tstep++){

      csSolve3.reset_dt(&csSolve3.lsol);

      simu3.tnow=simu3.tnow+simu3.dt;
      for(int ie=0; ie < simu3.macromesh.nbelems; ++ie){
        simu3.fd[ie].tnow=simu3.tnow;
      }
      for (int i=0; i<sizeCG; i++){
        csSolve3.lsol.sol[i] = wCG[i];
      }

      pb_pc.rhs_assembly(&pb_pc,&csSolve3);
      pb_pc.source_assembly(&csSolve3);
      pb_pc.bc_assembly(&csSolve3);

      //for (int i=0; i<size1varCG; i++){
      //  printf("i=%d, h=%.8e, u=%.8e, v=%.8e\n",i,csSolve3.lsol.rhs[3*i],csSolve3.lsol.rhs[3*i+1],csSolve3.lsol.rhs[3*i+2]);
      //}

      pb_pc.solvePC(&pb_pc,&simu3,csSolve3.lsol.sol,csSolve3.lsol.rhs);

      for (int i=0; i<sizeCG; i++){
        wCG[i] += csSolve3.lsol.sol[i];
      }
      //for (int i=0; i<size1varCG; i++){
      //  printf("i=%d, dh=%.8e, du=%.8e, dv=%.8e\n",i,csSolve3.lsol.sol[3*i],csSolve3.lsol.sol[3*i+1],csSolve3.lsol.sol[3*i+2]);
      //}

      int freq = (1 >= simu3.itermax_rk / 10)? 1 : simu3.itermax_rk / 10;
      if (tstep % freq == 0)
        printf("t=%f iter=%d/%d dt=%f\n", simu3.tnow, tstep, simu3.itermax_rk, simu3.dt);
    }

    PiInvertCgToDg(&csSolve3,model3.m,wCG,simu3.w);
    dd = L2error(&simu3);
    //for(int i=0;i<size1varDG;i++){
    //  printf("i=%d, h=%.8e, u=%.8e, v=%.8e, a=%.8e, b=%8.e, c=%8.e\n",i,simu3.w[i*6],simu3.w[i*6+1],simu3.w[i*6+2],simu3.w[i*6+3],simu3.w[i*6+4],simu3.w[i*6+5]);
    //}

    printf("erreur L2=%.12e\n", dd);
    dd = L2error_onefield(&simu3,0);
    real dd2 = L2error_onefield(&simu3,1);
    real dd3 = L2error_onefield(&simu3,2);
    printf("erreur L2 sur h=%.12e\n", dd);
    printf("erreur L2 sur u=%.12e\n", dd2);
    printf("erreur L2 sur v=%.12e\n", dd3);

    test = test && (dd<1.e-5);
    freeSimulation(&simu3);
  }

  printf("//////////////////////////////////////\n");
  printf("TEST FOUR: SW P ASF !\n");
  printf("//////////////////////////////////////\n");
  if(test4_ok==1){
    Model model4;

    model4.m = 6; 
    model4.NumFlux=ShallowWater_Rusanov_NumFlux;
    model4.InitData = TestSH_SteadyState_P_InitData;
    model4.ImposedData = TestSH_SteadyState_P_ImposedData;
    model4.BoundaryFlux = SteadyState_P_Rusanov_BoundaryFlux;
    model4.Source = ShallowWater_SteadyState_P_SourceTerm;

    int deg4[]={4, 4, 0};
    int raf4[]={4, 4, 1};

    assert(mesh.is2d);
    CheckMacroMesh(&mesh, deg4, raf4);
    Simulation simu4;
    EmptySimulation(&simu4);
    InitSimulation(&simu4, &mesh, deg4, raf4, &model4);

    ContinuousSolver csSolve4;
    int nbvar=3;
    int *listvar=calloc(nbvar,sizeof(int));
    listvar[0]=0;
    listvar[1]=1;
    listvar[2]=2;
    InitContinuousSolver(&csSolve4,&simu4,1,nbvar,listvar);

    real theta=0.5;
    simu4.theta=theta;
    simu4.dt = 0.001/4;
    simu4.vmax=_SPEED_WAVE;
    real tmax= 0.001;//1*simu4.dt;
    int itermax=tmax/simu4.dt;
    simu4.itermax_rk=itermax;
    
    csSolve4.lsol.solver_type=LU;//GMRES;
    csSolve4.lsol.tol=1.e-8;
    csSolve4.lsol.pc_type=NONE;//PHY_BASED_U2;
    csSolve4.lsol.iter_max=2000;
    csSolve4.lsol.restart_gmres=10;
    csSolve4.lsol.is_CG=true;

    csSolve4.bc_flux=Wave_BC_pressure_null;
    csSolve4.bc_assembly=BoundaryConditionFriedrichsAssembly;
    csSolve4.rhs_assembly=Source_Assembly;
    csSolve4.type_bc=2;

    int size1varDG = csSolve4.nb_dg_nodes;
    int size1varCG = csSolve4.nb_fe_nodes;
    int sizeDG = csSolve4.nb_dg_dof;
    int sizeCG = csSolve4.nb_fe_dof;

    real *wDG = calloc(sizeDG, sizeof(real));
    real *wCG = calloc(sizeCG, sizeof(real));
    real *resCG = calloc(sizeCG, sizeof(real));

    PB_PC pb_pc;    
    int mat2assemble[6] = {1, 1, 1, 1, 1, 1};
    Init_PBPC_SW_SchurVelocity_BCPressure(&simu4, &pb_pc, mat2assemble);
    Init_Parameters_PhyBasedPC(&pb_pc);
  
    real h=simu4.vmax*simu4.dt*simu4.theta;
    csSolve4.FluxMatrix = calloc(csSolve4.nb_phy_vars,sizeof(real));
    for (int i=0; i<csSolve4.nb_phy_vars; i++){
      csSolve4.FluxMatrix[i] = calloc(csSolve4.nb_phy_vars,sizeof(real));
    }
    for (int i=0; i<csSolve4.nb_phy_vars; i++){
      for (int j=0; j<csSolve4.nb_phy_vars; j++){
        csSolve4.FluxMatrix[i][j] = 0.0;
      }
    }

    csSolve4.FluxMatrix[0][1]=h;
    csSolve4.FluxMatrix[0][2]=h;
    csSolve4.FluxMatrix[1][0]=h;
    csSolve4.FluxMatrix[2][0]=h;
    csSolve4.lsol.mat_is_assembly=true;

    PiDgToCg(&csSolve4,model4.m,simu4.w,wCG);

    //for(int i=0;i<size1varDG;i++){
    //  printf("i=%d, h=%.8e, u=%.8e, v=%.8e, a=%.8e, b=%8.e, c=%8.e\n",i,simu4.w[i*6],simu4.w[i*6+1],simu4.w[i*6+2],simu4.w[i*6+3],simu4.w[i*6+4],simu4.w[i*6+5]);
    //}
    simu4.tnow=0;
    for(int ie=0; ie < simu4.macromesh.nbelems; ++ie){
      simu4.fd[ie].tnow=simu4.tnow;
    } 

    for(int tstep=0;tstep<simu4.itermax_rk;tstep++){

      csSolve4.reset_dt(&csSolve4.lsol);

      simu4.tnow=simu4.tnow+simu4.dt;
      for(int ie=0; ie < simu4.macromesh.nbelems; ++ie){
        simu4.fd[ie].tnow=simu4.tnow;
      }
      for (int i=0; i<sizeCG; i++){
        csSolve4.lsol.sol[i] = wCG[i];
      }

      pb_pc.rhs_assembly(&pb_pc,&csSolve4);
      pb_pc.source_assembly(&csSolve4);
      pb_pc.bc_assembly(&csSolve4);

      //for (int i=0; i<size1varCG; i++){
      //  printf("i=%d, h=%.8e, u=%.8e, v=%.8e\n",i,csSolve4.lsol.rhs[3*i],csSolve4.lsol.rhs[3*i+1],csSolve4.lsol.rhs[3*i+2]);
      //}

      pb_pc.solvePC(&pb_pc,&simu4,csSolve4.lsol.sol,csSolve4.lsol.rhs);

      for (int i=0; i<sizeCG; i++){
        wCG[i] += csSolve4.lsol.sol[i];
      }
      //for (int i=0; i<size1varCG; i++){
      //  printf("i=%d, dh=%.8e, du=%.8e, dv=%.8e\n",i,csSolve4.lsol.sol[3*i],csSolve4.lsol.sol[3*i+1],csSolve4.lsol.sol[3*i+2]);
      //}

      int freq = (1 >= simu4.itermax_rk / 10)? 1 : simu4.itermax_rk / 10;
      if (tstep % freq == 0)
        printf("t=%f iter=%d/%d dt=%f\n", simu4.tnow, tstep, simu4.itermax_rk, simu4.dt);
    }

    PiInvertCgToDg(&csSolve4,model4.m,wCG,simu4.w);
    dd = L2error(&simu4);
    //for(int i=0;i<size1varDG;i++){
    //  printf("i=%d, h=%.8e, u=%.8e, v=%.8e, a=%.8e, b=%8.e, c=%8.e\n",i,simu4.w[i*6],simu4.w[i*6+1],simu4.w[i*6+2],simu4.w[i*6+3],simu4.w[i*6+4],simu4.w[i*6+5]);
    //}

    printf("erreur L2=%.12e\n", dd);
    dd = L2error_onefield(&simu4,0);
    real dd2 = L2error_onefield(&simu4,1);
    real dd3 = L2error_onefield(&simu4,2);
    printf("erreur L2 sur h=%.12e\n", dd);
    printf("erreur L2 sur u=%.12e\n", dd2);
    printf("erreur L2 sur v=%.12e\n", dd3);

    test = test && (dd<1.e-5);
    freeSimulation(&simu4);
  }

#ifdef PARALUTION 
  paralution_end();
#endif 

  return test;
}


void TestSH_SteadyState_U_ImposedData(const real *x, const real t, real *w) {
  real alpha = 1.0; 
  real p = 3;
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
  real alpha = 1.0;
  real p = 3;
  real S_11, S_12,S_13, S_21, S_22, S_23, S_factor;
  real wexact[6];

  real h  = 1.0+alpha*(pow(x[0],p)-pow(x[0],2*p))*(pow(x[1],p)-pow(x[1],2*p));
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

  real h  = sqrt(2.0/g*(1.0+(x[0]-x[0]*x[0])*(x[1]-x[1]*x[1])));
  real hx = 0.5*(1.0-2.0*x[0])*(x[1]-x[1]*x[1])*sqrt(2.0/(g*(1.0+(x[0]-x[0]*x[0])*(x[1]-x[1]*x[1])))); 
  real hy = 0.5*(1.0-2.0*x[1])*(x[0]-x[0]*x[0])*sqrt(2.0/(g*(1.0+(x[0]-x[0]*x[0])*(x[1]-x[1]*x[1])))); 
  real u  = u01;
  real ux = 0.0;
  real uy = 0.0;
  real v  = u02;
  real vx = 0.0;
  real vy = 0.0;

  source[0]= h*(ux+vy) + (hx*u+hy*v);
  source[1]= h*(u*ux+v*uy) + g*h*hx;
  source[2]= h*(u*vx+v*vy) + g*h*hy;

  //S_h1=1.0/sqrt(2.0*g*(1.0+(x[0]-x[0]*x[0])*(x[1]-x[1]*x[1])));
  //S_h2=u01*(1-2.0*x[0])*(x[1]-x[1]*x[1])+u02*(1-2.0*x[1])*(x[0]-x[0]*x[0]);
  //
  //source[0]= S_h1*S_h2;
  //source[1]= (1-2.0*x[0])*(x[1]-x[1]*x[1]);//+u01*source[0];
  //source[2]= (1-2.0*x[1])*(x[0]-x[0]*x[0]);//+u02*source[0];
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


void SteadyStateOne_ImposedData(const real *xy, const real t, real *w) {

  real x=xy[0];
  real y=xy[1];

  w[0] = 10.0+exp(x)+exp(2*y);//+x*x+y*y*y;//10
  w[1] = x*y;
  w[2] = -y*y*0.5;

}


void SteadyStateOne_Source(const real *xy, const real t, const real *w, real *S){

  real x=xy[0];
  real y=xy[1];

  S[0] = 0;
  S[1] = exp(x);//2.0*x;//exp(x);
  S[2] = 2.0*exp(2.0*y);//3.0*y*y;//2*exp(2*y);

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



