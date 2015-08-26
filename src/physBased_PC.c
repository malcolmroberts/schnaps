#include "solverpoisson.h"
#include "solvercontinuous.h"
#include "geometry.h"
#include "quantities_vp.h"
#include "linear_solver.h"
#include "physBased_PC.h"
#include <unistd.h>

//void Wave_operators(void* pb_pc, int offset);
//void SW_operators(void* pb_pc, int offset);
//
//// \brief Copy given matrices inside pb_pc as initialization.
//void mat_assembly_PC(PB_PC* pb_pc, real DMat[4][4], real L1Mat[4][4], real L2Mat[4][4], real U1Mat[4][4], real U2Mat[4][4], real SchurMat[4][4][4]){
//  for (int i=0; i<4; i++){
//    for (int j=0; j<4; j++){
//      for (int k=0; k<4; k++){
//        if (SchurMat != NULL) pb_pc->Schur.diff_op2vec[k][i][j] = SchurMat[k][i][j];
//      }
//      //printf("i %d,j %d,D %f,L1 %f,L2 %f,U1 %f,U2 %f,Schur %f\n",i,j,DMat[i][j],L1Mat[i][j],L2Mat[i][j],U1Mat[i][j],U2Mat[i][j],SchurMat[0][i][j]);
//      if (DMat != NULL) pb_pc->D.diff_op[i][j]  = DMat[i][j];
//      if (L1Mat != NULL) pb_pc->L1.diff_op[i][j] = L1Mat[i][j];
//      if (L2Mat != NULL) pb_pc->L2.diff_op[i][j] = L2Mat[i][j];
//      if (U1Mat != NULL) pb_pc->U1.diff_op[i][j] = U1Mat[i][j];
//      if (U2Mat != NULL) pb_pc->U2.diff_op[i][j] = U2Mat[i][j];
//    }
//  }
//}
//
//void InitPhy_Wave(Simulation *simu, PB_PC* pb_pc, int* list_mat2assemble){
//
//  // system is linear
//  pb_pc->nonlinear=0;
//  pb_pc->list_mat2assemble = calloc(6, sizeof(int));
//  for (int i=0; i<6; i++) pb_pc->list_mat2assemble[i] = list_mat2assemble[i];
//  
//  // Initialization variables. Modify at will.
//  int nb_varD = 1;
//  int * listvarD = calloc(nb_varD,sizeof(int));
//  listvarD[0]=0;
//  int nb_varL1 = 1;
//  int * listvarL1 = calloc(nb_varL1,sizeof(int));
//  listvarL1[0]=0;
//  int nb_varL2 = 1;
//  int * listvarL2 = calloc(nb_varL2,sizeof(int));
//  listvarL2[0]=0;
//  int nb_varU1 = 1;
//  int * listvarU1 = calloc(nb_varU1,sizeof(int));
//  listvarU1[0]=1;
//  int nb_varU2 = 1;
//  int * listvarU2 = calloc(nb_varU2,sizeof(int));
//  listvarU2[0]=2;
//  int nb_varSchur = 2;
//  int * listvarSchur = calloc(nb_varSchur,sizeof(int));
//  listvarSchur[0]=1;
//  listvarSchur[1]=2;
//
//  // Initializing all solvers 
//  InitContinuousSolver(&pb_pc->D,simu,1,nb_varD,listvarD);
//  free(listvarD);
//  InitContinuousSolver(&pb_pc->L1,simu,1,nb_varL1,listvarL1);
//  free(listvarL1);
//  InitContinuousSolver(&pb_pc->L2,simu,1,nb_varL2,listvarL2);
//  free(listvarL2);
//  InitContinuousSolver(&pb_pc->U1,simu,1,nb_varU1,listvarU1);
//  free(listvarU1);
//  InitContinuousSolver(&pb_pc->U2,simu,1,nb_varU2,listvarU2);
//  free(listvarU2);
//  InitContinuousSolver(&pb_pc->Schur,simu,1,nb_varSchur,listvarSchur);
//  free(listvarSchur);
//
//  pb_pc->D.lsol.MatVecProduct=MatVect;
//  pb_pc->L1.lsol.MatVecProduct=MatVect;
//  pb_pc->L2.lsol.MatVecProduct=MatVect;
//  pb_pc->U1.lsol.MatVecProduct=MatVect;
//  pb_pc->U2.lsol.MatVecProduct=MatVect;
//  pb_pc->Schur.lsol.MatVecProduct=MatVect;
//  pb_pc->D.bc_assembly=ExactDirichletContinuousMatrix;
//  pb_pc->U1.bc_assembly=ExactDirichletContinuousMatrix;
//  pb_pc->U2.bc_assembly=ExactDirichletContinuousMatrix;
//  pb_pc->Schur.bc_assembly=ExactDirichletContinuousMatrix;
//
//  pb_pc->solver_prediction=LU;
//  pb_pc->solver_propagation=LU;
//  pb_pc->solver_correction=LU;
//
//  pb_pc->pc_prediction=NONE;
//  pb_pc->pc_propagation=NONE;
//  pb_pc->pc_correction=NONE;
//
//  pb_pc->tol_prediction=1.e-9;
//  pb_pc->tol_propagation=1.e-9;
//  pb_pc->tol_correction=1.e-9;
//
//  pb_pc->itermax_prediction=1000;
//  pb_pc->itermax_propagation=1000;
//  pb_pc->itermax_correction=1000;
//
//  pb_pc->restart_prediction=10;
//  pb_pc->restart_propagation=10;
//  pb_pc->restart_correction=10;
//
//  // Allocating all sub-equations' right-hand sides
//  pb_pc->rhs_prediction = calloc(pb_pc->D.lsol.neq,sizeof(real));
//  pb_pc->rhs_propagation = calloc(pb_pc->Schur.lsol.neq,sizeof(real));
//  pb_pc->rhs_correction = calloc(pb_pc->D.lsol.neq,sizeof(real));
//
//  // Associating problem matrices
//  pb_pc->mat_assembly = Wave_operators;
//  
//  // Since the problem is linear, problem's matrices are all independent on time.
//  // Thus, we build them immediately in the initialization
//  pb_pc->mat_assembly(pb_pc,0);
//  GenericOperator(pb_pc);
//
//}
//
//void InitPhy_SW(Simulation *simu, PB_PC* pb_pc, int* list_mat2assemble){
//
//  // system is non linear
//  pb_pc->nonlinear=1;
//  pb_pc->list_mat2assemble = calloc(6, sizeof(int));
//  for (int i=0; i<6; i++) pb_pc->list_mat2assemble[i] = list_mat2assemble[i];
//  
//  // Initialization variables. Modify at will.
//  int nb_varD = 1;
//  int * listvarD = calloc(nb_varD,sizeof(int));
//  listvarD[0]=0;
//  int nb_varL1 = 1;
//  int * listvarL1 = calloc(nb_varL1,sizeof(int));
//  listvarL1[0]=0;
//  int nb_varL2 = 1;
//  int * listvarL2 = calloc(nb_varL2,sizeof(int));
//  listvarL2[0]=0;
//  int nb_varU1 = 1;
//  int * listvarU1 = calloc(nb_varU1,sizeof(int));
//  listvarU1[0]=0;
//  int nb_varU2 = 1;
//  int * listvarU2 = calloc(nb_varU2,sizeof(int));
//  listvarU2[0]=1;
//  int nb_varSchur = 2;
//  int * listvarSchur = calloc(nb_varSchur,sizeof(int));
//  listvarSchur[0]=1;
//  listvarSchur[1]=2;
//
//  // Initializing all solvers 
//  InitContinuousSolver(&pb_pc->D,simu,1,nb_varD,listvarD);
//  free(listvarD);
//  InitContinuousSolver(&pb_pc->L1,simu,1,nb_varL1,listvarL1);
//  free(listvarL1);
//  InitContinuousSolver(&pb_pc->L2,simu,1,nb_varL2,listvarL2);
//  free(listvarL2);
//  InitContinuousSolver(&pb_pc->U1,simu,1,nb_varU1,listvarU1);
//  free(listvarU1);
//  InitContinuousSolver(&pb_pc->U2,simu,1,nb_varU2,listvarU2);
//  free(listvarU2);
//  InitContinuousSolver(&pb_pc->Schur,simu,1,nb_varSchur,listvarSchur);
//  free(listvarSchur);
//
//  pb_pc->D.lsol.MatVecProduct=MatVect;
//  pb_pc->L1.lsol.MatVecProduct=MatVect;
//  pb_pc->L2.lsol.MatVecProduct=MatVect;
//  pb_pc->U1.lsol.MatVecProduct=MatVect;
//  pb_pc->U2.lsol.MatVecProduct=MatVect;
//  pb_pc->Schur.lsol.MatVecProduct=MatVect;
//  pb_pc->D.bc_assembly=ExactDirichletContinuousMatrix;
//  pb_pc->U1.bc_assembly=ExactDirichletContinuousMatrix;
//  pb_pc->U2.bc_assembly=ExactDirichletContinuousMatrix;
//  pb_pc->Schur.bc_assembly=ExactDirichletContinuousMatrix;
//  // Allocating all sub-equations right-hand sides
//  pb_pc->rhs_prediction = malloc(pb_pc->D.lsol.neq*sizeof(real));
//  pb_pc->rhs_propagation = malloc(pb_pc->Schur.lsol.neq*sizeof(real));
//  pb_pc->rhs_correction = malloc(pb_pc->D.lsol.neq*sizeof(real));
//
//  // Associating problem matrices
//  pb_pc->mat_assembly = SW_operators;
//
//  //// Allocating all sub-equations right-hand sides
//  pb_pc->rhs_prediction = malloc(pb_pc->D.lsol.neq*sizeof(real));
//  pb_pc->rhs_propagation = malloc(pb_pc->Schur.lsol.neq*sizeof(real));
//  pb_pc->rhs_correction = malloc(pb_pc->D.lsol.neq*sizeof(real));
//
//}
//
//void InitParameters_PC(Simulation *simu, PB_PC* pb_pc){
//
//  pb_pc->solver_prediction=LU;//PAR_CG;
//  pb_pc->solver_propagation=LU;//PAR_CG;
//  pb_pc->solver_correction=LU;//LU;//PAR_CG;
//
//  pb_pc->pc_prediction=NONE;
//  pb_pc->pc_propagation=NONE;
//  pb_pc->pc_correction=NONE;
//
//  pb_pc->tol_prediction=1.e-9;
//  pb_pc->tol_propagation=1.e-9;
//  pb_pc->tol_correction=1.e-9;
//
//  pb_pc->itermax_prediction=100;
//  pb_pc->itermax_propagation=100;
//  pb_pc->itermax_correction=100;
//
//  pb_pc->restart_prediction=10;
//  pb_pc->restart_propagation=10;
//  pb_pc->restart_correction=10;
//
//}
//
//void Wave_operators(void* pb_pc, int offset){
//
//  PB_PC* PC = pb_pc;
//  real h=PC->D.simu->dt*PC->D.simu->theta*PC->D.simu->vmax;
//  real DMat[4][4]  = {{1,0,0,0},
//                      {0,0,0,0},
//                      {0,0,0,0},
//                      {0,0,0,0}};
//  real L1Mat[4][4] = {{0,h,0,0},
//                      {0,0,0,0},
//                      {0,0,0,0},
//                      {0,0,0,0}};
//  real L2Mat[4][4] = {{0,0,h,0},
//                      {0,0,0,0},
//                      {0,0,0,0},
//                      {0,0,0,0}};
//  real U1Mat[4][4] = {{0,h,0,0},
//                      {0,0,0,0},
//                      {0,0,0,0},
//                      {0,0,0,0}};
//  real U2Mat[4][4] = {{0,0,h,0},
//                      {0,0,0,0},
//                      {0,0,0,0},
//                      {0,0,0,0}};
//  real SchurMat[4][4][4] ={{{1.0,0,0,0},
//                            {0,h*h,0,0},
//                            {0,0,0,0},
//                            {0,0,0,0}},
//                           {{0,0,0,0},
//                            {0,0,h*h,0},
//                            {0,0,0,0},
//                            {0,0,0,0}},
//                           {{0,0,0,0},
//                            {0,0,0,0},
//                            {0,h*h,0,0},
//                            {0,0,0,0}},
//                           {{1.0,0,0,0},
//                            {0,0,0,0},
//                            {0,0,h*h,0},
//                            {0,0,0,0}}};
//
//  // Applying these matrices inside the different ContinuousSolver
//  mat_assembly_PC(PC,DMat,L1Mat,L2Mat,U1Mat,U2Mat,SchurMat);
//
//}
//
//
//void SW_operators(void* pb_pc, int offset){
//  PB_PC* PC = pb_pc;
//  real theta = PC->D.simu->theta;
//  real dt = PC->D.simu->dt;
//  real g = 10.;
//  real h = PC->D.simu->fd[0].wn[offset  ];
//  real u = PC->D.simu->fd[0].wn[offset+1];
//  real v = PC->D.simu->fd[0].wn[offset+2];
//  real DMat[4][4]  = {{1          ,0,0,0},
//                      {-theta*dt*u,0,0,0},
//                      {-theta*dt*v,0,0,0},
//                      {0          ,0,0,0}};
//  real L1Mat[4][4] = {{u,0,0,0},
//                      {-theta*dt*g*h,0,0,0},
//                      {0,0,0,0},
//                      {0,0,0,0}};
//  real L2Mat[4][4] = {{v,0,0,0},
//                      {0,0,0,0},
//                      {-theta*dt*g*h,0,0,0},
//                      {0,0,0,0}};
//  real U1Mat[4][4] = {{0,0,0,0},
//                      {-theta*dt*h,0,0,0},
//                      {0,0,0,0},
//                      {0,0,0,0}};
//  real U2Mat[4][4] = {{0,0,0,0},
//                      {0,0,0,0},
//                      {-theta*dt*h,0,0,0},
//                      {0,0,0,0}};
//  real SchurMat[4][4][4] ={{{h,theta*dt*h*u,theta*dt*h*v,0},
//                            {0,theta*theta*dt*dt*g*h*h,0,0},
//                            {0,0,0,0},
//                            {0,0,0,0}},
//                           {{0,0,0,0},
//                            {0,0,theta*theta*dt*dt*g*h*h,0},
//                            {0,0,0,0},
//                            {0,0,0,0}},
//                           {{0,0,0,0},
//                            {0,0,0,0},
//                            {0,theta*theta*dt*dt*g*h*h,0,0},
//                            {0,0,0,0}},
//                           {{h,theta*dt*h*u,theta*dt*h*v,0},
//                            {0,0,0,0},
//                            {0,0,theta*theta*dt*dt*g*h*h,0},
//                            {0,0,0,0}}};
//  mat_assembly_PC(PC,DMat,L1Mat,L2Mat,U1Mat,U2Mat,SchurMat);
//}

void InitPhy_Wave(Simulation *simu, PB_PC* pb_pc, int* list_mat2assemble){

  // system is linear
  pb_pc->nonlinear=0;


  pb_pc->list_mat2assemble = calloc(6, sizeof(int));
  for (int i=0; i<6; i++) pb_pc->list_mat2assemble[i] = list_mat2assemble[i];
  
  // Initialization variables. Modify at will.
  int nb_varD = 1;
  int * listvarD = calloc(nb_varD,sizeof(int));
  listvarD[0]=0;
  int nb_varL1 = 1;
  int * listvarL1 = calloc(nb_varL1,sizeof(int));
  listvarL1[0]=0;
  int nb_varL2 = 1;
  int * listvarL2 = calloc(nb_varL2,sizeof(int));
  listvarL2[0]=0;
  int nb_varU1 = 1;
  int * listvarU1 = calloc(nb_varU1,sizeof(int));
  listvarU1[0]=0;
  int nb_varU2 = 1;
  int * listvarU2 = calloc(nb_varU2,sizeof(int));
  listvarU2[0]=1;
  int nb_varSchur = 2;
  int * listvarSchur = calloc(nb_varSchur,sizeof(int));
  listvarSchur[0]=1;
  listvarSchur[1]=2;

  // Initializing all solvers 
  InitContinuousSolver(&pb_pc->D,simu,1,nb_varD,listvarD);
  free(listvarD);
  InitContinuousSolver(&pb_pc->L1,simu,1,nb_varL1,listvarL1);
  free(listvarL1);
  InitContinuousSolver(&pb_pc->L2,simu,1,nb_varL2,listvarL2);
  free(listvarL2);
  InitContinuousSolver(&pb_pc->U1,simu,1,nb_varU1,listvarU1);
  free(listvarU1);
  InitContinuousSolver(&pb_pc->U2,simu,1,nb_varU2,listvarU2);
  free(listvarU2);
  InitContinuousSolver(&pb_pc->Schur,simu,1,nb_varSchur,listvarSchur);
  free(listvarSchur);

  // Local operator matrices to build global Differential operators (Schur, Laplacian...)
  real h=simu->dt*simu->theta*simu->vmax;
  real DMat[4][4]  = {{1,0,0,0},
                      {0,0,0,0},
                      {0,0,0,0},
                      {0,0,0,0}};
  real L1Mat[4][4] = {{0,h,0,0},
                      {0,0,0,0},
                      {0,0,0,0},
                      {0,0,0,0}};
  real L2Mat[4][4] = {{0,0,h,0},
                      {0,0,0,0},
                      {0,0,0,0},
                      {0,0,0,0}};
  real U1Mat[4][4] = {{0,h,0,0},
                      {0,0,0,0},
                      {0,0,0,0},
                      {0,0,0,0}};
  real U2Mat[4][4] = {{0,0,h,0},
                      {0,0,0,0},
                      {0,0,0,0},
                      {0,0,0,0}};
  real SchurMat[4][4][4] ={{{1,0,0,0},
                            {0,h*h,0,0},
                            {0,0,0,0},
                            {0,0,0,0}},
                           {{0,0,0,0},
                            {0,0,h*h,0},
                            {0,0,0,0},
                            {0,0,0,0}},
                           {{0,0,0,0},
                            {0,0,0,0},
                            {0,h*h,0,0},
                            {0,0,0,0}},
                           {{1,0,0,0},
                            {0,0,0,0},
                            {0,0,h*h,0},
                            {0,0,0,0}}};

  // Applying these matrices inside the different ContinuousSolver
  InitMat_ContinuousSolver(pb_pc,DMat,L1Mat,L2Mat,U1Mat,U2Mat,SchurMat);

  // Assembling all operators' matrices
  GenericOperator(pb_pc);

  //// Allocating all sub-equations right-hand sides
  pb_pc->rhs_prediction = malloc(pb_pc->D.lsol.neq*sizeof(real));
  pb_pc->rhs_propagation = malloc(pb_pc->Schur.lsol.neq*sizeof(real));
  pb_pc->rhs_correction = malloc(pb_pc->D.lsol.neq*sizeof(real));

}

// \brief Copy given matrices inside pb_pc as initialization.
void InitMat_ContinuousSolver(PB_PC* pb_pc, real DMat[4][4], real L1Mat[4][4], real L2Mat[4][4], real U1Mat[4][4], real U2Mat[4][4], real SchurMat[4][4][4]){
  for (int i=0; i<4; i++){
    for (int j=0; j<4; j++){
      for (int k=0; k<4; k++){
        if (SchurMat != NULL) pb_pc->Schur.diff_op2vec[k][i][j] = SchurMat[k][i][j];
      }
      if (DMat != NULL) pb_pc->D.diff_op[i][j]  = DMat[i][j];
      if (L1Mat != NULL) pb_pc->L1.diff_op[i][j] = L1Mat[i][j];
      if (L2Mat != NULL) pb_pc->L2.diff_op[i][j] = L2Mat[i][j];
      if (U1Mat != NULL) pb_pc->U1.diff_op[i][j] = U1Mat[i][j];
      if (U2Mat != NULL) pb_pc->U2.diff_op[i][j] = U2Mat[i][j];
    }
  }
}

void solveIdentity(PB_PC* pb_pc, Simulation *simu, real* globalSol, real*globalRHS_DG){
  
  // 0)1) Reset everything (needed for time evolution)
  reset(pb_pc);


  // Parsing globalRHS (in DG) into a CG vector
  ContinuousSolver waveSolver;
  int nb_var = 3;
  int * listvarGlobal = calloc(nb_var, sizeof(int));
  listvarGlobal[0]=0;
  listvarGlobal[1]=1;
  listvarGlobal[2]=2;
  InitContinuousSolver(&waveSolver,simu,1,nb_var,listvarGlobal);
  free(listvarGlobal);
  real * globalRHS_CG = calloc(waveSolver.lsol.neq,sizeof(real));

  //for(int i=0;i<simu->wsize;i++){
  // globalSol[i]=globalRHS_DG[i];
  //}
  
  VectorDgToCg(&waveSolver, globalRHS_DG, globalRHS_CG);
  // Back from CG to DG
  VectorCgToDg(&waveSolver,globalRHS_CG,globalSol);

  /*real error=0;
  real norm=0;
  real normnum=0;
  for(int i=0;i<simu->wsize;i++){
    //printf("  i, dg, sol, diff %d %.7e %.7e %.7e \n",i,globalRHS_DG[i],globalSol[i],globalSol[i]-globalRHS_DG[i]);
    error=error+fabs(globalSol[i]-globalRHS_DG[i]);
    norm=norm+fabs(globalRHS_DG[i]);
    normnum=normnum+fabs(globalSol[i]);
    }
  //printf(" error cg to dg:%f %f %f %f\n",error,error/(norm+1.0),norm,normnum);
  //sleep(100);
  */
  freeContinuousSolver(&waveSolver);
  free(globalRHS_CG);
}

void solvePhy(PB_PC* pb_pc, Simulation *simu, real* globalSol, real*globalRHS_DG){
  
  // 0)1) Reset everything (needed for time evolution)
  reset(pb_pc);

  // Assembling all operators' matrices
  if(pb_pc->nonlinear == 1){
    GenericOperator(pb_pc);
  }

  // 0)2) Second, initialize all penalization right-hand-sides
  pb_pc->D.bc_assembly=ExactDirichletContinuousMatrix;
  pb_pc->Schur.bc_assembly=ExactDirichletContinuousMatrix;

  // Parsing globalRHS (in DG) into a CG vector
  ContinuousSolver waveSolver;
  int nb_var = 3;
  int * listvarGlobal = calloc(nb_var, sizeof(int));
  listvarGlobal[0]=0;
  listvarGlobal[1]=1;
  listvarGlobal[2]=2;
  InitContinuousSolver(&waveSolver,simu,1,nb_var,listvarGlobal);
  free(listvarGlobal);
  real * globalRHS_CG = calloc(waveSolver.lsol.neq,sizeof(real));
  VectorDgToCg(&waveSolver, globalRHS_DG, globalRHS_CG);

  // 1) PREDICTION STEP

  //pb_pc->D.lsol.solver_type=pb_pc->solver_prediction;
  //pb_pc->D.lsol.tol=pb_pc->tol_prediction;
  //pb_pc->D.lsol.pc_type=pb_pc->pc_prediction;
  //pb_pc->D.lsol.iter_max=pb_pc->itermax_prediction;
  //pb_pc->D.lsol.restart_gmres=pb_pc->restart_prediction;
  pb_pc->D.lsol.solver_type=LU;
  pb_pc->D.lsol.tol=1.e-9;
  pb_pc->D.lsol.pc_type=NONE;
  pb_pc->D.lsol.iter_max=50;
  pb_pc->D.lsol.MatVecProduct=MatVect;

  //printf("RHS assembly.....\n");
  for (int i=0;i<pb_pc->D.nb_fe_nodes;i++){
    pb_pc->D.lsol.rhs[i] = globalRHS_CG[i*3];
  }
  //printf("Solution...\n");
  pb_pc->D.bc_assembly(&pb_pc->D,&pb_pc->D.lsol);
  SolveLinearSolver(&pb_pc->D.lsol,simu);

  // 2) PROPAGATION STEP

  //pb_pc->Schur.lsol.solver_type=pb_pc->solver_propagation;
  //pb_pc->Schur.lsol.tol=pb_pc->tol_propagation;
  //pb_pc->Schur.lsol.pc_type=pb_pc->pc_propagation;
  //pb_pc->Schur.lsol.iter_max=pb_pc->itermax_propagation;
  //pb_pc->Schur.lsol.restart_gmres=pb_pc->restart_propagation;
  pb_pc->Schur.lsol.pc_type=NONE;
  pb_pc->Schur.lsol.iter_max=2000;
  pb_pc->Schur.lsol.solver_type=LU;
  pb_pc->Schur.lsol.tol=1.e-9;

  pb_pc->Schur.lsol.MatVecProduct=MatVect;
  pb_pc->L1.lsol.MatVecProduct=MatVect;
  pb_pc->L2.lsol.MatVecProduct=MatVect;
  // Parsing L1P, L2P into the "sol" of L1 and L2 (since it is unused).
  pb_pc->L1.lsol.MatVecProduct(&pb_pc->L1.lsol,pb_pc->D.lsol.sol,pb_pc->L1.lsol.sol);
  pb_pc->L2.lsol.MatVecProduct(&pb_pc->L2.lsol,pb_pc->D.lsol.sol,pb_pc->L2.lsol.sol);


  //printf("RHS assembly.....\n");
  for (int i=0;i<pb_pc->D.nb_fe_nodes;i++){
    pb_pc->Schur.lsol.rhs[i*2]   = globalRHS_CG[i*3+1] - pb_pc->L1.lsol.sol[i];
    pb_pc->Schur.lsol.rhs[i*2+1] = globalRHS_CG[i*3+2] - pb_pc->L2.lsol.sol[i];
  }

  //printf("Solution...\n");
  pb_pc->Schur.bc_assembly(&pb_pc->Schur,&pb_pc->Schur.lsol);
  SolveLinearSolver(&pb_pc->Schur.lsol,simu);


  // 3) CORRECTION STEP

  // Extracting both U1 and U2 from the previous solution of the propagation step
  real *solU1=calloc(pb_pc->U1.nb_fe_nodes, sizeof(real));
  real *solU2=calloc(pb_pc->U2.nb_fe_nodes, sizeof(real));
  extract2CGVectors(&pb_pc->U1,&pb_pc->U2,pb_pc->Schur.lsol.sol,solU1,solU2);

  //pb_pc->D.lsol.solver_type=pb_pc->solver_correction;
  //pb_pc->D.lsol.tol=pb_pc->tol_correction;
  //pb_pc->D.lsol.pc_type=pb_pc->pc_correction;
  //pb_pc->D.lsol.iter_max=pb_pc->itermax_correction;
  //pb_pc->D.lsol.restart_gmres=pb_pc->restart_correction;
  pb_pc->D.lsol.solver_type=LU;
  pb_pc->D.lsol.pc_type=NONE;
  pb_pc->D.lsol.iter_max=50;
  pb_pc->D.lsol.tol=1.e-9;

  pb_pc->D.lsol.MatVecProduct=MatVect;
  pb_pc->U1.lsol.MatVecProduct=MatVect;
  pb_pc->U2.lsol.MatVecProduct=MatVect;
  // Parsing U1u1, U2u2 into the "sol" of U1 and U2 (since it is unused).

  pb_pc->L1.lsol.MatVecProduct(&pb_pc->L1.lsol,solU1,pb_pc->L1.lsol.sol);
  pb_pc->L2.lsol.MatVecProduct(&pb_pc->L2.lsol,solU2,pb_pc->L2.lsol.sol);
  pb_pc->D.lsol.MatVecProduct(&pb_pc->D.lsol,pb_pc->D.lsol.sol,pb_pc->D.lsol.rhs);
  //pb_pc->U1.lsol.MatVecProduct(&pb_pc->U1.lsol,solU1,pb_pc->U1.lsol.sol);
  //pb_pc->U2.lsol.MatVecProduct(&pb_pc->U2.lsol,solU2,pb_pc->U2.lsol.sol);

  //printf("RHS assembly.....\n");
  for (int i=0;i<pb_pc->D.nb_fe_nodes;i++){
    pb_pc->D.lsol.rhs[i] += - pb_pc->L1.lsol.sol[i] - pb_pc->L2.lsol.sol[i];
  }

  //printf("Solution...\n");
  pb_pc->D.bc_assembly(&pb_pc->D,&pb_pc->D.lsol);
  SolveLinearSolver(&pb_pc->D.lsol,simu);

  // 4) OUTPUT STEP

  // Final concatenation
  cat2CGVectors(&pb_pc->D,&pb_pc->Schur,pb_pc->D.lsol.sol,pb_pc->Schur.lsol.sol,waveSolver.lsol.sol);

  // Back from CG to DG
  VectorCgToDg(&waveSolver,waveSolver.lsol.sol,globalSol);

  freeContinuousSolver(&waveSolver);
  free(globalRHS_CG);
  free(solU1);
  free(solU2);
}

void solvePhy_CG(PB_PC* pb_pc, Simulation *simu, real* globalSol, real*globalRHS){
  
  // 0)1) Reset everything (needed for time evolution)
  reset(pb_pc);

  // Assembling all operators' matrices
  if(pb_pc->nonlinear == 1){
    GenericOperator(pb_pc);
  }

  // 0)2) Second, initialize all penalization right-hand-sides
  pb_pc->D.bc_assembly=ExactDirichletContinuousMatrix;
  pb_pc->Schur.bc_assembly=ExactDirichletContinuousMatrix;

  // Parsing globalRHS (in DG) into a CG vector
  ContinuousSolver waveSolver;
  int nb_var = 3;
  int * listvarGlobal = calloc(nb_var, sizeof(int));
  listvarGlobal[0]=0;
  listvarGlobal[1]=1;
  listvarGlobal[2]=2;
  InitContinuousSolver(&waveSolver,simu,1,nb_var,listvarGlobal);
  free(listvarGlobal);

  // 1) PREDICTION STEP

  //pb_pc->D.lsol.solver_type=pb_pc->solver_prediction;
  //pb_pc->D.lsol.tol=pb_pc->tol_prediction;
  //pb_pc->D.lsol.pc_type=pb_pc->pc_prediction;
  //pb_pc->D.lsol.iter_max=pb_pc->itermax_prediction;
  //pb_pc->D.lsol.restart_gmres=pb_pc->restart_prediction;
  pb_pc->D.lsol.solver_type=LU;
  pb_pc->D.lsol.tol=1.e-9;
  pb_pc->D.lsol.pc_type=NONE;
  pb_pc->D.lsol.iter_max=50;
  pb_pc->D.lsol.MatVecProduct=MatVect;

  //printf("RHS assembly.....\n");
  for (int i=0;i<pb_pc->D.nb_fe_nodes;i++){
    pb_pc->D.lsol.rhs[i] = globalRHS[i*3];
  }
  //printf("Solution...\n");
  pb_pc->D.bc_assembly(&pb_pc->D,&pb_pc->D.lsol);
  SolveLinearSolver(&pb_pc->D.lsol,simu);

  // 2) PROPAGATION STEP

  //pb_pc->Schur.lsol.solver_type=pb_pc->solver_propagation;
  //pb_pc->Schur.lsol.tol=pb_pc->tol_propagation;
  //pb_pc->Schur.lsol.pc_type=pb_pc->pc_propagation;
  //pb_pc->Schur.lsol.iter_max=pb_pc->itermax_propagation;
  //pb_pc->Schur.lsol.restart_gmres=pb_pc->restart_propagation;
  pb_pc->Schur.lsol.pc_type=NONE;
  pb_pc->Schur.lsol.iter_max=2000;
  pb_pc->Schur.lsol.solver_type=LU;
  pb_pc->Schur.lsol.tol=1.e-9;

  pb_pc->Schur.lsol.MatVecProduct=MatVect;
  pb_pc->L1.lsol.MatVecProduct=MatVect;
  pb_pc->L2.lsol.MatVecProduct=MatVect;
  // Parsing L1P, L2P into the "sol" of L1 and L2 (since it is unused).
  pb_pc->L1.lsol.MatVecProduct(&pb_pc->L1.lsol,pb_pc->D.lsol.sol,pb_pc->L1.lsol.sol);
  pb_pc->L2.lsol.MatVecProduct(&pb_pc->L2.lsol,pb_pc->D.lsol.sol,pb_pc->L2.lsol.sol);


  //printf("RHS assembly.....\n");
  for (int i=0;i<pb_pc->D.nb_fe_nodes;i++){
    pb_pc->Schur.lsol.rhs[i*2]   = globalRHS[i*3+1] - pb_pc->L1.lsol.sol[i];
    pb_pc->Schur.lsol.rhs[i*2+1] = globalRHS[i*3+2] - pb_pc->L2.lsol.sol[i];
  }

  //printf("Solution...\n");
  pb_pc->Schur.bc_assembly(&pb_pc->Schur,&pb_pc->Schur.lsol);
  SolveLinearSolver(&pb_pc->Schur.lsol,simu);


  // 3) CORRECTION STEP

  // Extracting both U1 and U2 from the previous solution of the propagation step
  real *solU1=calloc(pb_pc->U1.nb_fe_nodes, sizeof(real));
  real *solU2=calloc(pb_pc->U2.nb_fe_nodes, sizeof(real));
  extract2CGVectors(&pb_pc->U1,&pb_pc->U2,pb_pc->Schur.lsol.sol,solU1,solU2);

  //pb_pc->D.lsol.solver_type=pb_pc->solver_correction;
  //pb_pc->D.lsol.tol=pb_pc->tol_correction;
  //pb_pc->D.lsol.pc_type=pb_pc->pc_correction;
  //pb_pc->D.lsol.iter_max=pb_pc->itermax_correction;
  //pb_pc->D.lsol.restart_gmres=pb_pc->restart_correction;
  pb_pc->D.lsol.solver_type=LU;
  pb_pc->D.lsol.pc_type=NONE;
  pb_pc->D.lsol.iter_max=50;
  pb_pc->D.lsol.tol=1.e-9;

  pb_pc->D.lsol.MatVecProduct=MatVect;
  pb_pc->U1.lsol.MatVecProduct=MatVect;
  pb_pc->U2.lsol.MatVecProduct=MatVect;
  // Parsing U1u1, U2u2 into the "sol" of U1 and U2 (since it is unused).

  pb_pc->L1.lsol.MatVecProduct(&pb_pc->L1.lsol,solU1,pb_pc->L1.lsol.sol);
  pb_pc->L2.lsol.MatVecProduct(&pb_pc->L2.lsol,solU2,pb_pc->L2.lsol.sol);
  pb_pc->D.lsol.MatVecProduct(&pb_pc->D.lsol,pb_pc->D.lsol.sol,pb_pc->D.lsol.rhs);
  //pb_pc->U1.lsol.MatVecProduct(&pb_pc->U1.lsol,solU1,pb_pc->U1.lsol.sol);
  //pb_pc->U2.lsol.MatVecProduct(&pb_pc->U2.lsol,solU2,pb_pc->U2.lsol.sol);

  //printf("RHS assembly.....\n");
  //for (int i=0;i<pb_pc->D.nb_fe_nodes;i++){
  //  pb_pc->D.lsol.rhs[i] += - pb_pc->L1.lsol.sol[i] - pb_pc->L2.lsol.sol[i];
  //}

  //printf("Solution...\n");
  pb_pc->D.bc_assembly(&pb_pc->D,&pb_pc->D.lsol);
  SolveLinearSolver(&pb_pc->D.lsol,simu);

  // 4) OUTPUT STEP

  // Final concatenation
  cat2CGVectors(&pb_pc->D,&pb_pc->Schur,pb_pc->D.lsol.sol,pb_pc->Schur.lsol.sol,waveSolver.lsol.sol);

  // Back from CG to DG
  for (int i=0;i<waveSolver.nb_fe_dof;i++){
    globalSol[i]=waveSolver.lsol.sol[i];
  }

  freeContinuousSolver(&waveSolver);
  free(solU1);
  free(solU2);
}

//void NewVectorDgToCg(ContinuousSolver * cs,real * rhsIn, real * rhsOut){
//  
//  field* f = &cs->simu->fd[0];
//  real* coeff = calloc(cs->nb_fe_nodes, sizeof(real));
//  real* rhsCopy = calloc(cs->nb_dg_nodes*f->model.m, sizeof(real));
//
//  for (int i=0;i<cs->nb_dg_nodes*f->model.m; i++) rhsCopy[i]=rhsIn[i];
//  
//  // right hand side assembly
//  for(int ino = 0; ino < cs->nb_fe_dof; ino++){
//    rhsOut[ino] = 0;
//  }
//
//  for(int ie = 0; ie < cs->nbel; ie++){  
//
//    int iemacro  = ie / (f->raf[0] * f->raf[1] * f->raf[2]);
//    int isubcell = ie % (f->raf[0] * f->raf[1] * f->raf[2]);
//    
// 
//    for(int iloc = 0; iloc < cs->nnodes; iloc++){
//      real wpg;
//      real xref[3];
//      int ilocmacro = iloc + isubcell * cs->nnodes;
//      ref_pg_vol(f->deg,f->raf,ilocmacro,xref,&wpg,NULL);
//      real dtau[3][3],codtau[3][3];
//      Ref2Phy(cs->simu->fd[iemacro].physnode,
//	      xref,NULL,0,NULL,
//	      dtau,codtau,NULL,NULL);
//      real det = dot_product(dtau[0], codtau[0]);	
//      int ino_dg = iloc + ie * cs->nnodes;
//      int ino_fe = cs->dg_to_fe_index[ino_dg];
//      coeff[ino_fe] += wpg * det; // We need to multiply by the member of dg node for 1 fe node.
//      for (int iv=0; iv<cs->nb_phy_vars;iv++){ 
//        int imem = f->varindex(f->deg,f->raf,f->model.m,
//		          ilocmacro,cs->list_of_var[iv]) ;
//
//        rhsOut[iv + cs->nb_phy_vars * ino_fe] += rhsCopy[imem] * wpg * det;//rhs[imem] * wpg * det; 
//      }
//    }
//  }
//  for (int i=0; i<cs->nb_fe_nodes;i++){
//    for (int iv=0; iv<cs->nb_phy_vars;iv++){
//      rhsOut[iv+cs->nb_phy_vars*i]/=coeff[i];//*coeff[i];
//    }
//  }
//  free(coeff);
//  free(rhsCopy);
//}

void VectorDgToCg(ContinuousSolver * cs,real * rhsIn, real * rhsOut){
  
  field* f = &cs->simu->fd[0];
  real* coeff = calloc(cs->nb_fe_nodes, sizeof(real));
  real* rhsCopy = calloc(cs->nb_dg_nodes*f->model.m, sizeof(real));

  for (int i=0;i<cs->nb_dg_nodes*f->model.m; i++) rhsCopy[i]=rhsIn[i];
  
  // Apply division by the discontinuous Garlerkin mass matrix
  int m = f->model.m;
  
  for(int ipg = 0; ipg < NPG(f->deg, f->raf); ipg++) {
    real dtau[3][3], codtau[3][3], xpgref[3], xphy[3], wpg;
    ref_pg_vol(f->deg, f->raf, ipg, xpgref, &wpg, NULL);
    Ref2Phy(f->physnode, // phys. nodes
      xpgref, // xref
      NULL, -1, // dpsiref, ifa
      xphy, dtau, // xphy, dtau
      codtau, NULL, NULL); // codtau, dpsi, vnds
    real det = dot_product(dtau[0], codtau[0]);
    for(int iv = 0; iv < f->model.m; iv++) {
      int imem = f->varindex(f->deg, f->raf, f->model.m, ipg, iv);
      rhsCopy[imem] /= (wpg * det);
    }
  }
  
  // right hand side assembly
  for(int ino = 0; ino < cs->nb_fe_dof; ino++){
    rhsOut[ino] = 0;
  }

  for(int ie = 0; ie < cs->nbel; ie++){  

    int iemacro  = ie / (f->raf[0] * f->raf[1] * f->raf[2]);
    int isubcell = ie % (f->raf[0] * f->raf[1] * f->raf[2]);
    
 
    for(int iloc = 0; iloc < cs->nnodes; iloc++){
      real wpg;
      real xref[3];
      int ilocmacro = iloc + isubcell * cs->nnodes;
      ref_pg_vol(f->deg,f->raf,ilocmacro,xref,&wpg,NULL);
      real dtau[3][3],codtau[3][3];
      Ref2Phy(cs->simu->fd[iemacro].physnode,
	      xref,NULL,0,NULL,
	      dtau,codtau,NULL,NULL);
      real det = dot_product(dtau[0], codtau[0]);	
      int ino_dg = iloc + ie * cs->nnodes;
      int ino_fe = cs->dg_to_fe_index[ino_dg];
      coeff[ino_fe] += wpg * det; // We need to multiply by the member of dg node for 1 fe node.
      for (int iv=0; iv<cs->nb_phy_vars;iv++){ 
        int imem = f->varindex(f->deg,f->raf,f->model.m,
		          ilocmacro,cs->list_of_var[iv]) ;

        rhsOut[iv + cs->nb_phy_vars * ino_fe] += rhsCopy[imem] * wpg * det;//rhs[imem] * wpg * det; 
      }
    }
  }
  free(coeff);
  free(rhsCopy);
}

void VectorCgToDg(ContinuousSolver * cs, real * rhsIn, real * rhsOut){

  field* f0 = &cs->simu->fd[0];

  //printf("Cg to Dg Copy...\n");
  
  for(int ie = 0; ie < cs->simu->macromesh.nbelems; ie++){  
    for(int ipg = 0;ipg < cs->npgmacrocell; ipg++){
      int ino_dg = ipg + ie * cs->npgmacrocell;
      int ino_fe = cs->dg_to_fe_index[ino_dg];
      for(int var =0; var < cs->nb_phy_vars; var++){ 
        int ipot_fe = cs->nb_phy_vars*ino_fe + cs->list_of_var[var];
        int ipot = f0->varindex(f0->deg,f0->raf,f0->model.m,
                                ipg,cs->list_of_var[var]);
        rhsOut[ipot]=rhsIn[ipot_fe];
      }
    }
  }

}

void GenericOperator(PB_PC* pb_pc){

  int nnodes = pb_pc->D.nnodes;
  field* f0 = &pb_pc->D.simu->fd[0];
  bool assembled = pb_pc->D.lsol.mat_is_assembly && 
                   pb_pc->L1.lsol.mat_is_assembly &&
                   pb_pc->L2.lsol.mat_is_assembly &&
                   pb_pc->U1.lsol.mat_is_assembly &&
                   pb_pc->U2.lsol.mat_is_assembly &&
                   pb_pc->Schur.lsol.mat_is_assembly;
  int nbel = pb_pc->D.nbel;
  int dg_to_fe_index[pb_pc->D.nb_dg_nodes];
  for (int i=0; i<pb_pc->D.nb_dg_nodes; i++){
    dg_to_fe_index[i]=pb_pc->D.dg_to_fe_index[i];
  }

  if(!assembled){
    for(int ie = 0; ie < nbel; ie++){  

        
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
          Ref2Phy(pb_pc->D.simu->fd[iemacro].physnode,
                  xref,dphiref_i,0,NULL,
                  dtau,codtau,dphi_i,NULL);
        
          real det = dot_product(dtau[0], codtau[0]);
          if (ilocmacro==ipgmacro){
            basisPhi_i[0]=1;
          }
          else
          {
            basisPhi_i[0]=0;
          }
          basisPhi_i[1]=dphi_i[0]/det;
          basisPhi_i[2]=dphi_i[1]/det;
          basisPhi_i[3]=dphi_i[2]/det;
          for(int jloc = 0; jloc < nnodes; jloc++){
            int jlocmacro = jloc + isubcell * nnodes;
            int jno_dg = jloc + ie * nnodes;
            int jno_fe = dg_to_fe_index[jno_dg];
            grad_psi_pg(f0->deg,f0->raf,jlocmacro,ipgmacro,dphiref_j);
            Ref2Phy(pb_pc->D.simu->fd[iemacro].physnode,
                    xref,dphiref_j,0,NULL,
                    dtau,codtau,dphi_j,NULL);
            if (jlocmacro==ipgmacro){
              basisPhi_j[0]=1;
            }
            else
            {
              basisPhi_j[0]=0;
            }
            basisPhi_j[1]=dphi_j[0]/det;
            basisPhi_j[2]=dphi_j[1]/det;
            basisPhi_j[3]=dphi_j[2]/det;
            real val;
            real res[4];
            ContinuousSolver * cs ;

            // Building D Matrix
            cs = &pb_pc->D;
            for (int i=0; i<4; i++){
              res[i]=0;
              for (int j=0; j<4; j++){
                res[i]+=basisPhi_j[j]*cs->diff_op[i][j];
              }
            }
            val = dot_product(basisPhi_i, res) * wpg * det  ;
            AddLinearSolver(&cs->lsol,ino_fe*cs->nb_phy_vars,jno_fe*cs->nb_phy_vars,val);

            // Building L1 Matrix
            cs = &pb_pc->L1;
            for (int i=0; i<4; i++){
              res[i]=0;
              for (int j=0; j<4; j++){
                res[i]+=basisPhi_j[j]*cs->diff_op[i][j];
              }
            }
            val = dot_product(basisPhi_i, res) * wpg * det  ;
            AddLinearSolver(&cs->lsol,ino_fe*cs->nb_phy_vars,jno_fe*cs->nb_phy_vars,val);

            // Building L2 Matrix
            cs = &pb_pc->L2;
            for (int i=0; i<4; i++){
              res[i]=0;
              for (int j=0; j<4; j++){
                res[i]+=basisPhi_j[j]*cs->diff_op[i][j];
              }
            }
            val = dot_product(basisPhi_i, res) * wpg * det  ;
            AddLinearSolver(&cs->lsol,ino_fe*cs->nb_phy_vars,jno_fe*cs->nb_phy_vars,val);

            // Building U1 Matrix
            cs = &pb_pc->U1;
            for (int i=0; i<4; i++){
              res[i]=0;
              for (int j=0; j<4; j++){
                res[i]+=basisPhi_j[j]*cs->diff_op[i][j];
              }
            }
            val = dot_product(basisPhi_i, res) * wpg * det  ;
            AddLinearSolver(&cs->lsol,ino_fe*cs->nb_phy_vars,jno_fe*cs->nb_phy_vars,val);

            // Building U2 Matrix
            cs = &pb_pc->U2;
            for (int i=0; i<4; i++){
              res[i]=0;
              for (int j=0; j<4; j++){
                res[i]+=basisPhi_j[j]*cs->diff_op[i][j];
              }
            }
            val = dot_product(basisPhi_i, res) * wpg * det  ;
            AddLinearSolver(&cs->lsol,ino_fe*cs->nb_phy_vars,jno_fe*cs->nb_phy_vars,val);


            // Building Schur Matrix
            cs = &pb_pc->Schur;
            for (int iv1=0; iv1<cs->nb_phy_vars; iv1++){
              for (int iv2=0; iv2<cs->nb_phy_vars; iv2++){
                for (int i=0; i<4; i++){
                  res[i]=0;
                  for (int j=0; j<4; j++){
                    res[i]+=basisPhi_j[j]*cs->diff_op2vec[2*iv1+iv2][i][j];
                  }
                }
                val = dot_product(basisPhi_i, res) * wpg * det  ;
                AddLinearSolver(&cs->lsol,ino_fe*cs->nb_phy_vars+iv1,jno_fe*cs->nb_phy_vars+iv2,val);
              } // end for iv1
            } // end for iv2

          } // end for jloc
        } // end for iloc
      } // end for ipg
    } // end for ie
  } // end if assembled
}


void freePB_PC(PB_PC* pb_pc){
  free(pb_pc->list_mat2assemble);
  freeContinuousSolver(&pb_pc->D);
  freeContinuousSolver(&pb_pc->L1);
  freeContinuousSolver(&pb_pc->L2);
  freeContinuousSolver(&pb_pc->U1);
  freeContinuousSolver(&pb_pc->U2);
  freeContinuousSolver(&pb_pc->Schur);
  free(pb_pc->rhs_prediction);
  free(pb_pc->rhs_propagation);
  free(pb_pc->rhs_correction);
  
}

void reset(PB_PC* pb_pc){
  
  for (int i=0;i<pb_pc->D.nb_fe_nodes; i++){
    // resetting sols
    pb_pc->D.lsol.sol[i]=0;
    pb_pc->L1.lsol.sol[i]=0;
    pb_pc->L2.lsol.sol[i]=0;
    pb_pc->U1.lsol.sol[i]=0;
    pb_pc->U2.lsol.sol[i]=0;
    pb_pc->Schur.lsol.sol[i*2]=0;
    pb_pc->Schur.lsol.sol[i*2+1]=0;
    // resetting rhss
    pb_pc->D.lsol.rhs[i]=0;
    pb_pc->L1.lsol.rhs[i]=0;
    pb_pc->L2.lsol.rhs[i]=0;
    pb_pc->U1.lsol.rhs[i]=0;
    pb_pc->U2.lsol.rhs[i]=0;
    pb_pc->Schur.lsol.rhs[i*2]=0;
    pb_pc->Schur.lsol.rhs[i*2+1]=0;
  }
}

void Wave_test(ContinuousSolver* cs,real theta, real dt){

  real h=cs->simu->vmax*dt*theta;
  real waveMat[9][4][4] ={{{1.0,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0}},
                          {{0,0,0,0},
                           {-h,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0}},
                          {{0,0,0,0},
                           {0,0,0,0},
                           {-h,0,0,0},
                           {0,0,0,0}},
                          {{0,0,0,0},
                           {-h,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0}},
                          {{1.0,0,0,0}, 
                           {0,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0}},
                          {{0,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0}},
                          {{0,0,0,0},
                           {0,0,0,0},
                           {-h,0,0,0},
                           {0,0,0,0}},
                          {{0,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0}},
                          {{1.0,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0}}};
  for (int i=0; i<9; i++){
    for (int j=0; j<4; j++){
      for (int k=0; k<4; k++){
        cs->diff_op3vec[i][j][k]=waveMat[i][j][k];
      }
    }
  }


}

void PiDgToCg(ContinuousSolver * cs,real * rhsIn, real * rhsOut){
  
 field* f = &cs->simu->fd[0];

 real ** Restriction;
 real **tabvar;
 real **tabvar_transform;
 real* coeff = calloc(cs->nb_fe_nodes, sizeof(real));

 Restriction = calloc(cs->nb_fe_nodes,sizeof(real*));
 for (int i=0;i<cs->nb_fe_nodes;i++){
   Restriction[i] = calloc(cs->nb_dg_nodes,sizeof(real));
 }
 
 tabvar = calloc(cs->nb_phy_vars,sizeof(real*));
 for (int i=0;i<cs->nb_phy_vars;i++){
   tabvar[i] = calloc(cs->nb_dg_nodes,sizeof(real));
 }

 tabvar_transform = calloc(cs->nb_phy_vars,sizeof(real*));
 for (int i=0;i<cs->nb_phy_vars;i++){
   tabvar_transform[i] = calloc(cs->nb_fe_nodes,sizeof(real));
 }


 for (int i=0;i<cs->nb_dg_nodes;i++){
   for (int iv=0;iv<cs->nb_phy_vars;iv++){
     tabvar[iv][i] = rhsIn[iv+i*cs->nb_phy_vars];
   }
 }

for (int ino_fe=0;ino_fe<cs->nb_fe_nodes;ino_fe++){
   for (int ino_dg=0;ino_dg<cs->nb_dg_nodes;ino_dg++){
     Restriction[ino_fe][ino_dg]=0.;
   }
 }


 
 for(int ie = 0; ie < cs->nbel; ie++){  

    int iemacro  = ie / (f->raf[0] * f->raf[1] * f->raf[2]);
    int isubcell = ie % (f->raf[0] * f->raf[1] * f->raf[2]);
    
    for(int iloc = 0; iloc < cs->nnodes; iloc++){
      real wpg;
      real xref[3];
      int ilocmacro = iloc + isubcell * cs->nnodes;
      ref_pg_vol(f->deg,f->raf,ilocmacro,xref,&wpg,NULL);
      real dtau[3][3],codtau[3][3];
      Ref2Phy(cs->simu->fd[iemacro].physnode,
	      xref,NULL,0,NULL,
	      dtau,codtau,NULL,NULL);
      real det = dot_product(dtau[0], codtau[0]);	
      int ino_dg = iloc + ie * cs->nnodes;
      int ino_fe = cs->dg_to_fe_index[ino_dg];
      Restriction[ino_fe][ino_dg]= wpg * det; // We need to multiply by the member of dg node for 1 fe node.
      coeff[ino_fe] += wpg * det; // We need to multiply by the member of dg node for 1 fe node.
    }
 }



 for (int ino_fe=0;ino_fe<cs->nb_fe_nodes;ino_fe++){
   for (int ino_dg=0;ino_dg<cs->nb_dg_nodes;ino_dg++){
     Restriction[ino_fe][ino_dg]/=coeff[ino_fe];
   }
 }


 for (int ino_fe=0;ino_fe<cs->nb_fe_nodes;ino_fe++){
   for(int iv=0;iv<cs->nb_phy_vars;iv++){
	 tabvar_transform[iv][ino_fe]=0.0;
   }
   for (int ino_dg=0;ino_dg<cs->nb_dg_nodes;ino_dg++){
     for(int iv=0;iv<cs->nb_phy_vars;iv++){
	 tabvar_transform[iv][ino_fe]+=Restriction[ino_fe][ino_dg]*tabvar[iv][ino_dg];
     }
   }
 }



 for (int i=0;i<cs->nb_fe_nodes;i++){
   for (int iv=0;iv<cs->nb_phy_vars;iv++){
     rhsOut[iv+i*cs->nb_phy_vars]=tabvar_transform[iv][i];
   }
 }

 free(coeff);

}

void PiInvertCgToDg(ContinuousSolver * cs,real * rhsIn, real * rhsOut){
  
 field* f = &cs->simu->fd[0];
 real ** Reconstruction;
 real **tabvar;
 real **tabvar_transform;
 real* coeff = calloc(cs->nb_dg_nodes, sizeof(real));

 Reconstruction = calloc(cs->nb_dg_nodes,sizeof(real*));
 for (int i=0;i<cs->nb_dg_nodes;i++){
   Reconstruction[i] = calloc(cs->nb_fe_nodes,sizeof(real));
 }
 
 tabvar = calloc(cs->nb_phy_vars,sizeof(real*));
 for (int i=0;i<cs->nb_phy_vars;i++){
   tabvar[i] = calloc(cs->nb_fe_nodes,sizeof(real));
 }

 tabvar_transform = calloc(cs->nb_phy_vars,sizeof(real*));
 for (int i=0;i<cs->nb_phy_vars;i++){
   tabvar_transform[i] = calloc(cs->nb_dg_nodes,sizeof(real));
 }



 for (int i=0;i<cs->nb_fe_nodes;i++){
   for (int iv=0;iv<cs->nb_phy_vars;iv++){
     tabvar[iv][i] = rhsIn[iv+i*cs->nb_phy_vars];
   }
 }

 for (int ino_dg=0;ino_dg<cs->nb_dg_nodes;ino_dg++){
   for (int ino_fe=0;ino_fe<cs->nb_fe_nodes;ino_fe++){
     Reconstruction[ino_dg][ino_fe]=0.;
   }
 }

for(int ie = 0; ie < cs->nbel; ie++){  

    int iemacro  = ie / (f->raf[0] * f->raf[1] * f->raf[2]);
    int isubcell = ie % (f->raf[0] * f->raf[1] * f->raf[2]);
    
    for(int iloc = 0; iloc < cs->nnodes; iloc++){
      real wpg;
      real xref[3];
      int ilocmacro = iloc + isubcell * cs->nnodes;
      ref_pg_vol(f->deg,f->raf,ilocmacro,xref,&wpg,NULL);
      real dtau[3][3],codtau[3][3];
      Ref2Phy(cs->simu->fd[iemacro].physnode,
	      xref,NULL,0,NULL,
	      dtau,codtau,NULL,NULL);
      real det = dot_product(dtau[0], codtau[0]);	
      int ino_dg = iloc + ie * cs->nnodes;
      int ino_fe = cs->dg_to_fe_index[ino_dg];
      Reconstruction[ino_dg][ino_fe]= wpg * det; // We need to multiply by the member of dg node for 1 fe node.
      coeff[ino_dg] = wpg * det; // We need to multiply by the member of dg node for 1 fe node.
    }
 }


 for (int ino_dg=0;ino_dg<cs->nb_dg_nodes;ino_dg++){
   for (int ino_fe=0;ino_fe<cs->nb_fe_nodes;ino_fe++){
     Reconstruction[ino_dg][ino_fe]/=coeff[ino_dg];
   }
 }
//////////////////////////////////////////////////////////////////////////////////
 real ** Restriction;
 real* coeff2 = calloc(cs->nb_fe_nodes, sizeof(real));

 Restriction = calloc(cs->nb_fe_nodes,sizeof(real*));
 for (int i=0;i<cs->nb_fe_nodes;i++){
   Restriction[i] = calloc(cs->nb_dg_nodes,sizeof(real));
 }

for (int ino_fe=0;ino_fe<cs->nb_fe_nodes;ino_fe++){
   for (int ino_dg=0;ino_dg<cs->nb_dg_nodes;ino_dg++){
     Restriction[ino_fe][ino_dg]=0.;
   }
 }


 
 for(int ie = 0; ie < cs->nbel; ie++){  

    int iemacro  = ie / (f->raf[0] * f->raf[1] * f->raf[2]);
    int isubcell = ie % (f->raf[0] * f->raf[1] * f->raf[2]);
    
    for(int iloc = 0; iloc < cs->nnodes; iloc++){
      real wpg;
      real xref[3];
      int ilocmacro = iloc + isubcell * cs->nnodes;
      ref_pg_vol(f->deg,f->raf,ilocmacro,xref,&wpg,NULL);
      real dtau[3][3],codtau[3][3];
      Ref2Phy(cs->simu->fd[iemacro].physnode,
	      xref,NULL,0,NULL,
	      dtau,codtau,NULL,NULL);
      real det = dot_product(dtau[0], codtau[0]);	
      int ino_dg = iloc + ie * cs->nnodes;
      int ino_fe = cs->dg_to_fe_index[ino_dg];
      Restriction[ino_fe][ino_dg]= wpg * det; // We need to multiply by the member of dg node for 1 fe node.
      coeff2[ino_fe] += wpg * det; // We need to multiply by the member of dg node for 1 fe node.
    }
 }



 for (int ino_fe=0;ino_fe<cs->nb_fe_nodes;ino_fe++){
   for (int ino_dg=0;ino_dg<cs->nb_dg_nodes;ino_dg++){
     Restriction[ino_fe][ino_dg]/=coeff2[ino_fe];
   }
 }

 real **pipi = calloc(cs->nb_dg_nodes,sizeof(real*));
 real **pipi2 = calloc(cs->nb_fe_nodes,sizeof(real*));
 for (int i=0;i<cs->nb_dg_nodes;i++){
   pipi[i] = calloc(cs->nb_dg_nodes,sizeof(real));
 }
 for (int i=0;i<cs->nb_fe_nodes;i++){
   pipi2[i] = calloc(cs->nb_fe_nodes,sizeof(real));
 }
 //for (int i=0; i<cs->nb_dg_nodes; i++){
 //  for (int j=0; j<cs->nb_dg_nodes; j++){
 //    pipi[i][j]=0.;
 //    for (int k=0; k<cs->nb_fe_nodes; k++){
 //      pipi[i][j] += Reconstruction[i][k]*Restriction[k][j];
 //    }
 //    if (j==cs->nb_dg_nodes-1){
 //      printf("%.12e;",pipi[i][j]);
 //    }
 //    else{
 //      printf("%.12e,",pipi[i][j]);
 //    }
 //  }
 //}
 //for (int i=0; i<cs->nb_fe_nodes; i++){
 //  for (int j=0; j<cs->nb_fe_nodes; j++){
 //    pipi2[i][j]=0.;
 //    for (int k=0; k<cs->nb_dg_nodes; k++){
 //      pipi2[i][j] += Restriction[i][k]*Reconstruction[k][j];
 //    }
 //    //if (j==cs->nb_fe_nodes-1){
 //    //  printf("%.12e;",pipi2[i][j]);
 //    //}
 //    //else{
 //    //  printf("%.12e,",pipi2[i][j]);
 //    //}
 //  }
 //}
 
//////////////////////////////////////////////////////////////////////////////////

 for (int ino_dg=0;ino_dg<cs->nb_dg_nodes;ino_dg++){
    for(int iv=0;iv<cs->nb_phy_vars;iv++){
	 tabvar_transform[iv][ino_dg]=0.0;
     }
     for (int ino_fe=0;ino_fe<cs->nb_fe_nodes;ino_fe++){
       for(int iv=0;iv<cs->nb_phy_vars;iv++){
	 tabvar_transform[iv][ino_dg]+=Reconstruction[ino_dg][ino_fe]*tabvar[iv][ino_fe];
     }
   }
 }


for (int i=0;i<cs->nb_dg_nodes;i++){
   for (int iv=0;iv<cs->nb_phy_vars;iv++){
     rhsOut[iv+i*cs->nb_phy_vars]=tabvar_transform[iv][i];
   }
 }

 free(coeff);

}
