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
//  pb_pc->D.bc_assembly=ExactDirichletContinuousMatrix_PC;
//  pb_pc->U1.bc_assembly=ExactDirichletContinuousMatrix_PC;
//  pb_pc->U2.bc_assembly=ExactDirichletContinuousMatrix_PC;
//  pb_pc->Schur.bc_assembly=ExactDirichletContinuousMatrix_PC;
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

void Init_PBPC_Wave_SchurVelocity_BCPressure(Simulation *simu, PB_PC* pb_pc, int* list_mat2assemble){
  real tab[4][4];
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
  real L1Mat[4][4] = {{0,0,0,0},
                      {-h,0,0,0},
                      {0,0,0,0},
                      {0,0,0,0}};
  real L2Mat[4][4] = {{0,0,0,0},
                      {0,0,0,0},
                      {-h,0,0,0},
                      {0,0,0,0}};
  real U1Mat[4][4] = {{0,0,0,0},
                      {-h,0,0,0},
                      {0,0,0,0},
                      {0,0,0,0}};
  real U2Mat[4][4] = {{0,0,0,0},
                      {0,0,0,0},
                      {-h,0,0,0},
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

  pb_pc->D.diff_op = calloc(pb_pc->D.nb_phy_vars*pb_pc->D.nb_phy_vars,sizeof(SDO));
  pb_pc->Schur.diff_op = calloc(pb_pc->Schur.nb_phy_vars*pb_pc->Schur.nb_phy_vars,sizeof(SDO));
  pb_pc->L1.diff_op = calloc(pb_pc->L1.nb_phy_vars*pb_pc->L1.nb_phy_vars,sizeof(SDO));
  pb_pc->L2.diff_op = calloc(pb_pc->L2.nb_phy_vars*pb_pc->L2.nb_phy_vars,sizeof(SDO));
  pb_pc->U1.diff_op = calloc(pb_pc->U1.nb_phy_vars*pb_pc->U1.nb_phy_vars,sizeof(SDO));
  pb_pc->U2.diff_op = calloc(pb_pc->U2.nb_phy_vars*pb_pc->U2.nb_phy_vars,sizeof(SDO));

  // Applying these matrices inside the different ContinuousSolver
  for (int i=0; i<4; i++){
    for (int j=0; j<4; j++){
      for (int k=0; k<4; k++){
        if (SchurMat != NULL) pb_pc->Schur.diff_op[k].DO[i][j] = SchurMat[k][i][j];
      }
      if (DMat != NULL) pb_pc->D.diff_op[0].DO[i][j]  = DMat[i][j];
      if (L1Mat != NULL) pb_pc->L1.diff_op[0].DO[i][j] = L1Mat[i][j];
      if (L2Mat != NULL) pb_pc->L2.diff_op[0].DO[i][j] = L2Mat[i][j];
      if (U1Mat != NULL) pb_pc->U1.diff_op[0].DO[i][j] = U1Mat[i][j];
      if (U2Mat != NULL) pb_pc->U2.diff_op[0].DO[i][j] = U2Mat[i][j];
    }
  }
  // Assembling all operators' matrices
  GenericOperator_PBPC_Velocity(pb_pc);

  // For the Schur the the condition is Neumann
  pb_pc->Schur.bc_flux=NULL;
  pb_pc->Schur.bc_assembly=NULL;
  
  pb_pc->D.bc_flux=RobinFlux_SchurPressure;
  pb_pc->D.bc_assembly=RobinBoundaryConditionAssembly;
  pb_pc->D.bc_assembly(&pb_pc->D);

  pb_pc->L1.bc_flux=BoundaryTerm_Xderivative;
  pb_pc->L1.bc_assembly=BoundaryConditionFriedrichsAssembly;
  pb_pc->U1.bc_flux=BoundaryTerm_Xderivative;
  pb_pc->U1.bc_assembly=BoundaryConditionFriedrichsAssembly;

  pb_pc->L2.bc_flux=BoundaryTerm_Yderivative;
  pb_pc->L2.bc_assembly=BoundaryConditionFriedrichsAssembly;
  pb_pc->U2.bc_flux=BoundaryTerm_Yderivative;
  pb_pc->U2.bc_assembly=BoundaryConditionFriedrichsAssembly;

 
  pb_pc->U2.bc_assembly(&pb_pc->U2);
  pb_pc->U1.bc_assembly(&pb_pc->U1);
  pb_pc->L2.bc_assembly(&pb_pc->L2);
  pb_pc->L1.bc_assembly(&pb_pc->L1);

  pb_pc->Schur.lsol.MatVecProduct=MatVect;
  pb_pc->L1.lsol.MatVecProduct=MatVect;
  pb_pc->L2.lsol.MatVecProduct=MatVect;
  pb_pc->D.lsol.MatVecProduct=MatVect;
  pb_pc->U1.lsol.MatVecProduct=MatVect;
  pb_pc->U2.lsol.MatVecProduct=MatVect;

  //// Allocating all sub-equations right-hand sides
  pb_pc->rhs_prediction = calloc(pb_pc->D.lsol.neq,sizeof(real));
  pb_pc->rhs_propagation = calloc(pb_pc->Schur.lsol.neq,sizeof(real));
  pb_pc->rhs_correction = calloc(pb_pc->D.lsol.neq,sizeof(real));

}

void Init_PBPC_Wave_SchurVelocity_BCVelocity(Simulation *simu, PB_PC* pb_pc, int* list_mat2assemble){
  real tab[4][4];
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
  real L1Mat[4][4] = {{0,0,0,0},
                      {-h,0,0,0},
                      {0,0,0,0},
                      {0,0,0,0}};
  real L2Mat[4][4] = {{0,0,0,0},
                      {0,0,0,0},
                      {-h,0,0,0},
                      {0,0,0,0}};
  real U1Mat[4][4] = {{0,0,0,0},
                      {-h,0,0,0},
                      {0,0,0,0},
                      {0,0,0,0}};
  real U2Mat[4][4] = {{0,0,0,0},
                      {0,0,0,0},
                      {-h,0,0,0},
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

  pb_pc->D.diff_op = calloc(pb_pc->D.nb_phy_vars*pb_pc->D.nb_phy_vars,sizeof(SDO));
  pb_pc->Schur.diff_op = calloc(pb_pc->Schur.nb_phy_vars*pb_pc->Schur.nb_phy_vars,sizeof(SDO));
  pb_pc->L1.diff_op = calloc(pb_pc->L1.nb_phy_vars*pb_pc->L1.nb_phy_vars,sizeof(SDO));
  pb_pc->L2.diff_op = calloc(pb_pc->L2.nb_phy_vars*pb_pc->L2.nb_phy_vars,sizeof(SDO));
  pb_pc->U1.diff_op = calloc(pb_pc->U1.nb_phy_vars*pb_pc->U1.nb_phy_vars,sizeof(SDO));
  pb_pc->U2.diff_op = calloc(pb_pc->U2.nb_phy_vars*pb_pc->U2.nb_phy_vars,sizeof(SDO));

  // Applying these matrices inside the different ContinuousSolver
  for (int i=0; i<4; i++){
    for (int j=0; j<4; j++){
      for (int k=0; k<4; k++){
        if (SchurMat != NULL) pb_pc->Schur.diff_op[k].DO[i][j] = SchurMat[k][i][j];
      }
      if (DMat != NULL) pb_pc->D.diff_op[0].DO[i][j]  = DMat[i][j];
      if (L1Mat != NULL) pb_pc->L1.diff_op[0].DO[i][j] = L1Mat[i][j];
      if (L2Mat != NULL) pb_pc->L2.diff_op[0].DO[i][j] = L2Mat[i][j];
      if (U1Mat != NULL) pb_pc->U1.diff_op[0].DO[i][j] = U1Mat[i][j];
      if (U2Mat != NULL) pb_pc->U2.diff_op[0].DO[i][j] = U2Mat[i][j];
    }
  }
  // Assembling all operators' matrices
  GenericOperator_PBPC_Velocity(pb_pc);

  // For the Schur the the condition is Neumann
  
  pb_pc->D.bc_flux=NULL;
  pb_pc->D.bc_assembly=NULL;
  pb_pc->Schur.bc_flux=NULL;
  pb_pc->Schur.bc_assembly=NULL;

  pb_pc->Schur.bc_flux=Dirichlet_Velocity;
  pb_pc->Schur.bc_assembly=BoundaryConditionFriedrichsAssembly;
  pb_pc->Schur.bc_assembly(&pb_pc->Schur);

  pb_pc->L1.bc_flux=BoundaryTerm_Xderivative;
  pb_pc->L1.bc_assembly=BoundaryConditionFriedrichsAssembly;
  pb_pc->U1.bc_flux=BoundaryTerm_Xderivative;
  pb_pc->U1.bc_assembly=BoundaryConditionFriedrichsAssembly;

  pb_pc->L2.bc_flux=BoundaryTerm_Yderivative;
  pb_pc->L2.bc_assembly=BoundaryConditionFriedrichsAssembly;
  pb_pc->U2.bc_flux=BoundaryTerm_Yderivative;
  pb_pc->U2.bc_assembly=BoundaryConditionFriedrichsAssembly;

 
  pb_pc->U2.bc_assembly(&pb_pc->U2);
  pb_pc->U1.bc_assembly(&pb_pc->U1);
  pb_pc->L2.bc_assembly(&pb_pc->L2);
  pb_pc->L1.bc_assembly(&pb_pc->L1);

  pb_pc->Schur.lsol.MatVecProduct=MatVect;
  pb_pc->L1.lsol.MatVecProduct=MatVect;
  pb_pc->L2.lsol.MatVecProduct=MatVect;
  pb_pc->D.lsol.MatVecProduct=MatVect;
  pb_pc->U1.lsol.MatVecProduct=MatVect;
  pb_pc->U2.lsol.MatVecProduct=MatVect;

  //// Allocating all sub-equations right-hand sides
  pb_pc->rhs_prediction = calloc(pb_pc->D.lsol.neq,sizeof(real));
  pb_pc->rhs_propagation = calloc(pb_pc->Schur.lsol.neq,sizeof(real));
  pb_pc->rhs_correction = calloc(pb_pc->D.lsol.neq,sizeof(real));

}


void Init_PBPC_Wave_SchurPressure_BCVelocity(Simulation *simu, PB_PC* pb_pc, int* list_mat2assemble){
 
  // system is linear
  pb_pc->nonlinear=0;


  pb_pc->list_mat2assemble = calloc(6, sizeof(int));
  for (int i=0; i<6; i++) pb_pc->list_mat2assemble[i] = list_mat2assemble[i];
  
  // Initialization variables. Modify at will.
  int nb_varD = 2;
  int * listvarD = calloc(nb_varD,sizeof(int));
  listvarD[0]=1;
  listvarD[1]=2;
  
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
  
  int nb_varSchur = 1;
  int * listvarSchur = calloc(nb_varSchur,sizeof(int));
  listvarSchur[0]=0;


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
  real SchurMat[4][4]  = {{1.0,0,0,0},
			  {0,h*h,0,0},
			  {0,0,h*h,0},
			  {0,0,0,0}};
  real L1Mat[4][4] = {{0,0,0,0},
                      {-h,0,0,0},
                      {0,0,0,0},
                      {0,0,0,0}};
  real L2Mat[4][4] = {{0,0,0,0},
                      {0,0,0,0},
                      {-h,0,0,0},
                      {0,0,0,0}};
  real U1Mat[4][4] = {{0,0,0,0},
                      {-h,0,0,0},
                      {0,0,0,0},
                      {0,0,0,0}};
  real U2Mat[4][4] = {{0,0,0,0},
                      {0,0,0,0},
                      {-h,0,0,0},
                      {0,0,0,0}};
  real DMat[4][4][4] ={{{1.0,0,0,0},
			{0,0,0,0},
			{0,0,0,0},
			{0,0,0,0}},
		       {{0.0,0,0,0},
			{0,0,0,0},
			{0,0,0,0},
			{0,0,0,0}},
		       {{0.0,0,0,0},
			{0,0,0,0},
			{0,0,0,0},
			{0,0,0,0}},
		       {{1.0,0,0,0},
			{0,0,0,0},
			{0,0,0,0},
			{0,0,0,0}}};

  pb_pc->D.diff_op = calloc(pb_pc->D.nb_phy_vars*pb_pc->D.nb_phy_vars,sizeof(SDO));
  pb_pc->Schur.diff_op = calloc(pb_pc->Schur.nb_phy_vars*pb_pc->Schur.nb_phy_vars,sizeof(SDO));
  pb_pc->L1.diff_op = calloc(pb_pc->L1.nb_phy_vars*pb_pc->L1.nb_phy_vars,sizeof(SDO));
  pb_pc->L2.diff_op = calloc(pb_pc->L2.nb_phy_vars*pb_pc->L2.nb_phy_vars,sizeof(SDO));
  pb_pc->U1.diff_op = calloc(pb_pc->U1.nb_phy_vars*pb_pc->U1.nb_phy_vars,sizeof(SDO));
  pb_pc->U2.diff_op = calloc(pb_pc->U2.nb_phy_vars*pb_pc->U2.nb_phy_vars,sizeof(SDO));
  
  // Applying these matrices inside the different ContinuousSolver
  for (int i=0; i<4; i++){
    for (int j=0; j<4; j++){
      for (int k=0; k<4; k++){
        if (DMat != NULL) pb_pc->D.diff_op[k].DO[i][j] = DMat[k][i][j];
      }
      if (SchurMat != NULL) pb_pc->Schur.diff_op[0].DO[i][j]  = SchurMat[i][j];
      if (L1Mat != NULL) pb_pc->L1.diff_op[0].DO[i][j] = L1Mat[i][j];
      if (L2Mat != NULL) pb_pc->L2.diff_op[0].DO[i][j] = L2Mat[i][j];
      if (U1Mat != NULL) pb_pc->U1.diff_op[0].DO[i][j] = U1Mat[i][j];
      if (U2Mat != NULL) pb_pc->U2.diff_op[0].DO[i][j] = U2Mat[i][j];
    }
  }

  // Assembling all operators' matrices
  GenericOperator_PBPC_Pressure(pb_pc);

  // For the Schur the the condition is Neumann
  
  pb_pc->D.bc_flux=NULL;
  pb_pc->Schur.bc_flux=NULL;
  pb_pc->D.bc_assembly=NULL;
  pb_pc->Schur.bc_assembly=NULL;

  pb_pc->D.bc_flux=Dirichlet_Velocity;
  pb_pc->D.bc_assembly=BoundaryConditionFriedrichsAssembly;
  pb_pc->D.bc_assembly(&pb_pc->D);

  pb_pc->L1.bc_flux=BoundaryTerm_Xderivative;
  pb_pc->L1.bc_assembly=BoundaryConditionFriedrichsAssembly;
  pb_pc->U1.bc_flux=BoundaryTerm_Xderivative;
  pb_pc->U1.bc_assembly=BoundaryConditionFriedrichsAssembly;

  pb_pc->L2.bc_flux=BoundaryTerm_Yderivative;
  pb_pc->L2.bc_assembly=BoundaryConditionFriedrichsAssembly;
  pb_pc->U2.bc_flux=BoundaryTerm_Yderivative;
  pb_pc->U2.bc_assembly=BoundaryConditionFriedrichsAssembly;

  pb_pc->U2.bc_assembly(&pb_pc->U2);
  pb_pc->U1.bc_assembly(&pb_pc->U1);
  pb_pc->L2.bc_assembly(&pb_pc->L2);
  pb_pc->L1.bc_assembly(&pb_pc->L1);

  pb_pc->Schur.lsol.MatVecProduct=MatVect;
  pb_pc->L1.lsol.MatVecProduct=MatVect;
  pb_pc->L2.lsol.MatVecProduct=MatVect;
  pb_pc->D.lsol.MatVecProduct=MatVect;
  pb_pc->U1.lsol.MatVecProduct=MatVect;
  pb_pc->U2.lsol.MatVecProduct=MatVect;

  //// Allocating all sub-equations right-hand sides
  pb_pc->rhs_prediction = calloc(pb_pc->D.lsol.neq,sizeof(real));
  pb_pc->rhs_propagation = calloc(pb_pc->Schur.lsol.neq,sizeof(real));
  pb_pc->rhs_correction = calloc(pb_pc->D.lsol.neq,sizeof(real));

}

void Init_PBPC_Wave_SchurPressure_BCPressure(Simulation *simu, PB_PC* pb_pc, int* list_mat2assemble){
 
  // system is linear
  pb_pc->nonlinear=0;


  pb_pc->list_mat2assemble = calloc(6, sizeof(int));
  for (int i=0; i<6; i++) pb_pc->list_mat2assemble[i] = list_mat2assemble[i];
  
  // Initialization variables. Modify at will.
  int nb_varD = 2;
  int * listvarD = calloc(nb_varD,sizeof(int));
  listvarD[0]=1;
  listvarD[1]=2;
  
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
  
  int nb_varSchur = 1;
  int * listvarSchur = calloc(nb_varSchur,sizeof(int));
  listvarSchur[0]=0;


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
  real SchurMat[4][4]  = {{1.0,0,0,0},
			  {0,h*h,0,0},
			  {0,0,h*h,0},
			  {0,0,0,0}};
  real L1Mat[4][4] = {{0,0,0,0},
                      {-h,0,0,0},
                      {0,0,0,0},
                      {0,0,0,0}};
  real L2Mat[4][4] = {{0,0,0,0},
                      {0,0,0,0},
                      {-h,0,0,0},
                      {0,0,0,0}};
  real U1Mat[4][4] = {{0,0,0,0},
                      {-h,0,0,0},
                      {0,0,0,0},
                      {0,0,0,0}};
  real U2Mat[4][4] = {{0,0,0,0},
                      {0,0,0,0},
                      {-h,0,0,0},
                      {0,0,0,0}};
  real DMat[4][4][4] ={{{1.0,0,0,0},
			{0,0,0,0},
			{0,0,0,0},
			{0,0,0,0}},
		       {{0.0,0,0,0},
			{0,0,0,0},
			{0,0,0,0},
			{0,0,0,0}},
		       {{0.0,0,0,0},
			{0,0,0,0},
			{0,0,0,0},
			{0,0,0,0}},
		       {{1.0,0,0,0},
			{0,0,0,0},
			{0,0,0,0},
			{0,0,0,0}}};

  pb_pc->D.diff_op = calloc(pb_pc->D.nb_phy_vars*pb_pc->D.nb_phy_vars,sizeof(SDO));
  pb_pc->Schur.diff_op = calloc(pb_pc->Schur.nb_phy_vars*pb_pc->Schur.nb_phy_vars,sizeof(SDO));
  pb_pc->L1.diff_op = calloc(pb_pc->L1.nb_phy_vars*pb_pc->L1.nb_phy_vars,sizeof(SDO));
  pb_pc->L2.diff_op = calloc(pb_pc->L2.nb_phy_vars*pb_pc->L2.nb_phy_vars,sizeof(SDO));
  pb_pc->U1.diff_op = calloc(pb_pc->U1.nb_phy_vars*pb_pc->U1.nb_phy_vars,sizeof(SDO));
  pb_pc->U2.diff_op = calloc(pb_pc->U2.nb_phy_vars*pb_pc->U2.nb_phy_vars,sizeof(SDO));
  
  // Applying these matrices inside the different ContinuousSolver
  for (int i=0; i<4; i++){
    for (int j=0; j<4; j++){
      for (int k=0; k<4; k++){
        if (DMat != NULL) pb_pc->D.diff_op[k].DO[i][j] = DMat[k][i][j];
      }
      if (SchurMat != NULL) pb_pc->Schur.diff_op[0].DO[i][j]  = SchurMat[i][j];
      if (L1Mat != NULL) pb_pc->L1.diff_op[0].DO[i][j] = L1Mat[i][j];
      if (L2Mat != NULL) pb_pc->L2.diff_op[0].DO[i][j] = L2Mat[i][j];
      if (U1Mat != NULL) pb_pc->U1.diff_op[0].DO[i][j] = U1Mat[i][j];
      if (U2Mat != NULL) pb_pc->U2.diff_op[0].DO[i][j] = U2Mat[i][j];
    }
  }
  
 // Assembling all operators' matrices
  GenericOperator_PBPC_Pressure(pb_pc);
  pb_pc->D.bc_flux=NULL;
  pb_pc->D.bc_assembly=NULL;
  
  pb_pc->Schur.bc_flux=RobinFlux_SchurPressure;
  pb_pc->Schur.bc_assembly=RobinBoundaryConditionAssembly;

  pb_pc->L1.bc_flux=BoundaryTerm_Xderivative;
  pb_pc->L1.bc_assembly=BoundaryConditionFriedrichsAssembly;
  pb_pc->U1.bc_flux=BoundaryTerm_Xderivative;
  pb_pc->U1.bc_assembly=BoundaryConditionFriedrichsAssembly;

  pb_pc->L2.bc_flux=BoundaryTerm_Yderivative;
  pb_pc->L2.bc_assembly=BoundaryConditionFriedrichsAssembly;
  pb_pc->U2.bc_flux=BoundaryTerm_Yderivative;
  pb_pc->U2.bc_assembly=BoundaryConditionFriedrichsAssembly;

  pb_pc->Schur.bc_assembly(&pb_pc->Schur);
  pb_pc->U2.bc_assembly(&pb_pc->U2);
  pb_pc->U1.bc_assembly(&pb_pc->U1);
  pb_pc->L2.bc_assembly(&pb_pc->L2);
  pb_pc->L1.bc_assembly(&pb_pc->L1);

  pb_pc->Schur.lsol.MatVecProduct=MatVect;
  pb_pc->L1.lsol.MatVecProduct=MatVect;
  pb_pc->L2.lsol.MatVecProduct=MatVect;
  pb_pc->D.lsol.MatVecProduct=MatVect;
  pb_pc->U1.lsol.MatVecProduct=MatVect;
  pb_pc->U2.lsol.MatVecProduct=MatVect;

  //// Allocating all sub-equations right-hand sides
  pb_pc->rhs_prediction = calloc(pb_pc->D.lsol.neq,sizeof(real));
  pb_pc->rhs_propagation = calloc(pb_pc->Schur.lsol.neq,sizeof(real));
  pb_pc->rhs_correction = calloc(pb_pc->D.lsol.neq,sizeof(real));
}




void Init_Parameters_PhyBasedPC(PB_PC* pb_pc){

  pb_pc->solver_prediction=LU;
  pb_pc->solver_propagation=LU;
  pb_pc->solver_correction=LU;

  pb_pc->pc_prediction=NONE;
  pb_pc->pc_propagation=NONE;
  pb_pc->pc_correction=NONE;

  pb_pc->tol_prediction=1.e-8;
  pb_pc->tol_propagation=1.e-8;
  pb_pc->tol_correction=1.e-8;

  pb_pc->itermax_prediction=1000;
  pb_pc->itermax_propagation=1000;
  pb_pc->itermax_correction=1000;
  
  pb_pc->restart_prediction=30;
  pb_pc->restart_propagation=30;
  pb_pc->restart_correction=30;
}


void PhyBased_PC_CG(PB_PC* pb_pc, Simulation *simu, real* globalSol, real*globalRHS){
  
  // 0)1) Reset everything (needed for time evolution)
  reset(pb_pc);

  /*for (int i=0;i<pb_pc->D.nb_fe_nodes;i++){
    printf(" 1) P[%d] = %8e, U[%d] = %8e, V[%d] = %8e\n", i, globalRHS[i*3], i, globalRHS[i*3+1], i, globalRHS[i*3+2]);
    }*/
  // Assembling all operators' matrices
  if(pb_pc->nonlinear == 1){
    GenericOperator_PBPC_Velocity(pb_pc);
  }

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

  pb_pc->D.lsol.solver_type=pb_pc->solver_prediction;
  pb_pc->D.lsol.tol=pb_pc->tol_prediction;
  pb_pc->D.lsol.pc_type=pb_pc->pc_prediction;
  pb_pc->D.lsol.iter_max=pb_pc->itermax_prediction;
  pb_pc->D.lsol.restart_gmres=pb_pc->restart_prediction;

  //printf("RHS assembly.....\n");
  for (int i=0;i<pb_pc->D.nb_fe_nodes;i++){
    pb_pc->D.lsol.rhs[i] = globalRHS[i*3];
  }
  //printf("Solution...\n");
  Advanced_SolveLinearSolver(&pb_pc->D.lsol,simu);
  /*for (int i=0;i<pb_pc->D.nb_fe_nodes;i++){
    printf("1) P[%d] = %8e\n", i, pb_pc->D.lsol.sol[i]);
    }*/
 
  // 2) PROPAGATION STEP

  pb_pc->Schur.lsol.solver_type=pb_pc->solver_propagation;
  pb_pc->Schur.lsol.tol=pb_pc->tol_propagation;
  pb_pc->Schur.lsol.pc_type=pb_pc->pc_propagation;
  pb_pc->Schur.lsol.iter_max=pb_pc->itermax_propagation;
  pb_pc->Schur.lsol.restart_gmres=pb_pc->restart_propagation;
  
  // Parsing L1P, L2P into the "sol" of L1 and L2 (since it is unused).
  pb_pc->L1.lsol.MatVecProduct(&pb_pc->L1.lsol,pb_pc->D.lsol.sol,pb_pc->L1.lsol.sol);
  pb_pc->L2.lsol.MatVecProduct(&pb_pc->L2.lsol,pb_pc->D.lsol.sol,pb_pc->L2.lsol.sol);

  //printf("RHS assembly.....\n");
  for (int i=0;i<pb_pc->D.nb_fe_nodes;i++){
    pb_pc->Schur.lsol.rhs[i*2]   = globalRHS[i*3+1] - pb_pc->L1.lsol.sol[i];
    pb_pc->Schur.lsol.rhs[i*2+1] = globalRHS[i*3+2] - pb_pc->L2.lsol.sol[i];
  }
  
  Advanced_SolveLinearSolver(&pb_pc->Schur.lsol,simu);
  /* for (int i=0;i<pb_pc->D.nb_fe_nodes;i++){
    printf("2) U[%d] = %8e, V[%d] = %8e\n", i, pb_pc->Schur.lsol.sol[i*2], i, pb_pc->Schur.lsol.sol[i*2+1]);
    }*/

  // 3) CORRECTION STEP
  // Extracting both U1 and U2 from the previous solution of the propagation step
  real *solU1=calloc(pb_pc->U1.nb_fe_nodes, sizeof(real));
  real *solU2=calloc(pb_pc->U2.nb_fe_nodes, sizeof(real));
  extract2CGVectors(&pb_pc->U1,&pb_pc->U2,pb_pc->Schur.lsol.sol,solU1,solU2);

  pb_pc->D.lsol.solver_type=pb_pc->solver_correction;
  pb_pc->D.lsol.tol=pb_pc->tol_correction;
  pb_pc->D.lsol.pc_type=pb_pc->pc_correction;
  pb_pc->D.lsol.iter_max=pb_pc->itermax_correction;
  pb_pc->D.lsol.restart_gmres=pb_pc->restart_correction;
  
  pb_pc->U1.lsol.MatVecProduct(&pb_pc->U1.lsol,solU1,pb_pc->U1.lsol.sol);
  pb_pc->U2.lsol.MatVecProduct(&pb_pc->U2.lsol,solU2,pb_pc->U2.lsol.sol);
  pb_pc->D.lsol.MatVecProduct(&pb_pc->D.lsol,pb_pc->D.lsol.sol,pb_pc->D.lsol.rhs);

  //printf("RHS assembly.....\n");
  for (int i=0;i<pb_pc->D.nb_fe_nodes;i++){
    pb_pc->D.lsol.rhs[i] += - pb_pc->U1.lsol.sol[i] - pb_pc->U2.lsol.sol[i];
  }
  //printf("Solution...\n");
  Advanced_SolveLinearSolver(&pb_pc->D.lsol,simu);

  /*for (int i=0;i<pb_pc->D.nb_fe_nodes;i++){
    printf("Prout Schur U : P[%d] = %8e, U[%d] = %8e, V[%d] = %8e\n", i, pb_pc->D.lsol.sol[i], i, pb_pc->Schur.lsol.sol[2*i], i, pb_pc->Schur.lsol.sol[2*i+1]);
  }*/
  
 
  // 4) OUTPUT STEP Final concatenation
  cat2CGVectors(&pb_pc->D,&pb_pc->Schur,pb_pc->D.lsol.sol,pb_pc->Schur.lsol.sol,globalSol);

  freeContinuousSolver(&waveSolver);
  free(solU1);
  free(solU2);
}


void PhyBased_PC_InvertSchur_CG(PB_PC* pb_pc, Simulation *simu, real* globalSol, real*globalRHS){
  reset(pb_pc);

  // Assembling all operators' matrices
  if(pb_pc->nonlinear == 1){
    GenericOperator_PBPC_Pressure(pb_pc);
  }

  /*for (int i=0;i<pb_pc->D.nb_fe_nodes;i++){
    printf(" 1) P[%d] = %8e, U[%d] = %8e, V[%d] = %8e\n", i, globalRHS[i*3], i, globalRHS[i*3+1], i, globalRHS[i*3+2]);
    }*/
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
  pb_pc->D.lsol.solver_type=pb_pc->solver_prediction;
  pb_pc->D.lsol.tol=pb_pc->tol_prediction;
  pb_pc->D.lsol.pc_type=pb_pc->pc_prediction;
  pb_pc->D.lsol.iter_max=pb_pc->itermax_prediction;
  pb_pc->D.lsol.restart_gmres=pb_pc->restart_prediction;

  //printf("RHS assembly.....\n");
  for (int i=0;i<pb_pc->D.nb_fe_nodes;i++){
    pb_pc->D.lsol.rhs[i*2]   = globalRHS[i*3+1];
    pb_pc->D.lsol.rhs[i*2+1] = globalRHS[i*3+2];
  }
  //printf("Solution prediction\n");
  Advanced_SolveLinearSolver(&pb_pc->D.lsol,simu);
  /*for (int i=0;i<pb_pc->D.nb_fe_nodes;i++){
    printf(" 1) U[%d] = %8e, V[%d] = %8e\n", i, pb_pc->D.lsol.sol[i*2], i, pb_pc->D.lsol.sol[i*2+1]);
    }*/

  // 2) PROPAGATION STEP
  pb_pc->Schur.lsol.solver_type=pb_pc->solver_propagation;
  pb_pc->Schur.lsol.tol=pb_pc->tol_propagation;
  pb_pc->Schur.lsol.pc_type=pb_pc->pc_propagation;
  pb_pc->Schur.lsol.iter_max=pb_pc->itermax_propagation;
  pb_pc->Schur.lsol.restart_gmres=pb_pc->restart_propagation;

  real *solU1=calloc(pb_pc->U1.nb_fe_nodes, sizeof(real));
  real *solU2=calloc(pb_pc->U2.nb_fe_nodes, sizeof(real));
  extract2CGVectors(&pb_pc->U1,&pb_pc->U2,pb_pc->D.lsol.sol,solU1,solU2);
  
  // Parsing L1P, L2P into the "sol" of L1 and L2 (since it is unused).
  pb_pc->L1.lsol.MatVecProduct(&pb_pc->L1.lsol,solU1,pb_pc->L1.lsol.sol);
  pb_pc->L2.lsol.MatVecProduct(&pb_pc->L2.lsol,solU2,pb_pc->L2.lsol.sol);

  //printf("RHS assembly.....\n");
  for (int i=0;i<pb_pc->D.nb_fe_nodes;i++){
    pb_pc->Schur.lsol.rhs[i]   = globalRHS[i*3] - pb_pc->L1.lsol.sol[i] - pb_pc->L2.lsol.sol[i];
  }
  
  //printf("Solution propagation\n");
  Advanced_SolveLinearSolver(&pb_pc->Schur.lsol,simu);
  /*for (int i=0;i<pb_pc->D.nb_fe_nodes;i++){
    printf(" 2) P[%d] = %8e\n", i, pb_pc->Schur.lsol.sol[i]);
    }*/

  
  // 3) CORRECTION STEP

  pb_pc->D.lsol.solver_type=pb_pc->solver_correction;
  pb_pc->D.lsol.tol=pb_pc->tol_correction;
  pb_pc->D.lsol.pc_type=pb_pc->pc_correction;
  pb_pc->D.lsol.iter_max=pb_pc->itermax_correction;
  pb_pc->D.lsol.restart_gmres=pb_pc->restart_correction;
  
  pb_pc->U1.lsol.MatVecProduct(&pb_pc->U1.lsol,pb_pc->Schur.lsol.sol,pb_pc->U1.lsol.sol);
  pb_pc->U2.lsol.MatVecProduct(&pb_pc->U2.lsol,pb_pc->Schur.lsol.sol,pb_pc->U2.lsol.sol);
  pb_pc->D.lsol.MatVecProduct(&pb_pc->D.lsol,pb_pc->D.lsol.sol,pb_pc->D.lsol.rhs);

  //printf("RHS assembly.....\n");
  for (int i=0;i<pb_pc->D.nb_fe_nodes;i++){
    pb_pc->D.lsol.rhs[i*2] = pb_pc->D.lsol.rhs[i*2]- pb_pc->U1.lsol.sol[i];
    pb_pc->D.lsol.rhs[i*2+1] =  pb_pc->D.lsol.rhs[i*2+1] - pb_pc->U2.lsol.sol[i];
  }

  //printf("Solution correction\n");
  Advanced_SolveLinearSolver(&pb_pc->D.lsol,simu);
  /*for (int i=0;i<pb_pc->D.nb_fe_nodes;i++){
    printf("Prout Schur P : P[%d] = %8e, U[%d] = %8e, V[%d] = %8e\n", i, pb_pc->Schur.lsol.sol[i], i, pb_pc->D.lsol.sol[2*i], i, pb_pc->D.lsol.sol[2*i+1]);
    }*/

  
  // 4) OUTPUT STEP Final concatenation
  cat2CGVectors(&pb_pc->D,&pb_pc->Schur,pb_pc->D.lsol.sol,pb_pc->Schur.lsol.sol,globalSol);

  freeContinuousSolver(&waveSolver);
  free(solU1);
  free(solU2);
}



void GenericOperator_PBPC_Velocity(PB_PC* pb_pc){

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
                res[i]+=basisPhi_j[j]*cs->diff_op[0].DO[i][j];
              }
            }
            val = dot_product(basisPhi_i, res) * wpg * det  ;
            AddLinearSolver(&cs->lsol,ino_fe*cs->nb_phy_vars,jno_fe*cs->nb_phy_vars,val);

            // Building L1 Matrix
            cs = &pb_pc->L1;
            for (int i=0; i<4; i++){
              res[i]=0;
              for (int j=0; j<4; j++){
                res[i]+=basisPhi_j[j]*cs->diff_op[0].DO[i][j];
              }
            }
            val = dot_product(basisPhi_i, res) * wpg * det  ;
            AddLinearSolver(&cs->lsol,ino_fe*cs->nb_phy_vars,jno_fe*cs->nb_phy_vars,val);

            // Building L2 Matrix
            cs = &pb_pc->L2;
            for (int i=0; i<4; i++){
              res[i]=0;
              for (int j=0; j<4; j++){
                res[i]+=basisPhi_j[j]*cs->diff_op[0].DO[i][j];
              }
            }
            val = dot_product(basisPhi_i, res) * wpg * det  ;
            AddLinearSolver(&cs->lsol,ino_fe*cs->nb_phy_vars,jno_fe*cs->nb_phy_vars,val);

            // Building U1 Matrix
            cs = &pb_pc->U1;
            for (int i=0; i<4; i++){
              res[i]=0;
              for (int j=0; j<4; j++){
                res[i]+=basisPhi_j[j]*cs->diff_op[0].DO[i][j];
              }
            }
            val = dot_product(basisPhi_i, res) * wpg * det  ;
            AddLinearSolver(&cs->lsol,ino_fe*cs->nb_phy_vars,jno_fe*cs->nb_phy_vars,val);

            // Building U2 Matrix
            cs = &pb_pc->U2;
            for (int i=0; i<4; i++){
              res[i]=0;
              for (int j=0; j<4; j++){
                res[i]+=basisPhi_j[j]*cs->diff_op[0].DO[i][j];
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
                    res[i]+=basisPhi_j[j]*cs->diff_op[2*iv1+iv2].DO[i][j];
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


void GenericOperator_PBPC_Pressure(PB_PC* pb_pc){

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
            cs = &pb_pc->Schur;
            for (int i=0; i<4; i++){
              res[i]=0;
              for (int j=0; j<4; j++){
                res[i]+=basisPhi_j[j]*cs->diff_op[0].DO[i][j];
              }
            }
            val = dot_product(basisPhi_i, res) * wpg * det  ;
            AddLinearSolver(&cs->lsol,ino_fe*cs->nb_phy_vars,jno_fe*cs->nb_phy_vars,val);

            // Building L1 Matrix
            cs = &pb_pc->L1;
            for (int i=0; i<4; i++){
              res[i]=0;
              for (int j=0; j<4; j++){
                res[i]+=basisPhi_j[j]*cs->diff_op[0].DO[i][j];
              }
            }
            val = dot_product(basisPhi_i, res) * wpg * det  ;
            AddLinearSolver(&cs->lsol,ino_fe*cs->nb_phy_vars,jno_fe*cs->nb_phy_vars,val);

            // Building L2 Matrix
            cs = &pb_pc->L2;
            for (int i=0; i<4; i++){
              res[i]=0;
              for (int j=0; j<4; j++){
                res[i]+=basisPhi_j[j]*cs->diff_op[0].DO[i][j];
              }
            }
            val = dot_product(basisPhi_i, res) * wpg * det  ;
            AddLinearSolver(&cs->lsol,ino_fe*cs->nb_phy_vars,jno_fe*cs->nb_phy_vars,val);

            // Building U1 Matrix
            cs = &pb_pc->U1;
            for (int i=0; i<4; i++){
              res[i]=0;
              for (int j=0; j<4; j++){
                res[i]+=basisPhi_j[j]*cs->diff_op[0].DO[i][j];
              }
            }
            val = dot_product(basisPhi_i, res) * wpg * det  ;
            AddLinearSolver(&cs->lsol,ino_fe*cs->nb_phy_vars,jno_fe*cs->nb_phy_vars,val);

            // Building U2 Matrix
            cs = &pb_pc->U2;
            for (int i=0; i<4; i++){
              res[i]=0;
              for (int j=0; j<4; j++){
                res[i]+=basisPhi_j[j]*cs->diff_op[0].DO[i][j];
              }
            }
            val = dot_product(basisPhi_i, res) * wpg * det  ;
            AddLinearSolver(&cs->lsol,ino_fe*cs->nb_phy_vars,jno_fe*cs->nb_phy_vars,val);


            // Building Schur Matrix
            cs = &pb_pc->D;
            for (int iv1=0; iv1<cs->nb_phy_vars; iv1++){
              for (int iv2=0; iv2<cs->nb_phy_vars; iv2++){
                for (int i=0; i<4; i++){
                  res[i]=0;
                  for (int j=0; j<4; j++){
                    res[i]+=basisPhi_j[j]*cs->diff_op[2*iv1+iv2].DO[i][j];
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
    for(int ivar=0;ivar<pb_pc->D.nb_phy_vars; ivar++){
      pb_pc->D.lsol.sol[i*pb_pc->D.nb_phy_vars+ivar]=0;
      pb_pc->D.lsol.rhs[i*pb_pc->D.nb_phy_vars+ivar]=0;
    }
     for(int ivar=0;ivar<pb_pc->L1.nb_phy_vars; ivar++){
      pb_pc->L1.lsol.sol[i*pb_pc->L1.nb_phy_vars+ivar]=0;
      pb_pc->L1.lsol.rhs[i*pb_pc->L1.nb_phy_vars+ivar]=0;
    }
     for(int ivar=0;ivar<pb_pc->L2.nb_phy_vars; ivar++){
      pb_pc->L2.lsol.sol[i*pb_pc->L2.nb_phy_vars+ivar]=0;
      pb_pc->L2.lsol.rhs[i*pb_pc->L2.nb_phy_vars+ivar]=0;
    }
     for(int ivar=0;ivar<pb_pc->U1.nb_phy_vars; ivar++){
      pb_pc->U1.lsol.sol[i*pb_pc->U1.nb_phy_vars+ivar]=0;
      pb_pc->U1.lsol.rhs[i*pb_pc->U1.nb_phy_vars+ivar]=0;
    }
     for(int ivar=0;ivar<pb_pc->U2.nb_phy_vars; ivar++){
      pb_pc->U2.lsol.sol[i*pb_pc->U2.nb_phy_vars+ivar]=0;
      pb_pc->U2.lsol.rhs[i*pb_pc->U2.nb_phy_vars+ivar]=0;
    }
     for(int ivar=0;ivar<pb_pc->Schur.nb_phy_vars; ivar++){
      pb_pc->Schur.lsol.sol[i*pb_pc->Schur.nb_phy_vars+ivar]=0;
      pb_pc->Schur.lsol.rhs[i*pb_pc->Schur.nb_phy_vars+ivar]=0;
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


void RobinFlux_SchurPressure(void * cs, real * xpg, real * w, real *vnorm, real * flux){
  ContinuousSolver * ps=cs;
  real lambda=-10000000;
  real Coef_diff=0;
  real p0=0;
  real Win[3];
    
  ps->simu->fd[0].model.ImposedData(xpg,ps->simu->dt,Win);
  p0=Win[0];

  Coef_diff=ps->simu->dt*ps->simu->vmax*ps->simu->theta;

  // u=- c dt theta nabla p
  // boudnary condition for wave lambda p-1/c (u,n)=g ==> * Coef_diff * (nabla p,n)+lambda p=g
  // ==> -(nabla p,n)= c* lambda / (Coeff*diff)-g *c /(Coeff*diff)
  // ==> -(Coeff*diff)**2 (nabla p,n)= c* lambda * (Coeff*diff)-g *c *(Coeff*diff)

  flux[0]=ps->simu->vmax * lambda* Coef_diff * w[0];//- ps->simu->vmax * Coef_diff * p0;
}

void Dirichlet_Velocity(void * cs, real * xpg, real * w, real *vnorm, real * flux){
  ContinuousSolver * ps=cs;
  real mu=1.0e14;//15;
  real Coef_diff=0;
  real u0_1 = 0;
  real u0_2 = 0;
  real Win[3];

  Coef_diff=ps->simu->dt*ps->simu->vmax*ps->simu->theta;

  flux[0]= mu * (w[0]*vnorm[0]*vnorm[0] + w[1]*vnorm[0]*vnorm[1])
   - mu * (u0_1*vnorm[0]*vnorm[0] + u0_2*vnorm[0]*vnorm[1]);
  flux[1]=  mu * (w[0]*vnorm[0]*vnorm[1] + w[1]*vnorm[1]*vnorm[1])
  -  mu * (u0_1*vnorm[0]*vnorm[1] + u0_2*vnorm[1]*vnorm[1]);
   
 
}


void BoundaryTerm_Xderivative(void * cs, real * xpg, real * w, real *vnorm, real * flux){
  ContinuousSolver * ps=cs;
  real lambda=0;
  real Coef_diff=0;
  real p0=0;

  Coef_diff=ps->simu->dt*ps->simu->vmax*ps->simu->theta;
  // term c dt theta n1
  flux[0]= Coef_diff * w[0] * vnorm[0];
}

void BoundaryTerm_Yderivative(void * cs, real * xpg, real * w, real *vnorm, real * flux){
  ContinuousSolver * ps=cs; 
  real lambda=0;
  real Coef_diff=0;
  real p0=0;

  Coef_diff=ps->simu->dt*ps->simu->vmax*ps->simu->theta;
  // term c dt theta n2
  flux[0]= Coef_diff * w[0] * vnorm[1];
}