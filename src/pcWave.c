#include "solverpoisson.h"
#include "solvercontinuous.h"
#include "geometry.h"
#include "quantities_vp.h"
#include "linear_solver.h"
#include "physBased_PC.h"
#include <unistd.h>

void Init_PBPC_Wave_SchurVelocity_BCPressure(Simulation *simu, PB_PC* pb_pc, int* list_mat2assemble){
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
  pb_pc->mat_assembly = GenericOperator_PBPC;
  pb_pc->mat_assembly(pb_pc, NULL);

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

  // Choose solver 
  pb_pc->solvePC = PhyBased_PC_CG;

}

void Init_PBPC_Wave_SchurVelocity_BCVelocity(Simulation *simu, PB_PC* pb_pc, int* list_mat2assemble){
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
  pb_pc->mat_assembly = GenericOperator_PBPC;
  pb_pc->mat_assembly(pb_pc, NULL);

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

  // Choose solver 
  pb_pc->solvePC = PhyBased_PC_CG;

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
  pb_pc->mat_assembly = GenericOperator_PBPC;
  pb_pc->mat_assembly(pb_pc, NULL);

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

  // Choose solver 
  pb_pc->solvePC = PhyBased_PC_InvertSchur_CG;

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
  pb_pc->mat_assembly = GenericOperator_PBPC;
  pb_pc->mat_assembly(pb_pc, NULL);

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

  // Choose solver 
  pb_pc->solvePC = PhyBased_PC_InvertSchur_CG;
}

