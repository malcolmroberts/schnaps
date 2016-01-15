#include "solverpoisson.h"
#include "solvercontinuous.h"
#include "geometry.h"
#include "quantities_vp.h"
#include "linear_solver.h"
#include "physBased_PC.h"
#include <unistd.h>
#include "pcSW.h"
#include "waterwave2d.h"

void Init_PBPC_SW_SchurVelocity_BCPressure(Simulation *simu, PB_PC* pb_pc, int* list_mat2assemble){
  // system is nonlinear
  pb_pc->nonlinear=1;


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

  pb_pc->D.diff_op = calloc(pb_pc->D.nb_phy_vars*pb_pc->D.nb_phy_vars,sizeof(SDO));
  pb_pc->Schur.diff_op = calloc(pb_pc->Schur.nb_phy_vars*pb_pc->Schur.nb_phy_vars,sizeof(SDO));
  pb_pc->L1.diff_op = calloc(pb_pc->L1.nb_phy_vars*pb_pc->L1.nb_phy_vars,sizeof(SDO));
  pb_pc->L2.diff_op = calloc(pb_pc->L2.nb_phy_vars*pb_pc->L2.nb_phy_vars,sizeof(SDO));
  pb_pc->U1.diff_op = calloc(pb_pc->U1.nb_phy_vars*pb_pc->U1.nb_phy_vars,sizeof(SDO));
  pb_pc->U2.diff_op = calloc(pb_pc->U2.nb_phy_vars*pb_pc->U2.nb_phy_vars,sizeof(SDO));

  // Assigning function to build all operators' matrices
  pb_pc->mat_assembly = GenericOperator_PBPC_NonLinear;
  pb_pc->loc_mat_assembly = Schur_ESF;
  pb_pc->rhs_assembly = GenericRHS_PBPC_NonLinear;
  pb_pc->loc_rhs_assembly = SW_RHS;
  pb_pc->bc_assembly = BoundaryConditionFriedrichsAssembly;
  pb_pc->loc_bc_assembly = Wave_BC_normalvelocity_null;
  pb_pc->source_assembly = Source_Assembly;

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

void Init_PBPC_SW_SchurVelocity_BCVelocity(Simulation *simu, PB_PC* pb_pc, int* list_mat2assemble){
  // system is nonlinear
  pb_pc->nonlinear=1;


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

  pb_pc->D.diff_op = calloc(pb_pc->D.nb_phy_vars*pb_pc->D.nb_phy_vars,sizeof(SDO));
  pb_pc->Schur.diff_op = calloc(pb_pc->Schur.nb_phy_vars*pb_pc->Schur.nb_phy_vars,sizeof(SDO));
  pb_pc->L1.diff_op = calloc(pb_pc->L1.nb_phy_vars*pb_pc->L1.nb_phy_vars,sizeof(SDO));
  pb_pc->L2.diff_op = calloc(pb_pc->L2.nb_phy_vars*pb_pc->L2.nb_phy_vars,sizeof(SDO));
  pb_pc->U1.diff_op = calloc(pb_pc->U1.nb_phy_vars*pb_pc->U1.nb_phy_vars,sizeof(SDO));
  pb_pc->U2.diff_op = calloc(pb_pc->U2.nb_phy_vars*pb_pc->U2.nb_phy_vars,sizeof(SDO));

  // Assigning function to build all operators' matrices
  pb_pc->mat_assembly = GenericOperator_PBPC_NonLinear;
  pb_pc->loc_mat_assembly = Schur_ESF;
  pb_pc->rhs_assembly = GenericRHS_PBPC_NonLinear;
  pb_pc->loc_rhs_assembly = SW_RHS;
  pb_pc->bc_assembly = BoundaryConditionFriedrichsAssembly;
  pb_pc->loc_bc_assembly = Wave_BC_normalvelocity_null;
  pb_pc->source_assembly = Source_Assembly;

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

void Schur_ESF(void* pb_pc, real* var){

  PB_PC* pc = (PB_PC*) pb_pc;

  real theta = pc->D.simu->theta;
  real dt = pc->D.simu->dt;
  real g = _GRAVITY;

  real h  = var[0];
  real hx = var[1];
  real hy = var[2];
  real hz = var[3];
  real u  = var[4];
  real ux = var[5];
  real uy = var[6];
  real uz = var[7];
  real v  = var[8];
  real vx = var[9];
  real vy = var[10];
  real vz = var[11];

  real DMat[4][4]  = {{1.0+theta*dt*(ux+vy),theta*dt*u,theta*dt*v,0},
                      {0,0,0,0},
                      {0,0,0,0},
                      {0,0,0,0}};
  real L1Mat[4][4] = {{u+theta*dt*g*hx+theta*dt*(u*ux+v*uy),theta*dt*g*h,0,0},
                      {0,0,0,0},
                      {0,0,0,0},
                      {0,0,0,0}};
  real L2Mat[4][4] = {{v+theta*dt*g*hy+theta*dt*(u*vx+v*vy),0,theta*dt*g*h,0},
                      {0,0,0,0},
                      {0,0,0,0},
                      {0,0,0,0}};
  real U1Mat[4][4] = {{theta*dt*hx,theta*dt*h,0,0},
                      {0,0,0,0},
                      {0,0,0,0},
                      {0,0,0,0}};
  real U2Mat[4][4] = {{theta*dt*hy,0,theta*dt*h,0},
                      {0,0,0,0},
                      {0,0,0,0},
                      {0,0,0,0}};
  real SchurMat[4][4][4] ={{{h + theta*dt*h*ux - theta*dt*(u*ux+v*uy)*theta*dt*hx 
                               , theta*dt*h*u  - theta*dt*(u*ux+v*uy)*theta*dt*h 
                               , theta*dt*h*v , 0},
                            {theta*theta*dt*dt*g*h*hx,theta*theta*dt*dt*g*h*h,0,0},
                            {0,0,0,0},
                            {0,0,0,0}},

                           {{theta*dt*h*uy - theta*dt*(u*ux+v*uy)*theta*dt*hy,
                             0,            - theta*dt*(u*ux+v*vy)*theta*dt*h ,0},
                            {theta*theta*dt*dt*g*h*hy,0,theta*theta*dt*dt*g*h*h,0},
                            {0,0,0,0},
                            {0,0,0,0}},

                           {{theta*dt*h*vx - theta*dt*(u*vx+v*vy)*theta*dt*hx,
                             0,            - theta*dt*(u*vx+v*vy)*theta*dt*h ,0},
                            {0,0,0,0},
                            {theta*theta*dt*dt*g*h*hx,theta*theta*dt*dt*g*h*h,0,0},
                            {0,0,0,0}},

                           {{h + theta*dt*h*vy - theta*dt*(u*vx+v*vy)*theta*dt*hy
                               , theta*dt*h*u
                               , theta*dt*h*v  - theta*dt*(u*vx+v*vy)*theta*dt*h , 0},
                            {0,0,0,0},
                            {theta*theta*dt*dt*g*h*hy,0,theta*theta*dt*dt*g*h*h,0},
                            {0,0,0,0}}};

  //real SchurMat[4][4][4] ={{{h + theta*dt*h*ux - (u+theta*dt*(u*ux+v*uy))*theta*dt*hx 
  //                             , theta*dt*h*u  - (u+theta*dt*(u*ux+v*uy))*theta*dt*h 
  //                             , theta*dt*h*v , 0},
  //                          {theta*theta*dt*dt*h*hx,theta*theta*dt*dt*h*h,0,0},
  //                          {0,0,0,0},
  //                          {0,0,0,0}},
  //                         {{theta*dt*h*uy - (u+theta*dt*(u*ux+v*uy))*theta*dt*hy,
  //                           0,            - (u+theta*dt*(u*ux+v*vy))*theta*dt*h ,0},
  //                          {theta*theta*dt*dt*h*hy,0,theta*theta*dt*dt*h*h,0},
  //                          {0,0,0,0},
  //                          {0,0,0,0}},
  //                         {{theta*dt*h*vx - (v+theta*dt*(u*vx+v*vy))*theta*dt*hx,
  //                           0,            - (v+theta*dt*(u*vx+v*vy))*theta*dt*h ,0},
  //                          {0,0,0,0},
  //                          {theta*theta*dt*dt*h*hx,theta*theta*dt*dt*h*h,0,0},
  //                          {0,0,0,0}},
  //                         {{h + theta*dt*h*vy - (v+theta*dt*(u*vx+v*vy))*theta*dt*hy
  //                             , theta*dt*h*u
  //                             , theta*dt*h*v  - (v+theta*dt*(u*vx+v*vy))*theta*dt*h , 0},
  //                          {0,0,0,0},
  //                          {theta*theta*dt*dt*h*hy,0,theta*theta*dt*dt*h*h,0},
  //                          {0,0,0,0}}};

  //real SchurMat[4][4][4] ={{{h + theta*dt*h*ux - (u+theta*dt*(u*ux+v*uy))*theta*dt*hx 
  //                             , theta*dt*h*u  - (u+theta*dt*(u*ux+v*uy))*theta*dt*h 
  //                             , theta*dt*h*v , 0},
  //                          {theta*theta*dt*dt*h*hx,theta*theta*dt*dt*h*h,0,0},
  //                          {0,0,0,0},
  //                          {0,0,0,0}},
  //                         {{theta*dt*h*uy,0,0,0},
  //                          {theta*theta*dt*dt*h*hy,0,theta*theta*dt*dt*h*h,0},
  //                          {0,0,0,0},
  //                          {0,0,0,0}},
  //                         {{theta*dt*h*vx,0,0,0},
  //                          {0,0,0,0},
  //                          {theta*theta*dt*dt*h*hx,theta*theta*dt*dt*h*h,0,0},
  //                          {0,0,0,0}},
  //                         {{h + theta*dt*h*vy - (v+theta*dt*(u*vx+v*vy))*theta*dt*hy
  //                             , theta*dt*h*u  - (v+theta*dt*(u*vx+v*vy))*theta*dt*h
  //                             , theta*dt*h*v , 0},
  //                          {0,0,0,0},
  //                          {theta*theta*dt*dt*h*hy,0,theta*theta*dt*dt*h*h,0},
  //                          {0,0,0,0}}};

  // Applying these matrices inside the different ContinuousSolver
  for (int i=0; i<4; i++){
    for (int j=0; j<4; j++){
      for (int k=0; k<4; k++){
        if (SchurMat != NULL) pc->Schur.diff_op[k].DO[i][j] = SchurMat[k][i][j];
      }
      if (DMat  != NULL)  pc->D.diff_op[0].DO[i][j] = DMat[i][j];
      if (L1Mat != NULL) pc->L1.diff_op[0].DO[i][j] = L1Mat[i][j];
      if (L2Mat != NULL) pc->L2.diff_op[0].DO[i][j] = L2Mat[i][j];
      if (U1Mat != NULL) pc->U1.diff_op[0].DO[i][j] = U1Mat[i][j];
      if (U2Mat != NULL) pc->U2.diff_op[0].DO[i][j] = U2Mat[i][j];
    }
  }
}

void Schur_ASF(void* pb_pc, real* var){

  PB_PC* pc = (PB_PC*) pb_pc;

  real theta = pc->D.simu->theta;
  real dt = pc->D.simu->dt;
  real g = _GRAVITY;

  real h  = var[0];
  real hx = var[1];
  real hy = var[2];
  real hz = var[3];
  real u  = var[4];
  real ux = var[5];
  real uy = var[6];
  real uz = var[7];
  real v  = var[8];
  real vx = var[9];
  real vy = var[10];
  real vz = var[11];

  real DMat[4][4]  = {{1.0+theta*dt*(ux+vy),theta*dt*u,theta*dt*v,0},
                      {0,0,0,0},
                      {0,0,0,0},
                      {0,0,0,0}};
  real L1Mat[4][4] = {{u+theta*dt*g*hx+theta*dt*(u*ux+v*uy),theta*dt*g*h,0,0},
                      {0,0,0,0},
                      {0,0,0,0},
                      {0,0,0,0}};
  real L2Mat[4][4] = {{v+theta*dt*g*hy+theta*dt*(u*vx+v*vy),0,theta*dt*g*h,0},
                      {0,0,0,0},
                      {0,0,0,0},
                      {0,0,0,0}};
  real U1Mat[4][4] = {{theta*dt*hx,theta*dt*h,0,0},
                      {0,0,0,0},
                      {0,0,0,0},
                      {0,0,0,0}};
  real U2Mat[4][4] = {{theta*dt*hy,0,theta*dt*h,0},
                      {0,0,0,0},
                      {0,0,0,0},
                      {0,0,0,0}};
  real SchurMat[4][4][4] ={{{h + theta*dt*h*ux
                               , theta*dt*h*u
                               , theta*dt*h*v , 0},
                            {theta*theta*dt*dt*h*hx,theta*theta*dt*dt*h*h,0,0},
                            {0,0,0,0},
                            {0,0,0,0}},
                           {{theta*dt*h*uy,0,0,0},
                            {theta*theta*dt*dt*h*hy,0,theta*theta*dt*dt*h*h,0},
                            {0,0,0,0},
                            {0,0,0,0}},
                           {{theta*dt*h*vx,0,0,0},
                            {0,0,0,0},
                            {theta*theta*dt*dt*h*hx,theta*theta*dt*dt*h*h,0,0},
                            {0,0,0,0}},
                           {{h + theta*dt*h*vy
                               , theta*dt*h*u
                               , theta*dt*h*v , 0},
                            {0,0,0,0},
                            {theta*theta*dt*dt*h*hy,0,theta*theta*dt*dt*h*h,0},
                            {0,0,0,0}}};

  // Applying these matrices inside the different ContinuousSolver
  for (int i=0; i<4; i++){
    for (int j=0; j<4; j++){
      for (int k=0; k<4; k++){
        if (SchurMat != NULL) pc->Schur.diff_op[k].DO[i][j] = SchurMat[k][i][j];
      }
      if (DMat  != NULL)  pc->D.diff_op[0].DO[i][j] = DMat[i][j];
      if (L1Mat != NULL) pc->L1.diff_op[0].DO[i][j] = L1Mat[i][j];
      if (L2Mat != NULL) pc->L2.diff_op[0].DO[i][j] = L2Mat[i][j];
      if (U1Mat != NULL) pc->U1.diff_op[0].DO[i][j] = U1Mat[i][j];
      if (U2Mat != NULL) pc->U2.diff_op[0].DO[i][j] = U2Mat[i][j];
    }
  }
}

// coeff = dt*det*wpg*basisPhi_i[0]
void SW_RHS(void* pb_pc, real* var, real coeff, real* loc_rhs){
  
  PB_PC* pc = (PB_PC*) pb_pc;

  real g = _GRAVITY;
  real h  = var[0];
  real hx = var[1];
  real hy = var[2];
  real hz = var[3];
  real u  = var[4];
  real ux = var[5];
  real uy = var[6];
  real uz = var[7];
  real v  = var[8];
  real vx = var[9];
  real vy = var[10];
  real vz = var[11];

  loc_rhs[0] = - (h*(ux+vy)+u*hx+v*hy+0.0*hz) * coeff;
  loc_rhs[1] = -   h * ( u*ux + v*uy + g*hx ) * coeff;
  loc_rhs[2] = -   h * ( u*vx + v*vy + g*hy ) * coeff;
  //printf("h=%.8e,hx=%.8e,hy=%.8e,u=%.8e,ux=%.8e,uy=%.8e,v=%.8e,vx=%.8e,vy=%.8e\n",h,hx,hy,u,ux,uy,v,vx,vy);
  //printf("loc_rhs = %.8e, %.8e, %.8e\n",loc_rhs[0], loc_rhs[1], loc_rhs[2]);
  
  
}
