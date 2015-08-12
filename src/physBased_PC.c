#include "solverpoisson.h"
#include "solvercontinuous.h"
#include "geometry.h"
#include "quantities_vp.h"
#include "linear_solver.h"
#include "physBased_PC.h"

void InitPhy_Wave(Simulation *simu, PB_PC* pb_pc, int* list_mat2assemble){
  
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
      //printf("i %d,j %d,D %f,L1 %f,L2 %f,U1 %f,U2 %f,Schur %f\n",i,j,DMat[i][j],L1Mat[i][j],L2Mat[i][j],U1Mat[i][j],U2Mat[i][j],SchurMat[0][i][j]);
      if (DMat != NULL) pb_pc->D.diff_op[i][j]  = DMat[i][j];
      if (L1Mat != NULL) pb_pc->L1.diff_op[i][j] = L1Mat[i][j];
      if (L2Mat != NULL) pb_pc->L2.diff_op[i][j] = L2Mat[i][j];
      if (U1Mat != NULL) pb_pc->U1.diff_op[i][j] = U1Mat[i][j];
      if (U2Mat != NULL) pb_pc->U2.diff_op[i][j] = U2Mat[i][j];
    }
  }
}

void solvePhy_wave(PB_PC* pb_pc, Simulation *simu, real* globalSol, real*globalRHS_DG){
  
  // 0)1) Reset everything (needed for time evolution)
  reset(pb_pc);

  // 0)2) Second, initialize all penalization right-hand-sides
  pb_pc->D.bc_assembly=ExactDirichletContinuousMatrix;
  pb_pc->U1.bc_assembly=ExactDirichletContinuousMatrix;
  pb_pc->U2.bc_assembly=ExactDirichletContinuousMatrix;
  pb_pc->Schur.bc_assembly=ExactDirichletContinuousMatrix;
  pb_pc->D.bc_assembly(&pb_pc->D,&pb_pc->D.lsol);
  pb_pc->U1.bc_assembly(&pb_pc->U1,&pb_pc->U1.lsol);
  pb_pc->U2.bc_assembly(&pb_pc->U2,&pb_pc->U2.lsol);
  pb_pc->Schur.bc_assembly(&pb_pc->Schur,&pb_pc->Schur.lsol);
  
  // Applying these right-hand-sides to the preconditioner's own right-hand-sides
  for (int i=0; i<pb_pc->D.nb_fe_nodes; i++){
    pb_pc->rhs_prediction[i]      = pb_pc->D.lsol.rhs[i];
    pb_pc->rhs_propagation[i*2]   = pb_pc->Schur.lsol.rhs[i*2];
    pb_pc->rhs_propagation[i*2+1] = pb_pc->Schur.lsol.rhs[i*2+1];
    pb_pc->rhs_correction[i]      = pb_pc->D.lsol.rhs[i];
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
  real * globalRHS_CG = calloc(waveSolver.lsol.neq,sizeof(real));
  VectorDgToCg(&waveSolver, globalRHS_DG, globalRHS_CG);

  // 1) PREDICTION STEP

  pb_pc->D.lsol.solver_type=LU;
  pb_pc->D.lsol.MatVecProduct=MatVect;

  //printf("RHS assembly.....\n");
  for (int i=0;i<pb_pc->D.nb_fe_nodes;i++){
    pb_pc->D.lsol.rhs[i] = pb_pc->rhs_prediction[i] + globalRHS_CG[i*3];
  }
  //printf("Solution...\n");
  SolveLinearSolver(&pb_pc->D.lsol,simu);

  // 2) PROPAGATION STEP

  pb_pc->Schur.lsol.solver_type=LU;
  pb_pc->Schur.lsol.MatVecProduct=MatVect;
  pb_pc->L1.lsol.MatVecProduct=MatVect;
  pb_pc->L2.lsol.MatVecProduct=MatVect;
  // Parsing L1P, L2P into the "sol" of L1 and L2 (since it is unused).
  pb_pc->L1.lsol.MatVecProduct(&pb_pc->L1.lsol,pb_pc->D.lsol.sol,pb_pc->L1.lsol.sol);
  pb_pc->L2.lsol.MatVecProduct(&pb_pc->L2.lsol,pb_pc->D.lsol.sol,pb_pc->L2.lsol.sol);


  //printf("RHS assembly.....\n");
  for (int i=0;i<pb_pc->D.nb_fe_nodes;i++){
    pb_pc->Schur.lsol.rhs[i*2]   = pb_pc->rhs_propagation[i*2] + globalRHS_CG[i*3+1] - pb_pc->L1.lsol.sol[i];
    pb_pc->Schur.lsol.rhs[i*2+1] = pb_pc->rhs_propagation[i*2+1] + globalRHS_CG[i*3+2] - pb_pc->L2.lsol.sol[i];
  }

  //printf("Solution...\n");
  SolveLinearSolver(&pb_pc->Schur.lsol,simu);


  // 3) CORRECTION STEP

  // Extracting both U1 and U2 from the previous solution of the propagation step
  real *solU1=calloc(pb_pc->U1.nb_fe_nodes, sizeof(real));
  real *solU2=calloc(pb_pc->U2.nb_fe_nodes, sizeof(real));
  extract2CGVectors(&pb_pc->U1,&pb_pc->U2,pb_pc->Schur.lsol.sol,solU1,solU2);

  pb_pc->D.lsol.solver_type=LU;
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
  SolveLinearSolver(&pb_pc->D.lsol,simu);

  // 4) OUTPUT STEP

  // Final concatenation
  cat2CGVectors(&pb_pc->D,&pb_pc->Schur,pb_pc->D.lsol.sol,pb_pc->Schur.lsol.sol,waveSolver.lsol.sol);

  // Back from CG to DG
  VectorCgToDg(&waveSolver,waveSolver.lsol.sol,globalSol);

  // TODO: implement free for solvercontinuous, linear_solver...
  freeContinuousSolver(&waveSolver,0);
  free(globalRHS_CG);
  free(solU1);
  free(solU2);
}


void physicPC_wave(Simulation *simu, real* globalSol, real* globalRHS){


  // Global Solver
  ContinuousSolver waveSolver;
  int nb_var = 3;
  int * listvarGlobal = malloc(nb_var * sizeof(int));
  listvarGlobal[0]=0;
  listvarGlobal[1]=1;
  listvarGlobal[2]=2;
  InitContinuousSolver(&waveSolver,simu,1,nb_var,listvarGlobal);


  ////////////////////////////////// Prediction pressure step.
  ContinuousSolver pressionSolver;
  nb_var = 1;
  int * listvar = malloc(nb_var * sizeof(int));
  listvar[0]=0;
  
  InitContinuousSolver(&pressionSolver,simu,1,nb_var,listvar);
  pressionSolver.lsol.solver_type=LU;

  real D[4][4] = {{1,0,0,0},
                  {0,0,0,0},
                  {0,0,0,0},
                  {0,0,0,0}};
  pressionSolver.matrix_assembly=GenericOperatorScalar_Continuous;//MatrixPoisson_Continuous;//MassMatrix;
  pressionSolver.bc_assembly=ExactDirichletContinuousMatrix;
  pressionSolver.postcomputation_assembly=NULL;//Computation_ElectricField_Poisson;
  pressionSolver.lsol.MatVecProduct=MatVect;

  //printf("Matrix assembly.....\n");
  pressionSolver.matrix_assembly(&pressionSolver,&pressionSolver.lsol);

  //printf("BC assembly.....\n");
  pressionSolver.bc_assembly(&pressionSolver,&pressionSolver.lsol);

  //printf("RHS assembly.....\n");
  real * rhsPression = calloc(pressionSolver.lsol.neq,sizeof(real));
  VectorDgToCg(&pressionSolver, globalRHS, rhsPression);
  for (int i=0;i<pressionSolver.nb_fe_nodes;i++){
    pressionSolver.lsol.rhs[i]+=rhsPression[i];
  }

  //printf("Solution...\n");
  SolveLinearSolver(&pressionSolver.lsol,simu);

  pressionSolver.lsol.solver_type=LU;
  pressionSolver.lsol.tol=1.e-11;


  //////////////////////////// Construction of the L operator
  
  real h1=simu->vmax*simu->dt*simu->theta;
  real DL1[4][4] = {{0,h1,0,0},
                   { 0,0,0,0},
                   { 0,0,0,0},
                   { 0,0,0,0}};
  real DL2[4][4] = {{0,0,h1,0},
                   { 0,0,0,0},
                   { 0,0,0,0},
                   { 0,0,0,0}};
  real DL3[4][4] = {{0,0,0,h1},
                   { 0,0,0,0},
                   { 0,0,0,0},
                   { 0,0,0,0}};
  ContinuousSolver L1;
  ContinuousSolver L2;
  nb_var = 1;
  int * listvarL = malloc(nb_var * sizeof(int));
  listvarL[0]=0;

  InitContinuousSolver(&L1,simu,1,nb_var,listvarL);
  InitContinuousSolver(&L2,simu,1,nb_var,listvarL);
  for (int i=0;i<4;i++){
    for (int j=0;j<4;j++){
      L1.diff_op[i][j]=DL1[i][j];
      L2.diff_op[i][j]=DL2[i][j];
    }
  }
  L1.matrix_assembly=GenericOperatorScalar_Continuous;
  L2.matrix_assembly=GenericOperatorScalar_Continuous;
  L1.matrix_assembly(&L1,&L1.lsol);
  L2.matrix_assembly(&L2,&L2.lsol);


  real *L1Psolved=calloc(pressionSolver.nb_fe_nodes, sizeof(real));
  real *L2Psolved=calloc(pressionSolver.nb_fe_nodes, sizeof(real));
  L1.lsol.MatVecProduct=MatVect;
  L1.lsol.MatVecProduct(&L1.lsol,pressionSolver.lsol.sol,L1Psolved);
  L2.lsol.MatVecProduct=MatVect;
  L2.lsol.MatVecProduct(&L2.lsol,pressionSolver.lsol.sol,L2Psolved);


  real *LPsolved=calloc(2*pressionSolver.nb_fe_nodes, sizeof(real));
  catGradients(&L1,&L2,L1Psolved,L2Psolved,LPsolved);

  //////////////////////////// Computation of the velocity
  ContinuousSolver velocitySolver;
  nb_var = 2;
  int * listvar2 = malloc(nb_var * sizeof(int));
  listvar2[0]=1;
  listvar2[1]=2;
  
  InitContinuousSolver(&velocitySolver,simu,1,nb_var,listvar2);
  velocitySolver.lsol.solver_type=LU;
  real h=h1;
  real PSchur[4][4][4] = {    {{1,0,0,0},
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
  for (int k=0;k<4;k++){
    for (int i=0;i<4;i++){
      for (int j=0;j<4;j++){
        velocitySolver.diff_op2vec[k][i][j]=PSchur[k][i][j];
      }
    }
  }

  velocitySolver.matrix_assembly=GenericOperator2Vec_Continuous;
  velocitySolver.bc_assembly=ExactDirichletContinuousMatrix;
  velocitySolver.postcomputation_assembly=NULL;//Computation_ElectricField_Poisson;

  //printf("Matrix assembly.....\n");
  velocitySolver.matrix_assembly(&velocitySolver,&velocitySolver.lsol);

  //printf("BC assembly.....\n");
  velocitySolver.bc_assembly(&velocitySolver,&velocitySolver.lsol);

  //printf("RHS assembly.....\n");
  real * rhsVelocity = calloc(velocitySolver.lsol.neq,sizeof(real));
  VectorDgToCg(&velocitySolver, globalRHS, rhsVelocity);
  for (int i=0;i<2*pressionSolver.nb_fe_nodes;i++){
    velocitySolver.lsol.rhs[i]+=rhsVelocity[i]-LPsolved[i];
  }

  //printf("Solution...\n");
  SolveLinearSolver(&velocitySolver.lsol,simu);


  // Contruction of the operator U
  L1.list_of_var[0]=0;
  L2.list_of_var[0]=1;
  real *solU1=calloc(L1.nb_fe_nodes, sizeof(real));
  real *solU2=calloc(L2.nb_fe_nodes, sizeof(real));
  extract2CGVectors(&L1,&L2,velocitySolver.lsol.sol,solU1,solU2);
  real *L1u1=calloc(L1.nb_fe_nodes, sizeof(real));
  real *L2u2=calloc(L2.nb_fe_nodes, sizeof(real));
  real *Mp=calloc(pressionSolver.nb_fe_nodes, sizeof(real));
  
  L1.lsol.MatVecProduct(&L1.lsol,solU1,L1u1);
  L2.lsol.MatVecProduct(&L2.lsol,solU2,L2u2);
  pressionSolver.lsol.MatVecProduct(&pressionSolver.lsol,pressionSolver.lsol.sol,Mp);
  // Correction step
  for (int i=0; i<pressionSolver.nb_fe_nodes; i++){
    pressionSolver.lsol.rhs[i] = Mp[i] - L1u1[i] - L2u2[i];
  }
  
  SolveLinearSolver(&pressionSolver.lsol,simu);
  
  // Final concatenation
  cat2CGVectors(&pressionSolver,&velocitySolver,pressionSolver.lsol.sol,velocitySolver.lsol.sol,waveSolver.lsol.sol);

  // Back from CG to DG
  VectorCgToDg(&waveSolver,waveSolver.lsol.sol,globalSol);

  
  free(LPsolved);
  free(L1Psolved);
  free(L2Psolved);
  free(solU1);
  free(solU2);
  free(L1u1);
  free(L2u2);
  free(Mp);

}

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
  freeContinuousSolver(&pb_pc->D,1);
  freeContinuousSolver(&pb_pc->L1,0);
  freeContinuousSolver(&pb_pc->L2,0);
  freeContinuousSolver(&pb_pc->U1,0);
  freeContinuousSolver(&pb_pc->U2,0);
  freeContinuousSolver(&pb_pc->Schur,0);
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
