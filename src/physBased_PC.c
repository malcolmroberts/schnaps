#include "solverpoisson.h"
#include "solvercontinuous.h"
#include "geometry.h"
#include "quantities_vp.h"
#include "linear_solver.h"
#include "physBased_PC.h"

//! class managing Physics Based Preconditioners
typedef struct PB_PC{

  int list_mat2assemble;
  ContinuousSolver* D;
  ContinuousSolver* Schur;
  ContinuousSolver* L1;
  ContinuousSolver* L2;
  ContinuousSolver* U1;
  ContinuousSolver* U2;

} PB_PC;

//void InitPhy_Wave(Simulation *simu, PB_PC* pb_pc){


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
  for (int i=0;i<4;i++){
    for (int j=0;j<4;j++){
      pressionSolver.diff_op[i][j]=D[i][j];
    }
  }
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

  pressionSolver.lsol.solver_type=GMRES;
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

void VectorDgToCg(ContinuousSolver * cs,real * rhs, real * rhs_out){
  
  field* f = &cs->simu->fd[0];
  real* coeff = calloc(cs->nb_fe_nodes, sizeof(real));
  real* rhsCopy = calloc(cs->nb_dg_nodes*f->model.m, sizeof(real));

  for (int i=0;i<cs->nb_dg_nodes*f->model.m; i++) rhsCopy[i]=rhs[i];
  
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
    rhs_out[ino] = 0;
  }

  for(int ie = 0; ie < cs->nbel; ie++){  

    int iemacro = ie / (f->raf[0] * f->raf[1] * f->raf[2]);
    int isubcell = ie % (f->raf[0] * f->raf[1] * f->raf[2]);
    
 
    for(int iloc = 0; iloc < cs->nnodes; iloc++){
      real wpg;
      real xref[3];
      //int ipgmacro = ipg + isubcell * nnodes;
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

        // TODO it is not ino_dg and ino_fe because it not depend of the number of the variable 
        rhs_out[iv + cs->nb_phy_vars * ino_fe] += rhsCopy[imem] * wpg * det;//rhs[imem] * wpg * det; 
      }
    }
  }
  //for (int i=0; i<cs->nb_fe_nodes;i++){
  //  for (int iv=0; iv<cs->nb_phy_vars;iv++){
  //    rhs_out[iv+cs->nb_phy_vars*i]/=coeff[i];//*coeff[i];
  //  }
  //}
  free(coeff);
  free(rhsCopy);
}

void VectorCgToDg(ContinuousSolver * cs, real * rhsIn, real * rhsOut){

  field* f0 = &cs->simu->fd[0];

  //printf("Cg to Dg Copy...\n");
  
  //printf("%d %d %d\n", cs->list_of_var[0],cs->list_of_var[1],cs->list_of_var[2]);
  // copy the potential at the right place
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





