#include "solverpoisson.h"
#include "solvercontinuous.h"
#include "geometry.h"
#include "quantities_vp.h"
#include "linear_solver.h"
#include "continuouspc.h"


void physicPC_wave(Simulation *simu, real* globalSol, real* globalRHS){


  // Prediction pressure step.
  ContinuousSolver pressionSolver;
  int nb_var = 1;
  int * listvar = malloc(nb_var * sizeof(int));
  listvar[0]=0;
  
  InitContinuousSolver(&pressionSolver,simu,1,nb_var,listvar);

  real D[4][4] = {{1,0,0,0},
                  {0,0,0,0},
                  {0,0,0,0},
                  {0,0,0,0}};
  for (int i=0;i<4;i++){
    for (int j=0;j<4;j++){
      pressionSolver.diff_op[i][j]=D[i][j];
    }
  }
  pressionSolver.matrix_assembly=GenericOperatorScalar_Continuous;//MassMatrix;
  pressionSolver.bc_assembly=NULL;// ExactDirichletContinuousMatrix;
  pressionSolver.postcomputation_assembly=NULL;//Computation_ElectricField_Poisson;

  printf("Matrix assembly.....\n");
  pressionSolver.matrix_assembly(&pressionSolver,&pressionSolver.lsol);

  printf("RHS assembly.....\n");
  VectorDgToCg(&pressionSolver,globalRHS);

  printf("Solution...\n");
  SolveLinearSolver(&pressionSolver.lsol);

  for (int i=0; i<pressionSolver.lsol.neq;i++){
    printf("Pouet %d, %f\n",i,pressionSolver.lsol.sol[i]);
  }

  // Construction of the L operator
  
  real DL1[4][4] = {{0,0,0,0},
                   {-1,0,0,0},
                   { 0,0,0,0},
                   { 0,0,0,0}};
  real DL2[4][4] = {{0,0,0,0},
                   { 0,0,0,0},
                   {-1,0,0,0},
                   { 0,0,0,0}};
  real DL3[4][4] = {{0,0,0,0},
                   { 0,0,0,0},
                   { 0,0,0,0},
                   {-1,0,0,0}};
  
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
  printf("Matrix assembly.....\n");
  L1.matrix_assembly(&L1,&L1.lsol);
  L2.matrix_assembly(&L2,&L2.lsol);
  // cs.lsol->MAtVECT L1Psolved, L2Psolved
  // Cat L1Psolved, L2Psolved
  // SUM WITH RU (plus bas)
  // Check
  //lsol->MatVecProduct=MatVect;
  //lsol->MatVecProduct(lsol,loc_x,loc_z);

  // Computation of the velocity
  ContinuousSolver velocitySolver;
  nb_var = 2;
  int * listvar2 = malloc(nb_var * sizeof(int));
  listvar2[0]=1;
  listvar2[1]=2;
  
  InitContinuousSolver(&velocitySolver,simu,1,nb_var,listvar2);
  real D2[4][4][4] = {    {{1,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0}},
                          {{0,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0}},
                          {{0,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0}},
                          {{1,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0}}};
                          
  for (int k=0;k<4;k++){
    for (int i=0;i<4;i++){
      for (int j=0;j<4;j++){
        velocitySolver.diff_op2vec[k][i][j]=D2[k][i][j];
      }
    }
  }

  velocitySolver.matrix_assembly=GenericOperator2Vec_Continuous;
  velocitySolver.bc_assembly=NULL;// ExactDirichletContinuousMatrix;
  velocitySolver.postcomputation_assembly=NULL;//Computation_ElectricField_Poisson;

  printf("Matrix assembly.....\n");
  velocitySolver.matrix_assembly(&velocitySolver,&velocitySolver.lsol);

  printf("RHS assembly.....\n");
  VectorDgToCg(&velocitySolver,globalRHS);

  printf("Solution...\n");
  SolveLinearSolver(&velocitySolver.lsol);
  for (int i=0; i<velocitySolver.lsol.neq;i++){
    printf("Pouet 2 %d, %f\n",i,velocitySolver.lsol.sol[i]);
  }
}

void VectorDgToCg(ContinuousSolver * cs,real * rhs){
  
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
    cs->lsol.rhs[ino] = 0;
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
        cs->lsol.rhs[iv + cs->nb_phy_vars * ino_fe] += rhsCopy[imem] * wpg * det; 
      }
    }
 
  }
  //for (int i=0; i<cs->nb_fe_nodes;i++){
  //  for (int iv=0; iv<cs->nb_phy_vars;iv++){
  //    cs->lsol.rhs[iv+cs->nb_phy_vars*i]/=coeff[i];//*coeff[i];
  //  }
  //}
  free(coeff);
  free(rhsCopy);
}

void VectorCgToDg(ContinuousSolver * cs,real * rhs){

  field* f0 = &cs->simu->fd[0];

  printf("Copy...\n");

  
  // copy the potential at the right place
  for(int var =0; var < cs->nb_phy_vars; var++){ 
    for(int ie = 0; ie < cs->simu->macromesh.nbelems; ie++){  
      for(int ipg = 0;ipg < cs->npgmacrocell; ipg++){
	int ino_dg = ipg + ie * cs->npgmacrocell;
	int ino_fe = cs->dg_to_fe_index[ino_dg];
	int ipot = f0->varindex(f0->deg,f0->raf,f0->model.m,
				ipg,cs->list_of_var[var]);
	rhs[ipot]=cs->lsol.sol[ino_fe];
	}
    }
  }

}





