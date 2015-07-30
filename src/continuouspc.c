#include "solverpoisson.h"
#include "solvercontinuous.h"
#include "geometry.h"
#include "quantities_vp.h"
#include "linear_solver.h"
#include "continuouspc.h"


void physicPC_wave(Simulation *simu, real* globalSol, real* globalRHS){

  ContinuousSolver pressionSolver;
  int nb_var = 1;
  int * listvar = malloc(nb_var * sizeof(int));
  listvar[0]=0;
  
  InitContinuousSolver(&pressionSolver,simu,1,nb_var,listvar);
  pressionSolver.matrix_assembly=MassMatrix;
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





  ContinuousSolver velocitySolver;
  nb_var = 2;
  int * listvar2 = malloc(nb_var * sizeof(int));
  listvar2[0]=1;
  listvar2[1]=2;
  
  InitContinuousSolver(&velocitySolver,simu,1,nb_var,listvar2);
  velocitySolver.matrix_assembly=MassMatrix;
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


void MassMatrix(void* ps, LinearSolver *lsol){

  ContinuousSolver *cs = ps;
   field* fg = &cs->simu->fd[0];

  if(!cs->lsol.mat_is_assembly){
    for(int ie = 0; ie < cs->nbel; ie++){  

      // local matrix 
      real aloc[cs->nnodes*cs->nb_phy_vars][cs->nnodes*cs->nb_phy_vars];
      for(int iloc = 0; iloc < cs->nnodes*cs->nb_phy_vars; iloc++){
	for(int jloc = 0; jloc < cs->nnodes*cs->nb_phy_vars; jloc++){
	  aloc[iloc][jloc] = 0;
	}
      }

      int iemacro = ie / (fg->raf[0] * fg->raf[1] * fg->raf[2]);
      int isubcell = ie % (fg->raf[0] * fg->raf[1] * fg->raf[2]);

  //    for(int ipg = 0;ipg < cs->nnodes; ipg++){
	  //real wpg;
	  //real xref[3];
	//int ipgmacro= ipg + isubcell * cs->nnodes;

	//ref_pg_vol(fg->deg,f0->raf,ipgmacro,xref,&wpg,NULL);

	for(int iloc = 0; iloc < cs->nnodes; iloc++){
	  real wpg;
	  real xref[3];
	  real dtau[3][3],codtau[3][3];
	  int ilocmacro = iloc + isubcell * cs->nnodes;
    ref_pg_vol(fg->deg,fg->raf,ilocmacro,xref,&wpg,NULL);
	  Ref2Phy(cs->simu->fd[iemacro].physnode,
		  xref,NULL,0,NULL,
		  dtau,codtau,NULL,NULL);
	  real det = dot_product(dtau[0], codtau[0]);
	  for (int iv=0; iv<cs->nb_phy_vars;iv++){
	    aloc[iloc*cs->nb_phy_vars+iv][iloc*cs->nb_phy_vars+iv] += wpg * det  ;
    }
//	}
      }


      for(int iloc = 0; iloc < cs->nnodes; iloc++){
	for(int jloc = 0; jloc < cs->nnodes; jloc++){
	  int ino_dg = iloc + ie * cs->nnodes;
	  int jno_dg = jloc + ie * cs->nnodes;
	  int ino_fe = cs->dg_to_fe_index[ino_dg];
	  int jno_fe = cs->dg_to_fe_index[jno_dg];
	  for (int iv=0;iv<cs->nb_phy_vars;iv++){
	    real val = aloc[iloc*cs->nb_phy_vars+iv][jloc*cs->nb_phy_vars+iv];
	    AddLinearSolver(&cs->lsol,ino_fe*cs->nb_phy_vars+iv,jno_fe*cs->nb_phy_vars+iv,val);
    }
	}
      }
   
    }
  } 
  
 }



