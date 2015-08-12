#include "solverpoisson.h"
#include "solvercontinuous.h"
#include "geometry.h"
#include "quantities_vp.h"
#include "linear_solver.h"


void Computation_ElectricField_Poisson(void * cs,LinearSolver* lsol){

  ContinuousSolver * ps=cs;
  
  ContinuousToDiscontinuous_Copy(ps,lsol);

  printf("Compute electric field...\n");

  
  for(int ie = 0; ie < ps->simu->macromesh.nbelems; ie++){
    ComputeElectricField(&ps->simu->fd[ie]);
  }
}


void RHSPoisson_Continuous(void * cs,LinearSolver* lsol){

  ContinuousSolver * ps=cs;
  
  field* f0 = &ps->simu->fd[0];

  real charge_average;
  charge_average=0;

  /*if(ps->type_bc == _Periodic_Poisson_BC){
    charge_average=Computation_charge_average(ps->simu,ps->w);
  }
  else {
    charge_average=0;
    }*/
  
  // right hand side assembly
  for(int ino = 0; ino < ps->nb_fe_dof; ino++){
    ps->lsol.rhs[ino] = 0;
  }

  real surf = 0;

  for(int ie = 0; ie < ps->nbel; ie++){  

    int iemacro = ie / (f0->raf[0] * f0->raf[1] * f0->raf[2]);
    int isubcell = ie % (f0->raf[0] * f0->raf[1] * f0->raf[2]);
    
    
 
    for(int iloc = 0; iloc < ps->nnodes; iloc++){
      real wpg;
      real xref[3];
      //int ipgmacro = ipg + isubcell * nnodes;
      int ilocmacro = iloc + isubcell * ps->nnodes;
      ref_pg_vol(f0->deg,f0->raf,ilocmacro,xref,&wpg,NULL);
      real dtau[3][3],codtau[3][3];
      Ref2Phy(ps->simu->fd[iemacro].physnode,
	      xref,NULL,0,NULL,
	      dtau,codtau,NULL,NULL);
      real det = dot_product(dtau[0], codtau[0]);	
      int ino_dg = iloc + ie * ps->nnodes;
      int ino_fe = ps->dg_to_fe_index[ino_dg];
      int imem = f0->varindex(f0->deg,f0->raf,f0->model.m,
			      ilocmacro,_INDEX_RHO) ;
	//+ iemacro * NPG(f0->deg,f0->raf) * f0->model.m ;
      real rho = ps->simu->fd[iemacro].wn[imem];
      ps->lsol.rhs[ino_fe] += (rho-charge_average)  * wpg * det ; // TODO: put the actual charge	
      surf += wpg * det ;  
    }
 
  }
}






