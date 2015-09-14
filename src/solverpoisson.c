#include "solverpoisson.h"
#include "solvercontinuous.h"
#include "geometry.h"
#include "quantities_vp.h"
#include "linear_solver.h"


void Computation_ElectricField_Poisson(void *cs, LinearSolver *lsol)
{
  ContinuousSolver *ps = cs;
  
  ContinuousToDiscontinuous_Copy(ps, lsol);

  printf("Compute electric field...\n");
  
  ComputeElectricField(ps->f);
}

void MatrixPoisson_Continuous(void *cs, LinearSolver *lsol)
{
  ContinuousSolver *ps = cs;

  field *f = ps->f;

  if(!ps->lsol.mat_is_assembly) {
    for(int ie = 0; ie < ps->nbel; ie++) {  

      // local matrix 
      real aloc[ps->nnodes][ps->nnodes];
      for(int iloc = 0; iloc < ps->nnodes; iloc++) {
	for(int jloc = 0; jloc < ps->nnodes; jloc++) {
	  aloc[iloc][jloc] = 0;
	}
      }

      int iemacro = ie / (f->raf[0] * f->raf[1] * f->raf[2]);
      int isubcell = ie % (f->raf[0] * f->raf[1] * f->raf[2]);

      /* real physnode[20][3]; */
      /* for(int ino = 0; ino < 20; ino++) { */
      /* 	int numnoe = ps->f->macromesh.elem2node[20 * iemacro + ino]; */
      /* 	for(int ii = 0; ii < 3; ii++) { */
      /* 	  physnode[ino][ii] = ps->f->macromesh.node[3 * numnoe + ii]; */
      /* 	} */
      /* } */
      //ref_pg_vol(ps->fd->interp_param+1,int ipg,
      // real* xpg,real* wpg,real* xpg_in);
      // grad_psi_pg(ps->fd->interp_param+1,ib,ipg,dphiref)
      // Ref2Phy(physnode,xref,dphiref,NULL,NULL,dtau,codtau,dphi,NULL);

      for(int ipg = 0;ipg < ps->nnodes; ipg++) {
	real wpg;
	real xref[3];
	int ipgmacro= ipg + isubcell * ps->nnodes;

	ref_pg_vol(f->deg, f->raf, ipgmacro, xref, &wpg, NULL);

	for(int iloc = 0; iloc < ps->nnodes; iloc++){
	  real dtau[3][3],codtau[3][3];
	  real dphiref_i[3],dphiref_j[3];
	  real dphi_i[3],dphi_j[3];
	  int ilocmacro = iloc + isubcell * ps->nnodes;
	  grad_psi_pg(f->deg, f->raf, ilocmacro, ipgmacro, dphiref_i);

	  real physnode[20][3];
	  for(int inoloc = 0; inoloc < 20; inoloc++) {
	    int ino = f->macromesh.elem2node[20 * iemacro + inoloc];
	    physnode[inoloc][0] = f->macromesh.node[3 * ino + 0];
	    physnode[inoloc][1] = f->macromesh.node[3 * ino + 1];
	    physnode[inoloc][2] = f->macromesh.node[3 * ino + 2];
	  }

	  Ref2Phy(physnode,
		  xref,dphiref_i,0,NULL,
		  dtau,codtau,dphi_i,NULL);
	  real det = dot_product(dtau[0], codtau[0]);
	  for(int jloc = 0; jloc < ps->nnodes; jloc++){
	    int jlocmacro = jloc + isubcell * ps->nnodes;
	    grad_psi_pg(f->deg, f->raf, jlocmacro, ipgmacro,dphiref_j);
	    Ref2Phy(physnode,
		    xref, dphiref_j, 0, NULL,
		    dtau, codtau, dphi_j, NULL);
	    aloc[iloc][jloc] += dot_product(dphi_i, dphi_j) * wpg / det  ;
	  }
	}
      }

      for(int iloc = 0; iloc < ps->nnodes; iloc++) {
	for(int jloc = 0; jloc < ps->nnodes; jloc++) {
	  int ino_dg = iloc + ie * ps->nnodes;
	  int jno_dg = jloc + ie * ps->nnodes;
	  int ino_fe = ps->dg_to_fe_index[ino_dg];
	  int jno_fe = ps->dg_to_fe_index[jno_dg];
	  real val = aloc[iloc][jloc];
	  AddLinearSolver(&ps->lsol, ino_fe, jno_fe, val);
	}
      }
   
    }
  } 
  
}

void RHSPoisson_Continuous(void *cs, LinearSolver* lsol)
{
  ContinuousSolver *ps = cs;
  field *f = ps->f;

  real charge_average;
  charge_average=0;

  if(ps->type_bc == _Periodic_Poisson_BC){
    //charge_average=Computation_charge_average(f,w);
  } else {
    charge_average=0;
  }
  
  // right hand side assembly
  for(int ino = 0; ino < ps->nb_fe_dof; ino++){
    ps->lsol.rhs[ino] = 0;
  }

  real surf = 0;

  for(int ie = 0; ie < ps->nbel; ie++){  

    int iemacro = ie / (f->raf[0] * f->raf[1] * f->raf[2]);
    int isubcell = ie % (f->raf[0] * f->raf[1] * f->raf[2]);
    
     for(int iloc = 0; iloc < ps->nnodes; iloc++){
      real wpg;
      real xref[3];
      //int ipgmacro = ipg + isubcell * nnodes;
      int ilocmacro = iloc + isubcell * ps->nnodes;
      ref_pg_vol(f->deg,f->raf,ilocmacro,xref,&wpg,NULL);
      real dtau[3][3],codtau[3][3];
      
      real physnode[20][3];
      for(int inoloc = 0; inoloc < 20; inoloc++) {
	int ino = f->macromesh.elem2node[20 * iemacro + inoloc];
	physnode[inoloc][0] = f->macromesh.node[3 * ino + 0];
	physnode[inoloc][1] = f->macromesh.node[3 * ino + 1];
	physnode[inoloc][2] = f->macromesh.node[3 * ino + 2];
      }

      Ref2Phy(physnode,
	      xref,NULL,0,NULL,
	      dtau,codtau,NULL,NULL);
      real det = dot_product(dtau[0], codtau[0]);	
      int ino_dg = iloc + ie * ps->nnodes;
      int ino_fe = ps->dg_to_fe_index[ino_dg];
      int imem = f->varindex(f->deg,f->raf,f->model.m,
			     iemacro, ilocmacro,_INDEX_RHO) ;
      //+ iemacro * NPG(f->deg,f->raf) * f->model.m ;
      real rho = ps->f[iemacro].wn[imem];
      ps->lsol.rhs[ino_fe] += (rho-charge_average)  * wpg * det ; // TODO: put the actual charge	
      surf += wpg * det ;  
    }
 
  }
}


