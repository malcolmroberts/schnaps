#include "solverpoisson.h"
#include "solvercontinuous.h"
#include "geometry.h"
#include "quantities_vp.h"
#include "linear_solver.h"


void VectorDgToCg(ContinuousSolver * ps,real * rhs){
  
  field* f0 = &ps->simu->fd[0];

  
  // right hand side assembly
  for(int ino = 0; ino < ps->nb_fe_dof; ino++){
    ps->lsol.rhs[ino] = 0;
  }

  real surf = 0;

  for(int ie = 0; ie < ps->nbel; ie++){  

    int iemacro = ie / (f0->raf[0] * f0->raf[1] * f0->raf[2]);
    int isubcell = ie % (f0->raf[0] * f0->raf[1] * f0->raf[2]);
    
    int iv =1; // TOOTOOTOTOTOTODO loop on variable
 
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
       int iv =1; // TOOTOOTOTOTOTODO loop on variable
      int imem = f0->varindex(f0->deg,f0->raf,nb_phy_vars,
			      ilocmacro,iv) ;

      // TODO it is not ino_dg and ino_fe because it not depend of the number of the variable 
      ps->lsol.rhs[iv + cs->nb_phy_vars * ino_fe] += rhs[imem]; // TODO: put the actual charge	
      surf += wpg * det ;  
    }
 
  }
}





