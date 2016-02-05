#include "collision.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "interpolation.h"
#include "geometry.h"
#include "skyline.h"
#include "quantities_vp.h"



/*void distribution_to_physic_entropy(field* f,real w,real *tw){
  *tw=log(w+1);
}

void physic_entropy_to_distribution(field* f,real w,real *tw){
  *tw=exp(w)-1;
  }*/


void Computation_charge_density(Simulation *simu){

  field * f=&simu->fd[0];
  
    for(int ipg=0;ipg<NPG(f->deg, f-> raf);ipg++){
      int imemc=f->varindex(f->deg, f->raf, f->model.m,ipg,_INDEX_RHO);
      simu->w[imemc]=0;
  
      for(int ielv=0;ielv<_NB_ELEM_V;ielv++){
	// loop on the local glops
	for(int iloc=0;iloc<_DEG_V+1;iloc++){
	  schnaps_real omega=wglop(_DEG_V,iloc);
	  schnaps_real vi=-_VMAX+ielv*_DV+_DV*glop(_DEG_V,iloc);
	  int ipgv=iloc+ielv*_DEG_V;
	  int imem=f->varindex(f->deg, f->raf, f->model.m,ipg,ipgv);
	  simu->w[imemc]+=omega*_DV*simu->w[imem];
	}
      }
    }
  
  
}


schnaps_real Computation_charge_average(Simulation *simu) {

  field * f=&simu->fd[0];
  schnaps_real average = 0;
  schnaps_real rho_imem = 0;
  schnaps_real size_domain = 0;


    // Loop on the glops (for numerical integration)
  const int npg = NPG(f->deg, f->raf);
    for(int ipg = 0; ipg < npg; ipg++) {
      int imem = f->varindex(f->deg, f->raf, f->model.m,
			     ipg, _INDEX_RHO);
	rho_imem = f->wn[imem];
      
      schnaps_real wpg, det;
      // Compute wpg, det, and the exact solution
      { 
	schnaps_real xphy[3], xpgref[3];
	schnaps_real dtau[3][3], codtau[3][3];
	// Get the coordinates of the Gauss point
	ref_pg_vol(f->deg, f->raf, ipg, xpgref, &wpg, NULL);
	schnaps_ref2phy(f->physnode, // phys. nodes
		xpgref, // xref
		NULL, -1, // dpsiref, ifa
		xphy, dtau, // xphy, dtau
		codtau, NULL, NULL); // codtau, dpsi, vnds
	det = dot_product(dtau[0], codtau[0]);
      }

        average += rho_imem * wpg * det;
	size_domain +=  wpg * det;

    }
    return average/size_domain;
}


void ComputeElectricField(field* f){

  int nraf[3] = {f->raf[0], 
		 f->raf[1],
		 f->raf[2]};
  
  int npg[3] = {f->deg[0] + 1, 
		f->deg[1] + 1,
		f->deg[2] + 1};
    
  int nbel =  nraf[0] * nraf[1] * nraf[2];
  int nnodes = npg[0] * npg[1] * npg[2] ;
  int npgmacrocell = nnodes * nraf[0] * nraf[1] * nraf[2];

  for (int ie = 0; ie < nbel; ie++){

    /* int iemacro = ie / (nraf[0] * nraf[1] * nraf[2]); */
    /* int isubcell = ie % (nraf[0] * nraf[1] * nraf[2]); */

    int iemacro = 0;
    int isubcell = ie; 

    // loop on the gauss points of the subcell
    for(int ipg = 0;ipg < nnodes; ipg++){
      //real wpg;
      schnaps_real xref[3];
      int ipgmacro= ipg + isubcell * nnodes;

      ref_pg_vol(f->deg,f->raf,ipgmacro,xref,NULL,NULL);
      int iex = f->varindex(f->deg,f->raf,f->model.m,
			    ipgmacro,_INDEX_EX);
      f->wn[iex] = 0;
      
      for(int ib=0; ib < nnodes; ib++){
	schnaps_real dtau[3][3],codtau[3][3];
	schnaps_real dphiref[3];
	schnaps_real dphi[3];
	int ibmacro = ib + isubcell * nnodes;
	grad_psi_pg(f->deg,f->raf,ibmacro,ipgmacro,dphiref);
	schnaps_ref2phy(f->physnode,xref,dphiref,0,NULL,
		  dtau,codtau,dphi,NULL);
	schnaps_real det = dot_product(dtau[0], codtau[0]);
	int ipot = f->varindex(f->deg,f->raf,f->model.m,
			   ibmacro,_INDEX_PHI);
	f->wn[iex] -= f->wn[ipot] * dphi[0] / det;
      }
    }
 

  }
}

/*void Compute_electric_field(field* f, real * w){

    
  int nraf[3] = {f->raf[0], 
		 f->raf[1],
		 f->raf[2]};
  
  int npg[3] = {f->deg[0] + 1, 
		f->deg[1] + 1,
		f->deg[2] + 1};

  int nnodes = npg[0] * npg[1] * npg[2] ;
 
  int npgmacrocell = nnodes * nraf[0] * nraf[1] * nraf[2];



    for(int ipg = 0;ipg < npgmacrocell; ipg++){
      //real wpg;
      real xref[3];

      ref_pg_vol(f->deg,f->raf,ipg,xref,NULL,NULL);
      int iex = f->varindex(f->deg,f->raf,f->model.m,
			    ipg,_INDEX_EX);
      w[iex] = 0;
      
      for(int ib=0; ib < npgmacrocell; ib++){
	real dtau[3][3],codtau[3][3];
	real dphiref[3];
	real dphi[3];
	grad_psi_pg(f->deg,f->raf,ib,ipg,dphiref);
	Ref2Phy(f->physnode,xref,dphiref,0,NULL,
		  dtau,codtau,dphi,NULL);
	real det = dot_product(dtau[0], codtau[0]);
	int ipot = f->varindex(f->deg,f->raf,f->model.m,
			   ib,_INDEX_PHI);
	w[iex] -= w[ipot] * dphi[0] / det;
      }
    }
  
    }*/




