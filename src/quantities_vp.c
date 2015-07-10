#include "collision.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "interpolation.h"
#include "geometry.h"
#include "skyline.h"
#include "quantities_vp.h"



void distribution_to_physic_entropy(field* f,real w,real *tw){
  *tw=log(w+1);
}

void physic_entropy_to_distribution(field* f,real w,real *tw){
  *tw=exp(w)-1;
}


void Computation_charge_density(field *f, real * w){
  
  for(int ie=0;ie<f->macromesh.nbelems;ie++){
    for(int ipg=0;ipg<NPG(f->interp_param+1);ipg++){
      int imemc=f->varindex(f->interp_param,ie,ipg,_INDEX_RHO);
      w[imemc]=0;
  
      for(int ielv=0;ielv<_NB_ELEM_V;ielv++){
	// loop on the local glops
	for(int iloc=0;iloc<_DEG_V+1;iloc++){
	  real omega=wglop(_DEG_V,iloc);
	  real vi=-_VMAX+ielv*_DV+_DV*glop(_DEG_V,iloc);
	  int ipgv=iloc+ielv*_DEG_V;
	  int imem=f->varindex(f->interp_param,ie,ipg,ipgv);
	  w[imemc]+=omega*_DV*w[imem];
	}
      }
    }
  }
  
}


real Computation_charge_average(field *f,real * w) {
  //int param[8] = {f->model.m, _DEGX, _DEGY, _DEGZ, _RAFX, _RAFY, _RAFZ, 0};
  real average = 0;
  real rho_imem = 0;
  real size_domain = 0;

  for (int ie = 0; ie < f->macromesh.nbelems; ie++) {
    // Get the physical nodes of element ie
    real physnode[20][3];
    for(int inoloc = 0; inoloc < 20; inoloc++) {
      int ino = f->macromesh.elem2node[20*ie+inoloc];
      physnode[inoloc][0] = f->macromesh.node[3 * ino + 0];
      physnode[inoloc][1] = f->macromesh.node[3 * ino + 1];
      physnode[inoloc][2] = f->macromesh.node[3 * ino + 2];
    }

    // Loop on the glops (for numerical integration)
    const int npg = NPG(f->interp_param + 1);
    for(int ipg = 0; ipg < npg; ipg++) {
	int imem = f->varindex(f->interp_param, ie, ipg, _INDEX_RHO);
	rho_imem = f->wn[imem];
      
      real wpg, det;
      // Compute wpg, det, and the exact solution
      { 
	real xphy[3], xpgref[3];
	real dtau[3][3], codtau[3][3];
	// Get the coordinates of the Gauss point
	ref_pg_vol(f->interp_param + 1, ipg, xpgref, &wpg, NULL);
	Ref2Phy(physnode, // phys. nodes
		xpgref, // xref
		NULL, -1, // dpsiref, ifa
		xphy, dtau, // xphy, dtau
		codtau, NULL, NULL); // codtau, dpsi, vnds
	det = dot_product(dtau[0], codtau[0]);
      }


        average += rho_imem * wpg * det;
	size_domain +=  wpg * det;

      
    }
  }
  return average/size_domain;
}


void ComputeElectricField(field* f){

  int nraf[3] = {f->interp_param[4], 
		 f->interp_param[5],
		 f->interp_param[6]};
  
  int npg[3] = {f->interp_param[1] + 1, 
		f->interp_param[2] + 1,
		f->interp_param[3] + 1};
    
  int nbel = f->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2];

  int nnodes = npg[0] * npg[1] * npg[2] ;
 
  int npgmacrocell = nnodes * nraf[0] * nraf[1] * nraf[2];


  for (int ie = 0; ie < nbel; ie++){

    // get the physical nodes of element ie

    int iemacro = ie / (nraf[0] * nraf[1] * nraf[2]);
    int isubcell = ie % (nraf[0] * nraf[1] * nraf[2]);

    real physnode[20][3];
    for(int ino = 0; ino < 20; ino++) {
      int numnoe = f->macromesh.elem2node[20 * iemacro + ino];
      for(int ii = 0; ii < 3; ii++) {
	physnode[ino][ii] = f->macromesh.node[3 * numnoe + ii];
      }
    }

    for(int ipg = 0;ipg < nnodes; ipg++){
      //real wpg;
      real xref[3];
      int ipgmacro= ipg + isubcell * nnodes;

      ref_pg_vol(f->interp_param+1,ipgmacro,xref,NULL,NULL);
      int iex = f->varindex(f->interp_param,iemacro,
			    ipgmacro,_INDEX_EX);
      f->wn[iex] = 0;
      
      for(int ib=0; ib < nnodes; ib++){
	real dtau[3][3],codtau[3][3];
	real dphiref[3];
	real dphi[3];
	int ibmacro = ib + isubcell * nnodes;
	grad_psi_pg(f->interp_param+1,ibmacro,ipgmacro,dphiref);
	Ref2Phy(physnode,xref,dphiref,0,NULL,
		  dtau,codtau,dphi,NULL);
	real det = dot_product(dtau[0], codtau[0]);
	int ipot = f->varindex(f->interp_param,iemacro,
			   ibmacro,_INDEX_PHI);
	f->wn[iex] -= f->wn[ipot] * dphi[0] / det;
      }
    }
 

  }
}

void Compute_electric_field(field* f, real * w){

  int nraf[3] = {f->interp_param[4], 
		 f->interp_param[5],
		 f->interp_param[6]};
  
  int npg[3] = {f->interp_param[1] + 1, 
		f->interp_param[2] + 1,
		f->interp_param[3] + 1};
    
  int nnodes = npg[0] * npg[1] * npg[2] ;
 
  int npgmacrocell = nnodes * nraf[0] * nraf[1] * nraf[2];


  for (int ie=0;ie<f->macromesh.nbelems;ie++){
    // get the physical nodes of element ie
    real physnode[20][3];
    for(int inoloc=0;inoloc<20;inoloc++){
      int ino=f->macromesh.elem2node[20*ie+inoloc];
      physnode[inoloc][0]=f->macromesh.node[3*ino+0];
      physnode[inoloc][1]=f->macromesh.node[3*ino+1];
      physnode[inoloc][2]=f->macromesh.node[3*ino+2];
    }

    for(int ipg = 0;ipg < npgmacrocell; ipg++){
      //real wpg;
      real xref[3];

      ref_pg_vol(f->interp_param+1,ipg,xref,NULL,NULL);
      int iex = f->varindex(f->interp_param,ie,
			    ipg,_INDEX_EX);
      w[iex] = 0;
      
      for(int ib=0; ib < npgmacrocell; ib++){
	real dtau[3][3],codtau[3][3];
	real dphiref[3];
	real dphi[3];
	grad_psi_pg(f->interp_param+1,ib,ipg,dphiref);
	Ref2Phy(physnode,xref,dphiref,0,NULL,
		  dtau,codtau,dphi,NULL);
	real det = dot_product(dtau[0], codtau[0]);
	int ipot = f->varindex(f->interp_param,ie,
			   ib,_INDEX_PHI);
	w[iex] -= w[ipot] * dphi[0] / det;
      }
    }
  }
  
}  






