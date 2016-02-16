#include "collision.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "interpolation.h"
#include "geometry.h"
#include "skyline.h"
#include "quantities_vp.h"
#include "collision.h"




void Computation_charge_density(Simulation *simu){

  KineticData *kd = &schnaps_kinetic_data;
  for(int ie = 0; ie < simu->macromesh.nbelems; ie++){
    field * f = simu->fd + ie; 
  
    for(int ipg=0;ipg<NPG(f->deg, f-> raf);ipg++){
      int imemc=f->varindex(f->deg, f->raf, f->model.m,ipg,kd->index_rho);
      f->wn[imemc]=0;
  
      for(int ielv=0;ielv<kd->nb_elem_v;ielv++){
	// loop on the local glops
	for(int iloc=0;iloc<kd->deg_v+1;iloc++){
	  schnaps_real omega=wglop(kd->deg_v,iloc);
	  schnaps_real vi=-kd->vmax+ielv*kd->dv+kd->dv*glop(kd->deg_v,iloc);
	  int ipgv=iloc+ielv*kd->deg_v;
	  int imem=f->varindex(f->deg, f->raf, f->model.m,ipg,ipgv);
	  f->wn[imemc]+=omega*kd->dv*simu->w[imem];
	}
      }
    }
  }
  
}

void Computation_Fluid_Quantities(Simulation *simu){

  KineticData *kd = &schnaps_kinetic_data;

  for(int ie = 0; ie < simu->macromesh.nbelems; ie++){
    field * f = simu->fd + ie; 
  
    for(int ipg=0;ipg<NPG(f->deg, f-> raf);ipg++){
      int imem_rho=f->varindex(f->deg, f->raf, f->model.m,ipg,kd->index_rho);
      int imem_U=f->varindex(f->deg, f->raf, f->model.m,ipg,kd->index_u);
      int imem_P=f->varindex(f->deg, f->raf, f->model.m,ipg,kd->index_P);
      int imem_T=f->varindex(f->deg, f->raf, f->model.m,ipg,kd->index_T);
      f->wn[imem_rho]=0;
      f->wn[imem_U]=0;
      f->wn[imem_T]=0;
      f->wn[imem_P]=0;

      schnaps_real rhoU=0,rho=0,U=0,tensor_P=0;
      for(int ielv=0;ielv<kd->nb_elem_v;ielv++){
	// loop on the local glops
	for(int iloc=0;iloc<kd->deg_v+1;iloc++){
	  schnaps_real omega=wglop(kd->deg_v,iloc);
	  schnaps_real vi=-kd->vmax+ielv*kd->dv+kd->dv*glop(kd->deg_v,iloc);
	  int ipgv=iloc+ielv*kd->deg_v;
	  int imem=f->varindex(f->deg, f->raf, f->model.m,ipg,ipgv);
	  rho+=omega*kd->dv*simu->w[imem];
	  rhoU+=omega*kd->dv*vi*simu->w[imem];
	  tensor_P+=omega*kd->dv*vi*vi*simu->w[imem];
	}
      }
    
      f->wn[imem_rho]=rho;
      f->wn[imem_U]=rhoU/rho;
      f->wn[imem_T]=tensor_P/rho;
      f->wn[imem_P]=(0.5*tensor_P-0.5*rhoU*rhoU/rho)*(kd->gamma-1);
           
    }
  }
  
}

schnaps_real Computation_Maxwellian(schnaps_real rho, schnaps_real U, schnaps_real T, schnaps_real v){ 

  schnaps_real maxw;
  schnaps_real my_pi= 4.0*atan(1.0);
  
  maxw= (rho/pow(2*my_pi*T,0.5))*exp(-pow(U-v,2.0)/(2.0*T));

  return maxw;

}
  
schnaps_real Computation_charge_average(Simulation *simu) {
   KineticData *kd = &schnaps_kinetic_data;
  field * f=&simu->fd[0];
  schnaps_real average = 0;
  schnaps_real rho_imem = 0;
  schnaps_real size_domain = 0;


    // Loop on the glops (for numerical integration)
  const int npg = NPG(f->deg, f->raf);
    for(int ipg = 0; ipg < npg; ipg++) {
      int imem = f->varindex(f->deg, f->raf, f->model.m,
			     ipg, kd->index_rho);
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
   KineticData *kd = &schnaps_kinetic_data;
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
			    ipgmacro,kd->index_ex);
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
			   ibmacro,kd->index_phi);	
	f->wn[iex] -= f->wn[ipot] * dphi[0] / det;
      }
    }
 

  }
}





void Collision_Source(Simulation *simu) {
  KineticData * kd=&schnaps_kinetic_data;
 

  Computation_Fluid_Quantities(simu);

  for(int ie = 0; ie < simu->macromesh.nbelems; ie++){
    field * f = simu->fd + ie;
    
    schnaps_real w_loc[f->model.m];
    schnaps_real source_loc[f->model.m];
    
    for(int ipg=0;ipg<NPG(f->deg, f-> raf);ipg++){
      for(int iv=0;iv<f->model.m;iv++){
	w_loc[iv] = simu->w[iv];
      }
      BGK_Source(w_loc, source_loc);  
      for(int iv=0;iv<f->model.m;iv++){
	int imemc=f->varindex(f->deg, f->raf, f->model.m, ipg, 0);
	simu->w[imemc] = simu->w[imemc]+ simu->dt*source_loc[iv];
      }

      
    }  
  }
};
