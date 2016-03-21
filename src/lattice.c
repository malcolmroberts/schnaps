#include "lattice.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "interpolation.h"
#include "geometry.h"
#include "skyline.h"



void Lattice_NumFlux(schnaps_real wL[],schnaps_real wR[],schnaps_real* vnorm,schnaps_real* flux){
 LatticeData * ld=&schnaps_lattice_data;
  for(int i = 0;i < ld->index_max_q+1;i++){
    schnaps_real vn=0;
    for(int dim = 0;dim < ld->d; dim++){
      vn += ld->q_tab[i][dim]*vnorm[dim];
    }
    schnaps_real vnp = vn>0 ? vn : 0;
    schnaps_real vnm = vn-vnp;
    flux[i] = vnp * wL[i] + vnm * wR[i];
    //flux[i]=0.0;
  }
  flux[ld->index_rho]=0; 
  flux[ld->index_ux]=0; 
  flux[ld->index_uy]=0;
  flux[ld->index_uz]=0;
  flux[ld->index_temp]=0; 
  flux[ld->index_p]=0; 
};





void Compute_distribution_eq(Simulation *simu, double * w_eq) {
  LatticeData * ld=&schnaps_lattice_data;

  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
    field * f = simu->fd+ie;
    
    for(int ipg=0;ipg<NPG(f->deg, f-> raf);ipg++){
      	int irho=f->varindex(f->deg, f->raf, f->model.m, ipg, ld->index_rho);
	int iux=f->varindex(f->deg, f->raf, f->model.m, ipg, ld->index_ux);
	int iuy=f->varindex(f->deg, f->raf, f->model.m, ipg, ld->index_uy);
	int it=f->varindex(f->deg, f->raf, f->model.m, ipg, ld->index_temp);

	double rho = f->wn[irho];
	double t = f->wn[it];
	double ux = f->wn[iux]/t;
	double uy = f->wn[iuy]/t;
	double u2=  ux*ux + uy*uy;
		
	for(int iv=0;iv<ld->index_max_q+1;iv++){
	  int ikin=f->varindex(f->deg, f->raf, f->model.m, ipg, iv);
	  
	  double uv = ux * ld->q_tab[iv][0] + uy * ld->q_tab[iv][1];
	  
	  w_eq[ikin]=ld->w_tab[iv]*rho*(1.0+uv+0.5*uv*uv-0.5*t*u2);	  
	}
    }
  }
}

void Compute_moments(Simulation *simu) {
  LatticeData * ld=&schnaps_lattice_data;

  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
    field * f = simu->fd+ie;
    
    for(int ipg=0;ipg<NPG(f->deg, f-> raf);ipg++){
      	int irho=f->varindex(f->deg, f->raf, f->model.m, ipg, ld->index_rho);
	int iux=f->varindex(f->deg, f->raf, f->model.m, ipg, ld->index_ux);
	int iuy=f->varindex(f->deg, f->raf, f->model.m, ipg, ld->index_uy);
	int iuz=f->varindex(f->deg, f->raf, f->model.m, ipg, ld->index_uz);
	int it=f->varindex(f->deg, f->raf, f->model.m, ipg, ld->index_temp);
	
	f->wn[irho]=0;
	f->wn[iux]=0;
	f->wn[iuy]=0;
	f->wn[iuz]=0;
	
	for(int iv=0;iv<ld->index_max_q+1;iv++){
	  int ikin=f->varindex(f->deg, f->raf, f->model.m, ipg, iv);
	  
	  f->wn[irho]+=f->wn[ikin];
	  f->wn[iux]+=f->wn[ikin]*ld->q_tab[iv][0];
	  f->wn[iuy]+=f->wn[ikin]*ld->q_tab[iv][1];
	}
        
	f->wn[iux]=f->wn[iux]/f->wn[irho];
	f->wn[iuy]=f->wn[iuy]/f->wn[irho];
    }
  }
}


void Compute_relaxation(Simulation *simu, schnaps_real * w_eq) {
  LatticeData * ld=&schnaps_lattice_data;
  double nu=0;
  nu=simu->dt/(ld->tau+0.5*simu->dt);
  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
    field * f = simu->fd+ie;
    
    for(int ipg=0;ipg<NPG(f->deg, f-> raf);ipg++){
	
	for(int iv=0;iv<ld->index_max_q+1;iv++){
	  int ikin=f->varindex(f->deg, f->raf, f->model.m, ipg, iv);
	  
	  f->wn[ikin]=f->wn[ikin]-nu*(f->wn[ikin]-w_eq[ikin]);	  

	}
    }
  }
}
