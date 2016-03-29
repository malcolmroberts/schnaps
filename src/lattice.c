#include "lattice.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "interpolation.h"
#include "geometry.h"
#include "skyline.h"

/**************************************************************************************/
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
    //printf("i :%i \t Flux :%e \n",i, flux[i]); 
  }
  flux[ld->index_rho]=0; 
  flux[ld->index_ux]=0; 
  flux[ld->index_uy]=0;
  flux[ld->index_uz]=0;
  flux[ld->index_temp]=0; 
  flux[ld->index_p]=0; 
}
/**************************************************************************************/
void Compute_distribution_eq(Simulation *simu, schnaps_real * w_eq) {
  LatticeData * ld=&schnaps_lattice_data;

  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
    field * f = simu->fd+ie;
    for(int ipg=0;ipg<NPG(f->deg, f-> raf);ipg++){
      	int irho=f->varindex(f->deg, f->raf, f->model.m, ipg, ld->index_rho);
        int iux=f->varindex(f->deg, f->raf, f->model.m, ipg, ld->index_ux);
        int iuy=f->varindex(f->deg, f->raf, f->model.m, ipg, ld->index_uy);
        int itemp=f->varindex(f->deg, f->raf, f->model.m, ipg, ld->index_temp);
        schnaps_real rho = f->wn[irho];
        schnaps_real temp = f->wn[itemp];
        schnaps_real  ux = f->wn[iux];
        schnaps_real uy = f->wn[iuy];
        schnaps_real uz= 0.0;
        schnaps_real p = 0.0;
        //schnaps_real uz = f->wn[iuz];
        //schnaps_real u2=  ux*ux + uy*uy;
        for(int iv=0;iv<ld->index_max_q+1;iv++){
          int ikin=f->varindex(f->deg, f->raf, f->model.m, ipg, iv);
          //schnaps_real uv = ux * ld->q_tab[iv][0] + uy * ld->q_tab[iv][1];
          //schnaps_real temp_val=ld->w_tab[iv]*rho*(1.0+uv+0.5*uv*uv-0.5*t*u2);	  
          schnaps_real temp_val=ld->feq(iv,ld, rho,ux,uy,uz,temp,p);
          w_eq[ikin]= temp_val;
        }
    } 
  }
}
//
/**************************************************************************************/
schnaps_real feq_isothermal_D2Q9(int i_node,void *lattice,schnaps_real rho,schnaps_real ux,schnaps_real uy,schnaps_real uz,schnaps_real temp,schnaps_real p){
    LatticeData *ld= lattice;
    schnaps_real u2= (ux * ux + uy * uy)/temp;
    schnaps_real uv= (ux *ld->q_tab[i_node][0] + uy * ld->q_tab[i_node][1])/temp;
    schnaps_real feq= ld->w_tab[i_node] * rho * (1.0+uv+0.5 * (uv * uv- u2));
    return feq;
}
/**************************************************************************************/
schnaps_real feq_isothermal_linearwave_D2Q9(int i_node,void *lattice,schnaps_real rho,schnaps_real ux,schnaps_real uy,schnaps_real uz,schnaps_real temp,schnaps_real p){
    LatticeData *ld= lattice;
    schnaps_real uv= (ux *ld->q_tab[i_node][0] + ld->q_tab[i_node][1]* uy)/temp;
    schnaps_real feq= ld->w_tab[i_node]* rho * (1.0+uv);
    return feq;
}
/**************************************************************************************/
void Compute_moments(Simulation *simu) {
  LatticeData * ld=&schnaps_lattice_data;
  //
  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
    field * f = simu->fd+ie;
    //
    for(int ipg=0;ipg<NPG(f->deg, f-> raf);ipg++){
      	int irho=f->varindex(f->deg, f->raf, f->model.m, ipg, ld->index_rho);
        int iux=f->varindex(f->deg, f->raf, f->model.m, ipg, ld->index_ux);
        int iuy=f->varindex(f->deg, f->raf, f->model.m, ipg, ld->index_uy);
        int iuz=f->varindex(f->deg, f->raf, f->model.m, ipg, ld->index_uz);
        int it=f->varindex(f->deg, f->raf, f->model.m, ipg, ld->index_temp);
      f->wn[irho]=0.0;
      f->wn[iux]=0.0;
      f->wn[iuy]=0.0;
      f->wn[iuz]=0.0;
      for(int iv=0;iv<ld->index_max_q+1;iv++){
        int ikin=f->varindex(f->deg, f->raf, f->model.m, ipg, iv);
        //
        f->wn[irho]+=f->wn[ikin];
        f->wn[iux]+=f->wn[ikin]*ld->q_tab[iv][0];
        f->wn[iuy]+=f->wn[ikin]*ld->q_tab[iv][1];
      };
	    f->wn[iux]=f->wn[iux]/f->wn[irho];
	    f->wn[iuy]=f->wn[iuy]/f->wn[irho];
    }
  }
}
/**************************************************************************************/

void Compute_relaxation(Simulation *simu, schnaps_real * w_eq) {
  LatticeData * ld=&schnaps_lattice_data;
  schnaps_real nu=0;
  nu=simu->dt/(ld->tau+0.5*simu->dt);
  //
  for(int ie = 0; ie < simu->macromesh.nbelems; ++ie) {
    field * f = simu->fd+ie;
    for(int ipg=0;ipg<NPG(f->deg, f-> raf);ipg++){
  	for(int iv=0;iv<ld->index_max_q+1;iv++){
      int ikin=f->varindex(f->deg, f->raf, f->model.m, ipg, iv);
      f->wn[ikin]=f->wn[ikin]-nu*(f->wn[ikin]-w_eq[ikin]);
    };
    };
  };
}
