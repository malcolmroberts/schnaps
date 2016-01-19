#include "collision.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "interpolation.h"
#include "geometry.h"
#include "skyline.h"
#include "quantities_vp.h"
#include "diagnostics_vp.h"

real L2VelError(field *f, real *x, real *w){

  real wex[_INDEX_MAX];
  real err2 = 0;
  real t = f->tnow;
  f->model.ImposedData(x, t, wex);
  // loop on the finite emlements
  for(int iel = 0; iel < _NB_ELEM_V; iel++){
    // loop on the local glops
    for(int iloc = 0; iloc < _DEG_V + 1; iloc++){
      real omega = wglop(_DEG_V, iloc);
      real vi = -_VMAX + iel*_DV + _DV * glop(_DEG_V, iloc);
      int ipg = iloc + iel * _DEG_V;
      err2 += omega * _DV * (w[ipg] - wex[ipg]) * (w[ipg] - wex[ipg]);
    }
  }
  
  return err2;
}

real L2_Kinetic_error(Simulation* simu){
  real error = 0;

  field *f =&simu->fd[0];

    // get the physical nodes of element ie

    // loop on the glops (for numerical integration)
  for(int ipg = 0; ipg < NPG(f->deg,f->raf); ipg++){
      real xpgref[3], xphy[3], wpg;
      real dtau[3][3], codtau[3][3];//,xpg[3];
      // get the coordinates of the Gauss point
      ref_pg_vol(f->deg,f->raf, ipg, xpgref, &wpg, NULL);
      Ref2Phy(f->physnode, // phys. nodes
	      xpgref,  // xref
	      NULL, -1, // dpsiref,ifa
	      xphy, dtau,  // xphy,dtau
	      codtau, NULL, NULL); // codtau,dpsi,vnds
      real det
	= dtau[0][0] * codtau[0][0]
	+ dtau[0][1] * codtau[0][1]
	+ dtau[0][2] * codtau[0][2]; 
      real w[f->model.m];
      for(int iv = 0;iv < f->model.m; iv++){
	int imem = f->varindex(f->deg,f->raf, f->model.m, ipg, iv);
	w[iv] = f->wn[imem];
      }
      // get the exact value
      error += L2VelError(f, xphy, w) * wpg * det;
    }
  
  return sqrt(error);
}

real local_kinetic_energy(field *f,real *x, real *w) {

  real wex[_INDEX_MAX];
  real l_ke=0;
  real t=f->tnow;
  
  // loop on the finite emlements
  for(int iel = 0; iel < _NB_ELEM_V; iel++){
    // loop on the local glops
    for(int iloc = 0; iloc < _DEG_V + 1; iloc++){
      real omega = wglop(_DEG_V, iloc);
      real vi = -_VMAX + iel * _DV + _DV * glop(_DEG_V, iloc);
      int ipg = iloc + iel * _DEG_V;
      l_ke += omega * _DV * w[ipg] * vi * vi ;
     }
  }
  return l_ke;
}


// TODO: do not store all diagnotics for all time, but instead just
// append to the output file.
void Energies(Simulation *simu, real *w, real k_energy, real e_energy, real t_energy,int first_diag) {
  
  k_energy = 0;
  e_energy = 0;
  t_energy = 0;

  field * f=&simu->fd[0];


    // loop on the glops (for numerical integration)
  for(int ipg = 0; ipg < NPG(f->deg,f->raf); ipg++){
      real xpgref[3], xphy[3], wpg;
      real dtau[3][3], codtau[3][3];//,xpg[3];
      // get the coordinates of the Gauss point
      ref_pg_vol(f->deg,f->raf, ipg, xpgref, &wpg, NULL);
      Ref2Phy(f->physnode, // phys. nodes
	      xpgref,  // xref
	      NULL,-1, // dpsiref,ifa
	      xphy,dtau,  // xphy,dtau
	      codtau,NULL,NULL); // codtau,dpsi,vnds
      real det
	= dtau[0][0] * codtau[0][0]
	+ dtau[0][1] * codtau[0][1]
	+ dtau[0][2] * codtau[0][2]; 
      real wn[f->model.m];
      for(int iv = 0; iv < _INDEX_MAX + 1; iv++){ 
	int imem = f->varindex(f->deg,f->raf, f->model.m, ipg, iv);
	wn[iv] = w[imem];
      }
      // get the exact value
      k_energy += local_kinetic_energy(f, xphy, wn) * wpg * det;
      e_energy += wn[_INDEX_EX] * wn[_INDEX_EX] * wpg * det;
  
    }
  
  
  t_energy = 0.5 * (e_energy + k_energy);
  
   simu->Diagnostics[simu->iter_time_rk + (first_diag-1) * simu->itermax_rk] = 0.5 * k_energy; 
   simu->Diagnostics[simu->iter_time_rk + (first_diag) * simu->itermax_rk] = 0.5 * e_energy; 
   simu->Diagnostics[simu->iter_time_rk + (first_diag+1) * simu->itermax_rk] = t_energy; 
}

void Charge_total(Simulation *simu, real *w, real t_charge,int first_diag) {
  
  t_charge=0;
  field *f=&simu->fd[0];

    // loop on the glops (for numerical integration)
  for(int ipg = 0; ipg < NPG(f->deg,f->raf); ipg++){
      real xpgref[3], xphy[3], wpg;
      real dtau[3][3], codtau[3][3];//,xpg[3];
      // get the coordinates of the Gauss point
      ref_pg_vol(f->deg,f->raf, ipg, xpgref, &wpg, NULL);
      Ref2Phy(f->physnode, // phys. nodes
	      xpgref,  // xref
	      NULL,-1, // dpsiref,ifa
	      xphy,dtau,  // xphy,dtau
	      codtau,NULL,NULL); // codtau,dpsi,vnds
      real det
	= dtau[0][0] * codtau[0][0]
	+ dtau[0][1] * codtau[0][1]
	+ dtau[0][2] * codtau[0][2]; 
      real wn[f->model.m];
      for(int iv = 0; iv < _INDEX_MAX + 1; iv++){ 
	int imem = f->varindex(f->deg,f->raf, f->model.m, ipg, iv);
	wn[iv] = w[imem];
      }
      t_charge += wn[_INDEX_RHO] * wpg * det;
    }
  

    simu->Diagnostics[simu->iter_time_rk + (first_diag+2) * simu->itermax_rk] = t_charge;
}


void Plot_Energies(Simulation *simu, real dt) { 
  int nb_diag = 0; 
  real e_energy = 0, k_energy = 0, t_energy = 0,t_charge=0; 
  FILE *Plot;
  Plot = fopen("Diagnostics.dat","w");

  for(int i = 1; i < simu->itermax_rk + 1; i++){
    simu->tnow = i*dt;
    k_energy = simu->Diagnostics[i]; 
    e_energy = simu->Diagnostics[i + simu->itermax_rk]; 
    t_energy = simu->Diagnostics[i + 2 * simu->itermax_rk]; 
    t_charge= simu->Diagnostics[i + 3 * simu->itermax_rk]; 
    fprintf(Plot, "%.11e %.11e %.11e %.11e %.15e\n", simu->tnow, k_energy, e_energy, t_energy,t_charge); 
  } 
  fclose(Plot); 
} 
