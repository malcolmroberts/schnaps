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
  }

  flux[ld->index_rho]=0; 
  flux[ld->index_ux]=0; 
  flux[ld->index_uy]=0;
  flux[ld->index_uz]=0;
  flux[ld->index_temp]=0; 
  flux[ld->index_p]=0; 

};





