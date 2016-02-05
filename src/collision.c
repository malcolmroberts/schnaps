#include "collision.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "interpolation.h"
#include "geometry.h"
#include "skyline.h"
#include "quantities_vp.h"



void VlasovP_Lagrangian_NumFlux(schnaps_real wL[],schnaps_real wR[],schnaps_real* vnorm,schnaps_real* flux){
  
  for(int i=0;i<_INDEX_MAX_KIN+1;i++){
    int j=i%_DEG_V; // local connectivity put in function
    int nel=i/_DEG_V; // element num (TODO : function)

    schnaps_real vn = (-_VMAX+nel*_DV +
		 _DV* glop(_DEG_V,j))*vnorm[0];
    
    schnaps_real vnp = vn>0 ? vn : 0;
    schnaps_real vnm = vn-vnp;

    flux[i] = vnp * wL[i] + vnm * wR[i];
  }
  // do not change the potential !
  // and the electric field
  flux[_INDEX_PHI]=0;  // flux for phi
  flux[_INDEX_EX]=0; // flux for E
  flux[_INDEX_RHO]=0; // flux for rho
  flux[_INDEX_VELOCITY]=0; // flux for u
  flux[_INDEX_PRESSURE]=0; // flux for p
  flux[_INDEX_TEMP]=0; // flux for e ou T

};


//! \brief compute compute the source term of the collision
//! model: electric force + true collisions
void VlasovP_Lagrangian_Source(const schnaps_real* x, const schnaps_real t, const schnaps_real* w, 
			       schnaps_real* source) {

  schnaps_real E=w[_INDEX_EX]; // electric field
  schnaps_real Md[_INDEX_MAX_KIN+1];
  schnaps_real db[_INDEX_MAX_KIN+1];
  for(int iv=0;iv<_INDEX_MAX_KIN+1;iv++){
    Md[iv]=0;
    db[iv]=0;
  }
  
  
  for(int iv=0;iv<_INDEX_MAX_KIN+1;iv++){
    source[iv]=0;
  }
  // no source on the potential for the moment
  source[_INDEX_PHI]=0;
  source[_INDEX_EX]=0;
  source[_INDEX_RHO]=0; //rho init
  source[_INDEX_VELOCITY]=0; // u init
  source[_INDEX_PRESSURE]=0; // p init
  source[_INDEX_TEMP]=0; 
  // loop on the finite emlements
  for(int iel=0;iel<_NB_ELEM_V;iel++){
    // loop on the local glops
    for(int kloc=0;kloc<_DEG_V+1;kloc++){
      schnaps_real omega=wglop(_DEG_V,kloc);
      int kpg=kloc+iel*_DEG_V;
      Md[kpg]+=omega*_DV;
      for(int iloc=0;iloc<_DEG_V+1;iloc++){
	int ipg=iloc+iel*_DEG_V;
	source[ipg]+=E*omega*w[kpg]*dlag(_DEG_V,iloc,kloc);
	if (iloc==kloc) db[ipg]+=E*omega*dlag(_DEG_V,iloc,kloc);
      }
    }
  }

  // upwinding
  if (E>0){
    source[_INDEX_MAX_KIN]-=E*w[_INDEX_MAX_KIN];
    db[_INDEX_MAX_KIN]-=E;
  }
  else {
    source[0]-=-E*w[0];
    db[0]-=-E;
  }

  for(int iv=0;iv<_INDEX_MAX_KIN+1;iv++){
    source[iv]/=Md[iv];
    //printf("%f ",source[iv]);
  }
  //printf("\n");
  //assert(1==2);
  

};



void VlasovP_Mass_modified(field *f,schnaps_real * w,void (*function)(field *f,schnaps_real w,schnaps_real *tw),schnaps_real* product){
  // give M^-1 * M_f(v)  
  schnaps_real tw;
  schnaps_real Mass[_INDEX_MAX_KIN+1];
  schnaps_real MassCollision[_INDEX_MAX];
  for(int iv=0;iv<_INDEX_MAX_KIN+1;iv++){
    Mass[iv]=0;

  }

  for(int iv=0;iv<_INDEX_MAX;iv++){
    MassCollision[iv]=0;
    product[iv]=0;
  }
  
  
  // loop on the finite emlements
  for(int iel=0;iel<_NB_ELEM_V;iel++){
    // loop on the local glops
    for(int kloc=0;kloc<_DEG_V+1;kloc++){
      schnaps_real omega=wglop(_DEG_V,kloc);
      int kpg=kloc+iel*_DEG_V;
      Mass[kpg]=omega*_DV;
      function(f,w[kpg],&tw);
      MassCollision[kpg]=omega*_DV*tw;
    }
  }

  for(int iv=0;iv<_INDEX_MAX_KIN+1;iv++){
    product[iv]=MassCollision[iv]/Mass[iv];

  }

  product[_INDEX_PHI]=w[_INDEX_PHI];
  product[_INDEX_EX]=w[_INDEX_EX];
  product[_INDEX_RHO]=w[_INDEX_RHO];
  product[_INDEX_VELOCITY]=w[_INDEX_VELOCITY];
  product[_INDEX_PRESSURE]=w[_INDEX_PRESSURE];
  product[_INDEX_TEMP]=w[_INDEX_TEMP];

};
