#include "collision.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "interpolation.h"
#include "geometry.h"
#include "skyline.h"
#include "quantities_collision.h"

//! \brief compute square of velocity L2 error
//! \param[in] w : values of f at glops
double L2VelError(double* x,double t,double *w){
  
  FILE * ver;
  ver = fopen( "vel_error.dat", "w" );

  double wex[_MV+6];
  double err2=0;
  CollisionImposedData(x, t,wex);
  // loop on the finite emlements
  for(int iel=0;iel<_NB_ELEM_V;iel++){
    // loop on the local glops
    for(int iloc=0;iloc<_DEG_V+1;iloc++){
      double omega=wglop(_DEG_V,iloc);
      double vi=-_VMAX+iel*_DV+_DV*glop(_DEG_V,iloc);
      int ipg=iloc+iel*_DEG_V;
      err2+=omega*_DV*(w[ipg]-wex[ipg])*(w[ipg]-wex[ipg]);
      fprintf(ver,"%f %f %f % f\n",vi,w[ipg],wex[ipg],w[ipg]-wex[ipg]);
    }
  }
  fclose(ver);
  return err2;
};


double Collision_ImposedKinetic_Data(double x[3],double t,double v){
  double f;
  double pi=4*atan(1.);
  f=exp(-(v-t)*(v-t))*exp(-36*((x[0]-v*t+0.5*t*t)-0.5)*((x[0]-v*t+0.5*t*t)-0.5));
  return f;
};

double L2_Kinetic_error(Field* f){

  double error=0;
  double error_space=0;
  double moy=0; // mean value
  double moy_space=0;


  for (int ie=0;ie<f->macromesh.nbelems;ie++){
    // get the physical nodes of element ie
    double physnode[20][3];
    for(int inoloc=0;inoloc<20;inoloc++){
      int ino=f->macromesh.elem2node[20*ie+inoloc];
      physnode[inoloc][0]=f->macromesh.node[3*ino+0];
      physnode[inoloc][1]=f->macromesh.node[3*ino+1];
      physnode[inoloc][2]=f->macromesh.node[3*ino+2];
    }

    // loop on the glops (for numerical integration)
    for(int ipg=0;ipg<NPG(f->interp_param+1);ipg++){
      double xpgref[3],xphy[3],wpg;
      double dtau[3][3],codtau[3][3];//,xpg[3];
      // get the coordinates of the Gauss point
      ref_pg_vol(f->interp_param+1,ipg,xpgref,&wpg,NULL);
      Ref2Phy(physnode, // phys. nodes
	      xpgref,  // xref
	      NULL,-1, // dpsiref,ifa
	      xphy,dtau,  // xphy,dtau
	      codtau,NULL,NULL); // codtau,dpsi,vnds
      double det
	= dtau[0][0] * codtau[0][0]
	+ dtau[0][1] * codtau[0][1]
	+ dtau[0][2] * codtau[0][2]; 
      double w[f->model.m];
      for(int iv=0;iv<f->model.m;iv++){
	int imem=f->varindex(f->interp_param,ie,ipg,iv);
	w[iv]=f->wn[imem];
      }
      // get the exact value
      error+=L2VelError(xphy,f->tnow,w)*wpg*det;
    }
  }

  return sqrt(error);
  //return sqrt(error)/sqrt(moy);
}


void Computation_charge_density(double* x,double t,double *w){
  
  w[_MV+2]=0.;
  
  for(int iel=0;iel<_NB_ELEM_V;iel++){
    // loop on the local glops
    for(int iloc=0;iloc<_DEG_V+1;iloc++){
      double omega=wglop(_DEG_V,iloc);
      double vi=-_VMAX+iel*_DV+_DV*glop(_DEG_V,iloc);
      int ipg=iloc+iel*_DEG_V;
      w[_MV+2]+=omega*_DV*w[ipg];
    }
  }
  


void Compute_electric_field(Field* f){


  for (int ie=0;ie<f->macromesh.nbelems;ie++){
    // get the physical nodes of element ie
    double physnode[20][3];
    for(int inoloc=0;inoloc<20;inoloc++){
      int ino=f->macromesh.elem2node[20*ie+inoloc];
      physnode[inoloc][0]=f->macromesh.node[3*ino+0];
      physnode[inoloc][1]=f->macromesh.node[3*ino+1];
      physnode[inoloc][2]=f->macromesh.node[3*ino+2];
    }

    // loop on the glops (for numerical integration)
    for(int ipg=0;ipg<NPG(f->interp_param+1);ipg++){
      double xpgref[3],xphy[3],wpg;
      double dtau[3][3],codtau[3][3];//,xpg[3];
      // get the coordinates of the Gauss point
      ref_pg_vol(f->interp_param+1,ipg,xpgref,&wpg,NULL);
      Ref2Phy(physnode, // phys. nodes
	      xpgref,  // xref
	      NULL,-1, // dpsiref,ifa
	      xphy,dtau,  // xphy,dtau
	      codtau,NULL,NULL); // codtau,dpsi,vnds
      double det
	= dtau[0][0] * codtau[0][0]
	+ dtau[0][1] * codtau[0][1]
	+ dtau[0][2] * codtau[0][2]; 
      double w[f->model.m];
      for(int iv=0;iv<f->model.m;iv++){
	int imem=f->varindex(f->interp_param,ie,ipg,iv);
	w[iv]=f->wn[imem];
      }
      // get the exact value
      error+=L2VelError(xphy,f->tnow,w)*wpg*det;
    }
  }

  return sqrt(error);
  //return sqrt(error)/sqrt(moy);
}

  
}
