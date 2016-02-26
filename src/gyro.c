//#include "collision.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "interpolation.h"
#include "geometry.h"
#include "simulation.h"
#include "gyro.h"
#include "quantities_vp.h"
#include "solverpoisson.h"


//centered flux num gyro
void GyroCenteredNumFlux(schnaps_real wL[],schnaps_real wR[],schnaps_real* vnorm,schnaps_real* flux){
  schnaps_real eps =0; //if not equal 0 => decentered flux
  KineticData *kd = &schnaps_kinetic_data;
  schnaps_real E_x =wL[kd->index_ex];
  schnaps_real E_y =wL[kd->index_ey];
  //real dv=kd->dv;
  //printf("test dv = %f \n",dv); 
  for(int i=0;i<kd->index_max_kin+1;i++){
    int j = i % kd->deg_v; // local connectivity put in function
    int nel=i/ kd->deg_v; // element num (TODO : function)  
    schnaps_real wm = (wL[i]+wR[i])/2;
    schnaps_real flux1 = E_y*wm;
    schnaps_real flux2 = -E_x*wm;
    schnaps_real v =- kd->vmax + nel* kd->dv +
      kd->dv * glop(kd->deg_v,j); // gauss_lob_point[j]

    //printf("v=%f\n",v);
    schnaps_real flux3 =v*wm;
  
    flux[i] = vnorm[0]*flux1+vnorm[1]*flux2+vnorm[2]*flux3-eps*(wR[i]-wL[i])/2;
    
  }
  flux[kd->index_phi] =0;
  flux[kd->index_rho] =0;
  flux[kd->index_ex]=0;
  flux[kd->index_ey]=0;
  flux[kd->index_ez]=0;
}

void GyroUpwindNumFlux(schnaps_real wL[],schnaps_real wR[],schnaps_real* vnorm,schnaps_real* flux){
  KineticData *kd = &schnaps_kinetic_data;
  schnaps_real E_x =wL[kd->index_ex];
  schnaps_real E_y =wL[kd->index_ey];
  
  for(int i=0;i<kd->index_max_kin+1;i++){
    schnaps_real vn = E_y * vnorm[0] - E_x *vnorm[1];    

    int j=i%kd->deg_v; // local connectivity put in function
    int nel=i/kd->deg_v; // element num (TODO : function)  
    schnaps_real v =-kd->vmax+ nel*kd->dv +
      kd->dv* glop(kd->deg_v,j); // gauss_lob_point[j]

    vn += v * vnorm[2];

    schnaps_real vnp = vn > 0 ? vn : 0;
    schnaps_real vnm = vn - vnp;
    //printf("v=%f\n",v);
    flux[i] = vnp * wL[i] + vnm * wR[i];
    //flux[i]=0;
  
  }
  flux[kd->index_phi] =0;
  flux[kd->index_rho] =0;
  flux[kd->index_ex]=0;
  flux[kd->index_ey]=0;
  flux[kd->index_ez]=0;
}

void GyroZeroNumFlux(schnaps_real wL[],schnaps_real wR[],schnaps_real* vnorm,schnaps_real* flux){
  
  KineticData *kd = &schnaps_kinetic_data;
  for(int i=0;i<kd->index_max_kin+1;i++){
    flux[i]=0;
  }
  flux[kd->index_phi] = 0;
  flux[kd->index_rho] =0;
  flux[kd->index_ex] = 0;
  flux[kd->index_ey] = 0;
  flux[kd->index_ez] = 0;
}

//! \brief compute square of velocity L2 error
//! \param[in] x  positions
//! \param[in] t  time
//! \param[in] w  values of f at glops
//! \returns the error
schnaps_real GyroL2VelError(schnaps_real* x,schnaps_real t,schnaps_real *w)
{
  KineticData *kd = &schnaps_kinetic_data;
  schnaps_real wex[kd->index_max];
  schnaps_real err2=0;
  GyroImposedData(x, t,wex);
  // loop on the finite emlements
  for(int iel=0;iel<kd->nb_elem_v;iel++){
    // loop on the local glops
    for(int iloc=0;iloc<kd->deg_v+1;iloc++){
      schnaps_real omega=wglop(kd->deg_v,iloc);
      schnaps_real vi=-kd->vmax+iel*kd->dv+kd->dv*glop(kd->deg_v,iloc);
      int ipg=iloc+iel*kd->deg_v;
      err2+=omega*kd->dv*(w[ipg]-wex[ipg])*(w[ipg]-wex[ipg]);
    }
  }
  return err2;
}

void GyroBoundaryFlux(schnaps_real x[3],schnaps_real t,schnaps_real wL[],schnaps_real* vnorm,
		      schnaps_real* flux)
{
  KineticData *kd = &schnaps_kinetic_data;
  schnaps_real wR[kd->index_max];
  GyroImposedData(x,t,wR);
  GyroUpwindNumFlux(wL,wR,vnorm,flux);
}

void GyroInitData(schnaps_real x[3],schnaps_real w[]){
  schnaps_real t=0;
  GyroImposedData(x,t,w);
}

void GyroImposedData(const schnaps_real x[3], const schnaps_real t, schnaps_real w[])
{
  KineticData *kd = &schnaps_kinetic_data;
  for(int i = 0; i <kd->index_max_kin + 1; i++){
    int j=i%kd->deg_v; // local connectivity put in function
    int nel=i/kd->deg_v; // element num (TODO : function)

    schnaps_real vi = (-kd->vmax+nel*kd->dv + kd->dv* glop(kd->deg_v,j));

    //printf("vi=%f\n",vi);
    //w[i]=cos(x[0]-vi*t);
    schnaps_real pi=4*atan(1.);
    schnaps_real xi = x[0]-t;
    //pour transport 1d
    // w[i] = cos(2*pi*xi);//*exp(-0.5*pow(xi-0.5,2));
    //w[i] = cos(2*pi*xi)*exp(-1/(1-pow(2*xi-1,2)));
    //w[i] = exp(-4*pow(xi-0.5,2));
    schnaps_real yi = x[1]+t;
    //pour transport 2D
    //w[i] = exp(-4*pow(yi-0.5,2))*exp(-4*pow(xi-0.5,2));
    //w[i] = cos(2*pi*yi);//*exp(-1/(1-pow(2*xi-1,2)));
    //w[i] = cos(2*pi*xi)*cos(2*pi*yi);//*exp(-0.5*pow(xi-0.5,2));
    //pour transport 3D
    //real zi=x[2]-vi*t;
    schnaps_real zi=x[2]-vi*t;
    //w[i] = exp(-4*pow(zi-0.5,2));//exp(-4*pow(yi-0.5,2))*exp(-4*pow(xi-0.5,2))*exp(-4*pow(zi-0.5,2));

    //w[i] = exp(-4*pow(zi-1.,2));// *Gyro_ImposedKinetic_Data(x,t,vi);
    //w[i] = cos(2 * pi * zi);
    w[i] = exp(-(vi - t) * (vi - t)/2);
  }
  // exact value of the potential
  // and electric field
  w[kd->index_phi]=0;
  w[kd->index_rho] =0;
  w[kd->index_ex]=0;//1;
  w[kd->index_ey]=0;//1;
  w[kd->index_ez]=1;

}



schnaps_real GyroImposedKineticData(const schnaps_real x[3], const schnaps_real t, schnaps_real v)
{
  schnaps_real f;
  f=exp(-pow((v-1.*t),2)/16.); //velocity transport, Ez=1
  //f=exp(-4*pow(xi-0.5,2))*exp(-pow((v-2.*t),2));
  return f;
}

schnaps_real GyroL2_Kinetic_error(field* f)
{
  schnaps_real error=0;
  schnaps_real error_space=0;
  schnaps_real moy=0; // mean value
  schnaps_real moy_space=0;


  // loop on the glops (for numerical integration)
  for(int ipg=0;ipg<NPG(f->deg, f->raf);ipg++){
    schnaps_real xpgref[3],xphy[3],wpg;
    schnaps_real dtau[3][3],codtau[3][3];//,xpg[3];
    // get the coordinates of the Gauss point
    ref_pg_vol(f->deg, f->raf, ipg, xpgref, &wpg, NULL);
    schnaps_ref2phy(f->physnode, // phys. nodes
		    xpgref,  // xref
		    NULL,-1, // dpsiref,ifa
		    xphy,dtau,  // xphy,dtau
		    codtau,NULL,NULL); // codtau,dpsi,vnds
    schnaps_real det
      = dtau[0][0] * codtau[0][0]
      + dtau[0][1] * codtau[0][1]
      + dtau[0][2] * codtau[0][2]; 
    schnaps_real w[f->model.m];
    for(int iv=0;iv<f->model.m;iv++){
      int imem=f->varindex(f->deg, f->raf, f->model.m,ipg,iv);
      w[iv]=f->wn[imem];
    }
    // get the exact value
    error+=GyroL2VelError(xphy,f->tnow,w)*wpg*det;
  }
  
  //moy=moy+weight*moy_space;

  return sqrt(error);
  //return sqrt(error)/sqrt(moy);
}

//! \brief compute compute the source term of the collision
//! model: electric force + true collisions
void GyroSource(const schnaps_real *x, const schnaps_real t, const schnaps_real *w, schnaps_real *source)
{
  KineticData *kd = &schnaps_kinetic_data;
  schnaps_real Ez=w[kd->index_ez]; // electric field
  schnaps_real Md[kd->mv];
  for(int iv=0;iv<kd->index_max_kin+1;iv++){
    Md[iv]=0;
    source[iv]=0;
  }
  // loop on the finite elements
  for(int iel=0;iel<kd->nb_elem_v;iel++){
    // loop on the local glops
    for(int kloc=0;kloc<kd->deg_v+1;kloc++){
      schnaps_real omega=wglop(kd->deg_v,kloc);
      int kpg=kloc+iel*kd->deg_v;
      Md[kpg]+=omega*kd->dv;
      for(int iloc=0;iloc<kd->deg_v+1;iloc++){
	int ipg=iloc+iel*kd->deg_v;
	source[ipg]+=Ez*omega*w[kpg]*dlag(kd->deg_v,iloc,kloc);
      }
    }
  }

  // upwinding !beta =1 (or 1/2?)
  if (Ez>0){
    source[kd->mv-1]+=-Ez*w[kd->mv-1]; 
    //source[0]+=Ez*w[0];
  }
  else {
    source[0]+=Ez*w[0];
    //source[kd->index_max_kin]+=-Ez*w[kd->index_max_kin];
  }

  for(int iv=0;iv<kd->index_max_kin+1;iv++){
    source[iv]/=Md[iv];
  }
  source[kd->index_phi]=0;
  source[kd->index_rho] =0;
  source[kd->index_ex]=0;
  source[kd->index_ey]=0;
  source[kd->index_ez]=0;
}

void Velocity_distribution_plot(schnaps_real *w)
{
  FILE *ver;
  ver = fopen( "fvel.dat", "w" );

  KineticData *kd = &schnaps_kinetic_data;
  // loop on the finite emlements
  for(int iel=0;iel<kd->nb_elem_v;iel++){
    // loop on the local glops
    for(int iloc=0;iloc<kd->deg_v+1;iloc++){
      schnaps_real omega=wglop(kd->deg_v,iloc);
      schnaps_real dv=kd->dv;
      schnaps_real vi=-kd->vmax+iel*dv+dv*glop(kd->deg_v,iloc);
      int ipg=iloc+iel*kd->deg_v;
      //err2+=omega*kd->dv*(w[ipg]-wex[ipg])*(w[ipg]-wex[ipg]);
      fprintf(ver,"%f %f\n",vi,w[ipg]);
      //printf("dv is  %f \n", dv);
      //printf("vi est %d %d %d %f\n", iel,iloc,ipg, vi);
    }
  }
  fclose(ver);
}

void UpdateGyroPoisson(void *si, schnaps_real *w) {
  KineticData *kd = &schnaps_kinetic_data;
  Simulation *simu = si;
  
  int type_bc = 1;
  
  Computation_charge_density(simu);
  static ContinuousSolver ps;
  static bool is_init = false;
  
  if (!is_init){
    is_init = true;
    int nb_var=1;
    int * listvar= malloc(nb_var * sizeof(int));
    listvar[0]=kd->index_phi;
    InitContinuousSolver(&ps,simu,1,nb_var,listvar);
    
    ps.matrix_assembly=ContinuousOperator_Poisson1D;
    ps.rhs_assembly=RHSPoisson_Continuous;
    ps.bc_assembly=Periodic_BoundaryCondition_Poisson1D;
    ps.postcomputation_assembly=Computation_ElectricField_Poisson;
    
    ps.lsol.solver_type = LU;
    ps.lsol.pc_type=NONE;
  }
  
  SolveContinuous2D(&ps);
  //freeContinuousSolver(&ps);
}


