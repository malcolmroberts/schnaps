#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "../test/test.h"
#include "collision.h"
#include "quantities_vp.h"
#include "diagnostics_vp.h"
#include "solverpoisson.h"


int TestCollision(void);

void Equilibrium_ImposedData(const schnaps_real x[3], const schnaps_real t,schnaps_real w[]);
void Equilibrium_InitData(schnaps_real x[3],schnaps_real w[]);
void Equilibrium_BoundaryFlux(schnaps_real x[3],schnaps_real t,schnaps_real wL[],schnaps_real* vnorm, schnaps_real* flux);

void Equilibrium2_ImposedData(const schnaps_real x[3], const schnaps_real t,schnaps_real w[]);
void Equilibrium2_InitData(schnaps_real x[3],schnaps_real w[]);
void Equilibrium2_BoundaryFlux(schnaps_real x[3],schnaps_real t,schnaps_real wL[],schnaps_real* vnorm, schnaps_real* flux);
void Equilibrium2_TotalSource(const schnaps_real* x, const schnaps_real t, const schnaps_real* w, schnaps_real* source);

void Collision_VlasovPoisson(void* field, schnaps_real *w);
void PlotVlasovPoisson(void* vf, schnaps_real * w);

int main(void) {
  
  // unit tests
    
  int resu=TestCollision();
	 
  if (resu) printf("landau test OK !\n");
  else printf("landau test failed !\n");

  return !resu;
} 


int TestCollision(void) {
  
  int test=0;
  int vec=1;
  schnaps_real k=0.5;
  schnaps_real pi=4.0*atan(1.0);

  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testcube.msh");
  schnaps_real A[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0,1}};
  schnaps_real x0[3] = {0, 0, 0};
  AffineMapMacroMesh(&mesh,A,x0);

  // try to detect a 1d mesh
  Detect1DMacroMesh(&mesh);
  bool is1d=mesh.is1d;
  assert(is1d);

  // mesh preparation
  //mesh.period[0]=1.0;
  BuildConnectivity(&mesh);

  
  Model model;
  schnaps_real degV=3;
  schnaps_real nbEV=20;
  KineticData * kd=&schnaps_kinetic_data;

  InitKineticData(kd,nbEV,degV);
  
  model.m=kd->index_max; // num of conservative variables f(vi) for each vi, phi, E, rho, u, p, e (ou T)
  /*model.NumFlux=VlasovP_Lagrangian_NumFlux;
  model.InitData =Equilibrium_InitData;
  model.ImposedData = Equilibrium_ImposedData;
  model.BoundaryFlux = Equilibrium_BoundaryFlux;
  model.Source = VlasovP_Lagrangian_Source;*/

  model.NumFlux=VlasovP_Lagrangian_NumFlux;
  model.InitData =Equilibrium2_InitData;
  model.ImposedData = Equilibrium2_ImposedData;
  model.BoundaryFlux = Equilibrium2_BoundaryFlux;
  model.Source = Equilibrium2_TotalSource;

  
  int deg[]={4, 0, 0};
  int raf[]={12, 1, 1};

  CheckMacroMesh(&mesh, deg, raf);
  Simulation simu;
  EmptySimulation(&simu);

  InitSimulation(&simu, &mesh, deg, raf, &model);
  simu.vmax = kd->vmax; // maximal wave speed
  simu.cfl=0.1;
  simu.nb_diags = 4;
  simu.pre_dtfields = Collision_VlasovPoisson;
  simu.post_dtfields=NULL;
  simu.update_after_rk = PlotVlasovPoisson;
 
  schnaps_real tmax = 0.1;
  RK2(&simu, tmax);

    // save the results and the error
  int iel = 2 * kd->nb_elem_v / 3;
  int iloc = kd->deg_v;
  printf("Trace vi=%f\n", -kd->vmax + iel * kd->dv + kd->dv * glop(kd->deg_v, iloc));
  PlotFields(iloc + iel * kd->deg_v, false, &simu, "sol","dgvisu_kin.msh");
  PlotFields(iloc + iel * kd->deg_v, true, &simu, "error","dgerror_kin.msh");
  
  Plot_Energies(&simu, simu.dt);

  schnaps_real dd_Kinetic = L2_Kinetic_error(&simu);
  
  printf("erreur kinetic L2=%.5e \n", dd_Kinetic);
  test= test && (dd_Kinetic < 1e-2);


 

  test= 1;

  return test; 

}

void Equilibrium_ImposedData(const schnaps_real x[3], const schnaps_real t, schnaps_real w[])
{
  //parameters of the case
  KineticData * kd=&schnaps_kinetic_data;
  schnaps_real k=0.5;
  schnaps_real eps = 0.001;
  schnaps_real my_pi= 4.0*atan(1.0);
  
  for(int i=0;i<kd->index_max_kin+1;i++){
    int j=i%kd->deg_v; // local connectivity put in function
    int nel=i/kd->deg_v; // element num (TODO : function)

    schnaps_real vi = (-kd->vmax+nel*kd->dv +
		 kd->dv* glop(kd->deg_v,j));
 
    w[i]=(1.0/sqrt(2.0*my_pi))*exp(-(vi*vi)/2.0);
  }
  // exact value of the potential
  // and electric field
  w[kd->index_phi]=0.0;
  w[kd->index_ex]=0.0;
  w[kd->index_rho]=1.0; //rho init
  w[kd->index_u]=0; // u init
  w[kd->index_P]=0.; // p init
  w[kd->index_T]=1.0; // e ou T init

};

void Equilibrium2_ImposedData(const schnaps_real x[3], const schnaps_real t, schnaps_real w[])
{
  //parameters of the case
  KineticData * kd=&schnaps_kinetic_data;
  schnaps_real my_pi= 4.0*atan(1.0);

  schnaps_real rho=4*my_pi*my_pi*(sin(2*my_pi*x[0]));
  
  for(int i=0;i<kd->index_max_kin+1;i++){
    int j=i%kd->deg_v; // local connectivity put in function
    int nel=i/kd->deg_v; // element num (TODO : function)

    schnaps_real vi = (-kd->vmax+nel*kd->dv +
		 kd->dv* glop(kd->deg_v,j));
 
    w[i]=(rho/sqrt(2.0*my_pi))*exp(-(vi*vi)/2.0);
  }
  // exact value of the potential
  // and electric field
  w[kd->index_phi]=sin(2*my_pi*x[0]);
  w[kd->index_ex]=-2*my_pi*cos(2*my_pi*x[0]);
  w[kd->index_ey]=0.0;
  w[kd->index_ez]=0.0;
  w[kd->index_rho]=rho;  //rho init
  w[kd->index_u]=0; // u init
  w[kd->index_P]=0.; // p init
  w[kd->index_T]=1.0; // e ou T init

};


void Equilibrium2_TotalSource(const schnaps_real* x, const schnaps_real t, const schnaps_real* w, 
			       schnaps_real* source) {
  KineticData * kd=&schnaps_kinetic_data;
  schnaps_real Transport_source[kd->index_max];
  schnaps_real M=0;
  schnaps_real my_pi= 4.0*atan(1.0);

  VlasovP_Lagrangian_Source(x,t,w,Transport_source);

   for(int iv=0;iv< kd->index_max_kin + 1;iv++){
      int j=iv%kd->deg_v; // local connectivity put in function
      int nel=iv/kd->deg_v; // element num (TODO : function)
   
      schnaps_real vn = (-kd->vmax+nel*kd->dv +
			 kd->dv* glop(kd->deg_v,j));
      
      M=(1.0/sqrt(2.0*my_pi))*exp(-(vn*vn)/2.0);	
      source[iv] = Transport_source[iv];
      source[iv] += (vn+(sin(2.0*my_pi*x[0]))*vn)*8.0*pow(my_pi,3.0)*cos(2.0*my_pi*x[0])*M;
   }
  source[kd->index_phi]=0.0;
  source[kd->index_ex]=0.0;
  source[kd->index_ey]=0.0;
  source[kd->index_ez]=0.0;
  source[kd->index_rho]=0.0;  //rho init
  source[kd->index_u]=0; // u init
  source[kd->index_P]=0.; // p init
  source[kd->index_T]=0.0;
};

void Equilibrium_InitData(schnaps_real x[3],schnaps_real w[]){

  schnaps_real t=0;
  Equilibrium_ImposedData(x,t,w);

};

void Equilibrium2_InitData(schnaps_real x[3],schnaps_real w[]){

  schnaps_real t=0;
  Equilibrium2_ImposedData(x,t,w);

};


void Equilibrium_BoundaryFlux(schnaps_real x[3],schnaps_real t,schnaps_real wL[],schnaps_real* vnorm,
				       schnaps_real* flux){
  KineticData * kd=&schnaps_kinetic_data;
  schnaps_real wR[kd->index_max];
  Equilibrium_ImposedData(x,t,wR);
  VlasovP_Lagrangian_NumFlux(wL,wR,vnorm,flux);
};

void Equilibrium2_BoundaryFlux(schnaps_real x[3],schnaps_real t,schnaps_real wL[],schnaps_real* vnorm,
				       schnaps_real* flux){
  KineticData * kd=&schnaps_kinetic_data;
  schnaps_real wR[kd->index_max];
  Equilibrium2_ImposedData(x,t,wR);
  VlasovP_Lagrangian_NumFlux(wL,wR,vnorm,flux);
};



void Collision_VlasovPoisson(void *si, schnaps_real *w) {
  Simulation *simu = si;
  KineticData * kd=&schnaps_kinetic_data;
  
  //Collision_Source(simu);
  
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
    ps.bc_assembly=ExactDirichletContinuousMatrix;
    ps.postcomputation_assembly=Computation_ElectricField_Poisson;
    
    ps.lsol.solver_type = LU;
    ps.lsol.pc_type=NONE;
  }
  
  SolveContinuous2D(&ps);
  //freeContinuousSolver(&ps);
}

void PlotVlasovPoisson(void *si, schnaps_real *w) {
  schnaps_real k_energy = 0, e_energy = 0, t_energy = 0, t_charge=0;
  
  Simulation *simu = si;
  
  Energies(simu, w, k_energy, e_energy, t_energy,1);
  Charge_total(simu,w,t_charge,4);
  si = simu; 
}







