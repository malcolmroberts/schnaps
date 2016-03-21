#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "../test/test.h"
#include "collision.h"
#include "quantities_vp.h"
#include "diagnostics_vp.h"
#include "solverpoisson.h"


int TestLattice_isothermal(void);

void Relaxation(void* simu);

void Moments(void * simu);


int main(void) {
  
  // unit tests
    
  int resu=TestLattice_isothermal();
	 
  if (resu) printf("landau test OK !\n");
  else printf("landau test failed !\n");

  return !resu;
} 


int TestLattice_isothermal(void) {
  
  int test=0;
  int vec=1;
  schnaps_real k=0.5;
  schnaps_real pi=4.0*atan(1.0);

  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testcube.msh");
  schnaps_real A[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0,1}};
  schnaps_real x0[3] = {0, 0, 0};
  AffineMapMacroMesh(&mesh,A,x0);

  bool is2d=mesh.is2d;
  assert(is2d);

  // mesh preparation
  mesh.period[0]=1.0;
  mesh.period[1]=1.0;
  BuildConnectivity(&mesh);

  
  Model model;
  LatticeData * ld=&schnaps_Lattice_data;

  InitLatticeData(2,9,0,1);
  
  model.m=ld->index_max; // num of conservative variables f(vi) for each vi, phi, E, rho, u, p, e (ou T)
  
  model.NumFlux=Lattice_NumFlux;

  model.InitData = DoubleShear_InitData;
  model.ImposedData = NULL;
  model.BoundaryFlux = NULL;
  model.Source = NULL;

  
  int deg[]={4, 0, 0};
  int raf[]={32, 1, 1};

  CheckMacroMesh(&mesh, deg, raf);
  Simulation simu;
  EmptySimulation(&simu);

  InitSimulation(&simu, &mesh, deg, raf, &model);
  simu.vmax = 2*ld->c; 
  simu.cfl=0.5;
  simu.nb_diags = 1;
  simu.pre_dtfields = Relaxation;
  simu.post_dtfiels = NULL;
  simu.update_after_rk = NULL;
 
  schnaps_real tmax = 0.1;

  RK1(&simu, tmax);
 

  PlotFields(index_rho, false, &simu, "sol","dgvisu_rho.msh");
  PlotFields(index_ux, false, &simu, "sol","dgvisu_ux.msh");
  PlotFields(index_uy, false, &simu, "sol","dgvisu_y.msh");

  test= 1;

  return test; 

}


void DoubleShear_InitData(const schnaps_real x[3], const schnaps_real t, schnaps_real w[])
{
  KineticData * kd=&schnaps_kinetic_data;
  schnaps_real my_pi= 4.0*atan(1.0);
  schnaps_real delta = 0.05, kappa=80.0;

  
  for(int i=0;i<kd->index_max_q;i++){
    w[i]=0;
  }

  w[kd->index_rho]=1.0;
  if(x[1]<0.5){
    w[kd->index_ux]=tanh(kappa*(x[1]-0.25));
  }
  else{
    w[kd->index_ux]=tanh(kappa*(0.75-x[1]));
  }    

  w[kd->index_rho] = 1.0
  w[kd->index_uy]= delta * sin(2.0 * my_pi*(x[0]+0.25));
  w[kd->index_uz]=0.0;
  w[kd->index_temp]=1./3.0; // p init
  w[kd->index_p]=1.0; // e ou T init
};


void Equilibrium_VelocityPerturbation_BoundaryFlux(schnaps_real x[3],schnaps_real t,schnaps_real wL[],schnaps_real* vnorm,
				       schnaps_real* flux){
  
  LatticeData * kd=&schnaps_lattice_data;

};


void Relaxation(void* simu){
  LatticeData * ld=&schnaps_lattice_data;
  double w_eq[ld->q];
  
  Compute_distribution_eq(simu,w_eq);
  Compute_relaxation(simu,w_eq);

}

void Moments(void* simu){
  
  Compute_moments(simu);

}
