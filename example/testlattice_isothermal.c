#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "../test/test.h"
#include "global.h"
#include "lattice.h"

int TestLattice_isothermal(void);

void Relaxation(void* s);
void Moments(void * s);

void DoubleShear_InitData(schnaps_real x[3], schnaps_real w[]);

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

  // mesh preparation
  mesh.period[0]=1.0;
  mesh.period[1]=1.0;
  BuildConnectivity(&mesh);

  
  Model model;
  LatticeData * ld=&schnaps_lattice_data;

  InitLatticeData(ld,2,9,0,1);
  
  model.m=ld->index_max; // num of conservative variables f(vi) for each vi, phi, E, rho, u, p, e (ou T)
  
  model.NumFlux=Lattice_NumFlux;

  model.InitData = DoubleShear_InitData;
  model.ImposedData = NULL;
  model.BoundaryFlux = NULL;
  model.Source = NULL;

  
  int deg[]={4, 4, 0};
  int raf[]={32, 32, 1};

  CheckMacroMesh(&mesh, deg, raf);
  Simulation simu;
  EmptySimulation(&simu);

  InitSimulation(&simu, &mesh, deg, raf, &model);
  simu.vmax = 2*ld->c; 
  simu.cfl=0.5;
  simu.nb_diags = 1;
  simu.pre_dtfields = Relaxation;
  simu.post_dtfields = NULL;
  simu.update_after_rk = NULL;
 
  schnaps_real tmax = 1.0;

  RK2(&simu, tmax);
 

  PlotFields(ld->index_rho, false, &simu, "sol","dgvisu_rho.msh");
  PlotFields(ld->index_ux, false, &simu, "sol","dgvisu_ux.msh");
  PlotFields(ld->index_uy, false, &simu, "sol","dgvisu_y.msh");

  test= 1;

  return test; 

}


void DoubleShear_InitData(schnaps_real x[3],schnaps_real w[])
{
  LatticeData * ld=&schnaps_lattice_data;
  schnaps_real my_pi= 4.0*atan(1.0);
  schnaps_real delta = 0.05, kappa=20.0;
  //
  ld->tau=1.0/1000.0;
  //
  
  for(int i=0;i<ld->index_max_q;i++){
    w[i]=0;
  }

  w[ld->index_rho]=1.0;
  if(x[1]<0.5){
    w[ld->index_ux]=tanh(kappa*(x[1]-0.25));
  }
  else{
    w[ld->index_ux]=tanh(kappa*(0.75-x[1]));
  }    

  w[ld->index_rho] = 1.0;
  w[ld->index_uy]= delta * sin(2.0 * my_pi*(x[0]+0.25));
  w[ld->index_uz]=0.0;
  w[ld->index_temp]=1./3.0; // p init
  w[ld->index_p]=1.0; // e ou T init
};


void Equilibrium_VelocityPerturbation_BoundaryFlux(schnaps_real x[3],schnaps_real t,schnaps_real wL[],schnaps_real* vnorm,
				       schnaps_real* flux){
  
  LatticeData * ld=&schnaps_lattice_data;

};


void Relaxation(void* s){
  Simulation * simu =s;
  LatticeData * ld=&schnaps_lattice_data;
  double w_eq[simu->wsize];
  
  Compute_distribution_eq(simu,w_eq);
  Compute_relaxation(simu,w_eq);

}

void Moments(void* s){
  Simulation * simu =s;
  
  Compute_moments(simu);

}
