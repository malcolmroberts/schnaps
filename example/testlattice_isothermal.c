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
void Vorticity_2D_computation(Simulation *simu, int ifield_ux, int ifield_uy);

void TestDerivative(Simulation *simu,int nbfield);


int main(void) {
  
  // unit tests
    
  int resu=TestLattice_isothermal();
	 
  if (resu) printf("lattice test OK !\n");
  else printf("lattice test failed !\n");

  return !resu;
} 
//
//

int TestLattice_isothermal(void) {
  
  int test=0;
  int vec=1;
  schnaps_real k=0.5;
  schnaps_real pi=4.0*atan(1.0);

  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testcube.msh");
  Detect2DMacroMesh(&mesh);
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
  ld->feq=&feq_isothermal_D2Q9;
  
  model.m=ld->index_max; // num of conservative variables f(vi) for each vi, phi, E, rho, u, p, e (ou T)
  
  model.NumFlux=Lattice_NumFlux;

  model.InitData = DoubleShear_InitData;
  model.ImposedData = NULL;
  model.BoundaryFlux = NULL;
  model.Source = NULL;

  
  int deg[]={4, 4, 0};
  int raf[]={16, 16, 1};

  CheckMacroMesh(&mesh, deg, raf);
  Simulation simu;
  EmptySimulation(&simu);

  InitSimulation(&simu, &mesh, deg, raf, &model);
  simu.vmax = 2*ld->c; 
  simu.cfl=1.0;
  simu.nb_diags = 1;
  simu.pre_dtfields = Relaxation;
  simu.post_dtfields = Moments;
  simu.update_after_rk = NULL;
 
  schnaps_real tmax = 0.00001;

  RK1(&simu, tmax);

  PlotFields(ld->index_rho, false, &simu, "rho","dgvisu_rho.msh");
  PlotFields(ld->index_ux, false, &simu, "ux","dgvisu_ux.msh");
  PlotFields(ld->index_uy, false, &simu, "uy","dgvisu_uy.msh");
  test= 1;
  return test; 
}


void DoubleShear_InitData(schnaps_real x[3],schnaps_real w[])
{
  LatticeData * ld=&schnaps_lattice_data;
  schnaps_real my_pi= 4.0*atan(1.0);
  schnaps_real delta = 0.05, kappa=20.0;
  //
  ld->tau=1.0/30000.0;
  //
  schnaps_real rho = 1.0;
  schnaps_real ux =0.0;
  schnaps_real uy =0.0;
  schnaps_real uz =0.0;
  schnaps_real temp =ld->c * ld-> c;
  schnaps_real p =1.0;
  //
  schnaps_real uref=0.1;
  if(x[1]<0.5){
    ux=uref * tanh(kappa*(x[1]-0.25));
  }
  else{
    ux=uref * tanh(kappa*(0.75-x[1]));
  }    
  uy = uref * delta * sin(2.0 * my_pi*(x[0]+0.25));
  //
  //
  w[ld->index_rho]=rho;
  w[ld->index_ux]  = ux;
  w[ld->index_uy]  = uy;
  w[ld->index_uz]  = uz;
  w[ld->index_temp] = temp;
  w[ld->index_p] = p;
  for(int i=0;i<ld->index_max_q+1;i++){
    w[i]= ld->feq(i,ld,rho,ux,uy,uz,temp,p);
  }
};

void Equilibrium_VelocityPerturbation_BoundaryFlux(schnaps_real x[3],schnaps_real t,schnaps_real wL[],schnaps_real* vnorm,
				       schnaps_real* flux){
  
  LatticeData * ld=&schnaps_lattice_data;

};


void Relaxation(void* s){
  Simulation * simu =s;
  LatticeData * ld=&schnaps_lattice_data;
  schnaps_real w_eq[simu->wsize];
  
  Compute_distribution_eq(simu,w_eq);
  Compute_relaxation(simu,w_eq);

}

void Moments(void* s){
  Simulation * simu =s;
  
  Compute_moments(simu);
}

//
void TestDerivative(Simulation *simu,int nbfield){
    field *f = simu->fd;
    int nb_dof=simu->wsize/f->model.m;
    int wdsize= 3 * nb_dof;
    schnaps_real *wd;
    printf(" WD size %i", wdsize);
    wd= calloc(wdsize,sizeof(schnaps_real));
    for (int i=0; i< wdsize; i++){
      wd[i]= 42.0;
    }
    Compute_derivative(simu,wd,nbfield);
    //
    //
    for (int i=0; i< nb_dof; i++){
      int i1= i;
      int i2= i+ nb_dof;
      int i3= i+ 2* nb_dof;
      printf("i %i \t gradient : %f \t %f \t %f\n", i, wd[i1],wd[i2],wd[i3]);
    }
    //
    if (wd != NULL){
      free(wd);
    }
}
//
void Vorticity_2D_computation(Simulation *simu, int ifield_ux, int ifield_uy){
    field *f = simu->fd;
    int nb_dof=simu->wsize/f->model.m;
    int wdsize= 3 * nb_dof;
    schnaps_real *temp_gradient;
    schnaps_real *vort;
    vort= calloc(nb_dof,sizeof(schnaps_real)); 
    temp_gradient= calloc(wdsize,sizeof(schnaps_real));
    //
    Compute_derivative(simu,temp_gradient,ifield_uy);
    //
    for (int i=0; i< nb_dof; i++){
      vort[i]= temp_gradient[i];
    };
    Compute_derivative(simu,temp_gradient,ifield_ux);
    for (int i=0; i< nb_dof; i++){
      vort[i]= vort[i]-temp_gradient[i+nb_dof];
    }; 
    //
/*    for (int i=0; i< nb_dof; i++){*/
/*      printf("i %i \t vort2d %f\n",i, vort[i]);*/
/*    };*/
    //
    if (temp_gradient != NULL){
      free(temp_gradient);
    }
    if (vort != NULL){
      free(vort);
    }
}
//
