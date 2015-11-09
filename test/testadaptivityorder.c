#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "test.h"
#include "collision.h"
#include "quantities_vp.h"
#include "solverpoisson.h"
#include "linear_solver.h"


void TestPoisson_ImposedData(const real x[3],const real t,real w[]);
void TestPoisson_InitData(real x[3],real w[]);
void TestPoisson_BoundaryFlux(real x[3],real t,real wL[],real* vnorm,
			      real* flux);

int main(void) 
{
  
  // unit tests
    
  int resu = Test_OrderAdaptivity();
	 
  if (resu) printf("1d poisson test OK !\n");
  else printf("1d poisson test failed !\n");

  return !resu;
}

int Test_OrderAdaptivity(void) 
{
  bool test = true;

#ifdef PARALUTION 
  paralution_begin();
#endif 

  // 2D meshes:
  // test/disque2d.msh
  // test/testdisque2d.msh
  // test/testmacromesh.msh
  // test/unit-cube.msh

  // 3D meshes"
  // test/testdisque.msh

  
  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testcube.msh");
  Detect1DMacroMesh(&mesh);
  bool is1d=mesh.is1d;
  assert(is1d);
  BuildConnectivity(&mesh);

  Model model;


  // num of conservative variables f(vi) for each vi, phi, E, rho, u,
  // p, e (ou T)
  model.m=_INDEX_MAX; 
  model.NumFlux = VlasovP_Lagrangian_NumFlux;
  model.Source = VlasovP_Lagrangian_Source;
  
  model.BoundaryFlux = TestPoisson_BoundaryFlux;
  model.InitData = TestPoisson_InitData;
  model.ImposedData = TestPoisson_ImposedData;
  model.Source = NULL;
 
  int deg[]={4, 0, 0};
  int raf[]={8, 1, 1};
    
  PrintMacroMesh(&mesh);
  //assert(1==2);
  //AffineMapMacroMesh(&(f.macromesh));
 
  CheckMacroMesh(&mesh, deg, raf);
  Simulation simu;

  InitSimulation(&simu, &mesh, deg, raf, &model);
  ContinuousSolver ps;
  
  int nb_var=1;
  int * listvar= malloc(nb_var * sizeof(int));
  listvar[0]=_INDEX_PHI;
  
  InitContinuousSolver(&ps,&simu,1,nb_var,listvar);

  ps.matrix_assembly=ContinuousOperator_Poisson1D;
  ps.rhs_assembly=RHSPoisson_Continuous;
  ps.bc_assembly= ExactDirichletContinuousMatrix;
  ps.postcomputation_assembly=Computation_ElectricField_Poisson;

  ps.lsol.solver_type=GMRES;
  ps.lsol.tol=1.0e-8;
  ps.lsol.pc_type=LO_POISSON;
  ps.lsol.iter_max=10000;
  ps.lsol.restart_gmres=30;
  ps.lsol.is_CG=true;

  SolveContinuous2D(&ps);


  real errl2 = L2error(&simu);

  printf("Erreur L2=%.12e\n",errl2);

  test = test && (errl2 < 2e-2);

  printf("Plot...\n");


  PlotFields(_INDEX_PHI, false, &simu, NULL, "dgvisu.msh");
  PlotFields(_INDEX_EX, false, &simu, NULL, "dgex.msh");

#ifdef PARALUTION 
  paralution_end();
#endif

  FreeMacroMesh(&mesh);




  return test;
}


void TestPoisson_ImposedData(const real x[3], const real t,real w[]){
  for(int i = 0; i < _INDEX_MAX_KIN + 1; i++){
    int j = i%_DEG_V; // local connectivity put in function
    int nel = i / _DEG_V; // element num (TODO : function)

    real vi = (-_VMAX + nel * _DV + _DV * glop(_DEG_V, j));

    w[i] = 1. / _VMAX;
  }
  // exact value of the potential
  // and electric field
  w[_INDEX_PHI] = x[0] * x[0] * x[0] * (1 - x[0]) * (1 - x[0]) * (1 - x[0]);
  w[_INDEX_EX] = 3.0 * x[0] * x[0] * (1 - x[0]) * (1 - x[0]) * (- 1. + 2. * x[0]);
  w[_INDEX_RHO] = 6*x[0]*(1 - x[0])*((1 - x[0])*(-1 + 2.*x[0]) - x[0]*(-1. +2.*x[0])-x[0]*(x[0]-1.0)); //rho init
  w[_INDEX_VELOCITY] = 0; // u init
  w[_INDEX_PRESSURE] = 0; // p init
  w[_INDEX_TEMP] = 0; // e ou T init

};

void TestPoisson_InitData(real x[3],real w[]){
  real t=0;
  TestPoisson_ImposedData(x,t,w);
};


void TestPoisson_BoundaryFlux(real x[3],real t,real wL[],real* vnorm,
				       real* flux){
  real wR[_MV+6];
  TestPoisson_ImposedData(x,t,wR);
  VlasovP_Lagrangian_NumFlux(wL,wR,vnorm,flux);
};


void restriction_Pq_P1(int q, real * Polynome_q, real * Polynome_1){
  real A[q+1][2];
  A[0][0]=1.0;
  A[0][1]=0.0;
  for(int ie = 1; ie<q; ie++){  
    A[ie][0]=0.0;
    A[ie][1]=0.0; 
  }
  A[q][0]=0.0;
  A[q][1]=1.0;

  ////// in general cas least square but in this case

  Polynome_1[0]=Polynome_q[0];
  Polynome_1[1]=Polynome_q[q];
  
}


