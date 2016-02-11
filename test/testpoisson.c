#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "test.h"
#include "collision.h"
#include "quantities_vp.h"
#include "solverpoisson.h"
#include "linear_solver.h"

int TestPoisson(void) ;

void TestPoisson_ImposedData(const schnaps_real x[3],const schnaps_real t,schnaps_real w[]);
void TestPoisson_InitData(schnaps_real x[3],schnaps_real w[]);
void TestPoisson_BoundaryFlux(schnaps_real x[3],schnaps_real t,schnaps_real wL[],schnaps_real* vnorm,
			      schnaps_real* flux);

int main(void) 
{
  
  // unit tests
    
  int resu = TestPoisson();
	 
  if (resu) printf("1d poisson test OK !\n");
  else printf("1d poisson test failed !\n");

  return !resu;
}

int TestPoisson(void) 
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
 
  int deg[]={3, 0, 0};
  int raf[]={30, 1, 1};
    
  PrintMacroMesh(&mesh);
  //assert(1==2);
  //AffineMapMacroMesh(&(f.macromesh));
 
  CheckMacroMesh(&mesh, deg, raf);
  Simulation simu;
  EmptySimulation(&simu);
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

#ifdef PARALUTION
  ps.lsol.solver_type =LU;
  ps.lsol.pc_type=NONE;
#else
  //ps.lsol.solver_type = GMRES;
  ps.lsol.solver_type = LU;
  ps.lsol.pc_type=NONE;
#endif 
  

  SolveContinuous2D(&ps);


  schnaps_real errl2 = L2error(&simu);

  printf("Erreur L2=%.12e\n",errl2);

  test = test && (errl2 < 2e-2);

  printf("Plot...\n");


  PlotFields(_INDEX_PHI, false, &simu, NULL, "dgvisu.msh");
  PlotFields(_INDEX_EX, false, &simu, NULL, "dgex.msh");

#ifdef PARALUTION 
  paralution_end();
#endif
  
  freeContinuousSolver(&ps);
  FreeMacroMesh(&mesh);




  return test;
}


void TestPoisson_ImposedData(const schnaps_real x[3], const schnaps_real t,schnaps_real w[]){
  for(int i = 0; i < _INDEX_MAX_KIN + 1; i++){
    int j = i%_DEG_V; // local connectivity put in function
    int nel = i / _DEG_V; // element num (TODO : function)

    schnaps_real vi = (-_VMAX + nel * _DV + _DV * glop(_DEG_V, j));

    w[i] = 1. / _VMAX;
  }
  // exact value of the potential
  // and electric field
  w[_INDEX_PHI] = x[0] * (1 - x[0]);
  w[_INDEX_EX] = - 1. + 2. * x[0];
  w[_INDEX_RHO] = 2.; //rho init
  w[_INDEX_VELOCITY] = 0; // u init
  w[_INDEX_PRESSURE] = 0; // p init
  w[_INDEX_TEMP] = 0; // e ou T init

};

void TestPoisson_InitData(schnaps_real x[3],schnaps_real w[]){

  schnaps_real t=0;
  TestPoisson_ImposedData(x,t,w);

};


void TestPoisson_BoundaryFlux(schnaps_real x[3],schnaps_real t,schnaps_real wL[],schnaps_real* vnorm,
				       schnaps_real* flux){
  schnaps_real wR[_MV+6];
  TestPoisson_ImposedData(x,t,wR);
  VlasovP_Lagrangian_NumFlux(wL,wR,vnorm,flux);
};


