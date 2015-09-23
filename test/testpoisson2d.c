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
    
  int resu = TestPoisson2d();
	 
  if (resu) printf("2d poisson test OK !\n");
  else printf("2d poisson test failed !\n");

  return !resu;
}

int TestPoisson2d(void) 
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
  ReadMacroMesh(&mesh,"../test/testdisque2d.msh");
  Detect2DMacroMesh(&mesh);
  bool is2d=mesh.is2d;
  assert(is2d);
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
 
  int deg[]={2, 2, 0};
  int raf[]={1, 1, 1};
   
 
  CheckMacroMesh(&mesh, deg, raf);
  Simulation simu;

  InitSimulation(&simu, &mesh, deg, raf, &model);

  ContinuousSolver ps;
  
  int nb_var=1;
  int * listvar= malloc(nb_var * sizeof(int));
  listvar[0]=_INDEX_PHI;
  
  InitContinuousSolver(&ps,&simu,1,nb_var,listvar);

  ps.diff_op=calloc(ps.nb_phy_vars*ps.nb_phy_vars,sizeof(SDO));
  real D[4][4] = {{0,0,0,0},
                 {0,1,0,0},
                 {0,0,1,0},
                 {0,0,0,0}};
  for (int i=0;i<4;i++){
    for (int j=0;j<4;j++){
      ps.diff_op[0].DO[i][j]=D[i][j];
    }
  }

  ps.matrix_assembly=GenericOperator_Continuous;
  ps.rhs_assembly=RHSPoisson_Continuous;
  ps.bc_assembly= ExactDirichletContinuousMatrix;
  ps.postcomputation_assembly=Computation_ElectricField_Poisson;

#ifdef PARALUTION
  ps.lsol.solver_type = PAR_LU;
  ps.lsol.pc_type=NONE;
#else
  ps.lsol.solver_type = LU;
  ps.lsol.pc_type=NONE;
#endif

  SolveContinuous2D(&ps);

  real errl2 = L2error(&simu);

  printf("Erreur L2=%f\n",errl2);

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

void TestPoisson_ImposedData(const real x[3], const real t, real w[])
{
  for(int i = 0; i < _INDEX_MAX; i++){
    w[i] = 0;
  }
  // exact value of the potential
  // and electric field
  w[_INDEX_PHI] = (x[0] * x[0] + x[1] * x[1])/4;
  w[_INDEX_EX] =  -x[0]/2;
  w[_INDEX_RHO] = -1; //rho init
  /* w[_INDEX_PHI] = x[0] ; */
  /* w[_INDEX_EX] =  -1; */
  /* w[_INDEX_RHO] = 0; //rho init */
}

void TestPoisson_InitData(real x[3], real w[])
{
  real t = 0;
  TestPoisson_ImposedData(x, t, w);
}

void TestPoisson_BoundaryFlux(real x[3], real t, real wL[], real *vnorm, 
			      real *flux)
{
  real wR[_INDEX_MAX];
  TestPoisson_ImposedData(x, t, wR);
  VlasovP_Lagrangian_NumFlux(wL, wR, vnorm, flux);
}


