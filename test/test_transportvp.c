#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "test.h"
#include "collision.h"
#include "quantities_vp.h"
#include "diagnostics_vp.h"
#include "solverpoisson.h"

void Test_TransportVP_ImposedData(const real *x, const real t, real *w);
void Test_TransportVP_InitData(real *x, real *w);
real TransportVP_ImposedKinetic_Data(const real *x, const real t, real v);
void Test_TransportVP_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
				   real *flux);

void UpdateVlasovPoisson(void *field, real *w);
void PlotVlasovPoisson(void *vf, real *w);

int main(void) {
  
  // unit tests
    
  int resu = Test_TransportVP();
	 
  if (resu) printf("poisson test OK !\n");
  else printf("poisson test failed !\n");

  return !resu;
} 

int Test_TransportVP(void) {

  bool test = true;

#ifdef PARALUTION 
  paralution_begin();
#endif 

  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testcube.msh");
  Detect1DMacroMesh(&mesh);
  
  bool is1d = mesh.is1d;
  assert(is1d);
  
  real A[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0,1}};
  real x0[3] = {0, 0, 0};
  AffineMapMacroMesh(&mesh,A,x0);

  BuildConnectivity(&mesh);

  Model model;
  
  model.m=_INDEX_MAX; // num of conservative variables f(vi) for each vi, phi, E, rho, u, p, e (ou T)
  model.NumFlux=VlasovP_Lagrangian_NumFlux;
  model.InitData = Test_TransportVP_InitData;
  model.ImposedData = Test_TransportVP_ImposedData;
  model.BoundaryFlux = Test_TransportVP_BoundaryFlux;
  model.Source = VlasovP_Lagrangian_Source;


  int deg[]={2, 0, 0};
  int raf[]={16, 1, 1};

  CheckMacroMesh(&mesh, deg, raf);
  Simulation simu;

  InitSimulation(&simu, &mesh, deg, raf, &model);
  simu.vmax = _VMAX; // maximal wave speed
  simu.cfl=0.2;
  simu.nb_diags = 3;
  simu.pre_dtfields = UpdateVlasovPoisson;
  simu.post_dtfields=NULL;
  simu.update_after_rk = PlotVlasovPoisson;


  real tmax = 0.03;
  RK4(&simu, tmax);

   // save the results and the error
  int iel = 2 * _NB_ELEM_V / 3;
  int iloc = _DEG_V;
  printf("Trace vi=%f\n", -_VMAX + iel * _DV + _DV * glop(_DEG_V, iloc));
  PlotFields(iloc + iel * _DEG_V, false, &simu, "sol","dgvisu_kin.msh");
  PlotFields(iloc + iel * _DEG_V, true, &simu, "error","dgerror_kin.msh");
  Plot_Energies(&simu, simu.dt);

  real dd_Kinetic = L2_Kinetic_error(&simu);
  
  printf("erreur kinetic L2=%lf\n", dd_Kinetic);
  test= test && (dd_Kinetic < 1e-2);

  FreeMacroMesh(&mesh);

  return test;
}

void Test_TransportVP_ImposedData(const real *x, const real t, real *w) {

  for(int i = 0; i <_INDEX_MAX_KIN + 1; i++) {
    int j = i % _DEG_V; // local connectivity put in function
    int nel = i / _DEG_V; // element num (TODO : function)

    real vi = (-_VMAX + nel * _DV + _DV * glop(_DEG_V, j));
 
    w[i] = TransportVP_ImposedKinetic_Data(x, t, vi);
  }
  // exact value of the potential and electric field
  w[_INDEX_PHI] = -x[0];
  w[_INDEX_EX] = 1;
  w[_INDEX_RHO] = 0.; //rho init
  w[_INDEX_VELOCITY] = 0; // u init
  w[_INDEX_PRESSURE] = 0; // p init
  w[_INDEX_TEMP] = 0; // e ou T init
}

void Test_TransportVP_InitData(real *x, real *w) {
  real t = 0;
  Test_TransportVP_ImposedData(x, t, w);
}

real TransportVP_ImposedKinetic_Data(const real *x, const real t, real v) {
  real f;
  real pi = 4.0 * atan(1.0);
  real xnew = 0, vnew = 0;
  f = exp(-(v - t) * (v - t)) *
    exp(-36 * ((x[0] - v * t + 0.5 * t * t) - 0.5)
	* ((x[0] - v * t + 0.5 * t * t) - 0.5));
  return f;
}

void Test_TransportVP_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
				       real *flux) {
  real wR[_INDEX_MAX];
  Test_TransportVP_ImposedData(x , t, wR);
  VlasovP_Lagrangian_NumFlux(wL, wR, vnorm, flux);
}

void UpdateVlasovPoisson(void *si, real *w) {
  Simulation *simu = si;
  
  int type_bc = 1;
  real bc_l = 1;
  real bc_r = 0;
    
  // Computation_charge_density(simu,simu->w);
  
  ContinuousSolver ps;
  
  int nb_var=1;
  int * listvar= malloc(nb_var * sizeof(int));
  listvar[0]=_INDEX_PHI;
  
  InitContinuousSolver(&ps,simu,1,nb_var,listvar);

  real D[4][4] = {{0,0,0,0},
                 {0,1,0,0},
                 {0,0,1,0},
                 {0,0,0,1}};
  for (int i=0;i<4;i++){
    for (int j=0;j<4;j++){
      ps.diff_op[i][j]=D[i][j];
    }
  }

  ps.matrix_assembly=GenericOperatorScalar_Continuous;
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

  //freeContinuousSolver(&ps);
  
}

void PlotVlasovPoisson(void *si, real *w) {
  real k_energy = 0, e_energy = 0, t_energy = 0;
  
  Simulation *simu = si;
  
  Energies(simu, w, k_energy, e_energy, t_energy,1);
  si = simu;
}
