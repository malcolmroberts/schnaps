#include "implicit.h"
#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "test.h"
#include "collision.h"
#include "solvercontinuous.h"
#include "waterwave2d.h"
#include "continuouspc.h"



void TestExtraction_ImposedData(const real *x, const real t, real *w);
void TestExtraction_InitData(real *x, real *w);
void TestExtraction_Source(const real *xy, const real t, const real *w, real *S);

int main(void) {
  
  // unit tests
    
  int resu = Test_Extraction();
	 
  if (resu) printf("wave periodic  test OK !\n");
  else printf("wave periodic test failed !\n");

  return !resu;
} 

int Test_Extraction(void) {

  bool test = true;
  real dd;


#ifdef PARALUTION 
  paralution_begin();
#endif 

  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testcube.msh");
  Detect2DMacroMesh(&mesh);
  
  real A[3][3] = {{_LENGTH_DOMAIN, 0, 0}, {0, _LENGTH_DOMAIN, 0}, {0, 0,1}};
  real x0[3] = {0, 0, 0};
  AffineMapMacroMesh(&mesh,A,x0);

  BuildConnectivity(&mesh);

  Model model;

  model.m = 3; 
  model.NumFlux=Wave_Upwind_NumFlux;
  model.InitData = TestExtraction_InitData;
  model.ImposedData = TestExtraction_ImposedData;
  model.BoundaryFlux = Wave_Upwind_BoundaryFlux;
  model.Source =NULL;// TestExtraction_Source;

  int deg[]={4, 4, 0};
  int raf[]={2, 2, 1};


  assert(mesh.is2d);

  CheckMacroMesh(&mesh, deg, raf);
  Simulation simu;

  InitSimulation(&simu, &mesh, deg, raf, &model);


  LinearSolver solver_implicit;
  LinearSolver solver_explicit;  

  real theta=0.5;
  simu.theta=theta;
  simu.dt=1;
  simu.vmax=_SPEED_WAVE;
  real tmax=1;
  
  int itermax=tmax/simu.dt+1;
  simu.itermax_rk=itermax;
  InitImplicitLinearSolver(&simu, &solver_implicit);
  InitImplicitLinearSolver(&simu, &solver_explicit);
  real *res = calloc(simu.wsize, sizeof(real));

  simu.tnow=0;
  for(int ie=0; ie < simu.macromesh.nbelems; ++ie){
    simu.fd[ie].tnow=simu.tnow;
  } 

  for(int tstep=0;tstep<simu.itermax_rk;tstep++){
  

    if(tstep==0){ 
      solver_implicit.mat_is_assembly=false;
      solver_explicit.mat_is_assembly=false;
    } 
    else 
      { 
      solver_implicit.mat_is_assembly=true;
      solver_explicit.mat_is_assembly=true;
    } 

    solver_implicit.rhs_is_assembly=false;
    solver_explicit.rhs_is_assembly=false;

    
    AssemblyImplicitLinearSolver(&simu, &solver_explicit,-(1.0-theta),simu.dt);
    simu.tnow=simu.tnow+simu.dt;
    for(int ie=0; ie < simu.macromesh.nbelems; ++ie){
      simu.fd[ie].tnow=simu.tnow;
    } 
    AssemblyImplicitLinearSolver(&simu, &solver_implicit,theta,simu.dt);
  

    MatVect(&solver_explicit, simu.w, res);

    for(int i=0;i<solver_implicit.neq;i++){
      solver_implicit.rhs[i]=-solver_explicit.rhs[i]+solver_implicit.rhs[i]+res[i];
    }
    physicPC_wave(&simu,simu.fd[0].wn,solver_implicit.rhs);
    //ContinuousSolver waveSolver;
    //int nb_var = 3;
    //int * listvarGlobal = malloc(nb_var * sizeof(int));
    //listvarGlobal[0]=0;
    //listvarGlobal[1]=1;
    //listvarGlobal[2]=2;
    //InitContinuousSolver(&waveSolver,&simu,1,nb_var,listvarGlobal);
    //VectorDgToCg(&waveSolver, simu.w, waveSolver.lsol.sol);
    //simu.dt=0;
    //physicPC_wave(&simu,simu.fd[0].wn,solver_implicit.rhs);
    //VectorCgToDg(&waveSolver, simu.w);
    //simu.dt=1;
    //SolveLinearSolver(&solver_implicit);

    //for(int i=0;i<solver_implicit.neq;i++){
    //  simu.w[i]=solver_implicit.sol[i];
    //}

  }
  
  
  // Apply division by the mass matrix
  //field *f = &simu.fd[0];
  //int m = f->model.m;
  //
  //for(int ipg = 0; ipg < NPG(f->deg, f->raf); ipg++) {
  //  real dtau[3][3], codtau[3][3], xpgref[3], xphy[3], wpg;
  //  ref_pg_vol(f->deg, f->raf, ipg, xpgref, &wpg, NULL);
  //  Ref2Phy(f->physnode, // phys. nodes
  //    xpgref, // xref
  //    NULL, -1, // dpsiref, ifa
  //    xphy, dtau, // xphy, dtau
  //    codtau, NULL, NULL); // codtau, dpsi, vnds
  //  real det = dot_product(dtau[0], codtau[0]);
  //  for(int iv = 0; iv < f->model.m; iv++) {
  //    int imem = f->varindex(f->deg, f->raf, f->model.m, ipg, iv);
  //    simu.fd[0].wn[imem] *= (wpg * det);
  //  }
  //}
  //physicPC_wave(&simu,simu.fd[0].dtwn,simu.fd[0].wn);
  //simu.fd[0].wn = simu.fd[0].dtwn;
 
  dd = L2error(&simu);

  printf("erreur L2=%.12e\n", dd);

  PlotFields(0,false, &simu, "p", "dgvisu_exp.msh");
  PlotFields(1,false, &simu, "u", "dgvisu_exu.msh");
  PlotFields(2,false, &simu, "v", "dgvisu_exv.msh");

#ifdef PARALUTION 
  paralution_end();
#endif 

  return test;
}



void TestExtraction_ImposedData(const real *xy, const real t, real *w) {

  real x=xy[0];
  real y=xy[1];


  //w[0] = x*(1-x)*y*(1-y)+1;
  //w[1] = 2*x*(1-x)*y*(1-y)+2;
  //w[2] = 3*x*(1-x)*y*(1-y)+3;

  w[0] = 10;
  w[1] = y*x+2;//x+2;// 
  w[2] = -y*y*0.5+5;//-y+3;//
  
  //w[0] = 1;
  //w[1] = x+2;// 
  //w[2] = -y+3;//

}

void TestExtraction_InitData(real *x, real *w) {
  real t = 0;
  TestExtraction_ImposedData(x, t, w);
}


void Wave_Upwind_BoundaryFlux(real *x, real t, real *wL, real *vnorm,
				       real *flux) {
  real wR[3];
  TestExtraction_ImposedData(x , t, wR);
  Wave_Upwind_NumFlux(wL, wR, vnorm, flux);
}

void TestExtraction_Source(const real *xy, const real t, const real *w, real *S){
  
  real x=xy[0];
  real y=xy[1];

  S[0] = 2*(1-2*x)*(y*(1-y))+3*(1-2*y)*(x*(1-x));
  S[1] = (1-2*x)*(y*(1-y));
  S[2] = (1-2*y)*(x*(1-x));

  S[0] *= _SPEED_WAVE;
  S[1] *= _SPEED_WAVE;
  S[2] *= _SPEED_WAVE;

}

