#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include "../test/test.h"
#include "collision.h"
#include "quantities_vp.h"
#include "solverpoisson.h"
#include "linear_solver.h"
#include "advanced_linear_solver.h"


void TestPoisson_ImposedData(const real x[3],const real t,real w[]);
void TestPoisson_InitData(real x[3],real w[]);
void TestPoisson_BoundaryFlux(real x[3],real t,real wL[],real* vnorm,
			      real* flux);

void Create_Polynome(ContinuousSolver *cs_HighOrder,int order, real * VecOut);

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
  
  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testcube.msh");
  Detect1DMacroMesh(&mesh);
  bool is1d=mesh.is1d;
  assert(is1d);
  BuildConnectivity(&mesh);

  Model model;
  model.m=_INDEX_MAX; 
  model.NumFlux = VlasovP_Lagrangian_NumFlux;
  model.Source = VlasovP_Lagrangian_Source;
  
  model.BoundaryFlux = TestPoisson_BoundaryFlux;
  model.InitData = TestPoisson_InitData;
  model.ImposedData = TestPoisson_ImposedData;
  model.Source = NULL;
 
  int deg[]={1, 0, 0};
  int raf[]={16, 1, 1};
    
  PrintMacroMesh(&mesh);
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
  
  ////////////////////////////////////////////////////
  /*int deg_ho[]={1, 0, 0};
  int raf_ho[]={40, 1, 1};
  Simulation simu3;
  EmptySimulation(&simu3);
  InitSimulation(&simu3, &mesh, deg_ho, raf_ho, &model);
  ContinuousSolver ps_ho;
  int nb_var3=1;
  int * listvar3= calloc(nb_var3,sizeof(int));
  listvar3[0]=0;
  InitContinuousSolver(&ps_ho,&simu3,1,nb_var3,listvar3);
  
  int deg_lo[]={1, 0, 0};
  int raf_lo[]={raf_ho[0], 1, 1};
  Simulation simu2;
  EmptySimulation(&simu2);
  InitSimulation(&simu2, &mesh, deg_lo, raf_lo, &model);
  ContinuousSolver ps_lo;
  int nb_var2=1;
  int * listvar22= calloc(nb_var2,sizeof(int));
  listvar22[0]=0;
  InitContinuousSolver(&ps_lo,&simu2,1,nb_var2,listvar22);
  
  real * Pol=NULL;
 
  real * Pol_ref=NULL;
  
  real * Pol_reduced=NULL;
  
  Pol=calloc(ps_ho.nb_fe_nodes,sizeof(real));
  Pol_ref=calloc(ps_ho.nb_fe_nodes,sizeof(real));
  
  Pol_reduced=calloc(ps_lo.nb_fe_nodes,sizeof(real));
  
  Create_Polynome(&ps_ho,deg_ho[0],Pol);
  Create_Polynome(&ps_ho,deg_ho[0],Pol_ref);
     
   Restriction1D_Pq_P1(&ps_lo,deg_ho[0],Pol,Pol_reduced);
   Interpolation1D_P1_Pq(&ps_lo,deg_ho[0],Pol_reduced,Pol);
   
   real error=0;
  for(int i=0;i<ps_ho.nb_fe_nodes;i++){
    error=error+(Pol[i]-Pol_ref[i])*(Pol[i]-Pol_ref[i]);
  }
  printf("error %.8e \n",sqrt(error));

  free(Pol);
  free(Pol_ref);
  free(Pol_reduced);*/
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





void Create_Polynome(ContinuousSolver *cs_HighOrder,int order, real * VecOut){
  field* f0 = &cs_HighOrder->simu->fd[0];
  real Polynome[order+1];
  real temp_poly[cs_HighOrder->nb_fe_nodes][order+1];

  Polynome[0]=4.0;
  Polynome[1]=-1.0;
  Polynome[2]=1.0;
  for(int i = 3; i< order+1; i++){
    Polynome[i]=0.0;
  }
  
  for(int ie = 0; ie< cs_HighOrder->nbel; ie++){
    
    int isubcell = ie % (f0->raf[0] * f0->raf[1] * f0->raf[2]);
    int iemacro = ie / (f0->raf[0] * f0->raf[1] * f0->raf[2]);
      
    real wpg;
    real xref[3], xref_begin[3], xref_end[3];
    real dtau[3][3],codtau[3][3];   
      
    for(int ipg = 0;ipg < cs_HighOrder->nnodes; ipg++){
	
      int ipgmacro= ipg + isubcell * cs_HighOrder->nnodes;
      
      ref_pg_vol(f0->deg,f0->raf,ipgmacro,xref,&wpg,NULL);
      real dtau[3][3],codtau[3][3];
      
      Ref2Phy(f0[iemacro].physnode,
	      xref,NULL,0,NULL,
	      dtau,codtau,NULL,NULL);
      
      int ino_dg = ipg + ie * cs_HighOrder->nnodes;
      int ino_fe = cs_HighOrder->dg_to_fe_index[ino_dg];
      for(int i = 0; i< order+1; i++){
	real deg=i;
	temp_poly[ino_fe][i]=pow(xref[0],deg)*Polynome[i];
      }
    }  
  }
  for(int iee = 0; iee<cs_HighOrder->nb_fe_nodes ; iee++){  
    for(int i = 0; i< order+1; i++){  
      VecOut[iee]+=temp_poly[iee][i];
    }
  }
}
