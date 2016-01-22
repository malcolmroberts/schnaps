#include "solverpoisson.h"
#include "solvercontinuous.h"
#include "geometry.h"
#include "quantities_vp.h"
#include "linear_solver.h"


void Computation_ElectricField_Poisson(void * cs){

  ContinuousSolver * ps=cs;
  
  ContinuousToDiscontinuous_Copy(ps);

  // printf("Compute electric field...\n");
  
  for(int ie = 0; ie < ps->simu->macromesh.nbelems; ie++){
    ComputeElectricField(&ps->simu->fd[ie]);
  }
}

void RHSPoisson_Continuous(void * cs){

  ContinuousSolver * ps=cs;
  
  field* f0 = &ps->simu->fd[0];
  
  // right hand side assembly
  for(int ino = 0; ino < ps->nb_fe_dof; ino++){
    ps->lsol.rhs[ino] = 0;
  }

  real surf = 0;

  for(int ie = 0; ie < ps->nbel; ie++){  

    int iemacro = ie / (f0->raf[0] * f0->raf[1] * f0->raf[2]);
    int isubcell = ie % (f0->raf[0] * f0->raf[1] * f0->raf[2]);  
 
    for(int iloc = 0; iloc < ps->nnodes; iloc++){
      real wpg;
      real xref[3];
      //int ipgmacro = ipg + isubcell * nnodes;
      int ilocmacro = iloc + isubcell * ps->nnodes;
      ref_pg_vol(f0->deg,f0->raf,ilocmacro,xref,&wpg,NULL);
      real dtau[3][3],codtau[3][3];
      Ref2Phy(ps->simu->fd[iemacro].physnode,
	      xref,NULL,0,NULL,
	      dtau,codtau,NULL,NULL);
      real det = dot_product(dtau[0], codtau[0]);	
      int ino_dg = iloc + ie * ps->nnodes;
      int ino_fe = ps->dg_to_fe_index[ino_dg];
      int imem = f0->varindex(f0->deg,f0->raf,f0->model.m,
			      ilocmacro,_INDEX_RHO) ;
	//+ iemacro * NPG(f0->deg,f0->raf) * f0->model.m ;
      real rho = ps->simu->fd[iemacro].wn[imem];
      ps->lsol.rhs[ino_fe] += rho * wpg * det ; 
      surf += wpg * det ;  
    }
 
  }
}


void RobinBoundaryConditionAssembly(void * cs){
  ContinuousSolver * ps=cs;
  m=ps->nb_phy_vars;
 
  if (!ps->lsol.mat_is_assembly){  

    for (int i=0; i<ps->simu->macromesh.nboundaryfaces;i++){
      int ifa = ps->simu->macromesh.boundaryface[i];
      int locfaL = ps->simu->macromesh.face2elem[4 * ifa + 1];
      int ie = ps->simu->macromesh.face2elem[4 * ifa ];
      int ieR = ps->simu->macromesh.face2elem[4 * ifa + 2];
      assert(ieR < 0);
      if (ieR<0){
	field *f = &ps->simu->fd[ie];
	
	for(int ipglf = 0;ipglf < NPGF(f->deg,f->raf,locfaL); ipglf++){
	  real xpgref[3], xpgref_in[3], wpg;


	  // Get the coordinates of the Gauss point and coordinates of a
	  // point slightly inside the opposite element in xref_in
	  int ipg = ref_pg_face(f->deg, f->raf, locfaL, ipglf, xpgref, &wpg, xpgref_in);
	  int ino_dg = ipg + ie * ps->npgmacrocell;
	  int ino_fe = ps->dg_to_fe_index[ino_dg];
	  // Normal vector at gauss point ipgL
	  real vnds[3], xpg[3];
	  {
	    real dtau[3][3], codtau[3][3];
	    Ref2Phy(f->physnode,
		    xpgref,
		    NULL, locfaL, // dpsiref, ifa
		    xpg, dtau,
		    codtau, NULL, vnds); // codtau, dpsi, vnds
	  }
	  // the boundary flux is an affine function
	  real flux0[m], wL[m];
	  for(int iv = 0; iv < m; iv++) {
	    wL[iv] = 0;
	  }
	  
	  ps->bc_flux(ps,xpg,wL,vnds,flux0);

	  for(int iv1 = 0; iv1 < m; iv1++) {
	    //int imem1 = fL->varindex(fL->deg, fL->raf,fL->model.m, ipgL, iv1) + offsetL;
	    int ipot_fe1 = ino_fe*ps->nb_phy_vars + iv1;
	    for(int iv = 0; iv < m; iv++) {
	      wL[iv] = (iv == iv1);
	    }
	    
	    real flux[m];
	    ps->bc_flux(ps,xpg,wL,vnds,flux);
	    
	    for(int iv2 = 0; iv2 < m; iv2++) {
	      // The basis functions is also the gauss point index
	      //int imem2 = fL->varindex(fL->deg, fL->raf,fL->model.m, ipgL, iv2) + offsetL;
	      int ipot_fe2 = ino_fe*ps->nb_phy_vars + iv2;
	      real val =  (flux[iv2]-flux0[iv2]) * wpg;
	      AddLinearSolver(&ps->lsol, ipot_fe2, ipot_fe1, val);
	    }
	  }	    
	}
      }
    }
 
    ps->lsol.mat_is_assembly=true;
  }
  
  //for(int var =0; var < ps->nb_phy_vars; var++){ 
  //for(int ie = 0; ie < ps->simu->macromesh.nbelems; ie++){  
  for (int i=0; i<ps->simu->macromesh.nboundaryfaces;i++){
    int ifa = ps->simu->macromesh.boundaryface[i];
    int locfaL = ps->simu->macromesh.face2elem[4 * ifa + 1];
    int ie = ps->simu->macromesh.face2elem[4 * ifa ];
    int ieR = ps->simu->macromesh.face2elem[4 * ifa + 2];
    assert(ieR < 0);
    if (ieR<0){
      field *f = &ps->simu->fd[ie];

      for(int ipglf = 0;ipglf < NPGF(f->deg,f->raf,locfaL); ipglf++){
	real xpgref[3], xpgref_in[3], wpg;
          
	// Get the coordinates of the Gauss point and coordinates of a
	// point slightly inside the opposite element in xref_in
	int ipg = ref_pg_face(f->deg, f->raf, locfaL, ipglf, xpgref, &wpg, xpgref_in);
	int ino_dg = ipg + ie * ps->npgmacrocell;
	int ino_fe = ps->dg_to_fe_index[ino_dg];
	// Normal vector at gauss point ipgL
	real vnds[3], xpg[3];
	{
	  real dtau[3][3], codtau[3][3];
	  Ref2Phy(f->physnode,
		  xpgref,
		  NULL, locfaL, // dpsiref, ifa
		  xpg, dtau,
		  codtau, NULL, vnds); // codtau, dpsi, vnds
	}
          
	// the boundary flux is an affine function
	real w0[f->model.m],flux0[f->model.m];
	for(int ivv=0; ivv < ps->nb_phy_vars; ivv++) w0[ivv]=0;
	ps->bc_flux(ps,xpg,w0,vnds,flux0);
	for(int var=0; var < ps->nb_phy_vars; var++){
	  int ipot = f->varindex(f->deg,f->raf,f->model.m,
				  ipg,ps->list_of_var[var]);
	  int ipot_fe = ino_fe*ps->nb_phy_vars + var;
	  real val = flux0[ps->list_of_var[var]] * wpg;
	  ps->lsol.rhs[ipot_fe] -= val;
	}
      }	    
    }
  }
}

void Periodic_BoundaryCondition_Poisson1D(void * cs){
  ContinuousSolver * ps=cs;
  field* f0 = &ps->simu->fd[0];

  real charge_average;
  charge_average=0;
  charge_average=Computation_charge_average(ps->simu);

  real surf = 0;

  for(int ie = 0; ie < ps->nbel; ie++){  

    int iemacro = ie / (f0->raf[0] * f0->raf[1] * f0->raf[2]);
    int isubcell = ie % (f0->raf[0] * f0->raf[1] * f0->raf[2]);  
 
    for(int iloc = 0; iloc < ps->nnodes; iloc++){
      real wpg;
      real xref[3];
      //int ipgmacro = ipg + isubcell * nnodes;
      int ilocmacro = iloc + isubcell * ps->nnodes;
      ref_pg_vol(f0->deg,f0->raf,ilocmacro,xref,&wpg,NULL);
      real dtau[3][3],codtau[3][3];
      Ref2Phy(ps->simu->fd[iemacro].physnode,
	      xref,NULL,0,NULL,
	      dtau,codtau,NULL,NULL);
      real det = dot_product(dtau[0], codtau[0]);	
      int ino_dg = iloc + ie * ps->nnodes;
      int ino_fe = ps->dg_to_fe_index[ino_dg];
      int imem = f0->varindex(f0->deg,f0->raf,f0->model.m,
			      ilocmacro,_INDEX_RHO) ;
	//+ iemacro * NPG(f0->deg,f0->raf) * f0->model.m ;
      real rho = ps->simu->fd[iemacro].wn[imem];
      ps->lsol.rhs[ino_fe] -= charge_average * wpg * det ; 
      surf += wpg * det ;  
    }
 
  }
  
  AddLinearSolver(&ps->lsol,0,0,1e20);
  AddLinearSolver(&ps->lsol,ps->lsol.neq-1,ps->lsol.neq-1,1e20);

  ps->lsol.rhs[0]=0;
  ps->lsol.rhs[ps->lsol.neq-1]=0;
}  

void ContinuousOperator_Poisson1D(void * cs){
  ContinuousSolver * ps=cs;  
  ps->diff_op=calloc(ps->nb_phy_vars*ps->nb_phy_vars,sizeof(SDO));
  real D[4][4] = {{0,0,0,0},
                 {0,1,0,0},
                 {0,0,0,0},
                 {0,0,0,0}};
  for (int i=0;i<4;i++){
    for (int j=0;j<4;j++){
      ps->diff_op[0].DO[i][j]=D[i][j];
    }
  }
  GenericOperator_Continuous(ps);
}

void ContinuousOperator_Poisson2D(void * cs){
  ContinuousSolver * ps=cs;
  ps->diff_op=calloc(ps->nb_phy_vars*ps->nb_phy_vars,sizeof(SDO));
  real D[4][4] = {{0,0,0,0},
                 {0,1,0,0},
                 {0,0,1,0},
                 {0,0,0,0}};
  for (int i=0;i<4;i++){
    for (int j=0;j<4;j++){
      ps->diff_op[0].DO[i][j]=D[i][j];
    }
  }
  GenericOperator_Continuous(ps);
}


void RobinFlux(void * cs, real * xpg, real * w, real *vnorm, real * flux){
  ContinuousSolver * ps=cs;
  real alpha=-10000000;
  real beta=0;
  real p0=0;
  real Win[3];
    
  ps->simu->fd[0].model.ImposedData(xpg,ps->simu->dt,Win);
  p0=Win[0];


  flux[0]= alpha * w[0]- beta * p0;
}





void LowerOrderPC_Poisson(LinearSolver * lsol, Simulation *simu, real* globalSol, real*globalRHS){
  Simulation simu2;
  EmptySimulation(&simu2);
  int deg_lo[]={1, simu->fd[0].deg[1], simu->fd[0].deg[2]};
  int raf_lo[]={simu->fd[0].raf[0], simu->fd[0].raf[1], simu->fd[0].raf[2]};
  int HighOrder=simu->fd[0].deg[0];
  
  InitSimulation(&simu2, &simu->macromesh, deg_lo, raf_lo, &simu->fd[0].model);
  ContinuousSolver ps;
  
  int nb_var=1;
  int * listvar= malloc(nb_var * sizeof(int));
  listvar[0]=_INDEX_PHI;
  InitContinuousSolver(&ps,&simu2,1,nb_var,listvar);
  
  ps.matrix_assembly=ContinuousOperator_Poisson1D;
  ps.rhs_assembly=NULL;

  /////////////////// Restriction /////////////////
  real * vector_reduced=NULL;
  vector_reduced=calloc(ps.nb_fe_nodes,sizeof(real));
  Restriction1D_Pq_P1(&ps,HighOrder,globalRHS,vector_reduced);

   for (int i=0;i<ps.nb_fe_nodes;i++){
     ps.lsol.rhs[i]=vector_reduced[i];
   }

  
   ///////////////////////////////////////////////
   
  ps.bc_assembly= ExactDirichletContinuousMatrix;
  ps.postcomputation_assembly=NULL;

  ps.lsol.solver_type=LU;
  ps.lsol.tol=1.0e-4;
  ps.lsol.pc_type=NONE;
  ps.lsol.iter_max=10000;
  ps.lsol.restart_gmres=30;
  ps.lsol.is_CG=true;

  ps.matrix_assembly(&ps);
  ps.bc_assembly(&ps);
  Advanced_SolveLinearSolver(&ps.lsol,ps.simu);
  
  ///////////////////////// Interpolation /////////////////
  for (int i=0;i<ps.nb_fe_nodes;i++){
    vector_reduced[i]=ps.lsol.sol[i];
  }
  Interpolation1D_P1_Pq(&ps,HighOrder,vector_reduced,globalSol);
  free(vector_reduced);
  //////////////////////////////////////////// /////////////////
}


void Interpolation1D_P1_Pq(ContinuousSolver *cs_LowOrder,int HighOrder, real * VecIn, real * VecOut){
  field* f0 = &cs_LowOrder->simu->fd[0];
  real Polynome_1[2];

  Simulation simu;
  EmptySimulation(&simu);
  int deg_ho[]={HighOrder, 0, 0};
  int raf_ho[]={f0->raf[0], f0->raf[1], f0->raf[2]};
  InitSimulation(&simu, &cs_LowOrder->simu->macromesh, deg_ho, raf_ho, &cs_LowOrder->simu->fd[0].model);
  
  ContinuousSolver cs_HighOrder;
  int nb_var=1;
  int * listvar= malloc(nb_var * sizeof(int));
  listvar[0]=_INDEX_PHI;
  InitContinuousSolver(&cs_HighOrder,&simu,1,nb_var,listvar);
 
      
  for(int ie = 0; ie< cs_HighOrder.nbel; ie++){
    Polynome_1[0]=VecIn[ie];
    Polynome_1[1]=VecIn[ie+cs_LowOrder->nnodes-1];
    
    int isubcell = ie % (f0->raf[0] * f0->raf[1] * f0->raf[2]);
    int iemacro = ie / (f0->raf[0] * f0->raf[1] * f0->raf[2]);

    real wpg;
    real xref[3], xref_begin[3], xref_end[3];
    real dtau[3][3],codtau[3][3];
    int ipgmacro_begin = isubcell * cs_HighOrder.nnodes;
    int ipgmacro_end = cs_HighOrder.nnodes - 1 + isubcell * cs_HighOrder.nnodes;

    ref_pg_vol(deg_ho,raf_ho,ipgmacro_begin,xref_begin,&wpg,NULL);
    Ref2Phy(f0[iemacro].physnode,
	       xref_begin,NULL,0,NULL,
	       dtau,codtau,NULL,NULL);

    ref_pg_vol(deg_ho,raf_ho,ipgmacro_end,xref_end,&wpg,NULL);
    Ref2Phy(f0[iemacro].physnode,
	       xref_end,NULL,0,NULL,
	       dtau,codtau,NULL,NULL);

    real h=xref_end[0]-xref_begin[0];    
    
    for(int ipg = 0;ipg < cs_HighOrder.nnodes; ipg++){
 
      int ino_dg= ipg + isubcell * cs_HighOrder.nnodes;
 
      ref_pg_vol(deg_ho,raf_ho,ino_dg,xref,&wpg,NULL);
      real dtau[3][3],codtau[3][3];    
      Ref2Phy(f0[iemacro].physnode,
	       xref,NULL,0,NULL,
	       dtau,codtau,NULL,NULL);
      
      int ino_fe = cs_HighOrder.dg_to_fe_index[ino_dg];
      VecOut[ino_fe]=((Polynome_1[1]-Polynome_1[0])/h)*(xref[0]-xref_begin[0])+Polynome_1[0];
    }
  }
}

void Restriction1D_Pq_P1(ContinuousSolver *cs_LowOrder,int HighOrder, real * VecIn, real * VecOut){
  field* f0 = &cs_LowOrder->simu->fd[0];

  Simulation simu;
  EmptySimulation(&simu);
  int deg_ho[]={HighOrder, 0, 0};
  int raf_ho[]={f0->raf[0], f0->raf[1], f0->raf[2]};
  InitSimulation(&simu, &cs_LowOrder->simu->macromesh, deg_ho, raf_ho, &cs_LowOrder->simu->fd[0].model);
  
  ContinuousSolver cs_HighOrder;
  int nb_var=1;
  int * listvar= malloc(nb_var * sizeof(int));
  listvar[0]=_INDEX_PHI;
  InitContinuousSolver(&cs_HighOrder,&simu,1,nb_var,listvar);

      
  for(int ie = 0; ie< cs_LowOrder->nbel; ie++){
    
    int isubcell = ie % (f0->raf[0] * f0->raf[1] * f0->raf[2]);
    int iemacro = ie / (f0->raf[0] * f0->raf[1] * f0->raf[2]);

    int ipg_LO_dg_begin = isubcell * cs_LowOrder->nnodes ;
    int ipg_LO_dg_end = cs_LowOrder->nnodes-1 + isubcell * cs_LowOrder->nnodes ;

    int ipg_HO_dg_begin = isubcell * cs_HighOrder.nnodes ;
    int ipg_HO_dg_end = cs_HighOrder.nnodes-1 + isubcell * cs_HighOrder.nnodes ;

    int ino_LO_fe_begin = cs_LowOrder->dg_to_fe_index[ipg_LO_dg_begin];
    int ino_LO_fe_end = cs_LowOrder->dg_to_fe_index[ipg_LO_dg_end];

    int ino_HO_fe_begin = cs_HighOrder.dg_to_fe_index[ipg_HO_dg_begin];
    int ino_HO_fe_end = cs_HighOrder.dg_to_fe_index[ipg_HO_dg_end];
    
    VecOut[ino_LO_fe_begin]=VecIn[ino_HO_fe_begin];
    VecOut[ino_LO_fe_end]=VecIn[ino_HO_fe_end];
    
  }
}
