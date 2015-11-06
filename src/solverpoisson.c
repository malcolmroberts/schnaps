#include "solverpoisson.h"
#include "solvercontinuous.h"
#include "geometry.h"
#include "quantities_vp.h"
#include "linear_solver.h"


void Computation_ElectricField_Poisson(void * cs){

  ContinuousSolver * ps=cs;
  
  ContinuousToDiscontinuous_Copy(ps);

  printf("Compute electric field...\n");

  
  for(int ie = 0; ie < ps->simu->macromesh.nbelems; ie++){
    ComputeElectricField(&ps->simu->fd[ie]);
  }
}

void RHSPoisson_Continuous(void * cs){

  ContinuousSolver * ps=cs;
  
  field* f0 = &ps->simu->fd[0];

  real charge_average;
  charge_average=0;

  /*if(ps->type_bc == _Periodic_Poisson_BC){
    charge_average=Computation_charge_average(ps->simu,ps->w);
  }
  else {
    charge_average=0;
    }*/
  
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
      ps->lsol.rhs[ino_fe] += (rho-charge_average)  * wpg * det ; // TODO: put the actual charge	
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
