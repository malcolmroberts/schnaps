#include "solvercontinuous.h"
#include "geometry.h"
#include "quantities_vp.h"
#include "linear_solver.h"
#include "macromesh.h"
#include "solverwave.h"



void BoundaryConditionFriedrichsAssembly(void * cs,LinearSolver* lsol){
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
	  
	  ps->bc_flux(ps,lsol,xpg,wL,vnds,flux0);

	  for(int iv1 = 0; iv1 < m; iv1++) {
	    //int imem1 = fL->varindex(fL->deg, fL->raf,fL->model.m, ipgL, iv1) + offsetL;
	    int ipot_fe1 = ino_fe*ps->nb_phy_vars + iv1;
	    for(int iv = 0; iv < m; iv++) {
	      wL[iv] = (iv == iv1);
	    }
	    
	    real flux[m];
	    ps->bc_flux(ps,lsol,xpg,wL,vnds,flux);
	    
	    for(int iv2 = 0; iv2 < m; iv2++) {
	      // The basis functions is also the gauss point index
	      //int imem2 = fL->varindex(fL->deg, fL->raf,fL->model.m, ipgL, iv2) + offsetL;
	      int ipot_fe2 = ino_fe*ps->nb_phy_vars + iv2;
	      real val =  (flux[iv2]-flux0[iv2]) * wpg;
	      AddLinearSolver(lsol, ipot_fe2, ipot_fe1, val);
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
	ps->bc_flux(ps,lsol,xpg,w0,vnds,flux0);
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



void Wave_test(ContinuousSolver* cs,real theta, real dt){

  real h=cs->simu->vmax*dt*theta;
  real tab[4][4];
  
  real waveMat[9][4][4] ={{{1.0,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0}},
                          {{0.0,0,0,0},
                           {-h,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0}},
                          {{0,0,0,0},
                           {0,0,0,0},
                           {-h,0,0,0},
                           {0,0,0,0}},
                          {{0,0,0,0},
                           {-h,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0}},
                          {{1.0,0,0,0}, 
                           {0,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0}},
                          {{0.0,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0}},
                          {{0,0,0,0},
                           {0,0,0,0},
                           {-h,0,0,0},
                           {0,0,0,0}},
                          {{0.0,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0}},
                          {{1.0,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0}}};
  
  cs->diff_op = calloc(cs->nb_phy_vars*cs->nb_phy_vars,sizeof(SDO));
  
  for (int i=0; i<cs->nb_phy_vars*cs->nb_phy_vars; i++){
    for (int j=0; j<4; j++){
      for (int k=0; k<4; k++){
        cs->diff_op[i].DO[j][k]=waveMat[i][j][k];
      }
    }
  }

  cs->FluxMatrix = calloc(cs->nb_phy_vars,sizeof(real));
  for (int i=0; i<cs->nb_phy_vars; i++){
    cs->FluxMatrix[i] = calloc(cs->nb_phy_vars,sizeof(real));
  }

  cs->FluxMatrix[0][1]=h;
  cs->FluxMatrix[0][2]=h;
  cs->FluxMatrix[1][0]=h;
  cs->FluxMatrix[2][0]=h;

  
}


void Wave_BC_normalvelocity_null(void * cs,LinearSolver* lsol, real * xpg, real * w, real *vnorm, real * flux){

  ContinuousSolver * ps=cs;
  real M[3][3];
  real BC[3][3];
  real lambda=0;
  real p0=0,u0_1=0,u0_2=0;

  real h=ps->FluxMatrix[0][1]/ps->simu->vmax;
  
  M[0][0]=lambda;
  M[0][1]=-(1.0/ps->simu->vmax)*vnorm[0];
  M[0][2]=-(1.0/ps->simu->vmax)*vnorm[1];
  M[1][0]=0;
  M[1][1]=0;
  M[1][2]=0;
  M[2][0]=0;
  M[2][1]=0;
  M[2][2]=0;

  BC[0][0]=M[0][0];
  BC[0][1]=(ps->FluxMatrix[0][1]*vnorm[0]+M[0][1]);
  BC[0][2]=(ps->FluxMatrix[0][2]*vnorm[1]+M[0][2]);
  BC[1][0]=(ps->FluxMatrix[1][0]*vnorm[0]+M[1][0]);
  BC[1][1]=M[1][1];
  BC[1][2]=M[1][2];
  BC[2][0]=(ps->FluxMatrix[2][0]*vnorm[1]+M[2][0]);
  BC[2][1]=M[2][1];
  BC[2][2]=M[2][2];

  flux[0]=(BC[0][0]*w[0]+BC[0][1]*w[1]+BC[0][2]*w[2])-h*(M[0][0]*p0+M[0][1]*u0_1+M[0][2]*u0_2);
  flux[1]=(BC[1][0]*w[0]+BC[1][1]*w[1]+BC[1][2]*w[2])-h*(M[1][0]*p0+M[1][1]*u0_1+M[1][2]*u0_2);
  flux[2]=(BC[2][0]*w[0]+BC[2][1]*w[1]+BC[2][2]*w[2])-h*(M[2][0]*p0+M[2][1]*u0_1+M[2][2]*u0_2);
}



void Wave_BC_pressure_imposed(void * cs,LinearSolver* lsol, real * xpg, real * w, real *vnorm, real * flux){

  ContinuousSolver * ps=cs;
  real M[3][3];
  real BC[3][3];
  real lambda=1.e+12;
  real p0=10,u0_1=0,u0_2=0;

  real h=ps->FluxMatrix[0][1];

  M[0][0]=lambda;
  M[0][1]=-(1.0/ps->simu->vmax)*vnorm[0];
  M[0][2]=-(1.0/ps->simu->vmax)*vnorm[1];
  M[1][0]=0;
  M[1][1]=0;
  M[1][2]=0;
  M[2][0]=0;
  M[2][1]=0;
  M[2][2]=0;

  BC[0][0]=M[0][0];
  BC[0][1]=(ps->FluxMatrix[0][1]*vnorm[0]+M[0][1]);
  BC[0][2]=(ps->FluxMatrix[0][2]*vnorm[1]+M[0][2]);
  BC[1][0]=(ps->FluxMatrix[1][0]*vnorm[0]+M[1][0]);
  BC[1][1]=M[1][1];
  BC[1][2]=M[1][2];
  BC[2][0]=(ps->FluxMatrix[2][0]*vnorm[1]+M[2][0]);
  BC[2][1]=M[2][1];
  BC[2][2]=M[2][2];

  flux[0]=(BC[0][0]*w[0]+BC[0][1]*w[1]+BC[0][2]*w[2])-h*(M[0][0]*p0+M[0][1]*u0_1+M[0][2]*u0_2);
  flux[1]=(BC[1][0]*w[0]+BC[1][1]*w[1]+BC[1][2]*w[2])-h*(M[1][0]*p0+M[1][1]*u0_1+M[1][2]*u0_2);
  flux[2]=(BC[2][0]*w[0]+BC[2][1]*w[1]+BC[2][2]*w[2])-h*(M[2][0]*p0+M[2][1]*u0_1+M[2][2]*u0_2);
}
