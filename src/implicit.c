#include "implicit.h"
#include <stdlib.h>


//void InternalCoupling(field *f,  LinearSolver *solver, int itest);
//void FluxCoupling(field *f,  LinearSolver *solver,int itest);

/* void InternalAssembly(field *f,  LinearSolver *solver,real theta, real dt); */
/* void FluxAssembly(field *f,  LinearSolver *solver,real theta, real dt); */
/* void InterfaceAssembly(field *f,  LinearSolver *solver,real theta, real dt); */
/* void SourceAssembly(field *f,  LinearSolver *solver,real theta, real dt); */
/* void MassAssembly(field *f,  LinearSolver *solver); */

/* void AssemblyImplicitLinearSolver(field *f, LinearSolver *solver,real theta, real dt); */


void InitImplicitLinearSolver(field *f, LinearSolver *solver)
{
  int neq = f->wsize;

  MatrixStorage ms = SKYLINE;
  Solver st = LU; 
  InitLinearSolver(solver,neq,&ms,&st);

  int itest = 1;

  for (int isky=0 ; isky < itest; isky++){
  
    InternalCoupling(f, solver, isky);
    FluxCoupling(f, solver, isky);
    InterfaceCoupling(f, solver, isky);

    if (isky == 0) AllocateLinearSolver(solver);
  }

  //DisplayLinearSolver(solver);
}

void AssemblyImplicitLinearSolver(field *f, LinearSolver *solver, real theta,
				  real tnow, real dt)
{
  if(solver->mat_is_assembly == false) {
    MassAssembly(f, solver);
    InternalAssembly(f, solver,theta,dt);
    FluxAssembly(f, solver,theta,dt);
    InterfaceAssembly(f, solver,theta,dt);
  }

  if(solver->rhs_is_assembly == false){
    for(int i=0;i<solver->neq;i++){
      solver->rhs[i]=0;
    }
    SourceAssembly(f, solver, theta, tnow, dt);
      
  }
  //DisplayLinearSolver(solver);
}

void ThetaTimeScheme(field *f, real tmax, real dt)
{
  LinearSolver solver_implicit;
  LinearSolver solver_explicit;  

  real theta=0.5;
  
  int itermax = tmax / dt+1;
  InitImplicitLinearSolver(f, &solver_implicit);
  InitImplicitLinearSolver(f, &solver_explicit);
  real *res = calloc(f->wsize, sizeof(real));

  real tnow = 0;

  for(int tstep=0;tstep<itermax;tstep++){
    if(tstep == 0) { 
      solver_implicit.mat_is_assembly=false;
      solver_explicit.mat_is_assembly=false;
    } else  { 
      solver_implicit.mat_is_assembly=true;
      solver_explicit.mat_is_assembly=true;
    } 

    solver_implicit.rhs_is_assembly=false;
    solver_explicit.rhs_is_assembly=false;

    
    AssemblyImplicitLinearSolver(f, &solver_explicit,-(1.0-theta), tnow, dt);
    tnow = tnow + dt;
    AssemblyImplicitLinearSolver(f, &solver_implicit, theta, tnow, dt);
  
    MatVect(&solver_explicit, f->wn, res);

    for(int i = 0; i < solver_implicit.neq; i++) {
      solver_implicit.rhs[i] = -solver_explicit.rhs[i] + solver_implicit.rhs[i]
	+ res[i];
    }
  
    SolveLinearSolver(&solver_implicit);

    // FIXME!!!!
    /* for(int i=0;i<solver_implicit.neq;i++){ */
    /*   f->w[i] = solver_implicit.sol[i]; */
    /* } */
    
    int freq = (1 >= itermax / 10)? 1 : itermax / 10;
    if (tstep % freq == 0)
      printf("t=%f iter=%d/%d dt=%f\n", tnow, tstep, itermax, dt);
  }
}

void InternalCoupling(field *f,  LinearSolver *solver, int isky){

  //for(int isky = 0; isky < itest; isky++){
  
  for(int ie = 0; ie < f->macromesh.nbelems; ie++){
    //    field *f = f->fd + ie; // &(f->fd[ie])
    int offsetw = f->wsize * ie;
    
    const int m = f->model.m;
    int deg[3] = {f->deg[0],
		  f->deg[1],
		  f->deg[2]};
    const int npg[3] = {deg[0] + 1,
			deg[1] + 1,
			deg[2] + 1};
    int raf[3] = {f->raf[0],
		   f->raf[1],
		   f->raf[2]};

    //const unsigned int sc_npg = npg[0] * npg[1] * npg[2];

    // Loop on the subcells
    for(int icL0 = 0; icL0 < raf[0]; icL0++) {
      for(int icL1 = 0; icL1 < raf[1]; icL1++) {
	for(int icL2 = 0; icL2 < raf[2]; icL2++) {

	  int icL[3] = {icL0, icL1, icL2};
	  // get the L subcell id
	  int ncL = icL[0] + raf[0] * (icL[1] + raf[1] * icL[2]);
	  // first glop index in the subcell
	  int offsetL = npg[0] * npg[1] * npg[2] * ncL;
	  // loop in the "cross" in the three directions
	  for(int dim0 = 0; dim0 < 3; dim0++) {
	    //for(int dim0 = 0; dim0 < 2; dim0++) {  // TODO : return to 3d !

	    // point p at which we compute the flux
	    for(int p0 = 0; p0 < npg[0]; p0++) {
	      for(int p1 = 0; p1 < npg[1]; p1++) {
		for(int p2 = 0; p2 < npg[2]; p2++) {

		  int p[3] = {p0, p1, p2};
		  int ipgL = offsetL + p[0] + npg[0] * (p[1] + npg[1] * p[2]);
		  for(int iv1 = 0; iv1 < m; iv1++) {
		    int imemL = f->varindex(deg, raf, m, ie, ipgL, iv1);// + offsetw;

		    int q[3] = {p[0], p[1], p[2]};
		    // loop on the direction dim0 on the "cross"
		    for(int iq = 0; iq < npg[dim0]; iq++) {
		      q[dim0] = (p[dim0] + iq) % npg[dim0];

		      int ipgR = offsetL + q[0] + npg[0] * (q[1] + npg[1] * q[2]);
		      for(int iv2 = 0; iv2 < m; iv2++) {
			int imemR = f->varindex(f->deg, f->raf, f->model.m,
						ie, ipgR, iv2);// + offsetw;
			if (isky ==0) IsNonZero(solver, imemL, imemR);
			if (isky ==1) AddLinearSolver(solver, imemL, imemR,1);
			// dtw[imems[temp]] += flux[iv] * wpgL;
		      }
		    } // iv1
		  } // iq
		} // p2
	      } // p1
	    } // p0

	    
	  } // dim loop


	} // icl2
      } //icl1
    } // icl0

  }
}

void FluxCoupling(field *f, LinearSolver *solver, int isky)
{
  for(int ie = 0; ie < f->macromesh.nbelems; ie++){
    //field *f = f->fd + ie; // &(f->fd[ie])
    //int offsetw = f->wsize * ie;
    const int m = f->model.m;

    const int raf[3] = {f->raf[0],
			 f->raf[1],
			 f->raf[2]};
    const int deg[3] = {f->deg[0],
			f->deg[1],
			f->deg[2]};
    const int npg[3] = {deg[0] + 1,
			deg[1] + 1,
			deg[2] + 1};

    // Loop on the subcells
    for(int icL0 = 0; icL0 < raf[0]; icL0++) {
      for(int icL1 = 0; icL1 < raf[1]; icL1++) {
	for(int icL2 = 0; icL2 < raf[2]; icL2++) {

	  int icL[3] = {icL0, icL1, icL2};

	  // Get the left subcell id
	  int ncL = icL[0] + raf[0] * (icL[1] + raf[1] * icL[2]);
	  // First glop index in the subcell
	  int offsetL = npg[0] * npg[1] * npg[2] * ncL;

	  // Sweeping subcell faces in the three directions
	  for(int dim0 = 0; dim0 < 3; dim0++) { 
	    
	    // Compute the subface flux only if we do not touch the
	    // subcell boundary along the current direction dim0
	    if (icL[dim0] != raf[dim0] - 1) {
	      int icR[3] = {icL[0], icL[1], icL[2]};
	      // The right cell index corresponds to an increment in
	      // the dim0 direction
	      icR[dim0]++;
	      int ncR = icR[0] + raf[0] * (icR[1] + raf[1] * icR[2]);
	      int offsetR = npg[0] * npg[1] * npg[2] * ncR;

	      // FIXME: write only write to L-values (and do both
	      // faces) to parallelise better.

	      const int altdim1[3] = {1, 0, 0};
	      const int altdim2[3] = {2, 2, 1};

	      // now loop on the left glops of the subface
	      //int dim1 = (dim0 + 1)%3, dim2 = (dim0+2)%3;
	      int dim1 = altdim1[dim0];
	      int dim2 = altdim2[dim0];
	      int iL[3];
	      iL[dim0] = deg[dim0];
	      for(iL[dim2] = 0; iL[dim2] < npg[dim2]; iL[dim2]++) {
		for(iL[dim1] = 0; iL[dim1] < npg[dim1]; iL[dim1]++) {
		  // find the right and left glops volume indices

		  int iR[3] = {iL[0], iL[1], iL[2]};
		  iR[dim0] = 0;

		  int ipgL = offsetL 
		    + iL[0] + (deg[0] + 1) * (iL[1] + (deg[1] + 1) * iL[2]);
		  int ipgR = offsetR 
		    + iR[0] + (deg[0] + 1) * (iR[1] + (deg[1] + 1) * iR[2]);
		  //printf("ipgL=%d ipgR=%d\n", ipgL, ipgR);

		  for(int iv1 = 0; iv1 < m; iv1++) {
		    int imemL = f->varindex(f->deg, f->raf, f->model.m, ie,
					    ipgL, iv1);
		    // + offsetw; 

		    // finally distribute the flux on the two sides
		    for(int iv2 = 0; iv2 < m; iv2++) {
		      int imemR = f->varindex(f->deg, f->raf, f->model.m, ie,
					      ipgR, iv2);
		      // + offsetw;
		      if (isky ==0) IsNonZero(solver, imemL, imemR);
		      if (isky ==1) {
			AddLinearSolver(solver, imemL, imemR,1);
			AddLinearSolver(solver, imemR, imemL,1);
		      }
		    }
		  }

		}  // face yhat loop
	      } // face xhat loop
	    } // endif internal face
	  } // dim loop
	} // subcell icl2 loop
      } // subcell icl1 loop
    } // subcell icl0 loop

  }
  

}

void InternalAssembly(field *f, LinearSolver *solver, real theta, real dt)
{

  for(int ie = 0; ie < f->macromesh.nbelems; ie++) {
    //    field *f = f->fd + ie;
    int offsetw = f->wsize * ie;
    
    const int m = f->model.m;

    int deg[3] = {f->deg[0],
		  f->deg[1],
		  f->deg[2]};
    const int npg[3] = {deg[0] + 1,
			deg[1] + 1,
			deg[2] + 1};
    int raf[3] = {f->raf[0],
		   f->raf[1],
		   f->raf[2]};
    
    real physnode[20][3];
    for(int inoloc = 0; inoloc < 20; inoloc++) {
      int ino = f->macromesh.elem2node[20*ie+inoloc];
      physnode[inoloc][0] = f->macromesh.node[3 * ino + 0];
      physnode[inoloc][1] = f->macromesh.node[3 * ino + 1];
      physnode[inoloc][2] = f->macromesh.node[3 * ino + 2];
    }
    
    const unsigned int sc_npg = npg[0] * npg[1] * npg[2];

    // Loop on the subcells
    for(int icL0 = 0; icL0 < raf[0]; icL0++) {
      for(int icL1 = 0; icL1 < raf[1]; icL1++) {
	for(int icL2 = 0; icL2 < raf[2]; icL2++) {
	  
	  int icL[3] = {icL0, icL1, icL2};
	  // get the L subcell id
	  int ncL = icL[0] + raf[0] * (icL[1] + raf[1] * icL[2]);
	  // first glop index in the subcell
	  int offsetL = npg[0] * npg[1] * npg[2] * ncL;
	  
	  // compute all of the xref for the subcell
	  real *xref0 = malloc(sc_npg * sizeof(real));
	  real *xref1 = malloc(sc_npg * sizeof(real));
	  real *xref2 = malloc(sc_npg * sizeof(real));
	  real *omega = malloc(sc_npg * sizeof(real));
	  int *imems = malloc(m * sc_npg * sizeof(int));
	  int pos = 0;
	  for(unsigned int p = 0; p < sc_npg; ++p) {
	    real xref[3];
	    real tomega;
	    
	    ref_pg_vol(f->deg, f->raf, offsetL + p, xref, &tomega, NULL);
	    xref0[p] = xref[0];
	    xref1[p] = xref[1];
	    xref2[p] = xref[2];
	    omega[p] = tomega;
	    
	    for(int im = 0; im < m; ++im) {
	      imems[pos++] = f->varindex(f->deg, f->raf, f->model.m, ie, offsetL
					 + p, im);// + offsetw;
	    }
	  }

	  // loop in the "cross" in the three directions
	  for(int dim0 = 0; dim0 < 3; dim0++) {
	    //for(int dim0 = 0; dim0 < 2; dim0++) {  // TODO : return to 3d !
	    // point p at which we compute the flux
	    
	    for(int p0 = 0; p0 < npg[0]; p0++) {
	      for(int p1 = 0; p1 < npg[1]; p1++) {
		for(int p2 = 0; p2 < npg[2]; p2++) {
		  real wL[m], flux[m];
		  int p[3] = {p0, p1, p2};
		  int ipgL = offsetL + p[0] + npg[0] * (p[1] + npg[1] * p[2]);
		  int q[3] = {p[0], p[1], p[2]};
		  // loop on the direction dim0 on the "cross"
		  for(int iq = 0; iq < npg[dim0]; iq++) {
		    q[dim0] = (p[dim0] + iq) % npg[dim0];
		    real dphiref[3] = {0, 0, 0};
		    // compute grad phi_q at glop p
		    dphiref[dim0] = dlag(deg[dim0], q[dim0], p[dim0]) 
		      * raf[dim0];
		    
		    real xrefL[3] = {xref0[ipgL - offsetL],
				     xref1[ipgL - offsetL],
				     xref2[ipgL - offsetL]};
		    real wpgL = omega[ipgL - offsetL];
		    /* real xrefL[3], wpgL; */
		    /* ref_pg_vol(f->interp_param+1,ipgL,xrefL, &wpgL, NULL); */

		    // mapping from the ref glop to the physical glop
		    real dtau[3][3], codtau[3][3], dphiL[3];
		    Ref2Phy(physnode,
			    xrefL,
			    dphiref, // dphiref
			    -1,  // ifa
			    NULL, // xphy
			    dtau,
			    codtau,
			    dphiL, // dphi
			    NULL);  // vnds


		    for(int iv1 = 0; iv1 < m; iv1++) {
		      int imemL = f->varindex(deg, raf, m,ie, ipgL, iv1);
		      //+ offsetw;
		      for(int iv = 0; iv < m; iv++) {
			wL[iv] = (iv == iv1);
		      }

		      f->model.NumFlux(wL, wL, dphiL, flux);

		      int ipgR = offsetL+q[0]+npg[0]*(q[1]+npg[1]*q[2]);
		      for(int iv2 = 0; iv2 < m; iv2++) {
			real val = theta * dt * flux[iv2] * wpgL;
			int imemR = f->varindex(f->deg, f->raf, f->model.m,
						ie, ipgR, iv2);// + offsetw;
			AddLinearSolver(solver, imemR, imemL,-val);
		      }
		    }
		  } // iq
		} // p2
	      } // p1
	    } // p0

	  } // dim loop

	  free(omega);
	  free(xref0);
	  free(xref1);
	  free(xref2);
	  free(imems);

	} // icl2
      } //icl1
    } // icl0

      
  }

}


void FluxAssembly(field *f, LinearSolver *solver, real theta, real dt)
{
  for(int ie = 0; ie < f->macromesh.nbelems; ie++){
    //    field *f = f->fd + ie; // &(f->fd[ie])
    int offsetw = f->wsize * ie;
    
    const int m = f->model.m;
    int raf[3] = {f->raf[0],
		   f->raf[1],
		   f->raf[2]};
    int deg[3] = {f->deg[0],
		  f->deg[1],
		  f->deg[2]};
    int npg[3] = {deg[0] + 1,
		  deg[1] + 1,
		  deg[2] + 1};

    real physnode[20][3];
    for(int inoloc = 0; inoloc < 20; inoloc++) {
      int ino = f->macromesh.elem2node[20*ie+inoloc];
      physnode[inoloc][0] = f->macromesh.node[3 * ino + 0];
      physnode[inoloc][1] = f->macromesh.node[3 * ino + 1];
      physnode[inoloc][2] = f->macromesh.node[3 * ino + 2];
    }
     
    // Loop on the subcells
    for(int icL0 = 0; icL0 < raf[0]; icL0++) {
      for(int icL1 = 0; icL1 < raf[1]; icL1++) {
	for(int icL2 = 0; icL2 < raf[2]; icL2++) {
	  
	  int icL[3] = {icL0, icL1, icL2};
	  
	  // Get the left subcell id
	  int ncL = icL[0] + raf[0] * (icL[1] + raf[1] * icL[2]);
	  // First glop index in the subcell
	  int offsetL = npg[0] * npg[1] * npg[2] * ncL;
	  
	  // Sweeping subcell faces in the three directions
	  for(int dim0 = 0; dim0 < 3; dim0++) { 
	    
	    // Compute the subface flux only if we do not touch the
	    // subcell boundary along the current direction dim0
	    if (icL[dim0] != raf[dim0] - 1) {
	      int icR[3] = {icL[0], icL[1], icL[2]};
	      // The right cell index corresponds to an increment in
	      // the dim0 direction
	      icR[dim0]++;
	      int ncR = icR[0] + raf[0] * (icR[1] + raf[1] * icR[2]);
	      int offsetR = npg[0] * npg[1] * npg[2] * ncR;
	      
	      // FIXME: write only write to L-values (and do both
	      // faces) to parallelise better.
	      
	      const int altdim1[3] = {1, 0, 0};
	      const int altdim2[3] = {2, 2, 1};
	      
	      // now loop on the left glops of the subface
	      //int dim1 = (dim0 + 1)%3, dim2 = (dim0+2)%3;
	      int dim1 = altdim1[dim0];
	      int dim2 = altdim2[dim0];
	      int iL[3];
	      iL[dim0] = deg[dim0];
	      for(iL[dim2] = 0; iL[dim2] < npg[dim2]; iL[dim2]++) {
		for(iL[dim1] = 0; iL[dim1] < npg[dim1]; iL[dim1]++) {
		  // find the right and left glops volume indices

		  int iR[3] = {iL[0], iL[1], iL[2]};
		  iR[dim0] = 0;
		  
		  int ipgL = offsetL 
		    + iL[0] + (deg[0] + 1) * (iL[1] + (deg[1] + 1) * iL[2]);
		  int ipgR = offsetR 
		    + iR[0] + (deg[0] + 1) * (iR[1] + (deg[1] + 1) * iR[2]);
		  //printf("ipgL=%d ipgR=%d\n", ipgL, ipgR);

		  // Compute the normal vector for integrating on the
		  // face
		  real vnds[3];
		  {
		    real xref[3], wpg3;
		    ref_pg_vol(f->deg, f->raf, ipgL, xref, &wpg3, NULL);
		    // mapping from the ref glop to the physical glop
		    real dtau[3][3], codtau[3][3];
		    Ref2Phy(physnode, // FIXME: L or R???
			    xref,
			    NULL, // dphiref
			    -1,  // ifa
			    NULL, // xphy
			    dtau,
			    codtau,
			    NULL, // dphi
			    NULL);  // vnds
		    // we compute ourself the normal vector because we
		    // have to take into account the subcell surface

		    real h1h2 = 1. / raf[dim1] / raf[dim2];
		    vnds[0] = codtau[0][dim0] * h1h2;
		    vnds[1] = codtau[1][dim0] * h1h2;
		    vnds[2] = codtau[2][dim0] * h1h2;
		  }


		  real wpg
		    = wglop(deg[dim1], iL[dim1])
		    * wglop(deg[dim2], iL[dim2]);
		  
		  
		  // numerical flux from the left and right state and
		  // normal vector
		  real wL[m], wR[m], flux[m];

		  for (int iv1 = 0; iv1 < m; iv1++){

			    

		    for(int iv = 0; iv < m; iv++) {
		      wL[iv] = (iv == iv1);
		      wR[iv] = 0;
		    }
		    int imem1 = f->varindex(f->deg, f->raf, f->model.m, ie,
					    ipgL, iv1);
		    //+offsetw;

		    f->model.NumFlux(wL, wR, vnds, flux);	

		    for(int iv2 = 0; iv2 < m; iv2++) {
		      int imem2 = f->varindex(f->deg, f->raf, f->model.m,
					      ie, ipgL, iv2); //+offsetw;
		      real val = theta * dt * flux[iv2] * wpg;		      
		      AddLinearSolver(solver, imem2, imem1, val);
		      
		      imem2 = f->varindex(f->deg, f->raf, f->model.m,
					  ie, ipgR, iv2);//+offsetw;
		      val = theta * dt * flux[iv2] * wpg;		      
		      AddLinearSolver(solver, imem2, imem1, -val);
		    }
		  
		    for(int iv = 0; iv < m; iv++) {
		      wL[iv] = 0;
		      wR[iv] = (iv == iv1);
		    }
		    imem1 = f->varindex(f->deg, f->raf, f->model.m,
					ie, ipgR, iv1); //+offsetw;


		    f->model.NumFlux(wL, wR, vnds, flux);

		    for(int iv2 = 0; iv2 < m; iv2++) {
		      int imem2 = f->varindex(f->deg, f->raf, f->model.m,
					      ie, ipgL, iv2) ; //+offsetw;
		      real val = theta * dt * flux[iv2] * wpg;
		      AddLinearSolver(solver, imem2, imem1, val);
		    
		      imem2 = f->varindex(f->deg, f->raf, f->model.m,
					  ie, ipgR, iv2); //+offsetw;
		      val = theta *dt * flux[iv2] * wpg;		    
		      AddLinearSolver(solver, imem2, imem1, -val);
		    }
		  }

		}  // face yhat loop
	      } // face xhat loop
	    } // endif internal face
	  } // dim loop
	} // subcell icl2 loop
      } // subcell icl1 loop
    } // subcell icl0 loop

  }
}

void MassAssembly(field *f, LinearSolver *solver)
{
  for(int ie = 0; ie < f->macromesh.nbelems; ie++) {
    //    field *f = f->fd + ie;
    int offsetw = f->wsize * ie;
    
    const int m = f->model.m;

    int deg[3] = {f->deg[0],
		  f->deg[1],
		  f->deg[2]};

    int raf[3] = {f->raf[0],
		   f->raf[1],
		   f->raf[2]};

    real physnode[20][3];
    for(int inoloc = 0; inoloc < 20; inoloc++) {
      int ino = f->macromesh.elem2node[20*ie+inoloc];
      physnode[inoloc][0] = f->macromesh.node[3 * ino + 0];
      physnode[inoloc][1] = f->macromesh.node[3 * ino + 1];
      physnode[inoloc][2] = f->macromesh.node[3 * ino + 2];
    }
    
    for(int ipg = 0; ipg < NPG(f->deg, f->raf); ipg++) {
      real dtau[3][3], codtau[3][3], xpgref[3], xphy[3], wpg;
      ref_pg_vol(f->deg, f->raf, ipg, xpgref, &wpg, NULL);
      Ref2Phy(physnode, // phys. nodes
	      xpgref, // xref
	      NULL, -1, // dpsiref, ifa
	      xphy, dtau, // xphy, dtau
	      codtau, NULL, NULL); // codtau, dpsi, vnds
      real det = dot_product(dtau[0], codtau[0]);
      for(int iv1 = 0; iv1 < m; iv1++) {
	int imem = f->varindex(deg, raf, m, ie, ipg, iv1); //+offsetw;
	real val = wpg * det;
	AddLinearSolver(solver, imem, imem,val);
      }
    }
  }
}

void SourceAssembly(field *f, LinearSolver *solver, real theta, real tnow,
		    real dt)
{
  if(f->model.Source != NULL) {
    for(int ie = 0; ie < f->macromesh.nbelems; ie++) {
      //field *f = f->fd + ie;
      int offsetw = f->wsize * ie;
    
      const int m = f->model.m;
    
      int deg[3] = {f->deg[0],
		    f->deg[1],
		    f->deg[2]};

      int raf[3] = {f->raf[0],
		     f->raf[1],
		     f->raf[2]};

      real physnode[20][3];
      for(int inoloc = 0; inoloc < 20; inoloc++) {
	int ino = f->macromesh.elem2node[20*ie+inoloc];
	physnode[inoloc][0] = f->macromesh.node[3 * ino + 0];
	physnode[inoloc][1] = f->macromesh.node[3 * ino + 1];
	physnode[inoloc][2] = f->macromesh.node[3 * ino + 2];
      }
      
      for(int ipg = 0; ipg < NPG(f->deg, f->raf); ipg++) {
	real dtau[3][3], codtau[3][3], xpgref[3], xphy[3], wpg;
	ref_pg_vol(f->deg, f->raf, ipg, xpgref, &wpg, NULL);
	Ref2Phy(physnode, // phys. nodes
		xpgref, // xref
		NULL, -1, // dpsiref, ifa
		xphy, dtau, // xphy, dtau
		codtau, NULL, NULL); // codtau, dpsi, vnds
	real det = dot_product(dtau[0], codtau[0]);
	real wL[m], source[m];
	/* for(int iv = 0; iv < m; ++iv){ */
	/* 	int imem = f->varindex(f->deg, f->raf, f->model.m, ipg, iv); */
	/* 	wL[iv] = w[imem]; */
	/* } */
	f->model.Source(xphy, tnow, wL, source);

	for(int iv1 = 0; iv1 < m; iv1++) {
	  int imem = f->varindex(deg, raf, m, ie, ipg, iv1); //+offsetw;
	  real val = theta * dt * source[iv1] * wpg * det;
	  solver->rhs[imem] += val;
	}
      }
    }
  }
  // assembly of the boundary terms

  int fsize =  f->wsize / f->macromesh.nbelems;

  for(int ifa = 0; ifa < f->macromesh.nbfaces; ifa++){
    int ieL = f->macromesh.face2elem[4 * ifa + 0];
    int locfaL = f->macromesh.face2elem[4 * ifa + 1];
    int ieR = f->macromesh.face2elem[4 * ifa + 2];
    //    field *fL = f->fd + ieL;
    //printf("iel=%d ier=%d\n",ieL,ieR);
    int offsetL = fsize * ieL;
    
    real physnodeL[20][3];
    for(int inoloc = 0; inoloc < 20; inoloc++) {
      int ino = f->macromesh.elem2node[20 * ieL + inoloc];
      physnodeL[inoloc][0] = f->macromesh.node[3 * ino + 0];
      physnodeL[inoloc][1] = f->macromesh.node[3 * ino + 1];
      physnodeL[inoloc][2] = f->macromesh.node[3 * ino + 2];
    }

    if (ieR < 0) {
 
      const unsigned int m = f->model.m;
      
      // Loop over the points on a single macro cell interface.
      for(int ipgfL = 0; ipgfL < NPGF(f->deg, f->raf, locfaL); ipgfL++) {
	
	real xpgref[3], xpgref_in[3], wpg;
	
	// Get the coordinates of the Gauss point and coordinates of a
	// point slightly inside the opposite element in xref_in
	int ipgL = ref_pg_face(f->deg, f->raf, locfaL, ipgfL, xpgref, &wpg,
			       xpgref_in);
	
	// Normal vector at gauss point ipgL
	real vnds[3], xpg[3];
	{
	  real dtau[3][3], codtau[3][3];
	  Ref2Phy(physnodeL,
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

	f->model.BoundaryFlux(xpg, tnow, wL, vnds, flux0);
	
	for(int iv2 = 0; iv2 < m; iv2++) {
	  int imem2 = f->varindex(f->deg, f->raf, f->model.m, ieL, ipgL, iv2);
	  //+ offsetL;
	  real val = theta * dt * flux0[iv2] * wpg;
	  solver->rhs[imem2] -= val;
	}
      }
    } // if ier < 0
  } // macroface loop
 
} // SourceAssembly


void InterfaceCoupling(field *f,  LinearSolver *solver, int itest)
{
  int fsize =  f->wsize / f->macromesh.nbelems;

  for(int ifa = 0; ifa < f->macromesh.nbfaces; ifa++){
    int ieL = f->macromesh.face2elem[4 * ifa + 0];
    int locfaL = f->macromesh.face2elem[4 * ifa + 1];
    int ieR = f->macromesh.face2elem[4 * ifa + 2];
    //field *fL = f->fd + ieL;
    //field *fR = NULL;
    int offsetR = -1;
    //printf("iel=%d ier=%d\n",ieL,ieR);
    int offsetL = fsize * ieL;
    if (ieR >= 0) {
      //fR = f->fd + ieR;
      offsetR = fsize * ieR;
    }

    const unsigned int m = f->model.m;

    // Loop over the points on a single macro cell interface.
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int ipgfL = 0; ipgfL < NPGF(f->deg, f->raf, locfaL); ipgfL++) {

      real xpgref[3], xpgref_in[3], wpg;
    
      // Get the coordinates of the Gauss point and coordinates of a
      // point slightly inside the opposite element in xref_in
      int ipgL = ref_pg_face(f->deg, f->raf, locfaL, ipgfL, xpgref, &wpg,
			     xpgref_in);


      real physnodeL[20][3];
      for(int inoloc = 0; inoloc < 20; inoloc++) {
	int ino = f->macromesh.elem2node[20 * ieL + inoloc];
	physnodeL[inoloc][0] = f->macromesh.node[3 * ino + 0];
	physnodeL[inoloc][1] = f->macromesh.node[3 * ino + 1];
	physnodeL[inoloc][2] = f->macromesh.node[3 * ino + 2];
      }

      
      // Normal vector at gauss point ipgL
      real vnds[3], xpg[3];
      {
	real dtau[3][3], codtau[3][3];
	Ref2Phy(physnodeL,
		xpgref,
		NULL, locfaL, // dpsiref, ifa
		xpg, dtau,
		codtau, NULL, vnds); // codtau, dpsi, vnds
      }
    
      if (ieR >= 0) {  // the right element exists


	real physnodeR[20][3];
	for(int inoloc = 0; inoloc < 20; inoloc++) {
	  int ino = f->macromesh.elem2node[20 * ieR + inoloc];
	  physnodeR[inoloc][0] = f->macromesh.node[3 * ino + 0];
	  physnodeR[inoloc][1] = f->macromesh.node[3 * ino + 1];
	  physnodeR[inoloc][2] = f->macromesh.node[3 * ino + 2];
	}

	real xrefL[3];
	{
	  real xpg_in[3];
	  Ref2Phy(physnodeL,
		  xpgref_in,
		  NULL, -1, // dpsiref, ifa
		  xpg_in, NULL,
		  NULL, NULL, NULL); // codtau, dpsi, vnds
	  PeriodicCorrection(xpg_in,f->period);
	  Phy2Ref(physnodeR, xpg_in, xrefL);
	}
      
	int ipgR = ref_ipg(f->deg, f->raf, xrefL);
	for (int iv1 = 0; iv1 < m; iv1++){
	  int imem1 = f->varindex(f->deg, f->raf, f->model.m, ieL, ipgL, iv1);
	  //+ offsetL;

	  for(int iv2 = 0; iv2 < m; iv2++) {
	    int imem2 = f->varindex(f->deg, f->raf, f->model.m, ieL, ipgL, iv2);
	    //+ offsetL;		  
	    IsNonZero(solver, imem2, imem1);
		      
	    imem2 = f->varindex(f->deg, f->raf, f->model.m, ieR, ipgR, iv2);
	      //+ offsetR;
	    IsNonZero(solver, imem2, imem1);
	  }
		  
	  imem1 = f->varindex(f->deg, f->raf, f->model.m, ieR, ipgR, iv1);
	  //+ offsetR;
	    
	  for(int iv2 = 0; iv2 < m; iv2++) {
	    int imem2 = f->varindex(f->deg, f->raf, f->model.m, ieL, ipgL, iv2);
	      //+ offsetL;
	    IsNonZero(solver, imem2, imem1);
		    
	    imem2 = f->varindex(f->deg, f->raf, f->model.m, ieR, ipgR, iv2);
	    //+ offsetR;		    
	    IsNonZero(solver, imem2, imem1);
	  }
	}

      } else { // The point is on the boundary.


	for(int iv1 = 0; iv1 < m; iv1++) {
	  int imem1 = f->varindex(f->deg, f->raf,f->model.m, ieL, ipgL, iv1);
	  //+ offsetL;

	  for(int iv2 = 0; iv2 < m; iv2++) {
	    int imem2 = f->varindex(f->deg, f->raf,f->model.m, ieL, ipgL, iv2);
	      //+ offsetL;
	    IsNonZero(solver, imem2, imem1);
	  }
	} // iv1

      } // else

  
    } // ipgfl

  } // macroface loop

}

void InterfaceAssembly(field *f,  LinearSolver *solver,real theta, real dt)
{
  int fsize =  f->wsize / f->macromesh.nbelems;

  for(int ifa = 0; ifa < f->macromesh.nbfaces; ifa++){
    int ieL = f->macromesh.face2elem[4 * ifa + 0];

    real physnodeL[20][3];
    for(int inoloc = 0; inoloc < 20; inoloc++) {
      int ino = f->macromesh.elem2node[20 * ieL + inoloc];
      physnodeL[inoloc][0] = f->macromesh.node[3 * ino + 0];
      physnodeL[inoloc][1] = f->macromesh.node[3 * ino + 1];
      physnodeL[inoloc][2] = f->macromesh.node[3 * ino + 2];
    }

    int locfaL = f->macromesh.face2elem[4 * ifa + 1];
    int ieR = f->macromesh.face2elem[4 * ifa + 2];
    //field *fL = f->fd + ieL;
    field *fR = NULL;
    int offsetR = -1;
    //printf("iel=%d ier=%d\n",ieL,ieR);
    int offsetL = fsize * ieL;
    if (ieR >= 0) {
      //      fR = f->fd + ieR;
      offsetR = fsize * ieR;
    }
  
    const unsigned int m = f->model.m;

    // Loop over the points on a single macro cell interface.
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int ipgfL = 0; ipgfL < NPGF(f->deg, f->raf, locfaL); ipgfL++) {

      real xpgref[3], xpgref_in[3], wpg;
    
      // Get the coordinates of the Gauss point and coordinates of a
      // point slightly inside the opposite element in xref_in
      int ipgL = ref_pg_face(f->deg, f->raf, locfaL, ipgfL, xpgref, &wpg,
			     xpgref_in);
    
      // Normal vector at gauss point ipgL
      real vnds[3], xpg[3];
      {
	real dtau[3][3], codtau[3][3];
	Ref2Phy(physnodeL,
		xpgref,
		NULL, locfaL, // dpsiref, ifa
		xpg, dtau,
		codtau, NULL, vnds); // codtau, dpsi, vnds
      }
    
      if (ieR >= 0) {  // the right element exists

	real physnodeR[20][3];
	for(int inoloc = 0; inoloc < 20; inoloc++) {
	  int ino = f->macromesh.elem2node[20 * ieR + inoloc];
	  physnodeR[inoloc][0] = f->macromesh.node[3 * ino + 0];
	  physnodeR[inoloc][1] = f->macromesh.node[3 * ino + 1];
	  physnodeR[inoloc][2] = f->macromesh.node[3 * ino + 2];
	}
    
	real xrefL[3];
	{
	  real xpg_in[3];
	  Ref2Phy(physnodeL,
		  xpgref_in,
		  NULL, -1, // dpsiref, ifa
		  xpg_in, NULL,
		  NULL, NULL, NULL); // codtau, dpsi, vnds
	  PeriodicCorrection(xpg_in,f->period);
	  Phy2Ref(physnodeR, xpg_in, xrefL);
	
	}
      
	int ipgR = ref_ipg(f->deg,f->raf, xrefL);

	real flux[m];
	real wL[m];
	real wR[m];

	for (int iv1 = 0; iv1 < m; iv1++){
	  for(int iv = 0; iv < m; iv++) {
	    wL[iv] = (iv == iv1);
	    wR[iv] = 0;
	  }
	  int imem1 = f->varindex(f->deg, f->raf, f->model.m, ieL, ipgL, iv1);
	  //+ offsetL;

	  // int_dL F(wL, wR, grad phi_ib)

	  f->model.NumFlux(wL, wR, vnds, flux);

	  // Add flux to both sides

	  for(int iv2 = 0; iv2 < m; iv2++) {
	    int imem2 = f->varindex(f->deg, f->raf, f->model.m, ieL, ipgL, iv2);
	    //+ offsetL;		  
	    real val = theta * dt * flux[iv2] * wpg;		      
	    AddLinearSolver(solver, imem2, imem1, val);
		      
	    imem2 = fR->varindex(fR->deg, fR->raf, f->model.m, ieR, ipgR, iv2);
	    //+ offsetR;
	    val = theta * dt * flux[iv2] * wpg;		      
	    AddLinearSolver(solver, imem2, imem1, -val);
	  }
		  
	  for(int iv = 0; iv < m; iv++) {
	    wL[iv] = 0;
	    wR[iv] = (iv == iv1);
	  }
	  imem1 = f->varindex(f->deg, f->raf, f->model.m, ieR, ipgR, iv1);
	    //+ offsetR;

	  f->model.NumFlux(wL, wR, vnds, flux);

	  for(int iv2 = 0; iv2 < m; iv2++) {
	    int imem2 = f->varindex(f->deg, f->raf, f->model.m, ieL, ipgL, iv2);
	    //+ offsetL;
	    real val = theta * dt * flux[iv2] * wpg;
	    AddLinearSolver(solver, imem2, imem1, val);
		    
	    imem2 = fR->varindex(fR->deg, fR->raf, f->model.m, ieR, ipgR, iv2);
	    //+ offsetR;	    
	    val = theta *dt * flux[iv2] * wpg;		    
	    AddLinearSolver(solver, imem2, imem1, -val);
	  }
	}

      } else { // The point is on the boundary.

	// the boundary flux is an affine function
	real flux0[m], wL[m];
	for(int iv = 0; iv < m; iv++) {
	  wL[iv] = 0;
	}
	f->model.BoundaryFlux(xpg, f->tnow, wL, vnds, flux0);

	for(int iv1 = 0; iv1 < m; iv1++) {
	  int imem1 = f->varindex(f->deg, f->raf,f->model.m, ieL, ipgL, iv1);
	  //+ offsetL;

	  for(int iv = 0; iv < m; iv++) {
	    wL[iv] = (iv == iv1);
	  }

	  real flux[m];
	  f->model.BoundaryFlux(xpg, f->tnow, wL, vnds, flux);

	  for(int iv2 = 0; iv2 < m; iv2++) {
	    // The basis functions is also the gauss point index
	    int imem2 = f->varindex(f->deg, f->raf, f->model.m, ieL, ipgL, iv2);
	    //+ offsetL;
	    real val = theta *dt * (flux[iv2]-flux0[iv2]) * wpg;
	    AddLinearSolver(solver, imem2, imem1, val);
	  }
	} // iv1

      } // else

  
    } // ipgfl

  } // macroface loop

}
