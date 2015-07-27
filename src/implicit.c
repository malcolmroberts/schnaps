#include "implicit.h"
#include <stdlib.h>


void InternalCoupling(Simulation *simu,  LinearSolver *solver, int itest);
void FluxCoupling(Simulation *simu,  LinearSolver *solver,int itest);

void InternalAssembly(Simulation *simu,  LinearSolver *solver);
void FluxAssembly(Simulation *simu,  LinearSolver *solver);
void SourceAssembly(Simulation *simu,  LinearSolver *solver);

void AssemblyImplicitLinearSolver(Simulation *simu, LinearSolver *solver);


void InitImplicitLinearSolver(Simulation *simu, LinearSolver *solver){

  int neq = simu->wsize;

  MatrixStorage ms = SKYLINE;
  Solver st = LU; 
  InitLinearSolver(solver,neq,&ms,&st);

  int itest = 2;

  for (int isky=0 ; isky < itest; isky++){
  
    InternalCoupling(simu, solver, isky);
    FluxCoupling(simu, solver, isky);

    if (isky == 0) AllocateLinearSolver(solver);
  }

  //DisplayLinearSolver(solver);

  

}

void AssemblyImplicitLinearSolver(Simulation *simu, LinearSolver *solver){

  
  InternalAssembly(simu, solver);
  FluxAssembly(simu, solver);

  //SourceAssembly(simu, solver);
  

}



void InternalCoupling(Simulation *simu,  LinearSolver *solver, int isky){

  //for(int isky = 0; isky < itest; isky++){
  
  for(int ie = 0; ie < simu->macromesh.nbelems; ie++){
    field *f = simu->fd + ie; // &(simu->fd[ie])

    const int m = f->model.m;
    int deg[3] = {f->deg[0],
		  f->deg[1],
		  f->deg[2]};
    const int npg[3] = {deg[0] + 1,
			deg[1] + 1,
			deg[2] + 1};
    int nraf[3] = {f->raf[0],
		   f->raf[1],
		   f->raf[2]};

    const unsigned int sc_npg = npg[0] * npg[1] * npg[2];

    // Loop on the subcells
    for(int icL0 = 0; icL0 < nraf[0]; icL0++) {
      for(int icL1 = 0; icL1 < nraf[1]; icL1++) {
	for(int icL2 = 0; icL2 < nraf[2]; icL2++) {

	  int icL[3] = {icL0, icL1, icL2};
	  // get the L subcell id
	  int ncL = icL[0] + nraf[0] * (icL[1] + nraf[1] * icL[2]);
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
		    int imemL = f->varindex(deg, nraf, m, ipgL, iv1);

		    int q[3] = {p[0], p[1], p[2]};
		    // loop on the direction dim0 on the "cross"
		    for(int iq = 0; iq < npg[dim0]; iq++) {
		      q[dim0] = (p[dim0] + iq) % npg[dim0];

		      int ipgR = offsetL + q[0] + npg[0] * (q[1] + npg[1] * q[2]);
		      for(int iv2 = 0; iv2 < m; iv2++) {
			int imemR = f->varindex(f->deg,f->raf,f->model.m, ipgR, iv2);
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


void FluxCoupling(Simulation *simu,  LinearSolver *solver, int isky){

  
  for(int ie = 0; ie < simu->macromesh.nbelems; ie++){
    field *f = simu->fd + ie; // &(simu->fd[ie])

    const int m = f->model.m;


    const int nraf[3] = {f->raf[0],
			 f->raf[1],
			 f->raf[2]};
    const int deg[3] = {f->deg[0],
			f->deg[1],
			f->deg[2]};
    const int npg[3] = {deg[0] + 1,
			deg[1] + 1,
			deg[2] + 1};

    // Loop on the subcells
    for(int icL0 = 0; icL0 < nraf[0]; icL0++) {
      for(int icL1 = 0; icL1 < nraf[1]; icL1++) {
	for(int icL2 = 0; icL2 < nraf[2]; icL2++) {

	  int icL[3] = {icL0, icL1, icL2};

	  // Get the left subcell id
	  int ncL = icL[0] + nraf[0] * (icL[1] + nraf[1] * icL[2]);
	  // First glop index in the subcell
	  int offsetL = npg[0] * npg[1] * npg[2] * ncL;

	  // Sweeping subcell faces in the three directions
	  for(int dim0 = 0; dim0 < 3; dim0++) { 
	    
	    // Compute the subface flux only if we do not touch the
	    // subcell boundary along the current direction dim0
	    if (icL[dim0] != nraf[dim0] - 1) {
	      int icR[3] = {icL[0], icL[1], icL[2]};
	      // The right cell index corresponds to an increment in
	      // the dim0 direction
	      icR[dim0]++;
	      int ncR = icR[0] + nraf[0] * (icR[1] + nraf[1] * icR[2]);
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
		    int imemL = f->varindex(f->deg, f->raf, f->model.m, ipgL, iv1); 

		    // finally distribute the flux on the two sides
		    for(int iv2 = 0; iv2 < m; iv2++) {
		      int imemR = f->varindex(f->deg, f->raf, f->model.m, ipgR, iv2);
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

void InternalAssembly(Simulation *simu,  LinearSolver *solver){

  for(int ie = 0; ie < simu->macromesh.nbelems; ie++){
    field *f = simu->fd + ie;
    
    
    const int m = f->model.m;

    int deg[3] = {f->deg[0],
		  f->deg[1],
		  f->deg[2]};
    const int npg[3] = {deg[0] + 1,
			deg[1] + 1,
			deg[2] + 1};
    int nraf[3] = {f->raf[0],
		   f->raf[1],
		   f->raf[2]};

    const unsigned int sc_npg = npg[0] * npg[1] * npg[2];


    // Loop on the subcells
    for(int icL0 = 0; icL0 < nraf[0]; icL0++) {
      for(int icL1 = 0; icL1 < nraf[1]; icL1++) {
	for(int icL2 = 0; icL2 < nraf[2]; icL2++) {
	  
	  int icL[3] = {icL0, icL1, icL2};
	  // get the L subcell id
	  int ncL = icL[0] + nraf[0] * (icL[1] + nraf[1] * icL[2]);
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
	      imems[pos++] = f->varindex(f->deg,f->raf,f->model.m, offsetL + p, im);
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
		      * nraf[dim0];
		    
		    real xrefL[3] = {xref0[ipgL - offsetL],
				     xref1[ipgL - offsetL],
				     xref2[ipgL - offsetL]};
		    real wpgL = omega[ipgL - offsetL];
		    /* real xrefL[3], wpgL; */
		    /* ref_pg_vol(f->interp_param+1,ipgL,xrefL, &wpgL, NULL); */

		    // mapping from the ref glop to the physical glop
		    real dtau[3][3], codtau[3][3], dphiL[3];
		    Ref2Phy(f->physnode,
			    xrefL,
			    dphiref, // dphiref
			    -1,  // ifa
			    NULL, // xphy
			    dtau,
			    codtau,
			    dphiL, // dphi
			    NULL);  // vnds


		    for(int iv1 = 0; iv1 < m; iv1++) {
		      int imemL = f->varindex(deg, nraf, m, ipgL, iv1);
		      for(int iv = 0; iv < m; iv++) {
			wL[iv] = (iv == iv1);
		      }

		      f->model.NumFlux(wL, wL, dphiL, flux);

		      int ipgR = offsetL+q[0]+npg[0]*(q[1]+npg[1]*q[2]);
		      for(int iv2 = 0; iv2 < m; iv2++) {
			real val = flux[iv2] * wpgL;
			int imemR = f->varindex(f->deg,f->raf,f->model.m, ipgR, iv2);
			AddLinearSolver(solver, imemL, imemR,val);
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


void FluxAssembly(Simulation *simu, LinearSolver *solver){


  for(int ie = 0; ie < simu->macromesh.nbelems; ie++){
    field *f = simu->fd + ie; // &(simu->fd[ie])

    const int m = f->model.m;

  

    int nraf[3] = {f->raf[0],
		   f->raf[1],
		   f->raf[2]};
    int deg[3] = {f->deg[0],
		  f->deg[1],
		  f->deg[2]};
    int npg[3] = {deg[0] + 1,
		  deg[1] + 1,
		  deg[2] + 1};

 
    // Loop on the subcells
    for(int icL0 = 0; icL0 < nraf[0]; icL0++) {
      for(int icL1 = 0; icL1 < nraf[1]; icL1++) {
	for(int icL2 = 0; icL2 < nraf[2]; icL2++) {
	  
	  int icL[3] = {icL0, icL1, icL2};
	  
	  // Get the left subcell id
	  int ncL = icL[0] + nraf[0] * (icL[1] + nraf[1] * icL[2]);
	  // First glop index in the subcell
	  int offsetL = npg[0] * npg[1] * npg[2] * ncL;
	  
	  // Sweeping subcell faces in the three directions
	  for(int dim0 = 0; dim0 < 3; dim0++) { 
	    
	    // Compute the subface flux only if we do not touch the
	    // subcell boundary along the current direction dim0
	    if (icL[dim0] != nraf[dim0] - 1) {
	      int icR[3] = {icL[0], icL[1], icL[2]};
	      // The right cell index corresponds to an increment in
	      // the dim0 direction
	      icR[dim0]++;
	      int ncR = icR[0] + nraf[0] * (icR[1] + nraf[1] * icR[2]);
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
		    Ref2Phy(f->physnode,
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

		    real h1h2 = 1. / nraf[dim1] / nraf[dim2];
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

		    f->model.NumFlux(wL, wR, vnds, flux);


		    int imem1 = f->varindex(f->deg, f->raf, f->model.m, ipgL, iv1);

		    for(int iv2 = 0; iv2 < m; iv2++) {
		      int imem2 = f->varindex(f->deg, f->raf, f->model.m, ipgL, iv2);		  
		      real val = flux[iv2] * wpg;		      
		      AddLinearSolver(solver, imem1, imem2, -val);
		      
		      imem2 = f->varindex(f->deg, f->raf, f->model.m, ipgR, iv2);
		      val = flux[iv2] * wpg;		      
		      AddLinearSolver(solver, imem1, imem2, val);
		    }
		  
		    for(int iv = 0; iv < m; iv++) {
		      wL[iv] = 0;
		      wR[iv] = (iv == iv1);
		    }


		    f->model.NumFlux(wL, wR, vnds, flux);

		    for(int iv2 = 0; iv2 < m; iv2++) {
		      int imem2 = f->varindex(f->deg, f->raf, f->model.m, ipgL, iv2);
		      real val = flux[iv2] * wpg;
		      AddLinearSolver(solver, imem1, imem2, -val);
		    
		      imem2 = f->varindex(f->deg, f->raf, f->model.m, ipgR, iv2);		    
		      val = flux[iv2] * wpg;		    
		      AddLinearSolver(solver, imem1, imem2, val);
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
