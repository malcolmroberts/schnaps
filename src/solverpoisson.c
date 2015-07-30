#include "solverpoisson.h"
#include "solvercontinuous.h"
#include "geometry.h"
#include "quantities_vp.h"
#include "linear_solver.h"


void Computation_ElectricField_Poisson(void * cs,LinearSolver* lsol){

  ContinuousSolver * ps=cs;
  
  ContinuousToDiscontinuous_Copy(ps,lsol);

  printf("Compute electric field...\n");

  
  for(int ie = 0; ie < ps->simu->macromesh.nbelems; ie++)
    ComputeElectricField(&ps->simu->fd[ie]);

}

void MatrixPoisson_Continuous(void * cs,LinearSolver* lsol){

  ContinuousSolver * ps=cs;

   field* f0 = &ps->simu->fd[0];

  if(!ps->lsol.mat_is_assembly){
    for(int ie = 0; ie < ps->nbel; ie++){  

      // local matrix 
      real aloc[ps->nnodes][ps->nnodes];
      for(int iloc = 0; iloc < ps->nnodes; iloc++){
	for(int jloc = 0; jloc < ps->nnodes; jloc++){
	  aloc[iloc][jloc] = 0;
	}
      }

      int iemacro = ie / (f0->raf[0] * f0->raf[1] * f0->raf[2]);
      int isubcell = ie % (f0->raf[0] * f0->raf[1] * f0->raf[2]);

      /* real physnode[20][3]; */
      /* for(int ino = 0; ino < 20; ino++) { */
      /* 	int numnoe = ps->simu->macromesh.elem2node[20 * iemacro + ino]; */
      /* 	for(int ii = 0; ii < 3; ii++) { */
      /* 	  physnode[ino][ii] = ps->simu->macromesh.node[3 * numnoe + ii]; */
      /* 	} */
      /* } */
      //ref_pg_vol(ps->fd->interp_param+1,int ipg,
      // real* xpg,real* wpg,real* xpg_in);
      // grad_psi_pg(ps->fd->interp_param+1,ib,ipg,dphiref)
      // Ref2Phy(physnode,xref,dphiref,NULL,NULL,dtau,codtau,dphi,NULL);


      for(int ipg = 0;ipg < ps->nnodes; ipg++){
	real wpg;
	real xref[3];
	int ipgmacro= ipg + isubcell * ps->nnodes;

	ref_pg_vol(f0->deg,f0->raf,ipgmacro,xref,&wpg,NULL);

	for(int iloc = 0; iloc < ps->nnodes; iloc++){
	  real dtau[3][3],codtau[3][3];
	  real dphiref_i[3],dphiref_j[3];
	  real dphi_i[3],dphi_j[3];
	  int ilocmacro = iloc + isubcell * ps->nnodes;
	  grad_psi_pg(f0->deg,f0->raf,ilocmacro,ipgmacro,dphiref_i);

	  Ref2Phy(ps->simu->fd[iemacro].physnode,
		  xref,dphiref_i,0,NULL,
		  dtau,codtau,dphi_i,NULL);
	  real det = dot_product(dtau[0], codtau[0]);
	  for(int jloc = 0; jloc < ps->nnodes; jloc++){
	    int jlocmacro = jloc + isubcell * ps->nnodes;
	    grad_psi_pg(f0->deg,f0->raf,jlocmacro,ipgmacro,dphiref_j);
	    Ref2Phy(ps->simu->fd[iemacro].physnode,
		    xref,dphiref_j,0,NULL,
		    dtau,codtau,dphi_j,NULL);
	    aloc[iloc][jloc] += dot_product(dphi_i, dphi_j) * wpg / det  ;
	  }
	}
      }


      for(int iloc = 0; iloc < ps->nnodes; iloc++){
	for(int jloc = 0; jloc < ps->nnodes; jloc++){
	  int ino_dg = iloc + ie * ps->nnodes;
	  int jno_dg = jloc + ie * ps->nnodes;
	  int ino_fe = ps->dg_to_fe_index[ino_dg];
	  int jno_fe = ps->dg_to_fe_index[jno_dg];
	  real val = aloc[iloc][jloc];
	  AddLinearSolver(&ps->lsol,ino_fe,jno_fe,val);
	}
      }
   
    }
  } 
  
 }

void RHSPoisson_Continuous(void * cs,LinearSolver* lsol){

  ContinuousSolver * ps=cs;
  
  field* f0 = &ps->simu->fd[0];

  real charge_average;
  charge_average=0;

  if(ps->type_bc == _Periodic_Poisson_BC){
    //charge_average=Computation_charge_average(f,w);
  }
  else {
    charge_average=0;
  }
  
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



void SolvePoisson1D(Simulation *simu,real * w,int type_bc, real bc_l, real bc_r,Solver solver_sys,PC precon){

  real charge_average;
  real *bounds = malloc(6 * sizeof(real));
  charge_average=0;

  field *f = &simu->fd[0];
  assert(simu->macromesh.nbelems == 1);

  if(type_bc == _Periodic_Poisson_BC){
    charge_average=Computation_charge_average(f,w);
    //printf(" chare average %e\n",charge_average);
    bc_l=0;
    bc_r=0;
  }
  else {
    charge_average=0;
  }
  
  // for the moment, works only for the 1d case
  assert(simu->macromesh.is1d);

  // assembly of the rigidity matrix

  LinearSolver sky;

  // number of equation of the Poisson solver
  // = number of nodes in the mesh
  int degx=f->interp.interp_param[1];
  int nelx=f->interp.interp_param[4];
  real xmin=simu->macromesh.xmin[0];
  real xmax=simu->macromesh.xmax[0];  //
  real dx=(xmax-xmin)/nelx;
  int neq=degx*nelx+1;
  
  // number of conservatives variables
  // = number of velocity glops + 1 (potential)
  int m = f->model.m;

  InitLinearSolver(&sky,neq,NULL,NULL); //InitSkyline(&sky, neq);
  
  sky.solver_type = solver_sys;
  sky.pc_type=precon;

  if(!sky.is_alloc){
  // compute the profile of the matrix
    for(int ie = 0; ie < nelx; ie++){
      for(int iloc = 0; iloc < degx + 1; iloc++){
	for(int jloc = 0; jloc < degx + 1; jloc++){
	  int ino = iloc + ie * degx;
	  int jno = jloc + ie * degx;
	  IsNonZero(&sky, ino, jno);
	}
      }
    }    
    AllocateLinearSolver(&sky);
  }

  if(!sky.mat_is_assembly){
  // local matrix (assuming identical finite elements)
    real aloc[degx+1][degx+1];
    for(int iloc=0;iloc<degx+1;iloc++){
      for(int jloc=0;jloc<degx+1;jloc++){
          aloc[iloc][jloc]=0;

      }
    }
    for(int ipg=0;ipg<degx+1;ipg++){
      real omega=wglop(degx,ipg);
      for(int iloc=0;iloc<degx+1;iloc++){
	for(int jloc=0;jloc<degx+1;jloc++){
	  real dxi=dlag(degx,iloc,ipg);
	  real dxj=dlag(degx,jloc,ipg);
	  aloc[iloc][jloc]+=dxi*dxj*omega/dx;
	}
      }
    }
    
    // assembly of the matrix
    for(int ie=0;ie<nelx;ie++){
      for(int iloc=0;iloc<degx+1;iloc++){
	for(int jloc=0;jloc<degx+1;jloc++){
	  int ino=iloc + ie * degx;
	  int jno=jloc + ie * degx;
	  real val = aloc[iloc][jloc];
	  AddLinearSolver(&sky,ino,jno,val);
	}
      }
    }
    
    // dirichlet boundary condition at the first and last location
    if(type_bc == _Dirichlet_Poisson_BC){
      AddLinearSolver(&sky,0,0,1e20);
      AddLinearSolver(&sky,neq-1,neq-1,1e20);
    }

    sky.mat_is_assembly=true;
  }

  //DisplayLinearSolver(&sky);
  //sleep(1000);
  
  // source assembly 
  real source[neq];
  for(int i=0;i<neq;i++){
    source[i]=0;
  }

  for(int ie=0;ie<nelx;ie++){
    for(int iloc=0;iloc<degx+1;iloc++){
      real omega=wglop(degx,iloc);
      int ino=iloc + ie * degx;  
      int imem=f->varindex(f->deg,f->raf,f->model.m,iloc+ie*(degx+1),_INDEX_RHO);
      real charge=w[imem];          
      source[ino]+= (charge-charge_average)*omega*dx;
    }
  }

  // Apply dirichlet Boundary condition
  
  source[0]=1e20*bc_l;
  source[neq-1]=1e20*bc_r;
  
  real solution[neq];
  sky.rhs=source;
  sky.sol=solution;
  SolveLinearSolver(&sky);


  // now put the solution at the right place
  for(int ie=0;ie<nelx;ie++){
    for(int ipg=0;ipg<degx+1;ipg++){
      // position in the continuous vector
      int ino=ipg + ie * degx;
      // position in the DG vector
      int imem=f->varindex(f->deg,f->raf,f->model.m,ipg+ie*(degx+1),_INDEX_PHI);
      w[imem]=solution[ino];
    }
  }
	
  FreeLinearSolver(&sky);

  //ComputeElectricField(f);
  Compute_electric_field(f,w);

}


