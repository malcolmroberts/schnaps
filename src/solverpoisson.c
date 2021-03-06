#include "solverpoisson.h"
#include "geometry.h"
#include "quantities_vp.h"
#include "linear_solver.h"

int CompareFatNode(const void* a,const void* b)
{
  FatNode* fna = (FatNode*) a;
  FatNode* fnb = (FatNode*) b;

  int r = fna->x_int[0]-fnb->x_int[0];
  if (r==0 || r==-1 || r==1)
    r = fna->x_int[1]-fnb->x_int[1];
  if (r==0 || r==-1 || r==1)
    r = fna->x_int[2]-fnb->x_int[2];
  if (r==0 || r==-1 || r==1)
    r=0;
  return r;

}

int BuildFatNodeList(field* f, FatNode* fn_list)
{
  int big_int = 1 << 28; // 2**28 = 268 435 456

  int *raf = f->interp_param + 4;
  int *deg = f->interp_param + 1;
  
  int nb_dg_nodes = NPG(raf, deg) * f->macromesh.nbelems;

  int ino=0;
  real* xmin=f->macromesh.xmin;
  real* xmax=f->macromesh.xmax;
  for(int ie = 0; ie < f->macromesh.nbelems; ie++) {
    
    real physnode[20][3];
    for(int inoloc = 0; inoloc < 20; inoloc++) {
      int ino = f->macromesh.elem2node[20 * ie + inoloc];
      physnode[inoloc][0] = f->macromesh.node[3 * ino + 0];
      physnode[inoloc][1] = f->macromesh.node[3 * ino + 1];
      physnode[inoloc][2] = f->macromesh.node[3 * ino + 2];
    }
    
    for(int ipg = 0; ipg < NPG(raf, deg); ipg++) {
      real xpg[3];
      real xref[3];
      int *raf = f->interp_param + 4;
      int *deg = f->interp_param + 1;
      ref_pg_vol(raf, deg, ipg, xref, NULL, NULL);
      Ref2Phy(physnode,
	      xref,
	      0, -1, // dphiref, ifa
              xpg, NULL,
	      NULL, NULL, NULL); // codtau, dphi, vnds

      fn_list[ino].x[0]=xpg[0];
      fn_list[ino].x[1]=xpg[1];
      fn_list[ino].x[2]=xpg[2];

      // convert points to "pixels" and round to nearest integer

      fn_list[ino].x_int[0]=(int) ((xpg[0]-xmin[0])/(xmax[0]-xmin[0]) * big_int)  ;
      fn_list[ino].x_int[1]=(int) ((xpg[1]-xmin[1])/(xmax[1]-xmin[1]) * big_int)  ;
      fn_list[ino].x_int[2]=(int) ((xpg[2]-xmin[2])/(xmax[2]-xmin[2]) * big_int)  ;
      
      fn_list[ino].dg_index = ino;

      ino++;
    }
  }

  assert(ino == nb_dg_nodes);

  qsort(fn_list, nb_dg_nodes, sizeof(FatNode),CompareFatNode);

  fn_list[0].fe_index=0;
  int fe_index=0;
  for(int ino=1;ino<nb_dg_nodes;ino++){
    if (CompareFatNode(fn_list+ino,fn_list+ino-1)!=0) fe_index++;
    fn_list[ino].fe_index=fe_index;
  }
  

  /* for(int ino=0;ino<nb_dg_nodes;ino++){ */
  /*   printf("ino=%d xyz= %f %f %f i_xyz=%d %d %d dg_index=%d fe_index=%d\n",ino, */
  /* 	   fn_list[ino].x[0],fn_list[ino].x[1],fn_list[ino].x[2], */
  /* 	   fn_list[ino].x_int[0],fn_list[ino].x_int[1],fn_list[ino].x_int[2], */
  /* 	   fn_list[ino].dg_index, */
  /* 	   fn_list[ino].fe_index); */
  /* } */

  

  return fe_index+1;

}

void InitPoissonSolver(PoissonSolver* ps, field* fd,int charge_index){

  ps->fd = fd;
  ps->charge_index = charge_index;

  int *raf = fd->interp_param + 4;
  int *deg = fd->interp_param + 1;
  
  ps->nb_dg_nodes = NPG(raf, deg) * fd->macromesh.nbelems;

  
  ps->fn_list = malloc(ps->nb_dg_nodes * sizeof(FatNode));
  assert(ps->fn_list);
  // paste the nodes of the DG mesh
  ps->nb_fe_nodes=BuildFatNodeList(fd,ps->fn_list);


  printf("nb dg nodes=%d ; nb fe nodes=%d\n",ps->nb_dg_nodes,ps->nb_fe_nodes);
  // build the connectivity array

  ps->dg_to_fe_index = malloc(ps->nb_dg_nodes * sizeof(int));
  assert(ps->dg_to_fe_index);

  for(int ino=0;ino<ps->nb_dg_nodes;ino++){
    /* printf("ino=%d idg=%d ife=%d\n",ino,ps->fn_list[ino].dg_index, */
    /* 	   ps->fn_list[ino].fe_index); */
    ps->dg_to_fe_index[ps->fn_list[ino].dg_index]=ps->fn_list[ino].fe_index;
    
  }
  /* for(int ino=0;ino<ps->nb_dg_nodes;ino++){ */
  /*   printf("idg=%d ife=%d\n",ino,ps->dg_to_fe_index[ino]); */
  /* } */


  // now construct the list of boundary nodes
  ps->is_boundary_node = malloc(ps->nb_fe_nodes * sizeof(int));
  assert(ps->is_boundary_node);
  for(int ino = 0; ino < ps->nb_fe_nodes; ino++) {
    ps->is_boundary_node[ino] = 0;
  }

  int nraf[3] = {ps->fd->interp_param[4], 
		 ps->fd->interp_param[5],
		 ps->fd->interp_param[6]};
  
  int npg[3] = {ps->fd->interp_param[1] + 1, 
		ps->fd->interp_param[2] + 1,
		ps->fd->interp_param[3] + 1};

  int npgcell = npg[0] * npg[1] * npg[2];
  int npgmacrocell = npgcell *  nraf[0] * nraf[1] * nraf[2];

  
  int nbel = ps->fd->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2];
  //printf("nraf=%d %d %d \n",nraf[0],nraf[1],nraf[2]);

  for (int ie = 0; ie < ps->fd->macromesh.nbelems; ie++) {
    int nbfa = 6;
    if (ps->fd->macromesh.is2d) nbfa = 4;
    for(int ifa = 0; ifa < nbfa; ifa++) {
      int ieR = ps->fd->macromesh.elem2elem[6*ie+ifa];
      if (ieR < 0) {
	int *raf = ps->fd->interp_param + 4;
	int *deg = ps->fd->interp_param + 1;
	for(int ipgf = 0; ipgf < NPGF(raf, deg, ifa); ipgf++) {
	  int* deg = ps->fd->interp_param + 1;
	  int* raf = ps->fd->interp_param + 4;
	  int ipg = ref_pg_face(raf, deg, ifa, ipgf, NULL, NULL, NULL);
	  int ino_dg = ipg + ie * npgmacrocell;
	  int ino_fe = ps->dg_to_fe_index[ino_dg];
	  /* printf("ie=%d ino_dg=%d ino_fe=%d boundary=%d\n", */
	  /* 	 ie,ino_dg,ino_fe,ps->is_boundary_node[ino_fe]); */
	  ps->is_boundary_node[ino_fe] = 1;

	}
      }
    }
  }

	
  InitLinearSolver(&ps->lsol,ps->nb_fe_nodes,NULL,NULL); //InitSkyline(&sky, neq);

  ps->lsol.rhs = malloc(ps->nb_fe_nodes * sizeof(real));
  assert(ps->lsol.rhs);

  ps->lsol.sol = malloc(ps->nb_fe_nodes * sizeof(real));
  assert(ps->lsol.sol);
  

}

void SolvePoisson1D(field *f,real * w,int type_bc, real bc_l, real bc_r,Solver solver_sys,PC precon){

  real charge_average;
  real *bounds = malloc(6 * sizeof(real));
  charge_average=0;

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
  assert(f->macromesh.is1d);

  // assembly of the rigidity matrix

  LinearSolver sky;

  // number of equation of the Poisson solver
  // = number of nodes in the mesh
  int degx=f->interp.interp_param[1];
  int nelx=f->interp.interp_param[4];
  real xmin=f->macromesh.xmin[0];
  real xmax=f->macromesh.xmax[0];  //
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

  if(!sky.is_assembly){
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

    sky.is_assembly=true;
  }

  //DisplayLinearSolver(&sky);
  //sleep(1000);
  
  // source assembly 
  real source[neq];
  for(int i=0;i<neq;i++){
    source[i]=0;
  }

  for(int ie = 0; ie < nelx; ie++) {
    MacroCell *mcell = f->mcell + ie;
    for(int iloc=0; iloc < degx + 1; iloc++) {
      real omega = wglop(degx, iloc);
      int ino = iloc + ie * degx;

      // FIXME: this used to have ie=0, how can that be correct?
      int imem = f->varindex(f->interp_param, iloc + ie * (degx + 1),
			     _INDEX_RHO);
      real charge = w[imem];          
      source[ino] += (charge - charge_average) * omega * dx;
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
  for(int ie = 0; ie < nelx; ie++) {
    for(int ipg=0;ipg<degx+1;ipg++) {
      // position in the continuous vector
      int ino = ipg + ie * degx;
      // position in the DG vector
      int imem = f->varindex(f->interp_param, ipg + ie * (degx + 1),_INDEX_PHI);
      w[imem]=solution[ino];
    }
  }
	
  FreeLinearSolver(&sky);

  //ComputeElectricField(f);
  Compute_electric_field(f,w);

}

void SolvePoisson2D(PoissonSolver* ps, int type_bc)
{
  real charge_average = 0.0;

  field* f=ps->fd;

  // value of the dirichlet bc;
  real bc=0;

  if(type_bc == _Periodic_Poisson_BC){
    //charge_average=Computation_charge_average(f,w);
  }
  else {
    charge_average=0;
  }
  
  printf("Init...\n");

  // number of equation of the Poisson solver
  int neq=ps->nb_fe_nodes;
  
  // number of conservatives variables
  // = number of velocity glops + 1 (potential)
  int m = f->model.m;

  int nraf[3] = {ps->fd->interp_param[4], 
		 ps->fd->interp_param[5],
		 ps->fd->interp_param[6]};

  real delta_x = 1. / nraf[0];
  real dv = pow(delta_x,2.);
  assert( nraf[0] == nraf[1]);
  assert( nraf[2] == 1);
  
  int npg[3] = {ps->fd->interp_param[1] + 1, 
		ps->fd->interp_param[2] + 1,
		ps->fd->interp_param[3] + 1};

  int nbel = ps->fd->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2];
  int nnodes = npg[0] * npg[1] * npg[2] ;
  int npgmacrocell = nnodes * nraf[0] * nraf[1] * nraf[2];
 
  printf("Allocation...\n");
  if(!ps->lsol.is_alloc){

    // compute the profile of the matrix
    for(int ie = 0; ie < nbel; ie++){
      for(int iloc = 0; iloc < nnodes; iloc++){
	for(int jloc = 0; jloc < nnodes; jloc++){
	  int ino_dg = iloc + ie * nnodes;
	  int jno_dg = jloc + ie * nnodes;
	  int ino_fe = ps->dg_to_fe_index[ino_dg];
	  int jno_fe = ps->dg_to_fe_index[jno_dg];
	  IsNonZero(&ps->lsol, ino_fe, jno_fe);
	}
      }
    }
    AllocateLinearSolver(&ps->lsol);
  } 
    
  
  printf("Assembly...\n");
  if(!ps->lsol.is_assembly){
    for(int ie = 0; ie < nbel; ie++){  

      // local matrix 
      real aloc[nnodes][nnodes];
      for(int iloc = 0; iloc < nnodes; iloc++){
	for(int jloc = 0; jloc < nnodes; jloc++){
	  aloc[iloc][jloc] = 0;
	}
      }

      int iemacro = ie / (nraf[0] * nraf[1] * nraf[2]);
      int isubcell = ie % (nraf[0] * nraf[1] * nraf[2]);

      real physnode[20][3];
      for(int ino = 0; ino < 20; ino++) {
	int numnoe = ps->fd->macromesh.elem2node[20 * iemacro + ino];
	for(int ii = 0; ii < 3; ii++) {
	  physnode[ino][ii] = ps->fd->macromesh.node[3 * numnoe + ii];
	}
      }
      //ref_pg_vol(ps->fd->interp_param+1,int ipg,
      // real* xpg,real* wpg,real* xpg_in);
      // grad_psi_pg(ps->fd->interp_param+1,ib,ipg,dphiref)
      // Ref2Phy(physnode,xref,dphiref,NULL,NULL,dtau,codtau,dphi,NULL);


      for(int ipg = 0;ipg < nnodes; ipg++){
	real wpg;
	real xref[3];
	int ipgmacro= ipg + isubcell * nnodes;

	int *raf = ps->fd->interp_param+4;
	int *deg = ps->fd->interp_param+1;
	
	ref_pg_vol(raf, deg, ipgmacro, xref, &wpg, NULL);

	for(int iloc = 0; iloc < nnodes; iloc++){
	  real dtau[3][3],codtau[3][3];
	  real dphiref_i[3],dphiref_j[3];
	  real dphi_i[3],dphi_j[3];
	  int ilocmacro = iloc + isubcell * nnodes;
	  grad_psi_pg(ps->fd->interp_param+1,ilocmacro,ipgmacro,dphiref_i);

	  Ref2Phy(physnode,xref,dphiref_i,0,NULL,
		  dtau,codtau,dphi_i,NULL);
	  real det = dot_product(dtau[0], codtau[0]);

	  for(int jloc = 0; jloc < nnodes; jloc++){
	    int jlocmacro = jloc + isubcell * nnodes;
	    grad_psi_pg(ps->fd->interp_param+1,jlocmacro,ipgmacro,dphiref_j);
	    Ref2Phy(physnode,xref,dphiref_j,0,NULL,
		    dtau,codtau,dphi_j,NULL);
	    aloc[iloc][jloc] += dot_product(dphi_i, dphi_j) * wpg / det  ;
	  }
	}
      }

      for(int iloc = 0; iloc < nnodes; iloc++){
	for(int jloc = 0; jloc < nnodes; jloc++){
	  int ino_dg = iloc + ie * nnodes;
	  int jno_dg = jloc + ie * nnodes;
	  int ino_fe = ps->dg_to_fe_index[ino_dg];
	  int jno_fe = ps->dg_to_fe_index[jno_dg];
	  real val = aloc[iloc][jloc];
	  AddLinearSolver(&ps->lsol,ino_fe,jno_fe,val);
	}
      }
   
    }

    for(int ino=0; ino<ps->nb_fe_nodes; ino++){
      real bigval = 1e20;
      if (ps->is_boundary_node[ino]) AddLinearSolver(&ps->lsol,ino,ino,bigval);
    }
    ps->lsol.is_assembly=true;
  } 
  

  printf("RHS assembly.....\n");

  // right hand side assembly
  for(int ino = 0; ino < ps->nb_fe_nodes; ino++){
    ps->lsol.rhs[ino] = 0;
  }

  real surf = 0;

  // FIXME: ie is the subcell, please rename to match rest of code
  // where ie is the macrocell.
  for(int ie = 0; ie < nbel; ie++){  
    int iemacro = ie / (nraf[0] * nraf[1] * nraf[2]);
    int isubcell = ie % (nraf[0] * nraf[1] * nraf[2]);

    MacroCell *mcell = f->mcell + iemacro;
    
    real physnode[20][3];
    for(int ino = 0; ino < 20; ino++) {
      int numnoe = ps->fd->macromesh.elem2node[20 * iemacro + ino];
      for(int ii = 0; ii < 3; ii++) {
	physnode[ino][ii] = ps->fd->macromesh.node[3 * numnoe + ii];
      }
    }
    
    for(int iloc = 0; iloc < nnodes; iloc++){
      real wpg;
      real xref[3];
      //int ipgmacro = ipg + isubcell * nnodes;
      int ilocmacro = iloc + isubcell * nnodes;
      int *raf = ps->fd->interp_param+4;
      int *deg = ps->fd->interp_param+1;
      ref_pg_vol(raf, deg, ilocmacro,xref,&wpg,NULL);
      real dtau[3][3],codtau[3][3];
      Ref2Phy(physnode,xref,NULL,0,NULL,
	      dtau,codtau,NULL,NULL);
      real det = dot_product(dtau[0], codtau[0]);	
      int ino_dg = iloc + ie * nnodes;
      int ino_fe = ps->dg_to_fe_index[ino_dg];
      int imem = ps->fd->varindex(ps->fd->interp_param, ilocmacro, _INDEX_RHO)
	+ mcell->woffset;
      real rho = ps->fd->wn[imem];
      ps->lsol.rhs[ino_fe] += rho  * wpg * det ; // TODO: put the actual charge	
      surf += wpg * det ;
    }
  }

  // apply non homogeneous dirichlet boundary conditions
  /* for(int ino=0; ino<ps->nb_fe_nodes; ino++){ */
  /*   real bigval = 1e20; */
  /*   if (ps->is_boundary_node[ino]) ps->rhs[ino]= 1 * bigval; */
  /* } */

  for(int ie = 0; ie < ps->fd->macromesh.nbelems; ie++){  
    MacroCell *mcell = f->mcell + ie;
    for(int ipg = 0;ipg < npgmacrocell; ipg++){
      int ino_dg = ipg + ie * npgmacrocell;
      int ino_fe = ps->dg_to_fe_index[ino_dg];
      int ipot = ps->fd->varindex(ps->fd->interp_param,
				  ipg,_INDEX_PHI) + mcell->woffset;
      if (ps->is_boundary_node[ino_fe]){
	real bigval = 1e20;
	ps->lsol.rhs[ino_fe] = ps->fd->wn[ipot] * bigval;
	//printf("ino_dg=%d ino_fe=%d ipot=%d\n",ino_dg,ino_fe,ipot);
      }
    }
  }
  printf("Solution...\n");

  SolveLinearSolver(&ps->lsol);

  printf("Copy...\n");

  // copy the potential at the right place 
  for(int ie = 0; ie < ps->fd->macromesh.nbelems; ie++){  
    MacroCell *mcell = f->mcell + ie;
    for(int ipg = 0;ipg < npgmacrocell; ipg++){
      int ino_dg = ipg + ie * npgmacrocell;
      int ino_fe = ps->dg_to_fe_index[ino_dg];
      int ipot = ps->fd->varindex(ps->fd->interp_param,
				  ipg,_INDEX_PHI) + mcell->woffset;
      ps->fd->wn[ipot]=ps->lsol.sol[ino_fe];
    }
  }

  printf("Compute electric field...\n");
  ComputeElectricField(ps->fd);

  printf("End SolvePoisson2D.\n");
}
