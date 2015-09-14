#include "solvercontinuous.h"
#include "geometry.h"
#include "quantities_vp.h"
#include "linear_solver.h"

int CompareFatNode(const void* a, const void* b)
{
  FatNode* fna = (FatNode*) a;
  FatNode* fnb = (FatNode*) b;

  int r = fna->x_int[0] - fnb->x_int[0];
  if (r==0 || r==-1 || r==1)
    r = fna->x_int[1] - fnb->x_int[1];
  if (r==0 || r==-1 || r==1)
    r = fna->x_int[2] - fnb->x_int[2];
  if (r==0 || r==-1 || r==1)
    r=0;
  return r;
}

int BuildFatNodeList(field *f, FatNode* fn_list)
{
  int big_int = 1 << 28; // 2**28 = 268 435 456

  int ino=0;
  real* xmin = f->macromesh.xmin;
  real* xmax = f->macromesh.xmax;
  int nb_dg_nodes;
  for(int ie = 0; ie < f->macromesh.nbelems; ie++) {
    
    nb_dg_nodes =  NPG(f->deg, f->raf) * f->macromesh.nbelems;

    for(int ipg = 0; ipg < NPG(f->deg, f->raf); ipg++) {
      real xpg[3];
      real xref[3];
      ref_pg_vol(f->deg, f->raf, ipg, xref, NULL, NULL);

      real physnode[20][3];
      for(int inoloc = 0; inoloc < 20; inoloc++) {
	int ino = f->macromesh.elem2node[20*ie+inoloc];
	physnode[inoloc][0] = f->macromesh.node[3 * ino + 0];
	physnode[inoloc][1] = f->macromesh.node[3 * ino + 1];
	physnode[inoloc][2] = f->macromesh.node[3 * ino + 2];
      }

      Ref2Phy(physnode,
	      xref,
	      0, -1, // dphiref, ifa
              xpg, NULL,
	      NULL, NULL, NULL); // codtau, dphi, vnds

      fn_list[ino].x[0]=xpg[0];
      fn_list[ino].x[1]=xpg[1];
      fn_list[ino].x[2]=xpg[2];

      // convert points to "pixels" and round to nearest integer

      fn_list[ino].x_int[0] = (int) ((xpg[0] - xmin[0]) / (xmax[0] - xmin[0])
				     * big_int);
      fn_list[ino].x_int[1] =(int) ((xpg[1] - xmin[1]) / (xmax[1] - xmin[1])
				    * big_int);
      fn_list[ino].x_int[2]=(int) ((xpg[2] - xmin[2]) / (xmax[2] - xmin[2])
				   * big_int);
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

void InitContinuousSolver(void *cs, field *f, int type_bc, int nb_phy_vars,
			  int *listvar)
{
  ContinuousSolver *ps = cs;

  ps->f = f;

  ps->type_bc = type_bc;
  ps->nb_phy_vars = nb_phy_vars;
  ps->list_of_var = listvar;

  ps->postcomputation_assembly = NULL;
  ps->matrix_assembly = NULL;
  ps->rhs_assembly = NULL;
 
  ps->nb_dg_nodes = NPG(f->deg, f->raf) * f->macromesh.nbelems;

  ps->nb_dg_dof= ps->nb_dg_nodes * ps->nb_phy_vars;
  
  ps->fn_list = malloc(ps->nb_dg_nodes * sizeof(FatNode));
  assert(ps->fn_list);
  // paste the nodes of the DG mesh
  ps->nb_fe_nodes=BuildFatNodeList(f, ps->fn_list);
  ps->nb_fe_dof= ps->nb_fe_nodes * ps->nb_phy_vars;

  printf("nb dg nodes=%d ; nb fe nodes=%d\n",ps->nb_dg_nodes,ps->nb_fe_nodes);
  // build the connectivity array

  ps->dg_to_fe_index = malloc(ps->nb_dg_nodes * sizeof(int));
  assert(ps->dg_to_fe_index);

  for(int ino=0;ino<ps->nb_dg_nodes;ino++){
    /* printf("ino=%d idg=%d ife=%d\n",ino,ps->fn_list[ino].dg_index, */
    /* 	   ps->fn_list[ino].fe_index); */
    ps->dg_to_fe_index[ps->fn_list[ino].dg_index]=ps->fn_list[ino].fe_index;
    
  }

  // now construct the list of boundary nodes
  ps->is_boundary_node = malloc(ps->nb_fe_nodes * sizeof(int));
  assert(ps->is_boundary_node);
  for(int ino = 0; ino < ps->nb_fe_nodes; ino++){
    ps->is_boundary_node[ino] = 0;
  }

  int nraf[3] = {f->raf[0],f->raf[1],f->raf[2]};
  
  int deg[3] = {f->deg[0],f->deg[1],f->deg[2]};

  int npg[3] = {deg[0] + 1, deg[1] + 1, deg[2] + 1};

  ps->nnodes = npg[0] * npg[1] * npg[2];
  ps->npgmacrocell = ps->nnodes *  nraf[0] * nraf[1] * nraf[2];

  
  ps->nbel = f->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2];

  for (int ie = 0; ie < f->macromesh.nbelems; ie++) {
    int nbfa = 6;
    if (f->macromesh.is2d)  nbfa = 4; 
    if (f->macromesh.is1d) nbfa = 2;
    for(int ii = 0; ii < nbfa; ii++) {
      int ifa = ii;
      if (f->macromesh.is1d) ifa = 2 * ii + 1;
      int ieR = f->macromesh.elem2elem[6*ie+ifa];
      if (ieR < 0) {
	for(int ipgf = 0; ipgf < NPGF(deg,nraf, ifa); ipgf++) {
	  printf("NPGF=%d ipgf=%d\n",NPGF(deg,nraf, ifa),ipgf);
	  int ipg = ref_pg_face(deg,nraf, ifa, ipgf,
				NULL, NULL, NULL);
	  int ino_dg = ipg + ie * ps->npgmacrocell;
	  int ino_fe = ps->dg_to_fe_index[ino_dg];
	  ps->is_boundary_node[ino_fe] = 1;
	  printf("ie=%d ino_dg=%d ino_fe=%d boundary=%d\n",
	  	 ie,ino_dg,ino_fe,ps->is_boundary_node[ino_fe]);

	}
      }
    }
  }

  int count_boundary = 0;

  for(int ino = 0; ino < ps->nb_fe_nodes; ino++){
    count_boundary += ps->is_boundary_node[ino];
  }

  printf("found %d boundary nodes (on %d fe nodes)\n",
	 count_boundary, ps->nb_fe_nodes);
	
  InitLinearSolver(&ps->lsol,ps->nb_fe_dof,NULL,NULL); //InitSkyline(&sky, neq);

  ps->lsol.rhs = malloc(ps->nb_fe_dof * sizeof(real));
  assert(ps->lsol.rhs);

  ps->lsol.sol = malloc(ps->nb_fe_dof * sizeof(real));
  assert(ps->lsol.sol);

  AllocateContinuousMatrix(ps,&ps->lsol);
}

void SolveContinuous2D(void *cs)
{
  ContinuousSolver *ps=cs;
  
  field *f = ps->f;

  printf("Init...\n");
  
  int nraf[3] = {f->raf[0],f->raf[1],f->raf[2]};
  int deg[3] = {f->deg[0],f->deg[1],f->deg[2]};

  
  // number of equation of the Poisson solver
  int neq=ps->nb_fe_dof; //(nodes)

  printf("Matrix assembly.....\n");
  if(ps->matrix_assembly != NULL){
    ps->matrix_assembly(ps,&ps->lsol);
  }
  else
    {
      printf("no matrix assembly ");
      exit(1);
    }
  
  printf("RHS assembly.....\n");
  if(ps->rhs_assembly != NULL){
    ps->rhs_assembly(ps,&ps->lsol);
  }
  else
    {
      printf("no rhs assembly ");
      exit(1);
    }

  printf("BC assembly.....\n");
  if(ps->bc_assembly != NULL){
    ps->bc_assembly(ps,&ps->lsol);
  }
  else
    {
      printf("no bc assembly ");
      exit(1);
    }
  
  printf("Solution...\n");

  SolveLinearSolver(&ps->lsol);

  
  printf("post computation assembly.....\n");
  if(ps->postcomputation_assembly != NULL){
    ps->postcomputation_assembly(ps,&ps->lsol);
  }
  printf("End SolvePoisson2D.\n");
}

void ContinuousToDiscontinuous_Copy(ContinuousSolver * cs,LinearSolver* lsol)
{
  field *f = cs->f;

  printf("Copy...\n");
  
  // copy the potential at the right place
  for(int var =0; var < cs->nb_phy_vars; var++){ 
    for(int ie = 0; ie < cs->f->macromesh.nbelems; ie++){  
      for(int ipg = 0;ipg < cs->npgmacrocell; ipg++){
	int ino_dg = ipg + ie * cs->npgmacrocell;
	int ino_fe = cs->dg_to_fe_index[ino_dg];
	int ipot = f->varindex(f->deg, f->raf, f->model.m,
			       ie, ipg, cs->list_of_var[var]);

	// FIXME!!!!!!!!
	//cs->f->wn[ipot]=cs->lsol.sol[ino_fe];
      }
    }
  }
}

void ExactDirichletContinuousMatrix(void * cs, LinearSolver* lsol)
{
  ContinuousSolver * ps=cs;
  
  field *f = ps->f;

  for(int ino=0; ino<ps->nb_fe_dof; ino++){
    real bigval = 1e20;
    if (ps->is_boundary_node[ino]) AddLinearSolver(&ps->lsol,ino,ino,bigval);
  }
  ps->lsol.mat_is_assembly=true;
  
  for(int var =0; var < ps->nb_phy_vars; var++){ 
    for(int ie = 0; ie < ps->f->macromesh.nbelems; ie++){  
      for(int ipg = 0;ipg < ps->npgmacrocell; ipg++){
	int ino_dg = ipg + ie * ps->npgmacrocell;
	int ino_fe = ps->dg_to_fe_index[ino_dg];
	int ipot = f->varindex(f->deg, f->raf, f->model.m,
			       ie, ipg, ps->list_of_var[var]);
	if (ps->is_boundary_node[ino_fe]){
	  real bigval = 1e20;

	  // FIXME!!!!!!
	  //ps->lsol.rhs[ino_fe] = ps->f->fd[ie].wn[ipot] * bigval;


	  //printf("ino_dg=%d ino_fe=%d ipot=%d\n",ino_dg,ino_fe,ipot);
	}
      }
    }
  }

}

void AllocateContinuousMatrix(void *cs,LinearSolver* lsol)
{
  ContinuousSolver * ps=cs;

  //static bool is_lu = false;

  field *f = ps->f;
  printf("Init...\n");

  // number of equation of the Poisson solver
  int neq=ps->nb_fe_dof; //(nodes)
  
  // number of conservatives variables
  // = number of velocity glops + 1 (potential)
  int m = ps->nb_phy_vars;

  int nraf[3] = {f->raf[0],f->raf[1],f->raf[2]};
  int deg[3] = {f->deg[0],f->deg[1],f->deg[2]};
  
  if (ps->f->macromesh.is2d) {
    assert( nraf[0] == nraf[1]);
    assert( nraf[2] == 1);
    assert(deg[2] == 0);
  }
  
  if (ps->f->macromesh.is1d) {
    assert( nraf[1] == 1);
    assert(deg[1] == 0);
    assert( nraf[2] == 1);
    assert(deg[2] == 0);
  }
 
  printf("Allocation...\n");
  if(!ps->lsol.is_alloc){

    // compute the profile of the matrix
    for(int ie = 0; ie < ps->nbel; ie++){
      for(int iloc = 0; iloc < ps->nnodes; iloc++){
	for(int jloc = 0; jloc < ps->nnodes; jloc++){
	  int ino_dg = iloc + ie * ps->nnodes;
	  int jno_dg = jloc + ie * ps->nnodes;
	  for(int var = 0; var < ps->nb_phy_vars; var++){
	    int ino_fe = ps->dg_to_fe_index[ino_dg]+var;
	    int jno_fe = ps->dg_to_fe_index[jno_dg]+var;
	    IsNonZero(&ps->lsol, ino_fe, jno_fe);
	  }
	}
      }
    }
    AllocateLinearSolver(&ps->lsol);
  } 
}


