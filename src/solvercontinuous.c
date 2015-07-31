#include "solvercontinuous.h"
#include "geometry.h"
#include "quantities_vp.h"
#include "linear_solver.h"


int CompareFatNode(const void* a,const void* b){

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

int BuildFatNodeList(Simulation *simu,FatNode* fn_list){


  int big_int = 1 << 28; // 2**28 = 268 435 456

  int ino=0;
  real* xmin=simu->macromesh.xmin;
  real* xmax=simu->macromesh.xmax;
  int nb_dg_nodes;
  for(int ie = 0; ie < simu->macromesh.nbelems; ie++) {
    
    field *f = &simu->fd[ie];

    nb_dg_nodes =  NPG(f->deg, f->raf) * simu->macromesh.nbelems;

    
    for(int ipg = 0; ipg < NPG(f->deg, f->raf); ipg++) {
      real xpg[3];
      real xref[3];
      ref_pg_vol(f->deg, f->raf, ipg, xref, NULL, NULL);
      Ref2Phy(f->physnode,
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

void InitContinuousSolver(void * cs, Simulation* simu,int type_bc,int nb_phy_vars, int * listvar){

  ContinuousSolver * ps=cs;

  ps->simu = simu;

  ps->type_bc=type_bc;
  ps->nb_phy_vars=nb_phy_vars;
  ps->list_of_var=listvar;

  ps->postcomputation_assembly=NULL;
  ps->matrix_assembly=NULL;
  ps->rhs_assembly=NULL;
 
  field *f0 = &simu->fd[0];

  ps->nb_dg_nodes = NPG(f0->deg, f0->raf)
    * simu->macromesh.nbelems;

  ps->nb_dg_dof= ps->nb_dg_nodes * ps->nb_phy_vars;

  
  ps->fn_list = malloc(ps->nb_dg_nodes * sizeof(FatNode));
  assert(ps->fn_list);
  // paste the nodes of the DG mesh
  ps->nb_fe_nodes=BuildFatNodeList(simu,ps->fn_list);
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

  int nraf[3] = {f0->raf[0],f0->raf[1],f0->raf[2]};
  
  int deg[3] = {f0->deg[0],f0->deg[1],f0->deg[2]};

  int npg[3] = {f0->npg[0],f0->npg[1],f0->npg[2]};

  ps->nnodes = npg[0] * npg[1] * npg[2];
  ps->npgmacrocell = ps->nnodes *  nraf[0] * nraf[1] * nraf[2];

  
  ps->nbel = simu->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2];

  for (int ie = 0; ie < simu->macromesh.nbelems; ie++) {
    int nbfa = 6;
    if (simu->macromesh.is2d)  nbfa = 4; 
    if (simu->macromesh.is1d) nbfa = 2;
    for(int ii = 0; ii < nbfa; ii++) {
      int ifa = ii;
      if (simu->macromesh.is1d) ifa = 2 * ii + 1;
      int ieR = simu->macromesh.elem2elem[6*ie+ifa];
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


void SolveContinuous2D(void* cs){

  ContinuousSolver * ps=cs;
  
  field* f0 = &ps->simu->fd[0];

  printf("Init...\n");
  
  int nraf[3] = {f0->raf[0],f0->raf[1],f0->raf[2]};
  int deg[3] = {f0->deg[0],f0->deg[1],f0->deg[2]};

  
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


void ContinuousToDiscontinuous_Copy(ContinuousSolver * cs,LinearSolver* lsol){
  
  field* f0 = &cs->simu->fd[0];

  printf("Copy...\n");

  
  // copy the potential at the right place
  for(int var =0; var < cs->nb_phy_vars; var++){ 
    for(int ie = 0; ie < cs->simu->macromesh.nbelems; ie++){  
      for(int ipg = 0;ipg < cs->npgmacrocell; ipg++){
	int ino_dg = ipg + ie * cs->npgmacrocell;
	int ino_fe = cs->dg_to_fe_index[ino_dg];
	int ipot = f0->varindex(f0->deg,f0->raf,f0->model.m,
				ipg,cs->list_of_var[var]);
	cs->simu->fd[ie].wn[ipot]=cs->lsol.sol[ino_fe];
	}
    }
  }

}


void ExactDirichletContinuousMatrix(void * cs,LinearSolver* lsol){
  ContinuousSolver * ps=cs;
  
  field* f0 = &ps->simu->fd[0];

  for(int ino=0; ino<ps->nb_fe_dof; ino++){
    real bigval = 1e20;
    if (ps->is_boundary_node[ino]) AddLinearSolver(&ps->lsol,ino,ino,bigval);
  }
  ps->lsol.mat_is_assembly=true;
  
  for(int var =0; var < ps->nb_phy_vars; var++){ 
    for(int ie = 0; ie < ps->simu->macromesh.nbelems; ie++){  
      for(int ipg = 0;ipg < ps->npgmacrocell; ipg++){
	int ino_dg = ipg + ie * ps->npgmacrocell;
	int ino_fe = ps->dg_to_fe_index[ino_dg];
	int ipot = f0->varindex(f0->deg,f0->raf,f0->model.m,
				ipg,ps->list_of_var[var]);
	if (ps->is_boundary_node[ino_fe]){
	  real bigval = 1e20;
	  ps->lsol.rhs[ino_fe] = ps->simu->fd[ie].wn[ipot] * bigval;
	  //printf("ino_dg=%d ino_fe=%d ipot=%d\n",ino_dg,ino_fe,ipot);
	}
      }
    }
  }

  
}

void AllocateContinuousMatrix(void *cs,LinearSolver* lsol){
    ContinuousSolver * ps=cs;

  //static bool is_lu = false;

  field* f0 = &ps->simu->fd[0];
  printf("Init...\n");

  // number of equation of the Poisson solver
  int neq=ps->nb_fe_dof; //(nodes)
  
  // number of conservatives variables
  // = number of velocity glops + 1 (potential)
  int m = ps->nb_phy_vars;

  int nraf[3] = {f0->raf[0],f0->raf[1],f0->raf[2]};
  int deg[3] = {f0->deg[0],f0->deg[1],f0->deg[2]};
  
  if (ps->simu->macromesh.is2d) {
    assert( nraf[0] == nraf[1]);
    assert( nraf[2] == 1);
    assert(deg[2] == 0);
  }
  
  if (ps->simu->macromesh.is1d) {
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


void GenericOperatorScalar_Continuous(void * cs,LinearSolver* lsol){

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

      for(int ipg = 0;ipg < ps->nnodes; ipg++){
	real wpg;
	real xref[3];
	int ipgmacro= ipg + isubcell * ps->nnodes;

	ref_pg_vol(f0->deg,f0->raf,ipgmacro,xref,&wpg,NULL);

	for(int iloc = 0; iloc < ps->nnodes; iloc++){
	  real dtau[3][3],codtau[3][3];
	  real dphiref_i[3],dphiref_j[3];
	  real dphi_i[3],dphi_j[3];
	  real basisPhi_i[4], basisPhi_j[4];
	  int ilocmacro = iloc + isubcell * ps->nnodes;
	  grad_psi_pg(f0->deg,f0->raf,ilocmacro,ipgmacro,dphiref_i);
	  Ref2Phy(ps->simu->fd[iemacro].physnode,
		  xref,dphiref_i,0,NULL,
		  dtau,codtau,dphi_i,NULL);

	  real det = dot_product(dtau[0], codtau[0]);
	  if (ilocmacro==ipgmacro){
	    basisPhi_i[0]=1;
    }
	  else
    {
      basisPhi_i[0]=0;
    }
    basisPhi_i[1]=dphi_i[0]/det;
    basisPhi_i[2]=dphi_i[1]/det;
    basisPhi_i[3]=dphi_i[2]/det;
	  for(int jloc = 0; jloc < ps->nnodes; jloc++){
	    int jlocmacro = jloc + isubcell * ps->nnodes;
	    grad_psi_pg(f0->deg,f0->raf,jlocmacro,ipgmacro,dphiref_j);
	    Ref2Phy(ps->simu->fd[iemacro].physnode,
		    xref,dphiref_j,0,NULL,
		    dtau,codtau,dphi_j,NULL);
	    if (jlocmacro==ipgmacro){
	      basisPhi_j[0]=1;
      }
	    else
      {
        basisPhi_j[0]=0;
      }
      basisPhi_j[1]=dphi_j[0]/det;
      basisPhi_j[2]=dphi_j[1]/det;
      basisPhi_j[3]=dphi_j[2]/det;
      real res[4] = {0, 0, 0, 0};
      for (int i=0; i<4; i++){
        for (int j=0; j<4; j++){
          res[i]+=basisPhi_j[j]*ps->diff_op[i][j];
        }
      }
	    aloc[iloc][jloc] += dot_product(basisPhi_i, res) * wpg * det  ;
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


void GenericOperator2Vec_Continuous(void * cs,LinearSolver* lsol){

  ContinuousSolver * ps=cs;

   field* f0 = &ps->simu->fd[0];
   if (ps->nb_phy_vars>2){
     printf("Not implemented for three variables.\n");
     exit;
   }

  if(!ps->lsol.mat_is_assembly){
    for(int ie = 0; ie < ps->nbel; ie++){  

      // local matrix 
      real aloc[ps->nnodes*ps->nb_phy_vars][ps->nnodes*ps->nb_phy_vars];
      for(int iloc = 0; iloc < ps->nnodes*ps->nb_phy_vars; iloc++){
	for(int jloc = 0; jloc < ps->nnodes*ps->nb_phy_vars; jloc++){
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
	  real basisPhi_i[4], basisPhi_j[4];
	  int ilocmacro = iloc + isubcell * ps->nnodes;
	  grad_psi_pg(f0->deg,f0->raf,ilocmacro,ipgmacro,dphiref_i);
	  Ref2Phy(ps->simu->fd[iemacro].physnode,
		  xref,dphiref_i,0,NULL,
		  dtau,codtau,dphi_i,NULL);

	  real det = dot_product(dtau[0], codtau[0]);
	  if (ilocmacro==ipgmacro){
	    basisPhi_i[0]=1;
    }
	  else
    {
      basisPhi_i[0]=0;
    }
    basisPhi_i[1]=dphi_i[0]/det;
    basisPhi_i[2]=dphi_i[1]/det;
    basisPhi_i[3]=dphi_i[2]/det;
	  for(int jloc = 0; jloc < ps->nnodes; jloc++){
	    int jlocmacro = jloc + isubcell * ps->nnodes;
	    grad_psi_pg(f0->deg,f0->raf,jlocmacro,ipgmacro,dphiref_j);
	    Ref2Phy(ps->simu->fd[iemacro].physnode,
		    xref,dphiref_j,0,NULL,
		    dtau,codtau,dphi_j,NULL);
	    if (jlocmacro==ipgmacro){
	      basisPhi_j[0]=1;
      }
	    else
      {
        basisPhi_j[0]=0;
      }
      basisPhi_j[1]=dphi_j[0]/det;
      basisPhi_j[2]=dphi_j[1]/det;
      basisPhi_j[3]=dphi_j[2]/det;
      for (int iv1=0; iv1<ps->nb_phy_vars; iv1++){
        for (int iv2=0; iv2<ps->nb_phy_vars; iv2++){
          real res[4] = {0, 0, 0, 0};
          for (int i=0; i<4; i++){
            for (int j=0; j<4; j++){
              res[i]+=basisPhi_j[j]*ps->diff_op2vec[2*iv1+iv2][i][j];
            }
          }
	        aloc[iv1+iloc*ps->nb_phy_vars][iv1+jloc*ps->nb_phy_vars] += dot_product(basisPhi_i, res) * wpg * det  ;
        }
	    }
    }
	}
      }


      for(int iloc = 0; iloc < ps->nnodes; iloc++){
	for(int jloc = 0; jloc < ps->nnodes; jloc++){
	  int ino_dg = iloc + ie * ps->nnodes;
	  int jno_dg = jloc + ie * ps->nnodes;
	  int ino_fe = ps->dg_to_fe_index[ino_dg];
	  int jno_fe = ps->dg_to_fe_index[jno_dg];
	  for (int iv=0;iv<ps->nb_phy_vars;iv++){
	    real val = aloc[iloc*ps->nb_phy_vars+iv][jloc*ps->nb_phy_vars+iv];
	    AddLinearSolver(&ps->lsol,ino_fe*ps->nb_phy_vars+iv,jno_fe*ps->nb_phy_vars+iv,val);
    }
	}
      }
   
    }
  } 
  
 }


void catCGVectors(void * cs,real *L1, real *L2, real *L){

  ContinuousSolver * ps = cs;
  int cc=0;
  for (int i=0; i<ps->nb_fe_nodes;i++)
  {
    L[cc] = L1[i];
    L[cc+1]=L2[i];
    cc+=2;
  }
}
