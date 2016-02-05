#include "linear_solver.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "skyline.h"
#include "skyline_spu.h"
#include "paralution_c.h"
#include "dpackfgmres.h"
#include "physBased_PC.h"

void InitLinearSolver(LinearSolver* lsol,int n,
		      MatrixStorage* matstor,
		      Solver* solvtyp){

  lsol->neq=n;
  lsol->storage_type = SKYLINE;
  lsol->solver_type = LU;
  lsol->pc_type = NONE;
  lsol->is_sym = false;
  lsol->is_init = false;
  lsol->is_alloc = false;
  lsol->rhs_is_assembly = false;
  lsol->mat_is_assembly = false;
  lsol->rhs=NULL;
  lsol->sol=NULL;
  lsol->MatVecProduct=NULL;
  lsol->tol=1.e-9;
  lsol->restart_gmres=1;
  lsol->iter_max=100;
  lsol->is_CG=false;

  if (matstor != NULL) lsol->storage_type = *matstor;
  if (solvtyp != NULL) lsol->solver_type = *solvtyp;

  if (matstor != NULL) {
    if (*matstor != SKYLINE_SPU){
      lsol->rhs=calloc(n,sizeof(schnaps_real));
      lsol->sol=calloc(n,sizeof(schnaps_real));
    }
  }
    
  switch(lsol->storage_type) {

    Skyline* sky;
    Skyline_SPU* sky_spu;

  case SKYLINE :
    sky = malloc(sizeof(Skyline));
    assert(sky);
    lsol->matrix = (void*) sky;
    InitSkyline(sky,n);
    lsol->is_init = true;
    break;

  case SKYLINE_SPU :
    sky_spu = malloc(sizeof(Skyline_SPU));
    assert(sky_spu);
    lsol->matrix = (void*) sky_spu;
    InitSkyline_SPU(sky_spu,n);
    lsol->rhs = sky_spu->rhs;
    lsol->sol = sky_spu->sol;
    lsol->is_init = true;
    break;

  case CSR :
    assert(1==2);
    break;

  default : 
    assert(1==2);

  }


}

void FreeLinearSolver(LinearSolver* lsol){

  assert(lsol->is_alloc);
  assert(lsol->rhs);
  assert(lsol->sol);
  // Free rhs
  if (lsol->rhs != NULL)
  {
    free(lsol->rhs);
  }
  // Free sol
  if (lsol->sol != NULL)
  {
    free(lsol->sol);
  }
  
  switch(lsol->storage_type) {

  case SKYLINE :
    FreeSkyline((Skyline*)lsol->matrix);
    // Free matrix
    if (lsol->matrix != NULL)
    {
      free(lsol->matrix);
    }
    break;

case SKYLINE_SPU :
    FreeSkyline_SPU((Skyline_SPU*)lsol->matrix);
    free(lsol->matrix);
    break;

  default : 
    assert(1==2);
   
  }

  lsol->is_alloc= false;
    

}


void IsNonZero(LinearSolver* lsol,int i,int j){

  assert(lsol->is_init);
  assert(!lsol->is_alloc);

  switch(lsol->storage_type) {

  case SKYLINE :
    SwitchOn((Skyline*)lsol->matrix,i,j);
    break;

  case SKYLINE_SPU :
    SwitchOn_SPU((Skyline_SPU*)lsol->matrix,i,j);
    break;

  default : 
    assert(1==2);
   
  }
  

} 

void AllocateLinearSolver(LinearSolver* lsol){

  assert(lsol->is_init);
  assert(!lsol->is_alloc);

  switch(lsol->storage_type) {

  case SKYLINE :
    AllocateSkyline((Skyline*)lsol->matrix);
    break;

  case SKYLINE_SPU :
    AllocateSkyline_SPU((Skyline_SPU*)lsol->matrix);
    break;

  default : 
    assert(1==2);
   
  }
  lsol->is_alloc=true;
}

void AddLinearSolver(LinearSolver* lsol,int i,int j,schnaps_real val){

  assert(lsol->is_init);
  assert(lsol->is_alloc);

  switch(lsol->storage_type) {

  case SKYLINE :
    AddSkyline((Skyline*)lsol->matrix,i,j,val);
    break;

  case SKYLINE_SPU :
    AddSkyline_SPU((Skyline_SPU*)lsol->matrix,i,j,val);
    break;

  default : 
    assert(1==2);
   
  }

}

void SetLinearSolver(LinearSolver* lsol,int i,int j,schnaps_real val){

  assert(lsol->is_init);
  assert(lsol->is_alloc);

  switch(lsol->storage_type) {

  case SKYLINE :
    SetSkyline((Skyline*)lsol->matrix,i,j,val);
    break;

  case SKYLINE_SPU :
    SetSkyline_SPU((Skyline_SPU*)lsol->matrix,i,j,val);
    break;

  default : 
    assert(1==2);
   
  }

} 

schnaps_real GetLinearSolver(LinearSolver* lsol,int i,int j){

  assert(lsol->is_init);
  assert(lsol->is_alloc);
  
  schnaps_real val;

  switch(lsol->storage_type) {

  case SKYLINE :
    val=GetSkyline((Skyline*)lsol->matrix,i,j);
    break;

  default : 
    assert(1==2);
   
  }
 
  return val;

} 


void DisplayLinearSolver(LinearSolver* lsol){

  assert(lsol->is_init);
  assert(lsol->is_alloc);
  
  switch(lsol->storage_type) {

  case SKYLINE :
    DisplaySkyline((Skyline*)lsol->matrix);
    break;

  case SKYLINE_SPU :
    DisplaySkyline_SPU((Skyline_SPU*)lsol->matrix);
    break;

  default : 
    assert(1==2);
  }

  printf("rhs=");
  for(int i=0; i<lsol->neq; i++){
    printf("%f ", lsol->rhs[i]);
  }

} 

void MatVect(void * system,schnaps_real x[],schnaps_real prod[]){
  int i,j;
  schnaps_real aij;
  LinearSolver* lsol=system;
  
  switch(lsol->storage_type) {

  case SKYLINE :

    //if (!((Skyline*)lsol->matrix)->is_lu){
     MatVectSkyline((Skyline*) lsol->matrix, x, prod);
     //}
     /*else
       {
	 for(i=0;i<lsol->neq;i++)
	   { 
	     prod[i]=0; 
	     for(j=0;j<lsol->neq;j++) { 
	  aij=GetLinearSolver(lsol,i,j); 
	  prod[i] += aij*x[j]; 
        } 
	} 
	}*/
     break;
    
  case SKYLINE_SPU :
    assert(1==2);
    break;

  default : 
    assert(1==2);
  }

  
}

void MatVectIn(void * system){
  int i,j;
  schnaps_real aij;
  LinearSolver* lsol=system;
  
  switch(lsol->storage_type) {

  case SKYLINE :
    assert(1==2); // not yet implemented and verified
    break;
    /* for(i=0;i<lsol->neq;i++) */
    /*   { */
    /* 	prod[i]=0; */
    /* 	for(j=0;j<lsol->neq;j++) { */
    /* 	  aij=GetLinearSolver(lsol,i,j); */
    /* 	  prod[i] += aij*x[j]; */
    /* 	} */
    /*   } */
    
  case SKYLINE_SPU :
    MatVectSkyline_SPU((Skyline_SPU*) lsol->matrix, NULL, NULL);
    break;

  default : 
    assert(1==2);
  }

  
}

void MatVect_SPU(void * system, starpu_data_handle_t sol_handle, starpu_data_handle_t rhs_handle){
  int i,j;
  schnaps_real aij;
  LinearSolver* lsol=system;
  
  switch(lsol->storage_type) {

  case SKYLINE :
    assert(1==2); // not yet implemented and verified
    break;
    /* for(i=0;i<lsol->neq;i++) */
    /*   { */
    /* 	prod[i]=0; */
    /* 	for(j=0;j<lsol->neq;j++) { */
    /* 	  aij=GetLinearSolver(lsol,i,j); */
    /* 	  prod[i] += aij*x[j]; */
    /* 	} */
    /*   } */
    
  case SKYLINE_SPU :
    MatVectSkyline_SPU((Skyline_SPU*) lsol->matrix, sol_handle, rhs_handle);
    break;

  default : 
    assert(1==2);
  }

  
}


void Vector_copy(schnaps_real x[],schnaps_real prod[],int N){
  int i;
 
    for(i=0;i<N;i++)
      {
	prod[i] = x[i];
      }      
}


schnaps_real Vector_norm2(schnaps_real x[],int  N){
  int i;
  schnaps_real norm=0;
 
  for(i=0;i<N;i++)
    {
      norm += x[i]*x[i];
    }
  norm=sqrt(norm);
  return norm;
}




schnaps_real Vector_prodot(schnaps_real x[],schnaps_real y[],int N){
  int i;
  schnaps_real prod;

  prod=0;
    for(i=0;i<N;i++)
      {
	prod += x[i]*y[i];
      }     
    return prod;
}



void LUDecompLinearSolver(LinearSolver* lsol){

  assert(lsol->is_init);
  assert(lsol->is_alloc);
  
  switch(lsol->storage_type) {

  case SKYLINE :
    FactoLU((Skyline*)lsol->matrix);
    break;

  case SKYLINE_SPU :
    FactoLU_SPU((Skyline_SPU*)lsol->matrix);
    break;

  default : 
    assert(1==2);
  }

}

void SolveLinearSolver(LinearSolver* lsol){
  
  assert(lsol->is_init);
  assert(lsol->is_alloc);
  assert(lsol->rhs);
  assert(lsol->sol);


  if(lsol->solver_type == LU) {
       Skyline* sky;
       Skyline_SPU* sky_spu;
       
       switch(lsol->storage_type) {
       case SKYLINE :
         sky=(Skyline*)lsol->matrix;
	 if (!sky->is_lu){
	   FactoLU(sky);
	  }
	 SolveSkyline(sky,lsol->rhs,lsol->sol);
	 break;
      
       case SKYLINE_SPU :
         sky_spu=(Skyline_SPU*)lsol->matrix;
	 if (!sky_spu->is_lu){
	   //UnRegisterSkyline_SPU(sky_spu);
	   //FactoLU_SPU(sky_spu);
	   printf("to do: insert factolu in the starpu system !!!!!\n");
	   //RegisterSkyline_SPU(sky_spu);
	  }
	 SolveSkyline_SPU(sky_spu);
	 break;
      
       default : 
	 assert(1==2);      
       }
  }
  else if(lsol->solver_type == GMRES) {
    GMRESSolver(lsol);
  }
  else {
#ifdef PARALUTION
    Solver_Paralution(lsol);
#else
    printf("paralution is not build this solver is not accessible.");
    assert(1==2);
#endif
  }  
   
}



void Solver_Paralution(LinearSolver* lsol){
  int * rows=NULL;
  int * cols=NULL;
  schnaps_real * coefs=NULL;
  
  schnaps_real * mat_coefs=NULL;
  schnaps_real * RHS=NULL;
  schnaps_real * Sol=NULL;
  char * solver;
  char * pc;
  char * storage;
  schnaps_real * residu=0; 
  int nnz=0,n=0,c=0;
  Skyline * mat;
  
  int basis_size_gmres=30, ILU_p=2,ILU_q=2;
  int* iter_final=0;
  int* ierr=0;
  int maxit=10000;
  schnaps_real norm_rhs=0;
  schnaps_real a_tol=0,r_tol=0,div_tol=1.e+8;

  storage="CSR";
  norm_rhs=Vector_norm2(lsol->rhs,lsol->neq);
  a_tol=lsol->tol*(1.0+1.e-20*norm_rhs);

  if(lsol->neq<61) {
    basis_size_gmres = (int) (lsol->neq/2);
    }
  else {
      basis_size_gmres = 30;
    }

  switch(lsol->solver_type){
  case PAR_CG :
    solver="CG";
    break;
  case PAR_GMRES :
    solver="GMRES";
    break;
  case PAR_FGMRES :
   solver="FGMRES";
    break;
  case PAR_BICGSTAB :
    solver="BiCGStab";
     break;
  case PAR_AMG :
    solver="AMG";
     break;
  case PAR_LU :
    solver="LU";
    storage="DENSE";
     break;
  case PAR_QR :
    solver="QR";
    storage="DENSE";
     break;
  default : 
    assert(1==2);   
  }


  switch(lsol->pc_type){
  case NONE :
    pc="None";
    break;   
  case PAR_JACOBI :
    pc="Jacobi";
    break;
  case PAR_ILU :
    pc="ILU";
    break;
  case PAR_MULTICOLOREDGS :
   pc="MultiColoredGS";
    break;
  case PAR_MULTICOLOREDSGS :
    pc="MultiColoredSGS";
     break;
  case PAR_MULTICOLOREDILU :
    pc="MultiColoredILU";
     break;
  case PAR_AMG_PC :
    pc="AMG";
     break;
  case PAR_ELIMI :
     pc="ELIMI";
     break;  
  default : 
    assert(1==2);   
  }


  
  n=lsol->neq;
  RHS = calloc(n,sizeof(schnaps_real));
  Sol = calloc(n,sizeof(schnaps_real));

  for(int i=0;i<n;i++){
    RHS[i] = (schnaps_real) lsol->rhs[i];
    Sol[i] = (schnaps_real )lsol->sol[i];   
  }



  
 switch(lsol->storage_type) {
  case SKYLINE :

    mat = lsol->matrix;
    
    nnz= mat->neq+2*mat->nmem;
  
    rows = (int*) malloc(nnz*sizeof(int)); 
    cols = (int*) malloc(nnz*sizeof(int));
    coefs = (schnaps_real*) malloc(nnz*sizeof(schnaps_real));
    assert(rows);
  
    for (int i=0;i< mat->neq; i++) {
      for (int j=0;j< mat->neq; j++) {
	if (j-i <= mat->prof[j] && i-j <= mat->prof[i]){
	  if (i==j){
	    coefs[c]=mat->vkgd[i];
	    rows[c]=i;
	    cols[c]=j;
	    c++;
	  }
	  else if ( j>i){
	    int k=mat->kld[j+1]-j+i;
	    coefs[c]=mat->vkgs[k];
	    rows[c]=i;
	    cols[c]=j;
	    c++; 
	  }
	  else {
	    int k=mat->kld[i+1]-i+j;
	    coefs[c]=mat->vkgi[k];
	    rows[c]=i;
	    cols[c]=j;
	    c++; 
	  }
	}
      }
    }    
    
    mat_coefs = malloc(nnz*sizeof(schnaps_real));
    for(int i=0;i<nnz;i++){
      mat_coefs[i] = (schnaps_real) coefs[i];
    }
    
#ifdef PARALUTION
    paralution_fortran_solve_coo(n,n,nnz,solver,storage,pc,storage,
    				 rows,cols,mat_coefs,RHS,a_tol,r_tol,div_tol,maxit,
    				 basis_size_gmres,ILU_p,ILU_q,Sol); 
#endif /* PARALUTION */
    break;

  case CSR :
    assert(1==2);
    break;

  default : 
    assert(1==2);
  }
  
  for(int i=0;i<n;i++){
    lsol->sol[i] = (schnaps_real) Sol[i];
  }
  
}



void GMRESSolver(LinearSolver* lsol){
  int revcom, colx, coly, colz, nbscal;
  int li_maxiter;
  int m,lwork,N;
  int * pt_m;
  int * pt_lwork;
  int * pt_Size;

  int irc[5+1];
  int icntl[8+1];
  int info[3+1];
  schnaps_real cntl[5+1];
  schnaps_real rinfo[2+1];
  schnaps_real sum,err,sum_rhs,lr_tol;
  schnaps_real * work;
  schnaps_real *loc_x;
  schnaps_real *loc_y;
  schnaps_real *loc_z;
  schnaps_real prodot=0.0;
  int res=0;
  int matvec=1, precondLeft=2, precondRight=3, dotProd=4;


  lsol->MatVecProduct=MatVect;
  
  res=init_dgmres(icntl,cntl);
  
  icntl[3]  = 6 ;           // output unit
  icntl[7]  = lsol->iter_max; // Maximum number of iterations
  icntl[4]  = 0; //!1            // preconditioner (1) = left preconditioner
  icntl[5]  = 1; ////3            // orthogonalization scheme
  icntl[6]  = 1; //1            // initial guess  (1) = user supplied guess
  icntl[8]  = 1; //1

  
  if ((lsol->pc_type == JACOBI) ||(lsol->pc_type == EXACT)){
    icntl[4] = 2;
  }
     
  cntl[1]  = lsol->tol; //       ! stopping tolerance
  cntl[2]  = 1.0;
  cntl[3]  = 1.0;
  cntl[4]  = 1.0;         
  cntl[5]  = 1.0;

  N = lsol->neq;

  if(lsol->restart_gmres == 1) {
    if(N<61) {
      m = (int) (N/2)+1;
    }
    else {
      m = 30;
    }
  }
  else {
      m = lsol->restart_gmres;
    }
  lwork = m*m + m*(N+5) + 5*N + m +1; //(+ one because  ) ??

  pt_m = &m;
  pt_Size = &N;
  pt_lwork = &lwork;

  work = calloc(lwork, sizeof(schnaps_real));
  loc_x = calloc(N, sizeof(schnaps_real));
  loc_y = calloc(N, sizeof(schnaps_real));
  loc_z = calloc(N, sizeof(schnaps_real));
  
  for(int ivec = 0; ivec < N; ivec++) {
    work[ivec+1]     = lsol->sol[ivec];                    
    work[N+ivec+1]    = lsol->rhs[ivec];
  }


  schnaps_real * Ax2=calloc(lsol->neq,sizeof(schnaps_real));
  lsol->MatVecProduct(lsol,lsol->sol,Ax2);
  schnaps_real errorb=0;
  schnaps_real errorb2=0;
   for(int i = 0; i < N; i++) {
     errorb=errorb+fabs((Ax2[i]-lsol->rhs[i])*(Ax2[i]-lsol->rhs[i]));                    
     errorb2=errorb2+fabs(lsol->sol[i]*lsol->sol[i]);                    
  }
  printf(" error gmres begin %.5e \n",sqrt(errorb)/(sqrt(errorb2)+1.0)); 

  //*****************************************
  //** Reverse communication implementation
  //*****************************************
  L10:    res=drive_dgmres(pt_Size,pt_Size,pt_m,pt_lwork,&work[1],&irc[1],&icntl[1],&cntl[1],&info[1],&rinfo[1]);

  revcom = irc[1];
  colx   = irc[2];
  coly   = irc[3];
  colz   = irc[4];
  nbscal = irc[5];

  for(int ivec = 0; ivec < N; ivec++) {
  
    loc_z[ivec]= work[colz+ivec];                    
    loc_x[ivec]= work[colx+ivec];
    loc_y[ivec]= work[coly+ivec];    
  }

  if (revcom == matvec) {                 // perform the matrix vector product
    // work(colz) <-- A * work(colx)   
    lsol->MatVecProduct(lsol,loc_x,loc_z);
    for(int ivec = 0; ivec < N; ivec++) {       
      work[colz+ivec]= loc_z[ivec];                    
      work[colx+ivec]= loc_x[ivec]; 
    }
     goto L10;
  }
  else if(revcom == precondLeft)  {                 // perform the matrix vector product
    // work(colz) <-- M-1 * work(colx)
    Vector_copy(loc_x,loc_z,N);
    for(int ivec = 0; ivec < N; ivec++) {
      work[colz+ivec]= loc_z[ivec];                    
      work[colx+ivec]= loc_x[ivec]; 
    }
     goto L10;
  }

  else if(revcom == precondRight)  {
    if(lsol->pc_type == EXACT){
      Exact_PC(lsol,loc_z,loc_x);
    }
    else if (lsol->pc_type == JACOBI){
      Jacobi_PC(lsol,loc_z,loc_x);
    }
    else {
      Vector_copy(loc_x,loc_z,N);
    }

    
    for(int ivec = 0; ivec < N; ivec++) {
      work[colz+ivec]= loc_z[ivec];                    
      work[colx+ivec]= loc_x[ivec]; 
    }
      goto L10;
  }

  else if(revcom == dotProd)  {// perform the matrix vector product
    // work(colz) <-- work(colx) work(coly)
    prodot=Vector_prodot(loc_x,loc_y,N);
    for(int ivec = 0; ivec < N; ivec++) {
      work[colx+ivec]= loc_x[ivec];
      work[coly+ivec]= loc_y[ivec]; 
    }
    work[colz]= prodot;  
    	 goto L10;
  }

  //******************************** end of GMRES reverse communication

  
  
  for(int ivec = 0; ivec < N; ivec++) {
    lsol->sol[ivec] = work[ivec+1];                    
  }

  schnaps_real * Ax=calloc(lsol->neq,sizeof(schnaps_real));
  lsol->MatVecProduct(lsol,lsol->sol,Ax);
  schnaps_real error=0;
  schnaps_real error2=0;
   for(int i = 0; i < N; i++) {
     error=error+fabs((Ax[i]-lsol->rhs[i])*(Ax[i]-lsol->rhs[i]));                    
     error2=error2+fabs(lsol->sol[i]*lsol->sol[i]);                    
  }
   printf(" error gmres end %.5e \n",sqrt(error)/(sqrt(error2)+1.0)); 
  
  
  free(work);
  free(loc_x);
  free(loc_y);
  free(loc_z);

}



void Jacobi_PC(LinearSolver *lsol, schnaps_real* sol, schnaps_real* rhs){

  for (int i=0;i<lsol->neq; i++){
    assert(GetLinearSolver(lsol,i,i)!=0);
    sol[i]=rhs[i]/GetLinearSolver(lsol,i,i);
  }
}

void Exact_PC(LinearSolver *lsol, schnaps_real* sol, schnaps_real* rhs){

  Skyline * mat;
  Skyline mat_copy;

  mat=(Skyline*)lsol->matrix; 
  
  mat_copy.is_alloc= mat->is_alloc;
  mat_copy.copy_is_alloc=mat->copy_is_alloc;
  mat_copy.is_sym=mat->is_sym;
  mat_copy.is_lu=mat->is_lu;

  mat_copy.neq=mat->neq;
  mat_copy.nmem=mat->nmem;

  mat_copy.vkgd=malloc(mat_copy.neq*sizeof(schnaps_real));
  mat_copy.vkgi=malloc(mat_copy.nmem*sizeof(schnaps_real));
  mat_copy.vkgs=malloc(mat_copy.nmem*sizeof(schnaps_real));
  

  mat_copy.prof=malloc(mat_copy.neq*sizeof(int));
  mat_copy.kld=malloc((mat_copy.neq+1)*sizeof(int));

  for (int i=0;i<mat->neq; i++){
    mat_copy.vkgd[i]=mat->vkgd[i];
  }
  for (int i=0;i<mat->nmem; i++){
    mat_copy.vkgi[i]=mat->vkgi[i];
    mat_copy.vkgs[i]=mat->vkgs[i];
  }

  for (int i=0;i<mat->neq+1; i++){
    mat_copy.kld[i]=mat->kld[i];
  }
  for (int i=0;i<mat->neq; i++){
    mat_copy.prof[i]=mat->prof[i];
  }
  
  FactoLU(&mat_copy);
  SolveSkyline(&mat_copy,rhs,sol);
  
}



