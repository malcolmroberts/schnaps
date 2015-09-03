#include "linear_solver.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "skyline.h"
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
  lsol->iter_max=10000;
  lsol->is_CG=false;

  if (matstor != NULL) lsol->storage_type = *matstor;
  if (solvtyp != NULL) lsol->solver_type = *solvtyp;

  lsol->rhs=calloc(n,sizeof(real));
  lsol->sol=calloc(n,sizeof(real));

  switch(lsol->storage_type) {

  Skyline* sky;

  case SKYLINE :
    sky = malloc(sizeof(Skyline));
    assert(sky);
    lsol->matrix = (void*) sky;
    InitSkyline(sky,n);
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

  default : 
    assert(1==2);
   
  }
  lsol->is_alloc=true;
}

void AddLinearSolver(LinearSolver* lsol,int i,int j,real val){

  assert(lsol->is_init);
  assert(lsol->is_alloc);

  switch(lsol->storage_type) {

  case SKYLINE :
    AddSkyline((Skyline*)lsol->matrix,i,j,val);
    break;

  default : 
    assert(1==2);
   
  }

}

void SetLinearSolver(LinearSolver* lsol,int i,int j,real val){

  assert(lsol->is_init);
  assert(lsol->is_alloc);

  switch(lsol->storage_type) {

  case SKYLINE :
    SetSkyline((Skyline*)lsol->matrix,i,j,val);
    break;

  default : 
    assert(1==2);
   
  }

} 

real GetLinearSolver(LinearSolver* lsol,int i,int j){

  assert(lsol->is_init);
  assert(lsol->is_alloc);
  
  real val;

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

  default : 
    assert(1==2);
  }

} 

void MatVect(void * system,real x[],real prod[]){
  int i,j;
  real aij;
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

  default : 
    assert(1==2);
  }

  
}

void MatVect_slow(void * system,real x[],real prod[]){
  int i,j;
  real aij;
  LinearSolver* lsol=system;
  
  for(i=0;i<lsol->neq;i++)
    { 
      prod[i]=0; 
      for(j=0;j<lsol->neq;j++) { 
	aij=GetLinearSolver(lsol,i,j); 
	prod[i] += aij*x[j]; 
      } 
    } 
   
  
}

void Vector_copy(real x[],real prod[],int N){
  int i;
 
    for(i=0;i<N;i++)
      {
	prod[i] = x[i];
      }      
}


real Vector_norm2(real x[],int  N){
  int i;
  real norm=0;
 
  for(i=0;i<N;i++)
    {
      norm += x[i]*x[i];
    }
  norm=sqrt(norm);
  return norm;
}




real Vector_prodot(real x[],real y[],int N){
  int i;
  real prod;

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

  default : 
    assert(1==2);
  }

}

void SolveLinearSolver(LinearSolver* lsol, Simulation* simu){
  
  assert(lsol->is_init);
  assert(lsol->is_alloc);
  assert(lsol->rhs);
  assert(lsol->sol);


  if(lsol->solver_type == LU) {
       Skyline* sky;
       switch(lsol->storage_type) {
       case SKYLINE :
         sky=(Skyline*)lsol->matrix;
	 if (!sky->is_lu){
	   FactoLU(sky);
	  }
	 SolveSkyline(sky,lsol->rhs,lsol->sol);
	 break;
      
       default : 
	 assert(1==2);      
       }
  }
  else if(lsol->solver_type == GMRES) {
    GMRESSolver(lsol,simu);
  }
  else {
#ifdef PARALUTION
    Solver_Paralution(lsol,simu);
#else
    printf("paralution is not build this solver is not accessible.");
    assert(1==2);
#endif
  }  
   
}



void InitJFLinearSolver(JFLinearSolver* lsol,int n,
		      Solver* solvtyp){

  lsol->neq=n;
  lsol->solver_type = GMRES;
  lsol->pc_type = NONE;
  lsol->rhs=NULL;
  lsol->sol=NULL;
  lsol->soln=NULL;
  lsol->MatVecProduct=NULL;
  lsol->NonlinearVector_computation=NULL;
  lsol->tol=1.e-6;
  lsol->restart_gmres=1;
  lsol->iter_max=10000;
  lsol->eps=0.000001;

  if (solvtyp != NULL) lsol->solver_type = *solvtyp;

  lsol->rhs=calloc(n,sizeof(real));
  lsol->sol=calloc(n,sizeof(real));
  lsol->soln=calloc(n,sizeof(real));

}

void FreeJFLinearSolver(JFLinearSolver* lsol){

  free(lsol->rhs);
  free(lsol->sol);
  free(lsol->soln);

}

void MatVecJacobianFree(Simulation *simu,void * system,real x[],real prod[]){
  int i,j;
  real aij;
  JFLinearSolver* lsol=system;
  real * U;
  real * Up;
  real * solnp;
  
  solnp=calloc(lsol->neq,sizeof(real));
  U=calloc(lsol->neq,sizeof(real));
  Up=calloc(lsol->neq,sizeof(real));

  
  for(i=0;i<lsol->neq;i++)
    {
	solnp[i]=lsol->soln[i]+lsol->eps*x[i];
    }
  
  lsol->NonlinearVector_computation(simu,system,lsol->soln,U);
  lsol->NonlinearVector_computation(simu,system,solnp,Up);
  
  for(i=0;i<lsol->neq;i++)
    {
      prod[i]=(Up[i]-U[i])/lsol->eps;
    }  

  free(solnp);
  free(U);
  free(Up);

  
}



void Solver_Paralution(LinearSolver* lsol, Simulation* simu){
  int * rows=NULL;
  int * cols=NULL;
  real * coefs=NULL;
  
  real * mat_coefs=NULL;
  real * RHS=NULL;
  real * Sol=NULL;
  char * solver;
  char * pc;
  char * storage;
  real * residu=0; 
  int nnz=0,n=0,c=0;
  Skyline * mat;
  
  int basis_size_gmres=30, ILU_p=2,ILU_q=2;
  int* iter_final=0;
  int* ierr=0;
  int maxit=100000;
  real norm_rhs=0;
  real a_tol=0,r_tol=0,div_tol=1.e+8;

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
  RHS = calloc(n,sizeof(real));
  Sol = calloc(n,sizeof(real));

  for(int i=0;i<n;i++){
    RHS[i] = (real) lsol->rhs[i];
    Sol[i] = (real )lsol->sol[i];   
  }



  
 switch(lsol->storage_type) {
  case SKYLINE :

    mat = lsol->matrix;
    
    nnz= mat->neq+2*mat->nmem;
  
    rows = (int*) malloc(nnz*sizeof(int)); 
    cols = (int*) malloc(nnz*sizeof(int));
    coefs = (real*) malloc(nnz*sizeof(real));
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
    
    mat_coefs = malloc(nnz*sizeof(real));
    for(int i=0;i<nnz;i++){
      mat_coefs[i] = (real) coefs[i];
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
    lsol->sol[i] = (real) Sol[i];
  }
  
}



void GMRESSolver(LinearSolver* lsol, Simulation* simu){
  int revcom, colx, coly, colz, nbscal;
  int li_maxiter;
  int m,lwork,N;
  int * pt_m;
  int * pt_lwork;
  int * pt_Size;

  int irc[5+1];
  int icntl[8+1];
  int info[3+1];
  real cntl[5+1];
  real rinfo[2+1];
  real sum,err,sum_rhs,lr_tol;
  real * work;
  real *loc_x;
  real *loc_y;
  real *loc_z;
  real prodot=0.0;
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

  PB_PC pb_pc;
  if((lsol->pc_type == PHY_BASED) || (lsol->pc_type == PHY_BASED_EXACT)){
     int mat2assemble[6] = {1, 1, 1, 1, 1, 1};
     //Init_PhyBasedPC_SchurVelocity_Wave(simu, &pb_pc, mat2assemble);
     Init_PhyBasedPC_SchurPressure_Wave(simu, &pb_pc, mat2assemble);
     Init_Parameters_PhyBasedPC(&pb_pc);
     icntl[4]  = 2;
  }
  else if ((lsol->pc_type == JACOBI) ||(lsol->pc_type == EXACT)){
    icntl[4] = 2;
  }
  if (lsol->pc_type == PHDF){
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

  work = calloc(lwork, sizeof(real));
  loc_x = calloc(N, sizeof(real));
  loc_y = calloc(N, sizeof(real));
  loc_z = calloc(N, sizeof(real));

  
  for(int ivec = 0; ivec < N; ivec++) {
    work[ivec+1]     = lsol->sol[ivec];                    
    work[N+ivec+1]    = lsol->rhs[ivec];
  }


  real * Ax2=calloc(lsol->neq,sizeof(real));
  lsol->MatVecProduct(lsol,lsol->sol,Ax2);
  real errorb=0;
  real errorb2=0;
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
    if(lsol->pc_type == PHY_BASED){    
       if (lsol->is_CG){
	 //PhyBased_PC_CG(&pb_pc,simu,loc_z,loc_x);
	 PhyBased_PC_InvertSchur_CG(&pb_pc,simu,loc_z,loc_x);
       }
       else {
        PhyBased_PC_DG(&pb_pc,simu,loc_z,loc_x);
        //solveIdentity(&pb_pc,simu,loc_z,loc_x);
      }
    }
    else if(lsol->pc_type == PHY_BASED_EXACT){
      if (lsol->is_CG){
	PhyBased_PC_Full(&pb_pc,simu,loc_z,loc_x);
      }
       else {
	 Vector_copy(loc_x,loc_z,N);
       }
    }
    else if(lsol->pc_type == EXACT){
      Exact_PC(lsol,simu,loc_z,loc_x);
    }
    else if (lsol->pc_type == JACOBI){
      Jacobi_PC(lsol,simu,loc_z,loc_x);
    }
    else if (lsol->pc_type == PHDF){
      PhyBasedPC_waveDF(lsol,simu,loc_z,loc_x);
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

  real * Ax=calloc(lsol->neq,sizeof(real));
  lsol->MatVecProduct(lsol,lsol->sol,Ax);
  real error=0;
  real error2=0;
   for(int i = 0; i < N; i++) {
     error=error+fabs((Ax[i]-lsol->rhs[i])*(Ax[i]-lsol->rhs[i]));                    
     error2=error2+fabs(lsol->sol[i]*lsol->sol[i]);                    
  }
   printf(" error gmres end %.5e \n",sqrt(error)/(sqrt(error2)+1.0)); 
  
  
  free(work);
  free(loc_x);
  free(loc_y);
  free(loc_z);
    if(lsol->pc_type == PHY_BASED){
      freePB_PC(&pb_pc);
    }

}




void SolveJFLinearSolver(JFLinearSolver* lsol,Simulation * simu){
  int revcom, colx, coly, colz, nbscal;
  int li_maxiter;
  int m,lwork,N;
  int * pt_m;
  int * pt_lwork;
  int * pt_Size;

  int irc[5+1];
  int icntl[8+1];
  int info[3+1];
  real cntl[5+1];
  real rinfo[2+1];
  real sum,err,sum_rhs,lr_tol;
  real * work;
  real *loc_x;
  real *loc_y;
  real *loc_z;
  real prodot=0.0;
  int res=0;
  int matvec=1, precondLeft=2, precondRight=3, dotProd=4;


  lsol->MatVecProduct=MatVecJacobianFree;

  res=init_dgmres(icntl,cntl);
  
  icntl[3]  = 6 ;           // output unit
  icntl[7]  = lsol->iter_max; // Maximum number of iterations
  icntl[4]  = 0; //!1            // preconditioner (1) = left preconditioner
  icntl[5]  = 1; ////3            // orthogonalization scheme
  icntl[6]  = 1; //1            // initial guess  (1) = user supplied guess
  icntl[8]  = 1; //1            
   
     
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
  lwork = m*m + m*(N+5) + 5*N + m + 1 +1; //(+ one because  ) ??

  pt_m = &m;
  pt_Size = &N;
  pt_lwork = &lwork;

  work = calloc(lwork, sizeof(real));
  loc_x = calloc(N, sizeof(real));
  loc_y = calloc(N, sizeof(real));
  loc_z = calloc(N, sizeof(real));
  
  for(int ivec = 0; ivec < N; ivec++) {
    work[ivec+1]     = lsol->sol[ivec];                    
    work[N+ivec+1]    = lsol->rhs[ivec];
  }

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
    lsol->MatVecProduct(simu,lsol,loc_x,loc_z);
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

  else if(revcom == precondRight)  {                 // perform the matrix vector product
    // work(colz) <-- M-1 * work(colx)  
    Vector_copy(loc_x,loc_z,N);
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


  free(work);
  free(loc_x);
  free(loc_y);
  free(loc_z);


}

void Jacobi_PC(LinearSolver *lsol, Simulation* simu, real* sol, real* rhs){

  for (int i=0;i<lsol->neq; i++){
    assert(GetLinearSolver(lsol,i,i)!=0);
    sol[i]=rhs[i]/GetLinearSolver(lsol,i,i);
  }
}

void Exact_PC(LinearSolver *lsol, Simulation* simu, real* sol, real* rhs){

  Skyline * mat;
  Skyline mat_copy;

  mat=(Skyline*)lsol->matrix; 
  
  mat_copy.is_alloc= mat->is_alloc;
  mat_copy.copy_is_alloc=mat->copy_is_alloc;
  mat_copy.is_sym=mat->is_sym;
  mat_copy.is_lu=mat->is_lu;

  mat_copy.neq=mat->neq;
  mat_copy.nmem=mat->nmem;

  mat_copy.vkgd=malloc(mat_copy.neq*sizeof(real));
  mat_copy.vkgi=malloc(mat_copy.nmem*sizeof(real));
  mat_copy.vkgs=malloc(mat_copy.nmem*sizeof(real));
  

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



void PhyBasedPC_waveDF(LinearSolver *lsol,Simulation * simu, real * Sol, real *RHS){


  int nwave;
  nwave=lsol->neq/2;
  LinearSolver Identity, schur,schur_approx,LDx,UDx;

  InitLinearSolver(&Identity,nwave,NULL,NULL);
  InitLinearSolver(&schur,nwave,NULL,NULL);
  InitLinearSolver(&schur_approx,nwave,NULL,NULL);
  InitLinearSolver(&LDx,nwave,NULL,NULL);
  InitLinearSolver(&UDx,nwave,NULL,NULL);

  real h=1.0/nwave;
  real dt=0,c=0;
  dt=simu->dt;
  c=simu->vmax;
  int approx_schur=0;
  
  ////////////// begin Identity ///////////
  for(int i=0;i<nwave;i++){   
     IsNonZero(&Identity,i,i); 
  }


  AllocateLinearSolver(&Identity);

  for(int i=0;i<nwave;i++){   
    AddLinearSolver(&Identity,i,i,1.0); 
  }
  ////////////// fin Identity ///////////

  ////////////// begin LDx ///////////
  for(int i=0;i<nwave;i++){
    if (i==0){ 
      IsNonZero(&LDx,0,0);
      IsNonZero(&LDx,0,1);
    }
    else if (i==nwave-1){ 
      IsNonZero(&LDx,nwave-1,nwave-1);
      IsNonZero(&LDx,nwave-1,nwave-2);
    } 
    else { 
      IsNonZero(&LDx,i,i);
      IsNonZero(&LDx,i,i-1);    
      IsNonZero(&LDx,i,i+1);
    } 
  }
  
  AllocateLinearSolver(&LDx);
  
  for(int i=0;i<nwave;i++){
    if (i==0){ 
      AddLinearSolver(&LDx,0,0,0.0);
      AddLinearSolver(&LDx,0,1,(c*dt)/(2.0*h));
    }
    else if (i==nwave-1){ 
      AddLinearSolver(&LDx,nwave-1,nwave-1,0.0);
      AddLinearSolver(&LDx,nwave-1,nwave-2,-(c*dt)/(2.0*h));
    } 
    else { 
      AddLinearSolver(&LDx,i,i,0.0);
      AddLinearSolver(&LDx,i,i-1,-(c*dt)/(2.0*h));    
      AddLinearSolver(&LDx,i,i+1,+(c*dt)/(2.0*h));
    } 
  }
  
  for(int i=0;i<nwave;i++){
    SetLinearSolver(&LDx,nwave-1,i,0.0);
    SetLinearSolver(&LDx,0,i,0.0);
  }
  ////////////// fin LDx ///////////


   ////////////// begin UDx ///////////
  for(int i=0;i<nwave;i++){
    if (i==0){ 
      IsNonZero(&UDx,0,0);
      IsNonZero(&UDx,0,1);
    }
    else if (i==nwave-1){ 
      IsNonZero(&UDx,nwave-1,nwave-1);
      IsNonZero(&UDx,nwave-1,nwave-2);
    } 
    else { 
      IsNonZero(&UDx,i,i);
      IsNonZero(&UDx,i,i-1);    
      IsNonZero(&UDx,i,i+1);
    } 
  }
  
  AllocateLinearSolver(&UDx);
  
  for(int i=0;i<nwave;i++){
    if (i==0){ 
      AddLinearSolver(&UDx,0,0,0.0);
      AddLinearSolver(&UDx,0,1,(c*dt)/(2.0*h));
    }
    else if (i==nwave-1){ 
      AddLinearSolver(&UDx,nwave-1,nwave-1,0.0);
      AddLinearSolver(&UDx,nwave-1,nwave-2,-(c*dt)/(2.0*h));
    } 
    else { 
      AddLinearSolver(&UDx,i,i,0.0);
      AddLinearSolver(&UDx,i,i-1,-(c*dt)/(2.0*h));    
      AddLinearSolver(&UDx,i,i+1,+(c*dt)/(2.0*h));
    } 
  }

  for(int i=0;i<nwave;i++){
    SetLinearSolver(&UDx,nwave-1,i,0.0);
    SetLinearSolver(&UDx,0,i,0.0);
    }
  ////////////// fin UDx ///////////
  
  
  ////////////// begin Schur_approx ///////////
  for(int i=0;i<nwave;i++){
    if (i==0){ 
      IsNonZero(&schur_approx,0,0);
      IsNonZero(&schur_approx,0,1);
    }
    else if (i==nwave-1){ 
      IsNonZero(&schur_approx,nwave-1,nwave-1);
      IsNonZero(&schur_approx,nwave-1,nwave-2);
    } 
    else { 
      IsNonZero(&schur_approx,i,i);
      IsNonZero(&schur_approx,i,i-1);    
      IsNonZero(&schur_approx,i,i+1);
    } 
  }
  
  AllocateLinearSolver(&schur_approx);
  
  for(int i=0;i<nwave;i++){
    if (i==0){ 
      AddLinearSolver(&schur_approx,0,0,0.0);
      AddLinearSolver(&schur_approx,0,1,(c*dt)/(2.0*h));
    }
    else if (i==nwave-1){ 
      AddLinearSolver(&schur_approx,nwave-1,nwave-1,0.0);
      AddLinearSolver(&schur_approx,nwave-1,nwave-2,-(c*dt)/(2.0*h));
    } 
    else { 
      AddLinearSolver(&schur_approx,i,i,1.0+2*((c*dt)/h)*((c*dt)/h));
      AddLinearSolver(&schur_approx,i,i-1,-((c*dt)/h)*((c*dt)/h));    
      AddLinearSolver(&schur_approx,i,i+1,-((c*dt)/h)*((c*dt)/h));
    } 
  }

  for(int i=0;i<nwave;i++){
    SetLinearSolver(&schur_approx,nwave-1,i,0.0);
    SetLinearSolver(&schur_approx,0,i,0.0);
  }
  SetLinearSolver(&schur_approx,nwave-1,nwave-1,1.0);
  SetLinearSolver(&schur_approx,0,0,1.0);


  
  ////////////// fin Schur_approx ///////////


  
   ////////////// begin Schur ///////////
  for(int i=0;i<nwave;i++){
    for(int j=0;j<nwave;j++){
      IsNonZero(&schur,i,j);    
    } 
  }
  
  AllocateLinearSolver(&schur);

  real ** A;
  real LMU=0.0;
  A=calloc(nwave,sizeof(real));
  for (int i=0;i<nwave;i++){
    A[i]=calloc(nwave,sizeof(real));
  }
  
  for (int i=0;i<nwave;i++){ 
    real InvertMii=1.0/GetLinearSolver(&Identity,i,i);
    for (int j=0;j<nwave;j++){
      A[i][j]=InvertMii*GetLinearSolver(&UDx,i,j);
    }
  }
  
  for (int i=0;i<nwave;i++){
    for (int j=0;j<nwave;j++){
      LMU=0.0;
      if(i>0 && i<nwave-1){
	for (int k=0;k<nwave;k++){
	  LMU+=GetLinearSolver(&LDx,i,k)*A[k][j];
	}
      }
      else { /// because by limit condition the last lign of Dx is equal to zero
	LMU=0;
      }
      real val=GetLinearSolver(&Identity,i,j)-LMU;
      SetLinearSolver(&schur,i,j,val);
    }
  }

  ////////////// fin Schur_approx ///////////
  

  //1 ) prediction
  for (int i=0;i<nwave;i++){
    Identity.rhs[i]   = RHS[i*2+1];
    schur.rhs[i] = RHS[i*2];
  }
  
  Identity.solver_type = LU;
  Identity.pc_type=JACOBI;
  Identity.iter_max=20000;
  Identity.tol=1.e-14;

  
   SolveLinearSolver(&Identity,simu);

  
  // 2) PROPAGATION STEP
  
  MatVect(&LDx,Identity.sol,LDx.sol); // Dx u 
  
  for (int i=0;i<nwave;i++){
    schur.rhs[i]   =  schur.rhs[i]- LDx.sol[i];
  }

  schur.solver_type = LU;
  schur.pc_type=NONE;//JACOBI;
  schur.iter_max=20000;
  schur.tol=1.e-14;

  if(approx_schur==1){
    for (int i=0;i<nwave;i++){
      schur_approx.rhs[i]   =  schur.rhs[i];
    }
    SolveLinearSolver(&schur_approx,simu);
    for (int i=0;i<nwave;i++){
      schur.sol[i]   =  schur_approx.sol[i];
    }
  }
  else {
    
    SolveLinearSolver(&schur,simu);

  }
  
  
  // 3) CORRECTION STEP
  
  MatVect(&UDx,schur.sol,UDx.sol); // Dx p

  
  //printf("RHS assembly.....\n");
  for (int i=0;i<nwave;i++){
    Identity.rhs[i] = Identity.sol[i]- UDx.sol[i];
  }

  //printf("Solution...\n");
   SolveLinearSolver(&Identity,simu);
 
  
  // 4) OUTPUT STEP
   for (int i=0;i<nwave;i++){
     Sol[2*i+1]   = Identity.rhs[i];
     Sol[2*i] = schur.rhs[i];
   }
 
}
