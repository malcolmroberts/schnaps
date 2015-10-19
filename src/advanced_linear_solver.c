#include "linear_solver.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "skyline.h"
#include "paralution_c.h"
#include "dpackfgmres.h"
#include "physBased_PC.h"
#include "advanced_linear_solver.h"
#include "linear_solver.h"

void Advanced_SolveLinearSolver(LinearSolver* lsol, Simulation* simu){
  
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
    Advanced_GMRESSolver(lsol,simu);
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




void Advanced_GMRESSolver(LinearSolver* lsol, Simulation* simu){
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
  if(lsol->pc_type == PHY_BASED_P1){
     int mat2assemble[6] = {1, 1, 1, 1, 1, 1};
     Init_PBPC_Wave_SchurPressure_BCVelocity(simu, &pb_pc, mat2assemble);
     Init_Parameters_PhyBasedPC(&pb_pc);
     icntl[4]  = 2;
  }
  if(lsol->pc_type == PHY_BASED_P2){
     int mat2assemble[6] = {1, 1, 1, 1, 1, 1};
     Init_PBPC_Wave_SchurPressure_BCPressure(simu, &pb_pc, mat2assemble);
     Init_Parameters_PhyBasedPC(&pb_pc);
     icntl[4]  = 2;
  }
   if(lsol->pc_type == PHY_BASED_U){
     int mat2assemble[6] = {1, 1, 1, 1, 1, 1};
     Init_PBPC_Wave_SchurVelocity(simu, &pb_pc, mat2assemble);
     Init_Parameters_PhyBasedPC(&pb_pc);
     icntl[4]  = 2;
  }
  else if ((lsol->pc_type == JACOBI) ||(lsol->pc_type == EXACT)){
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
    if(lsol->pc_type == PHY_BASED_P1 || lsol->pc_type == PHY_BASED_P2 ){
       if (lsol->is_CG){
	 //PhyBased_PC_CG(&pb_pc,simu,loc_z,loc_x);
	 PhyBased_PC_InvertSchur_CG(&pb_pc,simu,loc_z,loc_x);
       }
    }
    else if(lsol->pc_type == PHY_BASED_U){
	 PhyBased_PC_CG(&pb_pc,simu,loc_z,loc_x);
    }  
    else if(lsol->pc_type == EXACT){
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
    if(lsol->pc_type == PHY_BASED_P1 || lsol->pc_type == PHY_BASED_P2){
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





