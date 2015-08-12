#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "test.h"
#include "schnaps.h"

#define _NN 5

int main(void) {
  
  // unit tests
    
  int resu=TestLinearSolver();
	 

  if (resu) printf("Linear solver test OK !\n");
  else printf("Linear solver test failed !\n");

  return !resu;
} 




int TestLinearSolver(void){

  int test=1,test1=1,test2=1,test3=1,test4=1;
  Simulation simu;

  LinearSolver sky;

#ifdef PARALUTION 
  paralution_begin();
#endif 

  // preliminary work on the skyline struct
  // _NN is the size of the linear system to be solved
  InitLinearSolver(&sky,_NN,NULL,NULL);

  sky.solver_type = LU;
  sky.pc_type=NONE;

  real A[_NN][_NN];
  real vf[_NN],sol[_NN];


  A[0][0] = 0.2e1;
  A[0][1] = -0.1e1;
  A[0][2] = 0;
  A[0][3] = 0;
  A[0][4] = 0;
  A[1][0] = -0.1e1;
  A[1][1] = 0.2e1;
  A[1][2] = -0.1e1;
  A[1][3] = 0;
  A[1][4] = 0;
  A[2][0] = 0;
  A[2][1] = -0.1e1;
  A[2][2] = 0.2e1;
  A[2][3] = -0.1e1;
  A[2][4] = 0;
  A[3][0] = 0;
  A[3][1] = 0;
  A[3][2] = -0.1e1;
  A[3][3] = 0.2e1;
  A[3][4] = -0.1e1;
  A[4][0] = 0;
  A[4][1] = 0;
  A[4][2] = 0;
  A[4][3] = -0.1e1;
  A[4][4] = 0.2e1;
  vf[0] = 0;
  vf[1] = 0;
  vf[2] = 0;
  vf[3] = 0;
  vf[4] = 0.6e1;

  sol[0] = 0;
  sol[1] = 0;
  sol[2] = 0;
  sol[3] = 0;
  sol[4] = 0;



  // first mark the nonzero values in A
  for(int i=0;i<_NN;i++){
    for(int j=0;j<_NN;j++){
      if (A[i][j] != 0) IsNonZero(&sky,i,j);
      //if (i==j) SwitchOn(&sky,i,j);
    }
  }

  // once the nonzero positions are known allocate memory
  AllocateLinearSolver(&sky);

  // now set the nonzero terms
  for(int i=0;i<_NN;i++){
    for(int j=0;j<_NN;j++){
      if (A[i][j] != 0){
      	AddLinearSolver(&sky,i,j,A[i][j]);
      }
      /* if (i==j){ */
      /* 	SetLinearSolver(&sky,i,j,2); */
      /* } */
    }
  }



  // printf for checking...
  DisplayLinearSolver(&sky);


  // solve from a decomposed matrix
  // vf: rhs
  // sol: solution

  sky.rhs=vf;
  sky.sol=sol;
  SolveLinearSolver(&sky,&simu);


  // checking
  real verr=0;
  printf("sol of LU=");
  for(int i=0;i<_NN;i++){
    printf("%f ",sol[i]);
    verr+=fabs(sol[i]-i-1);
  }
  printf("\n");
  printf("\n");
  // deallocate memory
  FreeLinearSolver(&sky,0);
  

  test1= (verr<1e-10);


#ifdef PARALUTION
  // preliminary work on the skyline struct
  // _NN is the size of the linear system to be solved
  InitLinearSolver(&sky,_NN,NULL,NULL);

  sky.solver_type = PAR_GMRES;
  sky.pc_type=PAR_JACOBI;

  A[0][0] = 0.2e1;
  A[0][1] = -0.1e1;
  A[0][2] = 0;
  A[0][3] = 0;
  A[0][4] = 0;
  A[1][0] = -0.1e1;
  A[1][1] = 0.2e1;
  A[1][2] = -0.1e1;
  A[1][3] = 0;
  A[1][4] = 0;
  A[2][0] = 0;
  A[2][1] = -0.1e1;
  A[2][2] = 0.2e1;
  A[2][3] = -0.1e1;
  A[2][4] = 0;
  A[3][0] = 0;
  A[3][1] = 0;
  A[3][2] = -0.1e1;
  A[3][3] = 0.2e1;
  A[3][4] = -0.1e1;
  A[4][0] = 0;
  A[4][1] = 0;
  A[4][2] = 0;
  A[4][3] = -0.1e1;
  A[4][4] = 0.2e1;
  vf[0] = 0;
  vf[1] = 0;
  vf[2] = 0;
  vf[3] = 0;
  vf[4] = 0.6e1;

  sol[0] = 0;
  sol[1] = 0;
  sol[2] = 0;
  sol[3] = 0;
  sol[4] = 0;



  // first mark the nonzero values in A
  for(int i=0;i<_NN;i++){
    for(int j=0;j<_NN;j++){
      if (A[i][j] != 0) IsNonZero(&sky,i,j);
      //if (i==j) SwitchOn(&sky,i,j);
    }
  }

  // once the nonzero positions are known allocate memory
  AllocateLinearSolver(&sky);

  // now set the nonzero terms
  for(int i=0;i<_NN;i++){
    for(int j=0;j<_NN;j++){
      if (A[i][j] != 0){
      	AddLinearSolver(&sky,i,j,A[i][j]);
      }
      /* if (i==j){ */
      /* 	SetLinearSolver(&sky,i,j,2); */
      /* } */
    }
  }

  sky.rhs=vf;
  sky.sol=sol;
  SolveLinearSolver(&sky,&simu);


  // checking
  verr=0;
  printf("sol of paralution=");
  for(int i=0;i<_NN;i++){
    printf("%f ",sol[i]);
    verr+=fabs(sol[i]-i-1);
  }
  printf("\n");
  printf("\n");

  // deallocate memory
  FreeLinearSolver(&sky);
  

  test2= (verr<1e-6);
#endif

  InitLinearSolver(&sky,_NN,NULL,NULL);

  sky.solver_type = GMRES;
  sky.pc_type=NONE;
  
  // now test a symmetric matrix
  A[0][0] = 0.3e1;
  A[0][1] = -0.1e1;
  A[0][2] = 0;
  A[0][3] = 0;
  A[0][4] = -0.1e1;
  A[1][0] = -0.1e1;
  A[1][1] = 0.2e1;
  A[1][2] = -0.1e1;
  A[1][3] = 0;
  A[1][4] = 0;
  A[2][0] = 0;
  A[2][1] = -0.1e1;
  A[2][2] = 0.2e1;
  A[2][3] = -0.1e1;
  A[2][4] = 0;
  A[3][0] = 0;
  A[3][1] = 0;
  A[3][2] = -0.1e1;
  A[3][3] = 0.2e1;
  A[3][4] = -0.1e1;
  A[4][0] = -0.1e1;
  A[4][1] = 0;
  A[4][2] = 0;
  A[4][3] = -0.1e1;
  A[4][4] = 0.3e1;
  vf[0] = -0.4e1;
  vf[1] = 0;
  vf[2] = 0;
  vf[3] = 0;
  vf[4] = 0.10e2;

  sol[0] = 0;
  sol[1] = 0;
  sol[2] = 0;
  sol[3] = 0;
  sol[4] = 0;

  sky.is_sym=false;

  // first mark the nonzero values in A
  for(int i=0;i<_NN;i++){
    for(int j=0;j<_NN;j++){
      if (A[i][j] != 0) IsNonZero(&sky,i,j);
    }
  }

  // once the nonzero positions are known allocate memory
  AllocateLinearSolver(&sky);

  // now set the nonzero terms
  for(int i=0;i<_NN;i++){
    for(int j=0;j<_NN;j++){
      if (A[i][j] != 0){
      	AddLinearSolver(&sky,i,j,A[i][j]);
      }
    }
  }

  sky.rhs=vf;
  sky.sol=sol;
  SolveLinearSolver(&sky,&simu);


  // checking
  verr=0;
  printf("sol of gmres=");
  for(int i=0;i<_NN;i++){
    printf("%f ",sol[i]);
    verr+=fabs(sol[i]-i-1);
  }
  printf("\n");
  printf("\n");

  // deallocate memory
  FreeLinearSolver(&sky,0);

  test3 =  (verr<1e-6);


  int NPoisson=40;
  real h=1.0/NPoisson;
  

  InitLinearSolver(&sky,NPoisson,NULL,NULL);

  sky.solver_type = GMRES;
  sky.pc_type=NONE;
  sky.iter_max=80;
  sky.tol=1.e-7;
  
  for(int i=0;i<NPoisson;i++){
    if (i==0){ 
      IsNonZero(&sky,0,0);
      IsNonZero(&sky,0,1);
    } 
    else if (i==NPoisson-1){ 
      IsNonZero(&sky,NPoisson-1,NPoisson-1);
      IsNonZero(&sky,NPoisson-1,NPoisson-2);
    } 
    else { 
    IsNonZero(&sky,i,i);
    IsNonZero(&sky,i,i+1);
    IsNonZero(&sky,i,i-1);
    } 
  }
  // once the nonzero positions are known allocate memory
  AllocateLinearSolver(&sky);

  for(int i=0;i<NPoisson;i++){
    if (i==0){ 
      AddLinearSolver(&sky,0,0,2.0/(h*h));
      AddLinearSolver(&sky,0,1,-1.0/(h*h));
    } 
    else if (i==NPoisson-1){ 
      AddLinearSolver(&sky,NPoisson-1,NPoisson-1,2.0/(h*h));
      AddLinearSolver(&sky,NPoisson-1,NPoisson-2,-1.0/(h*h));
    } 
    else { 
      AddLinearSolver(&sky,i,i,2.0/(h*h));
      AddLinearSolver(&sky,i,i+1,-1.0/(h*h));
      AddLinearSolver(&sky,i,i-1,-1.0/(h*h));
    } 
    sky.sol[i]=0.0;
    sky.rhs[i]=2.0;
  }

  real bigval=1.e+15;
  AddLinearSolver(&sky,NPoisson-1,NPoisson-1,bigval);
  AddLinearSolver(&sky,0,0,bigval);
 
  SolveLinearSolver(&sky,&simu);


  // checking
  verr=0;
  printf("sol of laplacien with gmres=");
  for(int i=0;i<NPoisson;i++){
    printf("%f ",sky.sol[i]);
    verr+=fabs(sky.sol[i]-(i*h)*(1-i*h));
  }
  printf("\n");
  printf("\n");

  // deallocate memory
  FreeLinearSolver(&sky,0);

  test4 =  (verr<1e-5);



  if(test1==1 &&  test2==1 && test3==1 && test4==1) test=1;

#ifdef PARALUTION 
  paralution_end();
#endif 

  return test;

}
