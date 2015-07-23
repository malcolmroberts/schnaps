#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "test.h"
#include "schnaps.h"

#define _NN 3

int main(void) {
  
  // unit tests
    
  int resu=TestNonLinearSolver();
	 

  if (resu) printf("Linear solver test OK !\n");
  else printf("Linear solver test failed !\n");

  return !resu;
} 



void JFNonlinearVector_GX(void* system,field * f,real * solvector,real *nlvector){
  JFLinearSolver* lsol=system;

  nlvector[0] = 9*solvector[0]*solvector[0]+36*solvector[1]*solvector[1]+4*solvector[2]*solvector[2];
  nlvector[1] = solvector[0]*solvector[0]-2*solvector[1]*solvector[1]-20*solvector[2];
  nlvector[2] = solvector[0]*solvector[0]-solvector[1]*solvector[1]+solvector[2]*solvector[2];

}

void NonlinearVector_GX(void* system,field * f,real * solvector,real *nlvector){
  LinearSolver* lsol=system;

  nlvector[0] = 9*solvector[0]*solvector[0]+36*solvector[1]*solvector[1]+4*solvector[2]*solvector[2];
  nlvector[1] = solvector[0]*solvector[0]-2*solvector[1]*solvector[1]-20*solvector[2];
  nlvector[2] = solvector[0]*solvector[0]-solvector[1]*solvector[1]+solvector[2]*solvector[2];

}


int TestNonLinearSolver(void){

  int test=1,test1=1,test2=1,test3=1,test4=1;
  real A[_NN][_NN];
  real b[_NN],sold[_NN],soln[_NN],fsoln[_NN],rhs[_NN];
  double verr;

  LinearSolver sky;
  JFLinearSolver skyJF;
    Skyline * skymat;

  printf("/*****************************************************/ \n");
  printf("/*************** Newton + LU *************************/ \n");
  printf("/*****************************************************/ \n");
    
  InitLinearSolver(&sky,_NN,NULL,NULL);

  sky.solver_type = LU;
  sky.pc_type=NONE;
  sky.storage_type=SKYLINE;

  A[0][0] = 1.0e0;
  A[0][1] = 1.0e0;
  A[0][2] = 1.0e0;
  A[1][0] = 1.0e0;
  A[1][1] = 1.0e0;
  A[1][2] = 1.0e0;
  A[2][0] = 1.0e0;
  A[2][1] = 1.0e0;
  A[2][2] = 1.0e0;

  
  b[0] = 3.6e1;
  b[1] = 0;
  b[2] = 0;

  sold[0] = 0.0;
  sold[1] = 0.0;
  sold[2] = 0.0;

  rhs[0] = 0.0;
  rhs[1] = 0.0;
  rhs[2] = 0.0;
  
  soln[0] = 1;
  soln[1] = 1;
  soln[2] = 0;

  fsoln[0] = 0;
  fsoln[1] = 0;
  fsoln[2] = 0;
  
  sky.sol=sold;
  sky.rhs=rhs;

  for(int i=0;i<_NN;i++){
    for(int j=0;j<_NN;j++){
      if (A[i][j] != 0) IsNonZero(&sky,i,j);
      //if (i==j) SwitchOn(&sky,i,j);
    }
  }

  AllocateLinearSolver(&sky);
  
  for(int i=1;i<5;i++){

    NonlinearVector_GX(&sky,NULL,soln,fsoln);
    rhs[0] = b[0]-fsoln[0];
    rhs[1] = b[1]-fsoln[1];
    rhs[2] = b[2]-fsoln[2];
    

    A[0][0] = 1.8e1*soln[0];
    A[0][1] = 7.2e1*soln[1];
    A[0][2] = 0.8e1*soln[2];
    A[1][0] = 0.2e1*soln[0];
    A[1][1] = -0.4e1*soln[1];
    A[1][2] = -2.0e1;
    A[2][0] = 0.2e1*soln[0];
    A[2][1] = -0.2e1*soln[1];
    A[2][2] = 0.2e1*soln[2];

    
    // now set the nonzero terms
    for(int i=0;i<_NN;i++){
      for(int j=0;j<_NN;j++){
	if (A[i][j] != 0){
	  SetLinearSolver(&sky,i,j,A[i][j]);
	}
      }
    }

    skymat=(Skyline*)sky.matrix;
    skymat->is_lu=false;
  
    SolveLinearSolver(&sky);
    
    soln[0] = soln[0]+sky.sol[0];
    soln[1] = soln[1]+sky.sol[1];
    soln[2] = soln[2]+sky.sol[2];


    printf("linear solution : \n");
    printf(" Delta sol %f %f %f \n",sky.sol[0],sky.sol[1],sky.sol[2]);
    printf(" sol %f %f %f \n",soln[0],soln[1],soln[2]);
    printf("\n");
    printf("\n");

  }

  // checking
  verr=0;
  printf("nonlinear solution:");
  printf(" sol %f %f %f \n",soln[0],soln[1],soln[2]);
  verr=fabs(soln[0]-0.893628)+fabs(soln[1]-0.894527)+fabs(soln[2]+0.040083);
  printf("\n");
  printf("\n");

  // deallocate memory
  FreeLinearSolver(&sky);

  test1 =  (verr<1e-6);

  printf("/*****************************************************/ \n");
  printf("/*************** Newton + GMRES **********************/ \n");
  printf("/*****************************************************/ \n");

 InitLinearSolver(&sky,_NN,NULL,NULL);

  sky.solver_type = GMRES;
  sky.pc_type=NONE;
  sky.storage_type=SKYLINE;
  sky.restart_gmres=_NN;
  sky.tol=1.0e-8;
  sky.iter_max=20;

  A[0][0] = 1.0e0;
  A[0][1] = 1.0e0;
  A[0][2] = 1.0e0;
  A[1][0] = 1.0e0;
  A[1][1] = 1.0e0;
  A[1][2] = 1.0e0;
  A[2][0] = 1.0e0;
  A[2][1] = 1.0e0;
  A[2][2] = 1.0e0;

  
  b[0] = 3.6e1;
  b[1] = 0;
  b[2] = 0;

  sold[0] = 0.0;
  sold[1] = 0.0;
  sold[2] = 0.0;

  rhs[0] = 0.0;
  rhs[1] = 0.0;
  rhs[2] = 0.0;
  
  soln[0] = 1;
  soln[1] = 1;
  soln[2] = 0;

  fsoln[0] = 0;
  fsoln[1] = 0;
  fsoln[2] = 0;
  
  sky.sol=sold;
  sky.rhs=rhs;

  for(int i=0;i<_NN;i++){
    for(int j=0;j<_NN;j++){
      if (A[i][j] != 0) IsNonZero(&sky,i,j);
      //if (i==j) SwitchOn(&sky,i,j);
    }
  }

  AllocateLinearSolver(&sky);
  
  for(int i=1;i<5;i++){

    NonlinearVector_GX(&sky,NULL,soln,fsoln);
    rhs[0] = b[0]-fsoln[0];
    rhs[1] = b[1]-fsoln[1];
    rhs[2] = b[2]-fsoln[2];
    

    A[0][0] = 1.8e1*soln[0];
    A[0][1] = 7.2e1*soln[1];
    A[0][2] = 0.8e1*soln[2];
    A[1][0] = 0.2e1*soln[0];
    A[1][1] = -0.4e1*soln[1];
    A[1][2] = -2.0e1;
    A[2][0] = 0.2e1*soln[0];
    A[2][1] = -0.2e1*soln[1];
    A[2][2] = 0.2e1*soln[2];

    
    // now set the nonzero terms
    for(int i=0;i<_NN;i++){
      for(int j=0;j<_NN;j++){
	if (A[i][j] != 0){
	  SetLinearSolver(&sky,i,j,A[i][j]);
	}
      }
    }

    skymat=(Skyline*)sky.matrix;
    skymat->is_lu=false;
  
    SolveLinearSolver(&sky);
    
    soln[0] = soln[0]+sky.sol[0];
    soln[1] = soln[1]+sky.sol[1];
    soln[2] = soln[2]+sky.sol[2];


    printf("linear solution : \n");
    printf(" Delta sol %f %f %f \n",sky.sol[0],sky.sol[1],sky.sol[2]);
    printf(" sol %f %f %f \n",soln[0],soln[1],soln[2]);
    printf("\n");
    printf("\n");

  }

  // checking
  verr=0;
  printf("nonlinear solution:");
  printf(" sol %f %f %f \n",soln[0],soln[1],soln[2]);
  verr=fabs(soln[0]-0.893628)+fabs(soln[1]-0.894527)+fabs(soln[2]+0.040083);
  printf("\n");
  printf("\n");

  // deallocate memory
  FreeLinearSolver(&sky);

  test2 =  (verr<1e-6);


  printf("/*****************************************************/ \n");
  printf("/*************** Newton with Free Matrix + GMRES *****/ \n");
  printf("/*****************************************************/ \n");

  InitJFLinearSolver(&skyJF,_NN,NULL);

  skyJF.solver_type = GMRES;
  skyJF.pc_type=NONE;
  skyJF.eps=0.0001;
  skyJF.restart_gmres=_NN;
  skyJF.tol=1.0e-8;
  skyJF.iter_max=20;
  skyJF.NonlinearVector_computation=JFNonlinearVector_GX;
  
  
  b[0] = 3.6e1;
  b[1] = 0;
  b[2] = 0;

  sold[0] = 0.0;
  sold[1] = 0.0;
  sold[2] = 0.0;

  rhs[0] = 0.0;
  rhs[1] = 0.0;
  rhs[2] = 0.0;
  
  soln[0] = 1;
  soln[1] = 1;
  soln[2] = 0;

  fsoln[0] = 0;
  fsoln[1] = 0;
  fsoln[2] = 0;
  
  skyJF.sol=sold;
  skyJF.soln=soln;
  skyJF.rhs=rhs;
  
  for(int i=1;i<5;i++){
    JFNonlinearVector_GX(&skyJF,NULL,soln,fsoln);
    
    rhs[0] = b[0]-fsoln[0];
    rhs[1] = b[1]-fsoln[1];
    rhs[2] = b[2]-fsoln[2];
    
    SolveJFLinearSolver(&skyJF,NULL);
    
    soln[0] = soln[0]+skyJF.sol[0];
    soln[1] = soln[1]+skyJF.sol[1];
    soln[2] = soln[2]+skyJF.sol[2];

    printf("linear solution : \n");
    printf(" Delta sol %f %f %f \n",skyJF.sol[0],skyJF.sol[1],skyJF.sol[2]);
    printf(" sol %f %f %f \n",soln[0],soln[1],soln[2]);
    printf("\n");
    printf("\n");

  }

  // checking
  verr=0;
  printf("nonlinear solution:");
  printf(" sol %f %f %f \n",soln[0],soln[1],soln[2]);
  verr=fabs(soln[0]-0.893628)+fabs(soln[1]-0.894527)+fabs(soln[2]+0.040083);
  printf("\n");
  printf("\n");

  // deallocate memory
  FreeJFLinearSolver(&skyJF);

  test3 =  (verr<1e-6);
  

  if(test1==1 &&  test2==1 && test3==1) test=1;

  return test;

}