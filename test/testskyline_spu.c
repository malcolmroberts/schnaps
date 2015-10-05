#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "test.h"
#include "schnaps.h"

int main(void) {
  
  // unit tests
    
  int resu=TestSkyline_SPU();
	 

  if (resu) printf("Skyline_SPU test OK !\n");
  else printf("Skyline_SPU test failed !\n");

  return !resu;
} 




int TestSkyline_SPU(void){

  int test=1;

  Skyline_SPU sky;


#define _NN 5
  
  // preliminary work on the skyline struct
  // _NN is the size of the linear system to be solved
  InitSkyline_SPU(&sky,_NN);

  real A[_NN][_NN];
  real vf[_NN];
  real sol[_NN]={1,2,3,4,5};
  real vf2[_NN]={0,0,0,0,6};

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



  // first mark the nonzero values in A
  for(int i=0;i<_NN;i++){
    for(int j=0;j<_NN;j++){
      if (A[i][j] != 0) SwitchOn_SPU(&sky,i,j);
      //if (i==j) SwitchOn(&sky,i,j);
    }
  }

  // once the nonzero positions are known allocate memory
  AllocateSkyline_SPU(&sky);

  // now set the nonzero terms
  for(int i=0;i<_NN;i++){
    for(int j=0;j<_NN;j++){
      if (A[i][j] != 0){
      	SetSkyline_SPU(&sky,i,j,A[i][j]);
      }
      /* if (i==j){ */
      /* 	SetSkyline_SPU(&sky,i,j,2); */
      /* } */
    }
  }


  real vf3[_NN];

  // test the product non symmetric case
  for(int i=0; i < _NN; i++){
    vf3[i]=0;
    for(int j=0; j< _NN; j++){
      vf3[i] += A[i][j] * sol[j];
    }
  }

  
  for(int i=0; i < _NN; i++) vf2[i]=i*i;
  sky.sol = sol;
  sky.rhs = vf2;

  MatVectSkyline_SPU(&sky);
  UnRegisterSkyline_SPU(&sky);
  for(int i=0; i < _NN; i++) {
    printf("%d vf=%f vf3=%f\n",i,vf2[i],vf3[i]);
    test = test && fabs(vf2[i]-vf3[i]) < _SMALL;
  }
  

 
  
  // LU decomposition
  FactoLU_SPU(&sky);

  // printf for checking...
  DisplaySkyline_SPU(&sky);

  // solve from a decomposed matrix
  // vf: rhs
  // sol: solution
  sky.rhs = vf;
  sky.sol = sol;
  SolveSkyline_SPU(&sky);
  UnRegisterSkyline_SPU(&sky);


  // checking
  real verr=0;
  printf("sol=");
  for(int i=0;i<_NN;i++){
    printf("%f ",sol[i]);
    verr+=fabs(sol[i]-i-1);
  }
  printf("\n");

  // deallocate memory
  FreeSkyline_SPU(&sky);
  

  test= test && (verr < _SMALL);
  //assert(1==2);

  InitSkyline_SPU(&sky,_NN);
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

  sky.is_sym=true;

  // first mark the nonzero values in A
  for(int i=0;i<_NN;i++){
    for(int j=0;j<_NN;j++){
      if ((A[i][j] != 0) && (j >= i)) SwitchOn_SPU(&sky,i,j);
      //if (i==j) SwitchOn(&sky,i,j);
    }
  }

  // once the nonzero positions are known allocate memory
  AllocateSkyline_SPU(&sky);

  // now set the nonzero terms
  for(int i=0;i<_NN;i++){
    for(int j=0;j<_NN;j++){
      if ((A[i][j] != 0) && (j >= i)){
      	SetSkyline_SPU(&sky,i,j,A[i][j]);
      }
      /* if (i==j){ */
      /* 	SetSkyline_SPU(&sky,i,j,2); */
      /* } */
    }
  }

  // test the product symmetric case
  for(int i=0; i < _NN; i++){
    vf3[i]=0;
    for(int j=0; j< _NN; j++){
      vf3[i] += A[i][j] * sol[j];
    }
  }
  
  for(int i=0; i < _NN; i++) vf2[i]=10;
  sky.sol = sol;
  sky.rhs = vf2;

  MatVectSkyline_SPU(&sky);
  UnRegisterSkyline_SPU(&sky);
  for(int i=0; i < _NN; i++) {
    printf("%d vf=%f vf3=%f\n",i,vf2[i],vf3[i]);
    test = test && fabs(vf2[i]-vf3[i]) < _SMALL;
  }

 
  // LU decomposition
  FactoLU_SPU(&sky);

  // printf for checking...
  DisplaySkyline_SPU(&sky);

  // solve from a decomposed matrix
  // vf: rhs
  // sol: solution
  sky.rhs = vf;
  sky.sol = sol;
  SolveSkyline_SPU(&sky);
  UnRegisterSkyline_SPU(&sky);


  // checking
  verr=0;
  printf("sol=");
  for(int i=0;i<_NN;i++){
    printf("%f ",sol[i]);
    verr+=fabs(sol[i]-i-1);
  }
  printf("\nerror=%f small=%f test=%d \n",verr,_SMALL, (verr < _SMALL));

  // deallocate memory
  //FreeSkyline_SPU(&sky);

  test= test && (verr < _SMALL);


  return test;

}
