#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "test.h"
#include "schnaps.h"

int main(void) {
  
  // unit tests
    
  int resu=TestSkyline();
	 

  if (resu) printf("Skyline test OK !\n");
  else printf("Skyline test failed !\n");

  return !resu;
} 




int TestSkyline(void){

  int test=1;

  Skyline sky;


#define _NN 5
  
  InitSkyline(&sky,_NN);

  double A[_NN][_NN];
  double vf[_NN],sol[_NN];

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




  for(int i=0;i<_NN;i++){
    for(int j=0;j<_NN;j++){
      if (A[i][j] != 0) SwitchOn(&sky,i,j);
      //if (i==j) SwitchOn(&sky,i,j);
    }
  }

  AllocateSkyline(&sky);

  for(int i=0;i<_NN;i++){
    for(int j=0;j<_NN;j++){
      if (A[i][j] != 0){
      	SetSkyline(&sky,i,j,A[i][j]);
      }
      /* if (i==j){ */
      /* 	SetSkyline(&sky,i,j,2); */
      /* } */
    }
  }


  FactoLU(&sky);

  DisplaySkyline(&sky);

  SolveSkyline(&sky,vf,sol);


  double verr=0;
  printf("sol=");
  for(int i=0;i<_NN;i++){
    printf("%f ",sol[i]);
    verr+=fabs(sol[i]-i-1);
  }
  printf("\n");

  FreeSkyline(&sky);

  test= (verr<1e-10);

  return test;

}
