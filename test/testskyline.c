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

  A[0][0] = 2;
  A[0][1] = -0.5000000000e0;
  A[0][2] = 0;
  A[0][3] = 0;
  A[0][4] = -0.5000000000e0;
  A[1][0] = -1;
  A[1][1] = 2;
  A[1][2] = -1;
  A[1][3] = 0;
  A[1][4] = 0;
  A[2][0] = 0;
  A[2][1] = -1;
  A[2][2] = 2;
  A[2][3] = -1;
  A[2][4] = 0;
  A[3][0] = 0;
  A[3][1] = 0;
  A[3][2] = -1;
  A[3][3] = 2;
  A[3][4] = -1;
  A[4][0] = 0;
  A[4][1] = 0;
  A[4][2] = 0;
  A[4][3] = -1;
  A[4][4] = 2;


  for(int i=0;i<_NN;i++){
    for(int j=0;j<_NN;j++){
      if (A[i][j] != 0) SwitchOn(&sky,i,j);
    }
  }

  AllocateSkyline(&sky);

  for(int i=0;i<_NN;i++){
    for(int j=0;j<_NN;j++){
      if (A[i][j] != 0){
	SetSkyline(&sky,i,j,A[i][j]);
      }
    }
  }


  //FactoLU(&sky);

  DisplaySkyline(&sky);

  FreeSkyline(&sky);


  return test;

}
