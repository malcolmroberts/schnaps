#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <assert.h>
#include "suitesparse/klu.h"
#include "klu_csr.h"
#include "math.h"



void InitKLU(KLU* klu, int n){

  klu->is_alloc=false;
  klu->copy_is_alloc=false;
  klu->is_lu=false;

  klu->neq=n;
  
  klu_defaults (&Common) ;
 

  
  
}


void SwitchOnKLU(KLU* klu,int i,int j){


} 

void AllocateKLU(KLU* klu){

  assert(!klu->is_alloc);


  klu->is_alloc=true;


}

void AllocateCopyKLU(KLU* klu){

  assert(!klu->copy_is_alloc);

  klu->copy_is_alloc=true;


}

void AddKLU(KLU* klu,int i,int j,schnaps_real val){

  assert(klu->is_alloc);  


}

void SetKLU(KLU* klu,int i,int j,schnaps_real val){

  assert(klu->is_alloc);
  
  
} 

schnaps_real GetKLU(KLU* klu,int i,int j){



}

void DisplayKLU(KLU* klu){

  int n=klu->neq;

  //#define _FULL

#ifdef _FULL
  printf("\n");
  printf("\n");
  printf("\n");
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      printf("%.3e ", GetKLU(klu,i,j));
    }   
    printf("\n");
  }

#else
  printf("\n");
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      if (fabs(GetKLU(klu,i,j)) < 1e-8) {
	printf(" ");
      } else {
	printf("*");
      }
    }   
    printf("\n");
  }
#endif
  
}



void FactoKLU(KLU* klu){

  klu->is_lu=true;

}

void MatVectKLU(KLU * klu, schnaps_real * x, schnaps_real * prod) {

  //assert(!klu->is_lu);

}



void SolveKLU(KLU* klu,schnaps_real* vfg,schnaps_real* vu){
  assert(klu->is_lu);

}




void FreeKLU(KLU* klu){

  assert(klu->is_alloc);

 
}
