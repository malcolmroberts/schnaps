#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <assert.h>
#include "klu_csr.h"
#include "math.h"



void InitKLU(KLU* klu, int n){

  klu->is_alloc=false;
  klu->copy_is_alloc=false;
  klu->is_lu=false;

  klu->neq=n;
  
  klu_defaults (&klu->common) ;

  klu->T = cs_di_spalloc (0, 0, 1, 1, 1) ;
  
  klu->symbolic=NULL;
  
  klu->numeric=NULL;
}


void SwitchOnKLU(KLU* klu,int i,int j){

  assert(!klu->is_alloc);
  cs_di_entry(klu->T, i, j, 1);

} 

void AllocateKLU(KLU* klu){

  assert(!klu->is_alloc);

  klu->A = cs_di_compress(klu->T );
  cs_di_spfree(klu->T);

  klu->is_alloc=true;

  cs_di_dupl(klu->A);
  //cs_di_print(klu->A, 0);
  // check that the matrix is square and its size conforms to the prescribed one
  assert (klu->A->m == klu->neq);
  assert (klu->A->n == klu->neq);
  // symbolic factorization
  klu->symbolic = klu_analyze(klu->neq,klu->A->p,klu->A->i,&(klu->common));
  assert(klu->symbolic);
  //
}

void AllocateCopyKLU(KLU* klu){

  assert(!klu->copy_is_alloc);
  klu->Acopy = cs_di_add (klu->A, klu->A, 1, 0) ;

  klu->copy_is_alloc=true;
  //cs_di_print(klu->A, 0);


}

void AddKLU(KLU* klu,int i,int j,schnaps_real val){

  assert(klu->is_alloc);  

  for (int iloc=klu->A->p[j];iloc < klu->A->p[j+1];iloc++){
    if (klu->A->i[iloc]==i){
      klu->A->x[iloc]+=val;
      return;
    } 
  }
}

void SetKLU(KLU* klu,int i,int j,schnaps_real val){

  assert(klu->is_alloc);

  for (int iloc=klu->A->p[j];iloc < klu->A->p[j+1];iloc++){
    if (klu->A->i[iloc]==i){
      klu->A->x[iloc]=val;
      return;
    } 
  }
} 

schnaps_real GetKLU(KLU* klu,int i,int j){

  assert(klu->is_alloc);
  for (int iloc=klu->A->p[j];iloc < klu->A->p[j+1];iloc++){
    if (klu->A->i[iloc]==i){
      return klu->A->x[iloc];
    } 
  }
  return 0.0;
}

void DisplayKLU(KLU* klu){

  int n=klu->neq;
  printf(" Original ordering\n");
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
  printf(" KLU reordering\n");  
#ifdef _FULL
  printf("\n");
  printf("\n");
  printf("\n");
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      printf("%.3e ", GetKLU(klu,klu->symbolic->P[i],klu->symbolic->Q[j]));
    }   
    printf("\n");
  }

#else
  printf("\n");
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      if (fabs(GetKLU(klu,klu->symbolic->P[i],klu->symbolic->Q[j])) < 1e-8) {
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
  
  
  klu->numeric = klu_factor (klu->A->p, klu->A->i, klu->A->x, klu->symbolic, &(klu->common));
  assert(klu->numeric > 0);
  klu->is_lu=true;
}

void ReFactoKLU(KLU* klu){
  
  assert(klu->is_lu);
  int res = klu_refactor (klu->A->p, klu->A->i, klu->A->x, klu->symbolic,klu->numeric, &(klu->common));
  assert(res);

}

void MatVectKLU(KLU * klu, schnaps_real * x, schnaps_real * prod) {

  for (int i=0; i< klu->neq;i++){
    prod[i]=0.0;
  }
  for (int j=0; j< klu->neq;j++){
    for(int iloc=klu->A->p[j];iloc < klu->A->p[j+1];iloc++){
      prod[klu->A->i[iloc]]+= klu->A->x[iloc] * x[j]; 
    }
  }
}


void SolveKLU(KLU* klu,schnaps_real* rhs,schnaps_real* sol){
  assert(klu->is_lu);
  // TODO? :evaluate  the need for an in place version avoiding the copy
  for (int i=0;i < klu->neq;i++){
    sol[i]=rhs[i];
  }
    klu_solve (klu->symbolic, klu->numeric, klu->neq, 1, sol, &(klu->common));
}




void FreeKLU(KLU* klu){

  assert(klu->is_alloc);
  if(klu->symbolic){
    klu_free_symbolic(&(klu->symbolic),&(klu->common));
    klu->symbolic = NULL;
  }
  if (klu->is_lu){
    klu_free_numeric(&(klu->numeric),&(klu->common));
    klu->is_lu=false;
    klu->numeric = NULL;
  }
  //
  cs_di_spfree(klu->A);
  //
  if (klu->copy_is_alloc){
    free(klu->Acopy);
  }
}

