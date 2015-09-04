#ifdef _WITH_UMFPACK
#include "solverumfpack.h"
#endif
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "test.h"

int main(void) {
  // unit tests
  int resu=TestUmfPack();
  if (resu) printf("UmfPack test OK !\n");
  else printf("UmfPack test failed !\n");
  return !resu;
} 

int TestUmfPack()
{
  int test = 1;
#ifdef _WITH_UMFPACK
  // reference element
  smalltestumfpack();
#endif
  return test;
}
