#ifndef _INTERFACE_H
#define _INTERFACE_H

#include "field.h"


typedef struct Interface{

  field *fL;
  field *fR;

  int locfaL,locfaR;

  int npgL, npgR;

  int* vol_indexL;
  int * vol_indexR;

  int wsizeL,wsizeR;

  real *wL;
  real *wR;


} Interface;


void ExtractInterface(Interface* inter, int side);

int VarindexFace(int npg, int m, int ipgf, int iv); 

#endif
