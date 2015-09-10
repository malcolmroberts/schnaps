#ifndef _INTERFACE_H
#define _INTERFACE_H

#include "field.h"


typedef struct Interface{

  field *fL;
  field *fR;

  int npgL, npgR;
  

  int wsizeL,wsizeR;

  real *wL;
  real *wR;


} Interface;


void ExtractInterface(Interface* inter, int side);


#endif
