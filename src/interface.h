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

  real *vnds;
  real *xpg;


} Interface;


void ExtractInterface(Interface* inter, int side);


void InterfaceExplicitFlux(Interface* inter, int side);
void ComputeFluxes(Interface* inter);

int VarindexFace(int npg, int m, int ipgf, int iv); 

#endif
