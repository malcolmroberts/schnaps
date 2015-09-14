#include "interface.h"
#include <assert.h>
#include <stdio.h>

int VarindexFace(int npg, int m, int ipgf, int iv){

  return ipgf * m + iv;

}



void ExtractInterface(Interface* inter, int side){

  field *fd;
  int locfa,npgf;
  real *wf;

  if (side == 0) {
    fd = inter->fL;
    locfa = inter->locfaL;
    npgf = inter->npgL;
    wf = inter->wL;
  }

  if (side == 1) {
    fd = inter->fR;
    locfa = inter->locfaR;
    npgf = inter->npgR;
    wf = inter->wR;
 }
  
  assert(side ==1 || side == 0);

  if (fd !=NULL){
    for(int ipgf = 0; ipgf < npgf; ipgf++){
      int ipgv = ref_pg_face(fd->deg, fd->raf, locfa, ipgf,
			     NULL, NULL, NULL);
      for(int iv=0; iv < fd->model.m; iv++){
	int imem = fd->varindex(fd->deg, fd->raf, fd->model.m,
				ipgv, iv);
	int imemf = VarindexFace(npgf, fd->model.m, ipgf, iv);
	wf[imemf] = fd->wn[imem];
      }
    }
  }
    

}

