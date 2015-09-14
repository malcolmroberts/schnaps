#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "test.h"
#include "schnaps.h"

int main(void) 
{
  // Unit tests
  int resu=TestPIC();
  if (resu) printf("PIC test OK !\n");
  else printf("PIC test failed !\n");
  return !resu;
} 

void Maxwell2DConstInitData(real *x, real *w) 
{
  w[0]=0;
  w[1]=0;
  w[2]=1;
  w[3]=0;
  w[4]=0;
  w[5]=0;
  w[6]=0;
}

// some unit tests of the macromesh code
int TestPIC()
{
  MacroMesh m;

  bool test=true;

  field f;
  init_empty_field(&f);

  f.deg[0] = 3;  // x direction degree
  f.deg[1] = 3;  // y direction degree
  f.deg[2] = 0;  // z direction degree
  f.raf[0] = 4;  // x direction refinement
  f.raf[1] = 4;  // y direction refinement
  f.raf[2] = 1;  // z direction refinement

  // test gmsh file reading
  ReadMacroMesh(&f.macromesh, "../test/testmacromesh.msh");
  BuildConnectivity(&f.macromesh);
  CheckMacroMesh(&f.macromesh, f.deg, f.raf);
  //PrintMacroMesh(&m);

  PIC pic;

  InitPIC(&pic,1); 
  CreateParticles(&pic,&f.macromesh);
  PlotParticles(&pic,&f.macromesh);

  f.model.m = 7; // num of conservative variables

  /* f.model.NumFlux = Maxwell2DNumFlux; */
  /* f.model.BoundaryFlux = Maxwell2DBoundaryFlux; */
  f.model.InitData = Maxwell2DConstInitData;
  /* f.model.ImposedData = Maxwell2DImposedData; */
  f.varindex = GenericVarindex;
    

  Initfield(&f);

  // place the particle at (0,1,0) and v=(1,0,0)
  pic.xv[0]=0;
  pic.xv[1]=1;
  pic.xv[2]=0.5;
  real xref[3];
  pic.cell_id[0]=NumElemFromPoint(&f.macromesh,pic.xv,xref);
  pic.xv[0]=xref[0];  
  pic.xv[1]=xref[1];  
  pic.xv[2]=xref[2];  
  pic.xv[3]=1;
  pic.xv[4]=0;
  pic.xv[5]=0;

  f.pic = &pic;

  real final_pos_phy[3]={0,-1,0.5};
  real final_pos[3];
  int final_cell=NumElemFromPoint(&f.macromesh,
				  final_pos_phy,final_pos);

  pic.dt=0.001;
  for(int iter=0;iter<3141;iter++){
    PushParticles(&f, &pic);
  }
  PlotParticles(&pic,&f.macromesh);
 
  printf("Dist=%f\n",Dist(pic.xv,final_pos));

  test = test && (Dist(pic.xv,final_pos) < 1e-3);

  test = test && (final_cell == pic.cell_id[0]);

  AccumulateParticles(&f,f.wn);

  return test;
}
