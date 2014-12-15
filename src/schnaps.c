/* -*- eval: (folding-mode); -*- */
#include "schnaps.h"
#include <stdio.h>
#include <assert.h>
#include <string.h>

int main(int argc, char *argv[]){

  // default value
  char castest[20] = "transport";

  if(argc == 2){
    strcpy(castest, argv[1]);
  }
  else if(argc != 1){
    printf("Bad number argument !\n");
    assert(1==2);
  }

  if(strcmp(castest, "transport") == 0){

    Field f;
    f.model.m = 1; // only one conservative variable
    f.model.NumFlux = TransNumFlux2d;
    f.model.BoundaryFlux = TransBoundaryFlux2d;
    f.model.InitData = TransInitData2d;
    f.model.ImposedData = TransImposedData2d;
    f.varindex = GenericVarindex;

    f.interp.interp_param[0] = 1;  // _M
    f.interp.interp_param[1] = 3;  // x direction degree
    f.interp.interp_param[2] = 3;  // y direction degree
    f.interp.interp_param[3] = 0;  // z direction degree
    f.interp.interp_param[4] = 4;  // x direction refinement
    f.interp.interp_param[5] = 4;  // y direction refinement
    f.interp.interp_param[6] = 1;  // z direction refinement

    // Read the gmsh file
    ReadMacroMesh(&(f.macromesh), "disque.msh");
    //ReadMacroMesh(&(f.macromesh), "geo/cube.msh");
    // Try to detect a 2d mesh
    bool is2d = Detect2DMacroMesh(&(f.macromesh));
    assert(is2d);
    //PrintMacroMesh(&(f.macromesh));

    // Mesh preparation
    BuildConnectivity(&(f.macromesh));

    // AffineMapMacroMesh(&(f.macromesh));

    // Prepare the initial fields
    InitField(&f);
    f.is2d = true;

    // prudence...
    CheckMacroMesh(&(f.macromesh), f.interp.interp_param + 1);

    printf("cfl param =%f\n", f.hmin);

    // Apply the DG scheme time integration by RK2 scheme up to final
    // time = 1.
    RK2(&f, 1.0);

    // Save the results and the error
    PlotField(0, false, &f, "V", "dgvisu.msh");
    PlotField(0, true, &f, "Err", "dgerror.msh");

    double dd = L2error(&f);

    printf("erreur L2=%f\n", dd);
    return 0;
  }
  else if(strcmp(castest, "mhd") == 0){
    Field f;
    f.model.m=9;
    f.model.NumFlux=MHDNumFlux;
    f.model.BoundaryFlux=MHDBoundaryFlux;
    f.model.InitData=MHDInitData;
    f.model.ImposedData=MHDImposedData;
    f.varindex=GenericVarindex;


    f.interp.interp_param[0]=9;  // _M
    f.interp.interp_param[1]=1;  // x direction degree
    f.interp.interp_param[2]=1;  // y direction degree
    f.interp.interp_param[3]=0;  // z direction degree
    f.interp.interp_param[4]=2;  // x direction refinement
    f.interp.interp_param[5]=1000;  // y direction refinement
    f.interp.interp_param[6]=1;  // z direction refinement


    // read the gmsh file
    ReadMacroMesh(&(f.macromesh),"geo/cartesianGrid2d.msh");
    // try to detect a 2d mesh
    bool is2d=Detect2DMacroMesh(&(f.macromesh));
    assert(is2d);

    // mesh preparation
    BuildConnectivity(&(f.macromesh));

    //AffineMapMacroMesh(&(f.macromesh));

    // prepare the initial fields
    InitField(&f);
    f.is2d=true;

    // prudence...
    CheckMacroMesh(&(f.macromesh),f.interp.interp_param+1);

    printf("cfl param =%f\n",f.hmin);

    PlotField(0,(1==0),&f, "Rho", "dginit.msh");

    // apply the DG scheme
    // time integration by RK2 scheme
    // up to final time = 1.
    RK2(&f,1.0);

    Gnuplot(&f,1,0.0,"data1D.dat");
    // save the results and the error
    PlotField(0,(1==0),&f, "Rho", "dgvisu.msh");

    return 0;
  }
  else{
    printf("Bad argument !\n");
    assert(1==2);
  }
};
