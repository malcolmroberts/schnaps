#include "test.h"
#include "schnaps.h"
#include<stdio.h>
#include <assert.h>
#include <math.h>

int TestCache()
{
  int test = true;

  int deg[] = {2, 2, 0};
  int raf[] = {3, 3, 1};
  
  Model model;
  //model.cfl = 0.05;

  real cfl = 0.05;
  
  int mx = 2;
  int my = 1;

  /*   field f; */
  /*   init_empty_field(&f); */
  
  //f.varindex = GenericVarindex;
  model.vlasov_mz = 1;
  model.NumFlux = vlaTransNumFlux2d;
  model.vlasov_mx = mx;
  model.vlasov_my = my;
  
  model.BoundaryFlux = vlaTransBoundaryFlux2d;
  model.InitData = vlaTransInitData2d;
  model.ImposedData = vlaTransImposedData2d;

  real vmax = 1.0;
  Simulation simu;
  EmptySimulation(&simu);
  simu.vmax = vmax;
  simu.cfl = cfl;
  
  model.m = model.vlasov_mx * model.vlasov_my * model.vlasov_mz;

  MacroMesh mesh;
  ReadMacroMesh(&mesh, "../test/testmacromesh.msh");

  set_vlasov_params(&model); 

  // Try to detect a 2d mesh
  Detect2DMacroMesh(&mesh);
  assert(mesh.is2d);
  BuildConnectivity(&mesh);
  CheckMacroMesh(&mesh, deg, raf);

  InitSimulation(&simu, &mesh, deg, raf, &model);

  int ic[3]={0,0,0};
  int ix[3]={0,1,0};
  int ipg,ic0[3],ix0[3];

  for(int ii=0;ii<3;ii++){
    ic0[ii]=ic[ii];
    ix0[ii]=ix[ii];
  }

  xyz_to_ipg(raf,deg,ic,ix,&ipg);
  ipg_to_xyz(raf,deg,ic,ix,&ipg);

  for(int ii=0;ii<3;ii++){
    test = test && (ic0[ii]==ic[ii]);
    test = test && (ix0[ii]==ix[ii]);
  }

  for(int ii=0;ii<3;ii++){
    ic[ii]=raf[ii]-1;
    ix[ii]=deg[ii];
  }

  int worksize;  // number of gauss points in the macrocell

  xyz_to_ipg(raf,deg,ic,ix,&worksize);
  worksize++;
  
  printf("worksize=%d\n",worksize);

  int groupsize; // number of gauss points in a subcell

  
  int cache_size_in = (deg[0] + 3) * (deg[1] + 3) * (deg[2] + 3) * model.m;
  printf("cache_size_in=%d\n", cache_size_in);
  int cache_size_out = (deg[0] + 1) * (deg[1] + 1) * (deg[2] + 1) * model.m;
  printf("cache_size_out=%d\n", cache_size_out);
    
  int prec = 4;
  printf("memory usage=%d\n",(cache_size_out + cache_size_in) * prec);

  int group_id = 0;
  assert(group_id < worksize / groupsize);
  // kernel simulation
  for(int local_id = 0;  local_id < groupsize; ++local_id) {
    int global_id = local_id + group_id * groupsize;

    // get the central cell indices
    int ic0[3],ix0[3];
    ipg_to_xyz(raf, deg, ic0, ix0, &global_id);

    // prefetching
    int nfetch=(model.m * cache_size_in) / cache_size_out;
    printf("nfetch=%d\n",nfetch);
    
    for(int ifetch=0; ifetch < nfetch; ++ifetch) {

      // (local_id,ifetch)->(ic[3],ix[3],iw)-> imem (varindex)
      int cache_in=local_id + ifetch * groupsize;
      int iw = local_id % model.m;
      int ipg_in= cache_in / nfetch ;

      // (local_id,ifetch)-> imem2 (varindex_local)
      
      // wnloc[imem2]=wn[imem];
    }
    // (local_id,ifetch)->(ic[3],ix[3],iw)-> imem (varindex)
    // (local_id,ifetch)-> imem2 (varindex_local)
    
    // wnloc[imem2]=wn[imem]; if imem2 < cache_size_in;     
    

    // barrier

    // work on wnloc and dtwnloc ..........


    // barrier

    // back to global mem
    for(int ifetch = 0; ifetch < model.m; ++ifetch) {
      // (local_id,ifetch)->(ic[3],ix[3],iw)-> imem (varindex)
      // (local_id,ifetch)-> imem2 (varindex_local2)
      // dtwn[imem]=dtwnloc[imem2]
    }
    
  }

  return test;
}
  
int main(void) {
  // Unit tests
  int resu = TestCache();

  if (resu)
    printf("field test OK !\n");
  else
    printf("field test failed !\n");
  return !resu;
} 
