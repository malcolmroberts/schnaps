#include "schnaps.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "../test/test.h"
#include "global.h"
#include "lattice.h"

int TestLattice_isothermal(void);

void Relaxation(void* s);
void Moments(void * s);
void Dummy_InitData(schnaps_real x[3],schnaps_real w[]);
void LinearWave_InitData(schnaps_real x[3],schnaps_real w[]);
void DoubleShear_InitData(schnaps_real x[3], schnaps_real w[]);
void Vorticity_2D_computation(Simulation *simu, int ifield_ux, int ifield_uy);
void Plot_Intermediate_fields(void *s,schnaps_real *w);
void TestDerivative(Simulation *simu,int nbfield);
void Compute_and_dump_Vorticity_2d(Simulation *simu,char *filename);
void PlotGenericFieldAsciiSparse(Simulation* simu,schnaps_real *w_in,char* fieldname,char *filename);
void PlotGenericFieldBinSparse(Simulation* simu,schnaps_real *w_in,char* fieldname,char *filename);
int main(void) {
  
  // unit tests
    
  int resu=TestLattice_isothermal();
	 
  if (resu) printf("lattice test OK !\n");
  else printf("lattice test failed !\n");

  return !resu;
} 
//
//

int TestLattice_isothermal(void) {
  
  int test=0;
  int vec=1;
  schnaps_real k=0.5;
  schnaps_real pi=4.0*atan(1.0);

  MacroMesh mesh;
  ReadMacroMesh(&mesh,"../test/testcube.msh");
  Detect2DMacroMesh(&mesh);
  //
  schnaps_real A[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  schnaps_real x0[3] = {0, 0, 0};
  AffineMapMacroMesh(&mesh,A,x0);

  // mesh preparation
  mesh.period[0]=1.0;
  mesh.period[1]=1.0;
  BuildConnectivity(&mesh);

  
  Model model;
  LatticeData * ld=&schnaps_lattice_data;  
  ld->feq=&feq_isothermal_D2Q9;
  //schnaps_real csound= 1.0/sqrt(3.0);
  schnaps_real csound= 1.0;
  InitLatticeData(ld,2,9,0,csound);
  
  model.m=ld->index_max; // num of conservative variables f(vi) for each vi, phi, E, rho, u, p, e (ou T)
  //
  model.NumFlux=Lattice_NumFlux;
  //
  //model.InitData = Dummy_InitData;
  model.InitData = DoubleShear_InitData;
  //model.InitData = LinearWave_InitData;

  model.ImposedData = NULL;
  model.BoundaryFlux = NULL;
  model.Source = NULL;

  
  int deg[]={4, 4, 0};
  int raf[]={32,32, 1};

  CheckMacroMesh(&mesh, deg, raf);
  Simulation simu;
  EmptySimulation(&simu);

  InitSimulation(&simu, &mesh, deg, raf, &model);
  //Moments(&simu);
  //simu.vmax = 2*ld->c; 
  simu.vmax = ld->c; 
  simu.cfl=1.0;
  simu.nb_diags = 1;
  simu.pre_dtfields = Relaxation;
  simu.post_dtfields = Moments;
  simu.update_after_rk = Plot_Intermediate_fields;
  //simu.update_after_rk = NULL;
  schnaps_real tmax = 25.0;
  //
  RK2(&simu, tmax);
  //
/*  Compute_and_dump_Vorticity_2d(&simu,"dgvisu_vort.msh");*/
/*  PlotFieldsAsciiSparse(ld->index_rho, false, &simu, "rho","dgvisu_rho.msh");*/
/*  PlotFieldsAsciiSparse(ld->index_ux, false, &simu, "ux","dgvisu_ux.msh");*/
/*  PlotFieldsAsciiSparse(ld->index_uy, false, &simu, "uy","dgvisu_uy.msh");*/
  //
  Moments(&simu);
  PlotFieldsBinSparse(ld->index_rho, false, &simu, "rho","dgvisu_rho_fin.msh");
/*  PlotFields(ld->index_rho, false, &simu, "rho","dgvisu_rho_test.msh");*/
  Compute_and_dump_Vorticity_2d(&simu,"dgvisu_vort_fin.msh");
/*  PlotFieldsBinSparse(1, false, &simu, "f1","dgvisu_f1_bin.msh");*/
/*  PlotFieldsBinSparse(2, false, &simu, "f2","dgvisu_f1_bin.msh");*/
  PlotFieldsBinSparse(ld->index_ux, false, &simu, "ux","dgvisu_ux_fin.msh");
  PlotFieldsBinSparse(ld->index_uy, false, &simu, "uy","dgvisu_uy_fin.msh");
  //
  test= 1;
  return test; 
}
void Dummy_InitData(schnaps_real x[3],schnaps_real w[])
{
  LatticeData * ld=&schnaps_lattice_data;
  schnaps_real my_pi= 4.0*atan(1.0);
  //
  ld->tau=1.0/2000.0;
  //
  schnaps_real uref=1.0;
  schnaps_real lam= 20.0;
  schnaps_real amp= 0.2;
  //
  schnaps_real tx= x[0]-0.5;
  schnaps_real ty= x[1]-0.5;
  schnaps_real rho = 1.0 + amp*exp(-lam*(tx*tx + ty*ty));
  schnaps_real ux =uref;
  schnaps_real uy =uref;
  schnaps_real uz =0.0;
  schnaps_real temp =(ld->c * ld-> c)/3.0;
  schnaps_real p =1.0;
  //
  //
  w[ld->index_rho]=rho;
  w[ld->index_ux]  = ux;
  w[ld->index_uy]  = uy;
  w[ld->index_uz]  = uz;
  w[ld->index_temp] = temp;
  w[ld->index_p] = p;
  for(int i=0;i<ld->index_max_q+1;i++){
    w[i]= ld->feq(i,ld,rho,ux,uy,uz,temp,p);
  }
};


void DoubleShear_InitData(schnaps_real x[3],schnaps_real w[])
{
  LatticeData * ld=&schnaps_lattice_data;
  schnaps_real my_pi= 4.0*atan(1.0);
  schnaps_real delta = 0.05, kappa=20.0;
  //
  ld->tau=1.0/10000.0;
  //
  schnaps_real rho = 1.0;
  schnaps_real ux =0.0;
  schnaps_real uy =0.0;
  schnaps_real uz =0.0;
  schnaps_real temp =(ld->c * ld->c)/3.0;
  schnaps_real p =1.0;
  //
  //printf(" c0 :%f \t temp:%f\n",ld->c, temp);
  schnaps_real uref=0.04;
  if(x[1]<0.5){
    ux=uref * tanh(kappa*(x[1]-0.25));
  }
  else{
    ux=uref * tanh(kappa*(0.75-x[1]));
  }    
  uy = uref * delta * sin(2.0 * my_pi*(x[0]+0.25));
  //
  w[ld->index_rho]=rho;
  w[ld->index_ux]  = ux;
  w[ld->index_uy]  = uy;
  w[ld->index_uz]  = uz;
  w[ld->index_temp] = temp;
  w[ld->index_p] = p;
  for(int i=0;i<ld->index_max_q+1;i++){
    w[i]= ld->feq(i,ld,rho,ux,uy,uz,temp,p);
  }
};
void LinearWave_InitData(schnaps_real x[3],schnaps_real w[])
{
  LatticeData * ld=&schnaps_lattice_data;
  schnaps_real my_pi= 4.0*atan(1.0);
  schnaps_real kx = 1.0, ky= 1.0;
  schnaps_real uref=0.1;
  schnaps_real pert=0.1;
  //
  ld->tau=1.0/10000.0;
  //
  schnaps_real phix= 2.0 * my_pi * kx * x[0];
  schnaps_real phiy= 2.0 * my_pi * ky * x[1];
  //
  schnaps_real rho = 1.0 + pert * cos(phix);
  schnaps_real ux =uref * pert * sin(phiy);
  schnaps_real uy =0.0;
  schnaps_real uz =0.0;
  schnaps_real temp =(ld->c * ld-> c)/3.0;
  schnaps_real p =1.0;
  //
  w[ld->index_rho]=rho;
  w[ld->index_ux]  = ux;
  w[ld->index_uy]  = uy;
  w[ld->index_uz]  = uz;
  w[ld->index_temp] = temp;
  w[ld->index_p] = p;
  for(int i=0;i<ld->index_max_q+1;i++){
    w[i]= ld->feq(i,ld,rho,ux,uy,uz,temp,p);
  }
};

void Equilibrium_VelocityPerturbation_BoundaryFlux(schnaps_real x[3],schnaps_real t,schnaps_real wL[],schnaps_real* vnorm,
				       schnaps_real* flux){
  
  LatticeData * ld=&schnaps_lattice_data;

};


void Relaxation(void* s){
  Simulation * simu =s;
  LatticeData * ld=&schnaps_lattice_data;
  schnaps_real w_eq[simu->wsize];
  Compute_distribution_eq(simu,w_eq);
  Compute_relaxation(simu,w_eq);
}

void Moments(void* s){
  Simulation * simu =s;
  
  Compute_moments(simu);
}
/**************************************************************************************/
void Plot_Intermediate_fields(void *s,schnaps_real *w){
  Simulation *simu=s;
  LatticeData *ld=&schnaps_lattice_data;
  int diagperiod=500;
  int it=simu->iter_time_rk;
  if ((it)%diagperiod ==0){
    printf("Dumping fields at it=%i\n",it);
  char rho_filename[sizeof "dgvisu_rho_00000000.msh"];
  sprintf(rho_filename,"dgvisu_rho_%08d.msh",it);
  char ux_filename[sizeof "dgvisu_ux_00000000.gmsh"];
  sprintf(ux_filename,"dgvisu_ux_%08d.msh",it);
  char uy_filename[sizeof "dgvisu_uy_00000000.gmsh"];
  sprintf(uy_filename,"dgvisu_uy_%08d.msh",it);
  char vort_filename[sizeof "dgvisu_vort_00000000.msh"];
  sprintf(vort_filename,"dgvisu_vort_%08d.msh",it);
  //
/*  PlotFieldsAsciiSparse(ld->index_rho, false, simu, "rho",rho_filename);*/
/*  PlotFieldsAsciiSparse(ld->index_ux, false, simu, "ux",ux_filename);*/
/*  PlotFieldsAsciiSparse(ld->index_uy, false, simu, "uy",uy_filename);*/
  //
  PlotFieldsBinSparse(ld->index_rho, false, simu, "rho",rho_filename);
  PlotFieldsBinSparse(ld->index_ux, false, simu, "ux",ux_filename);
  PlotFieldsBinSparse(ld->index_uy, false, simu, "uy",uy_filename);
  Compute_and_dump_Vorticity_2d(simu,vort_filename);
  };
};
/**************************************************************************************/
void TestDerivative(Simulation *simu,int nbfield){
    field *f = simu->fd;
    int nb_dof=simu->wsize/f->model.m;
    int wdsize= 3 * nb_dof;
    schnaps_real *wd;
    printf(" WD size %i", wdsize);
    wd= calloc(wdsize,sizeof(schnaps_real));
    for (int i=0; i< wdsize; i++){
      wd[i]= 42.0;
    }
    Compute_derivative(simu,wd,nbfield);
    //
    //
    for (int i=0; i< nb_dof; i++){
      int i1= i;
      int i2= i+ nb_dof;
      int i3= i+ 2* nb_dof;
      printf("i %i \t gradient : %f \t %f \t %f\n", i, wd[i1],wd[i2],wd[i3]);
    }
    //
    if (wd != NULL){
      free(wd);
    }
}
/**************************************************************************************/
void Vorticity_2D_computation(Simulation *simu, int ifield_ux, int ifield_uy){
    field *f = simu->fd;
    int nb_dof=simu->wsize/f->model.m;
    int wdsize= 3 * nb_dof;
    schnaps_real *temp_gradient;
    schnaps_real *vort;
    vort= calloc(nb_dof,sizeof(schnaps_real)); 
    temp_gradient= calloc(wdsize,sizeof(schnaps_real));
    //
    Compute_derivative(simu,temp_gradient,ifield_uy);
    //
    for (int i=0; i< nb_dof; i++){
      vort[i]= temp_gradient[i];
    };
    Compute_derivative(simu,temp_gradient,ifield_ux);
    for (int i=0; i< nb_dof; i++){
      vort[i]= vort[i]-temp_gradient[i+nb_dof];
    }; 
    //
    if (temp_gradient != NULL){
      free(temp_gradient);
    }
    if (vort != NULL){
      free(vort);
    }
}
/**************************************************************************************/
void Compute_and_dump_Vorticity_2d(Simulation *simu,char *filename){
    LatticeData *ld=&schnaps_lattice_data;
    field *f = simu->fd;
    int ifield_ux=ld->index_ux;
    int ifield_uy=ld->index_uy;
    int nb_dof=simu->wsize/f->model.m;
    int wdsize= 3 * nb_dof;
    schnaps_real *temp_gradient;
    schnaps_real *vort;
    vort= calloc(nb_dof,sizeof(schnaps_real)); 
    temp_gradient= calloc(wdsize,sizeof(schnaps_real));
    //
    Compute_derivative(simu,temp_gradient,ifield_uy);
    //
    for (int i=0; i< nb_dof; i++){
      vort[i]= temp_gradient[i];
    };
    Compute_derivative(simu,temp_gradient,ifield_ux);
    schnaps_real max_vort=-10000.0;
    for (int i=0; i< nb_dof; i++){
      vort[i]= vort[i]-temp_gradient[i+nb_dof];
      if (vort[i] > max_vort){
        max_vort=vort[i];
      }
    };
    printf("Max vort %f\n",max_vort); 
    //
    if (temp_gradient != NULL){
      free(temp_gradient);
    }
    //
    //PlotGenericFieldAsciiSparse(simu,vort,"vorticity",filename);
    PlotGenericFieldBinSparse(simu,vort,"vorticity",filename);
    //
    if (vort != NULL){
      free(vort);
    }
}
/**************************************************************************************/
////// * test routine for field plot with sparser interpolation
//////////////////////////////////////
void PlotGenericFieldAsciiSparse(Simulation* simu,schnaps_real *w_in,char* fieldname,char *filename) {
  schnaps_real hexa64ref[3 * 64] = {
    0, 0, 3, 3, 0, 3, 3, 3, 3, 0, 3, 3, 0, 0, 0, 3, 0, 0, 3, 3, 0, 0, 3, 0,
    1, 0, 3, 2, 0, 3, 0, 1, 3, 0, 2, 3, 0, 0, 2, 0, 0, 1, 3, 1, 3, 3, 2, 3,
    3, 0, 2, 3, 0, 1, 2, 3, 3, 1, 3, 3, 3, 3, 2, 3, 3, 1, 0, 3, 2, 0, 3, 1,
    1, 0, 0, 2, 0, 0, 0, 1, 0, 0, 2, 0, 3, 1, 0, 3, 2, 0, 2, 3, 0, 1, 3, 0,
    1, 1, 3, 1, 2, 3, 2, 2, 3, 2, 1, 3, 1, 0, 2, 2, 0, 2, 2, 0, 1, 1, 0, 1,
    0, 1, 2, 0, 1, 1, 0, 2, 1, 0, 2, 2, 3, 1, 2, 3, 2, 2, 3, 2, 1, 3, 1, 1,
    2, 3, 2, 1, 3, 2, 1, 3, 1, 2, 3, 1, 1, 1, 0, 2, 1, 0, 2, 2, 0, 1, 2, 0,
    1, 1, 2, 2, 1, 2, 2, 2, 2, 1, 2, 2, 1, 1, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1};
  for(int i = 0; i < 3 * 64; ++i)
    hexa64ref[i] /= 3.0;

  FILE * gmshfile;
  gmshfile = fopen(filename, "w" );

  //int param[8] = {f->model.m, _DEGX, _DEGY, _DEGZ, _RAFX, _RAFY, _RAFZ, 0};
  int nraf[3] = {simu->fd[0].raf[0],
    simu->fd[0].raf[1],
    simu->fd[0].raf[2]};
  int deg[3] = {simu->fd[0].deg[0],
    simu->fd[0].deg[1],
    simu->fd[0].deg[2]};
  const int npg[3] = {deg[0] + 1,
    deg[1] + 1,
    deg[2] + 1};
  const unsigned int sc_npg = npg[0] * npg[1] * npg[2];
  //
  // Refinement size in each direction
  schnaps_real hh[3] = {1.0 / nraf[0], 1.0 / nraf[1], 1.0 / nraf[2]};
  // Header
  fprintf(gmshfile, "$MeshFormat\n2.2 0 %d\n", (int) sizeof(schnaps_real));
  //
  int nb_plotnodes = simu->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2] * 64;
  fprintf(gmshfile, "$EndMeshFormat\n$Nodes\n%d\n", nb_plotnodes);
  //
  schnaps_real *value = calloc(nb_plotnodes , sizeof(schnaps_real));
  assert(value);
  int nodecount = 0;
  // Nodes
  //int npgv = NPG(deg, nraf);
  for(int i = 0; i < simu->macromesh.nbelems; i++) {
    // Get the nodes of element L
    field *f = simu->fd + i;
    // Loop on the macro elem subcells
    int icL[3];
    // Loop on the subcells
    for(icL[0] = 0; icL[0] < nraf[0]; icL[0]++) {
    for(icL[1] = 0; icL[1] < nraf[1]; icL[1]++) {
    for(icL[2] = 0; icL[2] < nraf[2]; icL[2]++) {
      for(int ino = 0; ino < 64; ino++) {
        schnaps_real Xr[3] = { hh[0] * (icL[0] + hexa64ref[3 * ino + 0]),
             hh[1] * (icL[1] + hexa64ref[3 * ino + 1]),
             hh[2] * (icL[2] + hexa64ref[3 * ino + 2]) };
        for(int ii = 0; ii < 3; ii++) {
          assert(Xr[ii] < 1 +  1e-10);
          assert(Xr[ii] > -1e-10);
        }
        schnaps_real Xphy[3];
        schnaps_ref2phy(simu->fd[i].physnode, Xr, NULL, -1, Xphy, NULL,  NULL, NULL, NULL);

        schnaps_real Xplot[3] = {Xphy[0], Xphy[1], Xphy[2]};
        schnaps_real testpsi = 0;
        ////////////////////////////////////////
        int jcL[3];
        int jcLmin[3]={0,0,0};
        int jcLmax[3]={nraf[0],nraf[1],nraf[2]};
        // restrict projection to neighbouring subcells
        for (int k=0; k<3;k++){
          if (icL[k] > 0){
            jcLmin[k] = icL[k]-1;
          }
          if (icL[k] < nraf[k]-1){
            jcLmax[k] = icL[k]+1;
          }
        };
        // loop on neighbouring subcells
        for (jcL[0]= jcLmin[0]; jcL[0] < jcLmax[0];jcL[0]++){
        for (jcL[1]= jcLmin[1]; jcL[1] < jcLmax[1];jcL[1]++){
        for (jcL[2]= jcLmin[2]; jcL[2] < jcLmax[2];jcL[2]++){
          // index cell
          int nc = jcL[0] + nraf[0] * (jcL[1] + nraf[1] * jcL[2]);
          // first glop index in the subcell
          int first_gp_cell = (deg[0]+1) * (deg[1]+1) * (deg[2]+1) * nc;
          for(int ipg=0;ipg<sc_npg;ipg++){
            int index_glob_igp=f->varindex(f->deg,f->raf,f->model.m, first_gp_cell + ipg,0);
              schnaps_real psi;
              psi_ref_subcell(f->deg, f->raf, icL, index_glob_igp, Xr, &psi, NULL);
              testpsi += psi;
              value[nodecount] += psi * w_in[index_glob_igp];
          }; // end loop subcell gauss points
        }; //end loop neighbour subcell 2
        };//end loop neighbour subcell 1
        }; //end loop neighbour subcell 0
        assert(fabs(testpsi-1) < _SMALL);
      ///////////////////////////////////////
      // Compare with an exact solution
      nodecount++;
      fprintf(gmshfile, "%d %f %f %f\n", nodecount,
        Xplot[0], Xplot[1], Xplot[2]);
      } // end loop Hex64 fe nodes
    } // end loop subcell 2
    } // end loop subcell 1
    } // end loop subcell 0
  } // end loop macrocells

  fprintf(gmshfile, "$EndNodes\n");
  // Elements
  fprintf(gmshfile, "$Elements\n");
  fprintf(gmshfile, "%d\n",
  simu->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2]);
  int elm_type = 92;
  int num_tags = 0;
  for(int i = 0; i < simu->macromesh.nbelems; i++) {
    // Loop on the subcells
    int icL[3];
    for(icL[0] = 0; icL[0] < nraf[0]; icL[0]++) {
    for(icL[1] = 0; icL[1] < nraf[1]; icL[1]++) {
  	for(icL[2] = 0; icL[2] < nraf[2]; icL[2]++) {
    // Get the subcell id
    int ncL = icL[0] + nraf[0] * (icL[1] + nraf[1] * icL[2]);
    // Global subcell id
    int numelem = ncL + i * nraf[0] * nraf[1] * nraf[2] + 1;
    fprintf(gmshfile, "%d ", numelem);
    fprintf(gmshfile, "%d ", elm_type);
    fprintf(gmshfile, "%d ", num_tags);
    for(int ii = 0; ii < 64; ii++) {
      int numnoe = 64 * (i * nraf[0] * nraf[1] * nraf[2] + ncL) + ii  + 1;
      fprintf(gmshfile, "%d ", numnoe);
    }
    fprintf(gmshfile, "\n");
    }
    }
    }
  } 
  fprintf(gmshfile, "$EndElements\n");
  // Now display data
  fprintf(gmshfile, "$NodeData\n");
  fprintf(gmshfile, "1\n");
  if(fieldname == NULL)
    fprintf(gmshfile, "\"field : anonymous\n");
  else
    fprintf(gmshfile, "\"field: %s\"\n", fieldname);
  //
  schnaps_real t = 0;
  fprintf(gmshfile, "1\n%f\n3\n0\n1\n", t);
  fprintf(gmshfile, "%d\n", nb_plotnodes);
  //
  for(int ino = 0; ino < nb_plotnodes; ino++) {
    fprintf(gmshfile, "%d %f\n", ino + 1, value[ino]);
  }
  fprintf(gmshfile, "\n$EndNodeData\n");
  //
  fclose(gmshfile);
  free(value);
}
void PlotGenericFieldBinSparse(Simulation* simu,schnaps_real *w_in,char* fieldname,char *filename){
  schnaps_real hexa64ref[3 * 64] = {
    0, 0, 3, 3, 0, 3, 3, 3, 3, 0, 3, 3, 0, 0, 0, 3, 0, 0, 3, 3, 0, 0, 3, 0,
    1, 0, 3, 2, 0, 3, 0, 1, 3, 0, 2, 3, 0, 0, 2, 0, 0, 1, 3, 1, 3, 3, 2, 3,
    3, 0, 2, 3, 0, 1, 2, 3, 3, 1, 3, 3, 3, 3, 2, 3, 3, 1, 0, 3, 2, 0, 3, 1,
    1, 0, 0, 2, 0, 0, 0, 1, 0, 0, 2, 0, 3, 1, 0, 3, 2, 0, 2, 3, 0, 1, 3, 0,
    1, 1, 3, 1, 2, 3, 2, 2, 3, 2, 1, 3, 1, 0, 2, 2, 0, 2, 2, 0, 1, 1, 0, 1,
    0, 1, 2, 0, 1, 1, 0, 2, 1, 0, 2, 2, 3, 1, 2, 3, 2, 2, 3, 2, 1, 3, 1, 1,
    2, 3, 2, 1, 3, 2, 1, 3, 1, 2, 3, 1, 1, 1, 0, 2, 1, 0, 2, 2, 0, 1, 2, 0,
    1, 1, 2, 2, 1, 2, 2, 2, 2, 1, 2, 2, 1, 1, 1, 2, 1, 1, 2, 2, 1, 1, 2, 1};
  for(int i = 0; i < 3 * 64; ++i)
    hexa64ref[i] /= 3.0;
  //
  FILE * gmshfile;
  gmshfile = fopen(filename, "w" );
  //int param[8] = {f->model.m, _DEGX, _DEGY, _DEGZ, _RAFX, _RAFY, _RAFZ, 0};
  int nraf[3] = {simu->fd[0].raf[0],
		 simu->fd[0].raf[1],
		 simu->fd[0].raf[2]};
  int deg[3] = {simu->fd[0].deg[0],
		simu->fd[0].deg[1],
		simu->fd[0].deg[2]};
		//
    const int npg[3] = {deg[0] + 1,
		      deg[1] + 1,
		      deg[2] + 1};
    const unsigned int sc_npg = npg[0] * npg[1] * npg[2];

  // Refinement size in each direction
  schnaps_real hh[3] = {1.0 / nraf[0], 1.0 / nraf[1], 1.0 / nraf[2]};
  // Header
  fprintf(gmshfile, "$MeshFormat\n2.2 1 %d\n", (int) sizeof(schnaps_real));
  int one = 1;
  fwrite(&one, sizeof(int), 1, gmshfile);
  fprintf(gmshfile, "$EndMeshFormat\n");
  //  End Header
  int nb_plotnodes = simu->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2] * 64;
  fprintf(gmshfile, "$Nodes\n%d\n", nb_plotnodes);
  //
  schnaps_real *value = calloc(nb_plotnodes , sizeof(schnaps_real));
  assert(value);
  int nodecount = 0;
  // Nodes
  for(int i = 0; i < simu->macromesh.nbelems; i++) {
    // Get the nodes of element L
    field *f = simu->fd + i;
    // Loop on the macro elem subcells
    int icL[3];
    // Loop on the subcells
    for(icL[0] = 0; icL[0] < nraf[0]; icL[0]++) {
    for(icL[1] = 0; icL[1] < nraf[1]; icL[1]++) {
    for(icL[2] = 0; icL[2] < nraf[2]; icL[2]++) {
    //
      for(int ino = 0; ino < 64; ino++) {
        schnaps_real Xr[3] = { hh[0] * (icL[0] + hexa64ref[3 * ino + 0]),
             hh[1] * (icL[1] + hexa64ref[3 * ino + 1]),
             hh[2] * (icL[2] + hexa64ref[3 * ino + 2]) };
        for(int ii = 0; ii < 3; ii++) {
          assert(Xr[ii] < 1 +  1e-10);
          assert(Xr[ii] > -1e-10);
        };
        schnaps_real Xphy[3];
        schnaps_ref2phy(simu->fd[i].physnode, Xr, NULL, -1, Xphy, NULL,  NULL, NULL, NULL);
        double Xplot[3] ={ (double)Xphy[0],  (double)Xphy[1],  (double)Xphy[2]};
        schnaps_real testpsi = 0;
        int nb_nonzero_psi=0;
        int jcL[3];
        int jcLmin[3]={0,0,0};
        int jcLmax[3]={nraf[0],nraf[1],nraf[2]};
        // restrict projectio to neighbouring subcells
        for (int k=0; k<3;k++){
          if (icL[k] > 0){
            jcLmin[k] = icL[k]-1;
          }
          if (icL[k] < nraf[k]-1){
            jcLmax[k] = icL[k]+1;
          }
        };
        // loop on neighbouring subcells
        for (jcL[0]= jcLmin[0]; jcL[0] < jcLmax[0];jcL[0]++){
        for (jcL[1]= jcLmin[1]; jcL[1] < jcLmax[1];jcL[1]++){
        for (jcL[2]= jcLmin[2]; jcL[2] < jcLmax[2];jcL[2]++){
          // index cell
          int nc = jcL[0] + nraf[0] * (jcL[1] + nraf[1] * jcL[2]);
          // first glop index in the subcell
          int first_gp_cell = (deg[0]+1) * (deg[1]+1) * (deg[2]+1) * nc;
          for(int ipg=0;ipg<sc_npg;ipg++){
            int index_glob_igp=f->varindex(f->deg,f->raf,f->model.m, first_gp_cell + ipg,0);
            schnaps_real psi;
            psi_ref_subcell(f->deg, f->raf, icL, index_glob_igp, Xr, &psi, NULL);
            testpsi += psi;
            value[nodecount] += psi * w_in[index_glob_igp];
          };
        };
        };
        };
        assert(fabs(testpsi-1) < _SMALL);
        nodecount++;
        fwrite(&nodecount, sizeof(int), 1, gmshfile);
        fwrite(Xplot, sizeof(double), 3, gmshfile);
      }
    };// end for icl[0]
    };// end for icl[1]
    };// end for icl[2]
  }// end for i macro elements
  fprintf(gmshfile, "$EndNodes\n");
  // Elements
  fprintf(gmshfile, "$Elements\n");
  int nb_elements=simu->macromesh.nbelems * nraf[0] * nraf[1] * nraf[2]; 
  fprintf(gmshfile, "%d\n",nb_elements);
  int elm_type = 92; //
  int num_tags = 0;
  int num_elm_follow= nb_elements;
  int elem_header[3] = {elm_type, num_elm_follow, num_tags};
  fwrite(elem_header, sizeof(int), 3, gmshfile);
  for(int i = 0; i < simu->macromesh.nbelems; i++) {
      // Loop on the subcells
      int icL[3];
      for(icL[0] = 0; icL[0] < nraf[0]; icL[0]++) {
      for(icL[1] = 0; icL[1] < nraf[1]; icL[1]++) {
      for(icL[2] = 0; icL[2] < nraf[2]; icL[2]++) {
      // Get the subcell id
      int ncL = icL[0] + nraf[0] * (icL[1] + nraf[1] * icL[2]);
      // Global subcell id
      int numelem = ncL + i * nraf[0] * nraf[1] * nraf[2] + 1;
      int elm_data_size=1+num_tags+64; 
      int elem_data[elm_data_size];
      elem_data[0]= numelem;
        for(int ii = 0; ii < 64; ii++) {
          int numnoe = 64 * (i * nraf[0] * nraf[1] * nraf[2] + ncL) + ii  + 1;
          elem_data[ii+1]= numnoe;
        };
      fwrite(elem_data, sizeof(int), elm_data_size, gmshfile);
      };
      };
      };
  };
  fprintf(gmshfile, "$\nEndElements\n");
  // Node Data
  fprintf(gmshfile, "$NodeData\n");
  fprintf(gmshfile, "1\n");
  if(fieldname == NULL)
    fprintf(gmshfile, "field : anonymous");
  else
  fprintf(gmshfile, "\"field: %s\"\n", fieldname);
  schnaps_real t = 0;
  fprintf(gmshfile, "1\n%f\n3\n0\n1\n", t);
  fprintf(gmshfile, "%d\n", nb_plotnodes);

  for(int ino = 1; ino < nb_plotnodes+1; ino++) {
    fwrite(&(ino),sizeof(int),1,gmshfile);
    fwrite(&value[ino-1],sizeof(double),1,gmshfile);
  };
  fprintf(gmshfile, "\n$EndNodeData\n");

  fclose(gmshfile);
  free(value);
} 
/////////////////////////////////////********
