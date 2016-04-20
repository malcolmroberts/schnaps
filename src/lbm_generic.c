#include "lbm_generic.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include <math.h>
#include "implicit.h"
//
/*****************************************************************************************/
/*****************************************************************************************/
// basic moments routines
/*****************************************************************************************/
/*****************************************************************************************/
void NewMPolyDescriptor(MPolyDescriptor *mpd, int ndim, int nc, char * name){
  DestroyMPolyDescriptor(mpd);
  mpd->ndim=ndim;
  mpd->nc=nc;
  mpd->c= (schnaps_real*) calloc(nc,sizeof(schnaps_real));
  mpd->e = (int **) calloc(mpd->nc,sizeof(int *));
  for (int i=0;i< mpd->nc;i++){
    mpd->e[i]= (int *) calloc(ndim,sizeof(int));
  };
  for  (int i=0;i< mpd->nc;i++){
    for (int k=0;k< mpd->ndim;k++){
      mpd->e[i][k]=1;
  };
  };
  if (name !=NULL){
  mpd->name=strdup(name);
  }
}
/*****************************************************************************************/
void SetMPolyDescriptor(MPolyDescriptor *mpd, schnaps_real c_in[mpd->nc],int e_in[mpd->nc][mpd->ndim]){
  for (int i=0;i< mpd->nc;i++){
    mpd->c[i]=c_in[i];
    for (int k=0; k< mpd->ndim;k++){
      mpd->e[i][k]=e_in[i][k];
    };
  };
  //
  char tmp_poly[LBM_MAX_POLYSTRSIZE];
  tmp_poly[0]='\0';
  for (int i=0;i<mpd->nc;i++){
    char tmp_str[LBM_MAX_POLYSTRSIZE];
    tmp_str[0]='\0';
    sprintf(tmp_str,"%f",mpd->c[i]);
    bool has_nonnull_exp=false;
    for (int k=0;k< mpd->ndim;k++){
      if (mpd->e[i][k] >0){
        has_nonnull_exp=true;
      } //endif
    } // k
    if (has_nonnull_exp){
      strcat(tmp_str," * ");
      for (int k=0;k< mpd->ndim;k++){
        if (mpd->e[i][k]>0){
        char varname[2]={LBM_POLYVARNAMES[k],'\0'};
        strcat(tmp_str,varname);
        char tmp_exp[LBM_MAX_POLYSTRSIZE];
        sprintf(tmp_exp,"^%i",mpd->e[i][k]);
        strcat(tmp_str,tmp_exp);
        } //endif
      } //k
    }// endif
    int l=strlen(tmp_poly)+strlen(tmp_str)+1;
    if (l < LBM_MAX_POLYSTRSIZE){
    strcat(tmp_poly,tmp_str);
    }
    l=strlen(tmp_poly)+4;
    if (l < LBM_MAX_POLYSTRSIZE){
    if (i< mpd->nc-1){
      strcat(tmp_poly," + ");
    }
    }
  } //i
  int len=strlen(tmp_poly);
  mpd->poly_string=(char*) malloc((len+1));
  strcpy(mpd->poly_string,tmp_poly);
}
/*****************************************************************************************/
void NewMPolyDescriptorFromModel(MPolyDescriptor *mpd,const MPolyData *model, char *name){
  if (name==NULL){
    NewMPolyDescriptor(mpd,model->ndim,model->nc,model->s);
  }
  else{
  NewMPolyDescriptor(mpd,model->ndim,model->nc,name);
  }
  SetMPolyDescriptor(mpd,model->c,model->e);
}

/*****************************************************************************************/
void DisplayMPolyDescriptor(MPolyDescriptor *mpd){
  if (mpd){
  if (mpd->poly_string){
  char *temp_str= calloc(LBM_MAX_POLYSTRSIZE+LBM_MAX_MACROVARNAMESIZE+4,sizeof(char));
  strcpy(temp_str,mpd->name);
  strcat(temp_str," : ");
  strcat(temp_str,mpd->poly_string);
  puts(temp_str);
  free(temp_str);
  }
  }
}
/*****************************************************************************************/
schnaps_real EvaluatePoly(MPolyDescriptor *mpd,schnaps_real *V){
  schnaps_real result=0.0;
  for (int i=0;i<mpd->nc;i++){
    schnaps_real temp_val=mpd->c[i];
    for (int k=0;k<mpd->ndim;k++){
    for (int q=0;q< mpd->e[i][k];q++){
      temp_val=temp_val*V[k];
    } //q
    } //k
    result=result+temp_val;
  }//i
  return result;
}
/*****************************************************************************************/
void DestroyMPolyDescriptor(MPolyDescriptor *mpd){
  if (mpd->c) free(mpd->c);
  if( mpd->e){
  for (int i=0;i< mpd->nc;i++){
    free(mpd->e[i]);
  };
  free(mpd->e);
  }
  if (mpd->name) free(mpd->name);
  if (mpd->poly_string) free(mpd->poly_string);
  mpd->ndim=0;
  mpd->nc=0;
}
//
/*****************************************************************************************/
/*****************************************************************************************/
// LBModelDescriptor routines
/*****************************************************************************************/
/*****************************************************************************************/
void NewLBModelDescriptor(LBModelDescriptor *lb,int d, int nb_macro,int q){
  DestroyLBModelDescriptor(lb);
  assert(d>0);
  lb->d=d;
  assert(q>0);
  lb->q=q;
  assert(nb_macro>0);
  lb->nb_macro=nb_macro;
  //
  lb->macro_names= (char**) calloc(lb->nb_macro,sizeof(char*));
  for (int i=0;i< lb->nb_macro;i++){
    lb->macro_names[i]= (char*) calloc(LBM_MAX_MACROVARNAMESIZE,sizeof(char));
  }
  //
  lb->vi=(schnaps_real**) calloc(lb->q,sizeof(schnaps_real*));
  //
  for (int i=0;i< lb->q;i++){
    lb->vi[i]= (schnaps_real*) calloc(lb->d,sizeof(schnaps_real));
  };
  //
  lb->Moments= (MPolyDescriptor*) calloc(lb->q,sizeof(MPolyDescriptor));
  for (int i=0;i< lb->q;i++){
    lb->Moments[i]=MPolyDescriptor_NULL;
  };
  //
  lb->inode_min= (int *) calloc(lb->q,sizeof(int));
  lb->inode_max= (int *) calloc(lb->q,sizeof(int));
  //
  lb->M = (schnaps_real**) calloc(lb->q,sizeof(schnaps_real*));
  for (int i=0;i< lb->q;i++){
    lb->M[i]= (schnaps_real*) calloc(lb->q,sizeof(schnaps_real));
  };
  //
  lb->s = (schnaps_real*) calloc(lb->q,sizeof(schnaps_real));
  //
/*  lb->is_relaxed= (bool*) calloc(lb->q,sizeof(schnaps_real));*/
/*  //*/
/*  lb->is_conserved= (bool*) calloc(lb->q,sizeof(schnaps_real));*/
  //
  
}
/*****************************************************************************************/
/*****************************************************************************************/
void DestroyLBModelDescriptor(LBModelDescriptor *lb){
  if (lb->macro_names){
    for (int i=0;i< lb->nb_macro;i++){
      free(lb->macro_names[i]);
    }
    free(lb->macro_names);
  }
  if (lb->vi){
    for (int i=0; i<lb->q ;i++){
      free(lb->vi[i]);
    }
    free(lb->vi);
  }
  if( lb->Moments){
    for (int i=0;i<lb->q;i++){
      DestroyMPolyDescriptor(&(lb->Moments[i]));
    }
    free(lb->Moments);
  }
  if (lb->M){
    for (int i=0;i<lb->q;i++){
      free(lb->M[i]);
    };
    free(lb->M);
  }
  if (lb->inode_min) free(lb->inode_min);
  if (lb->inode_max) free(lb->inode_max);
  if (lb->s) free(lb->s);
/*  if (lb->is_relaxed) free(lb->is_relaxed);*/
/*  if (lb->is_conserved) free(lb->is_conserved);*/
  lb->macro_names=NULL;
  lb->vi=NULL;
  lb->Moments=NULL;
  lb->M=NULL;
  lb->s=NULL;
/*  lb->is_conserved=NULL;*/
/*  lb->is_relaxed=NULL;*/
  lb->d=0;
  lb->q=0;
  lb->nb_macro=0;
/*  lb->nc=0;*/
/*  lb->nr=0;*/
  lb->f_to_macro=NULL;
  lb->macro_to_f=NULL;
  lb->feq=NULL;
  lb->meq=NULL;
}
/******************************************************************************************/
void ComputeLBModelDescriptorMomentMatrix(LBModelDescriptor *lb){
    for (int i=0;i< lb->q;i++){
      for (int j=lb->inode_min[i];j<=lb->inode_max[i];j++){
        lb->M[i][j]=EvaluatePoly(&(lb->Moments[i]),&(lb->vi[j][0]));
      }
    }
}
/*******************************************************************************************/
void DisplayLBModelDescriptorMomentMatrix(LBModelDescriptor *lb){
  for (int i=0;i< lb->q;i++){
    DisplayMPolyDescriptor(&(lb->Moments[i]));
  }
  for (int i=0;i< lb->q;i++){
    for (int j=0;j< lb->q;j++){
      printf("%f\t",lb->M[i][j]);
    }
    printf("\n");
  }
}
/*****************************************************************************************/
void LBM_Dummy_InitMacroData(schnaps_real x[3],schnaps_real w[]){
  LatticeBoltzmannSimData *lsd=&schnaps_lbm_simdata;
  for (int i=0;i<lsd->lb_model->nb_macro;i++){
    w[i]=0.0;
  }
}
//
void LBM_Dummy_InitData(schnaps_real x[3],schnaps_real w[]){
  LatticeBoltzmannSimData *lsd=&schnaps_lbm_simdata;
  for (int i=0;i< lsd->lb_model->q;i++){
    w[i]=0.0;
  }
}
/*****************************************************************************************/
void LBM_Dummy_InitData_OneNode(schnaps_real x[3],schnaps_real w[]){
  LatticeBoltzmannSimData *lsd=&schnaps_lbm_simdata;
    w[0]=0.0;
}
/*****************************************************************************************/
/*****************************************************************************************/
// Lattice Boltmann Simulation objects routines
/*****************************************************************************************/
/*****************************************************************************************/
void InitLBMSimulation( LBMSimulation *lbsimu,LatticeBoltzmannSimData *lsd, MacroMesh *mesh, int deg[3], int raf[3]){
  lbsimu->d=lsd->lb_model->d;
  lbsimu->nb_macro_fields=lsd->lb_model->nb_macro;
  lbsimu->q=lsd->lb_model->q;
  //
  lbsimu->macro_model.m=lsd->lb_model->nb_macro;
  lbsimu->macro_model.NumFlux=NULL;
  if (lbsimu->macro_model.InitData==NULL){
  lbsimu->macro_model.InitData = LBM_Dummy_InitMacroData;
  }
  //
  EmptySimulation(&(lbsimu->macro_simu));
  InitSimulation(&(lbsimu->macro_simu), mesh, deg, raf, &(lbsimu->macro_model));
  lbsimu->macro_simu.pre_dtfields = NULL;
  lbsimu->macro_simu.post_dtfields = NULL;
  lbsimu->macro_simu.update_after_rk = NULL;
  //
  lbsimu->micro_model.m=lsd->lb_model->q;
  lbsimu->micro_model.NumFlux=NULL;
  lbsimu->micro_model.InitData = LBM_Dummy_InitData;
  lbsimu->micro_model.ImposedData = NULL;
  lbsimu->micro_model.BoundaryFlux = NULL;
  lbsimu->micro_model.Source = NULL;
  EmptySimulation(&(lbsimu->micro_simu));
  InitSimulation(&(lbsimu->micro_simu), mesh, deg, raf, &(lbsimu->micro_model));
  lbsimu->micro_simu.pre_dtfields = LBM_pre_dtfields_wrapper;
  lbsimu->micro_simu.post_dtfields = LBM_post_dtfields_wrapper;
  lbsimu->micro_simu.update_after_rk = LBM_update_after_rk_wrapper;
  // set initial data to equilibrium one computed with initial macro data
  LB_Relaxation_bgk_f_full(lbsimu);
  lsd->current_lb_sim=lbsimu;
  //
  lbsimu->pre_advec=NULL;
  lbsimu->post_advec=NULL;
  lbsimu->post_tstep=NULL;
  //
  lbsimu->diag_2d_period=0.0;
  lbsimu->collect_diags=NULL;
}
void FreeBMSimulation(LBMSimulation *lbsimu){
  freeSimulation(&(lbsimu->macro_simu));
  freeSimulation(&(lbsimu->micro_simu));
}
//******************************************************************************//

void LB_Relaxation_bgk_f( void *lbs){
  LBMSimulation *lbsimu=lbs;
  LatticeBoltzmannSimData *lsd=&schnaps_lbm_simdata;
  Simulation *macsimu=&(lbsimu->macro_simu);
  Simulation *micsimu=&(lbsimu->micro_simu);
  assert(lsd->lb_model->feq);
  schnaps_real rate=lsd->lb_model->s[0]; // TODO use something simpler ?
    //
    for(int ie = 0; ie < macsimu->macromesh.nbelems; ++ie) {
      field * f = macsimu->fd+ie;
      field * fmic=micsimu->fd+ie;
      for(int ipg=0;ipg<NPG(f->deg, f-> raf);ipg++){
        // recover macro quantities;
        schnaps_real wmac[lsd->lb_model->nb_macro];
        for(int iv=0;iv<lsd->lb_model->nb_macro;iv++){
            int imac=f->varindex(f->deg, f->raf, f->model.m, ipg, iv);
            wmac[iv]=f->wn[imac];
        } //iv macro quantities
        // now compute equilibria and relax for all micro quantities
        for(int inode=0;inode<lsd->lb_model->q;inode++){
          schnaps_real feq=lsd->lb_model->feq(inode,lsd->lb_model->nb_macro,wmac);
            int imic=fmic->varindex(fmic->deg, fmic->raf, fmic->model.m, ipg, inode);
            fmic->wn[imic]=fmic->wn[imic]-rate*(fmic->wn[imic]-feq);
        }
      }//ipg
      //
    }//ie
}
void LB_Relaxation_bgk_f_full( void *lbs){
  LBMSimulation *lbsimu=lbs;
  LatticeBoltzmannSimData *lsd=&schnaps_lbm_simdata;
  assert(lsd->lb_model->feq);
  Simulation *macsimu=&(lbsimu->macro_simu);
  Simulation *micsimu=&(lbsimu->micro_simu);
    //
    for(int ie = 0; ie < macsimu->macromesh.nbelems; ++ie) {
      field * f = macsimu->fd+ie;
      field * fmic=micsimu->fd+ie;
      for(int ipg=0;ipg<NPG(f->deg, f-> raf);ipg++){
        // recover macro quantities;
        schnaps_real wmac[lsd->lb_model->nb_macro];
        for(int iv=0;iv<lsd->lb_model->nb_macro;iv++){
            int imac=f->varindex(f->deg, f->raf, f->model.m, ipg, iv);
            wmac[iv]=f->wn[imac];
        } //iv macro quantities
        // now compute equilibria and relax for all micro quantities
        for(int inode=0;inode<lsd->lb_model->q;inode++){
            int imic=fmic->varindex(fmic->deg, fmic->raf, fmic->model.m, ipg, inode);
            fmic->wn[imic]=lsd->lb_model->feq(inode,lsd->lb_model->nb_macro,wmac);
        }
      }//ipg
      //
    }//ie
}
/*//void LB_Relaxation_Moments( LBMSimulation *lbsimu);*/
void LB_ComputeMacroFromMicro(void *lbs){
  LBMSimulation *lbsimu=lbs;
  LatticeBoltzmannSimData *lsd=&schnaps_lbm_simdata;
  Simulation *macsimu=&(lbsimu->macro_simu);
  Simulation *micsimu=&(lbsimu->micro_simu);
  for(int ie = 0; ie < macsimu->macromesh.nbelems; ++ie) {
    field * f = macsimu->fd+ie;
    field * fmic=micsimu->fd+ie;
    for(int ipg=0;ipg<NPG(f->deg, f-> raf);ipg++){
      schnaps_real wmac[lsd->lb_model->nb_macro];
      schnaps_real wmic[lsd->lb_model->q];
      for(int inode=0;inode<lsd->lb_model->q;inode++){
      int imic=fmic->varindex(fmic->deg, fmic->raf, fmic->model.m, ipg, inode);
      wmic[inode]=fmic->wn[imic];
      } //inode
      lsd->lb_model->f_to_macro(wmic,wmac);
      for(int iv=0;iv<lsd->lb_model->nb_macro;iv++){
          int imac=f->varindex(f->deg, f->raf, f->model.m, ipg, iv);
          f->wn[imac]=wmac[iv];
      };
    }//ipg
  } //ie
}
/***********************************************************************************/
// wrappers for compatibility with RK schemes
void LBM_pre_dtfields_wrapper( void * simu){
  LatticeBoltzmannSimData *lsd=&schnaps_lbm_simdata;
  LBMSimulation *lbsimu=lsd->current_lb_sim;
  printf("pre dt wrapper");
  // update time data 
  lbsimu->tnow=lbsimu->micro_simu.tnow;
  lbsimu->iter_time=lbsimu->micro_simu.iter_time_rk;
  lbsimu->itermax=lbsimu->micro_simu.itermax_rk;
  // redirect call to global function
  lbsimu->pre_advec(lbsimu);
}
void LBM_post_dtfields_wrapper( void * simu){
  LatticeBoltzmannSimData *lsd=&schnaps_lbm_simdata;
  LBMSimulation *lbsimu=lsd->current_lb_sim;
  printf("post dt wrapper");
  // update time data 
  lbsimu->tnow=lbsimu->micro_simu.tnow;
  lbsimu->iter_time=lbsimu->micro_simu.iter_time_rk;
  lbsimu->itermax=lbsimu->micro_simu.itermax_rk;
  // redirect call to global function
  lbsimu->post_advec(lbsimu);
}
void LBM_update_after_rk_wrapper( void *simu,schnaps_real *w){
  LatticeBoltzmannSimData *lsd=&schnaps_lbm_simdata;
  LBMSimulation *lbsimu=lsd->current_lb_sim;
  printf("update_after_rk dt wrapper");
  // update time data
  lbsimu->tnow=lbsimu->micro_simu.tnow;
  lbsimu->iter_time=lbsimu->micro_simu.iter_time_rk;
  lbsimu->itermax=lbsimu->micro_simu.itermax_rk;
  // redirect call to global function
  lbsimu->post_tstep(lbsimu,w);
}
//******************************************************************************//
//******************************************************************************//
//******************************************************************************//
// flux functions for LB models
//******************************************************************************//
//******************************************************************************//
void LBM_OneLatticeNumFlux(schnaps_real *wL, schnaps_real *wR, schnaps_real *vnorm, schnaps_real *flux){
  LatticeBoltzmannSimData *lsd=&schnaps_lbm_simdata;
  for(int i = 0;i < lsd->lb_model->q;i++){
    schnaps_real vn=0;
    for(int dim = 0;dim < lsd->lb_model->d; dim++){
      vn += lsd->lb_model->vi[i][dim]*vnorm[dim];
    }
    schnaps_real vnp = vn>0 ? vn : 0;
    schnaps_real vnm = vn-vnp;
    flux[i] = vnp * wL[i] + vnm * wR[i];
  }
}
void LBM_OneNodeNumFlux(schnaps_real *wL, schnaps_real *wR, schnaps_real *vnorm, schnaps_real *flux){
  LatticeBoltzmannSimData *lsd=&schnaps_lbm_simdata;
  int inode=lsd->current_node_index;
  schnaps_real vn=0;
  for(int dim = 0;dim < lsd->lb_model->d; dim++){
    vn += lsd->lb_model->vi[inode][dim]*vnorm[dim];
  }
  schnaps_real vnp = vn>0 ? vn : 0;
  schnaps_real vnm = vn-vnp;
  flux[0] = vnp * wL[0] + vnm * wR[0];
}
/***********************************************************************************/
// Time schemes
/*************************************************************************************/
void LBMThetaTimeScheme(LBMSimulation *lbsimu,schnaps_real theta, schnaps_real tmax, schnaps_real dt){
  LatticeBoltzmannSimData *lsd=&schnaps_lbm_simdata;
  int nb_nodes=lsd->lb_model->q; // velocity nodes on the lattice(s)
  int ig_glob=0,ig_node=0;
  field *f_glob,*f_node;
  Simulation *micsimu=&(lbsimu->micro_simu);
  Simulation *macsimu=&(lbsimu->macro_simu);
  // 
  int nraf[3] = {micsimu->fd[0].raf[0],
		 micsimu->fd[0].raf[1],
		 micsimu->fd[0].raf[2]};
  int deg[3] = {micsimu->fd[0].deg[0],
		micsimu->fd[0].deg[1],
		micsimu->fd[0].deg[2]};  
  int nb_ipg = NPG(deg, nraf);
  //
  Simulation *simu_advec;
  simu_advec=malloc(sizeof(Simulation));
  EmptySimulation(simu_advec);
  InitSimulation(simu_advec, &(micsimu->macromesh), deg, nraf, &(lbsimu->model_advec));
  //
  simu_advec->vmax=lbsimu->vmax;
  simu_advec->cfl=0.0;
  simu_advec->nb_diags=0;
  simu_advec->pre_dtfields=NULL;
  simu_advec->post_dtfields=NULL;
  simu_advec->update_after_rk=NULL;
  //
  LinearSolver solver_imp[nb_nodes]; 
  LinearSolver solver_exp[nb_nodes];
  schnaps_real *res_advec=  calloc(simu_advec->wsize,sizeof(schnaps_real*)); 
  //
  lbsimu->dt=dt;
  int itermax=(int) (tmax/dt)+1;
  printf("Called with tmax=%f dt=%f Nb iterations:%i \n",tmax,dt,itermax);
  lbsimu->itermax=itermax;
  // important sync micro simu params for diag compatibility with RK which uses only the micro simu);
  lbsimu->micro_simu.itermax_rk=itermax; 
  lbsimu->micro_simu.dt=dt;
  int freq = (1 >= lbsimu->itermax / 10)? 1 : lbsimu->itermax / 10;
  lbsimu->tnow=0.0;
  simu_advec->dt=dt;
  simu_advec->itermax_rk=itermax;
  simu_advec->tnow=0.0;
  for(int ie=0; ie < micsimu->macromesh.nbelems; ++ie){
    printf("ie:%i",ie);
    micsimu->fd[ie].tnow=lbsimu->tnow;
    macsimu->fd[ie].tnow=lbsimu->tnow;
    simu_advec->fd[ie].tnow=simu_advec->tnow;
  }
  // Diagnostics (this should be elsewhere, some timetraces  module ?
  int mac_size_diags = macsimu->nb_diags * itermax;
  int mic_size_diags = micsimu->nb_diags * itermax;
  if(macsimu->nb_diags != 0) {
     macsimu->Diagnostics = malloc(mac_size_diags * sizeof(schnaps_real));
  };
  if(micsimu->nb_diags != 0) {
     micsimu->Diagnostics = malloc(mic_size_diags * sizeof(schnaps_real));
  };
  //
  time_t t_start,t_end; // time measurements for op factorization
  t_start=time(NULL);
  //  Solvers Init/Assembly
  printf("Sparse Linear Solvers init");
  for (int isim=0;isim < nb_nodes;isim++){
    lsd->current_node_index=isim;
    //
    InitImplicitLinearSolver(simu_advec, &solver_imp[isim]);
    InitImplicitLinearSolver(simu_advec, &solver_exp[isim]);
    //
  }
  // End Operators init
  //
  // Time loop start
  for (int iter=0;iter < itermax;iter++){
    //
    if (iter%freq == 0){
      printf(" iter %i/%i t=%f\n",iter,itermax,lbsimu->tnow);
    }
    //
    lbsimu->iter_time=iter;
    macsimu->iter_time_rk=iter;
    micsimu->iter_time_rk=iter;
    //
    if (iter==0){
      t_start=time(NULL);
    }
    //
    if (lbsimu->pre_advec !=NULL){
      lbsimu->pre_advec(lbsimu);
    };
    // now loop on velocity nodes
    for(int isim=0;isim < nb_nodes;isim++){
      lsd->current_node_index=isim;
      // dispatch main w to per node w's
      for(int ie = 0; ie < simu_advec->macromesh.nbelems; ie++){
        f_glob=micsimu->fd+ie;
        for(int ipg=0; ipg < nb_ipg;ipg++){
          f_node= simu_advec->fd+ie;
          ig_glob=f_glob->varindex(f_glob->deg,f_glob->raf,f_glob->model.m,ipg,isim);
          ig_node=f_node->varindex(f_node->deg,f_node->raf,1,ipg,0);
          f_node->wn[ig_node]= f_glob->wn[ig_glob];
        //
        }; // ipg end loop glops
      }; // ie end loop macroelements
      // end of data dispatch
      // Solvers assembly and factorization if necessary
      solver_imp[isim].rhs_is_assembly=false;
      solver_exp[isim].rhs_is_assembly=false;
      if (iter==0){
          solver_imp[isim].mat_is_assembly=false;
          solver_exp[isim].mat_is_assembly=false;
      }
      else
      {
          solver_imp[isim].mat_is_assembly=true;
          solver_exp[isim].mat_is_assembly=true;
      };
      //
      simu_advec->tnow=lbsimu->tnow;
      for(int ie=0; ie < simu_advec->macromesh.nbelems; ++ie){
        simu_advec->fd[ie].tnow=simu_advec->tnow;
      }
      AssemblyImplicitLinearSolver(simu_advec, &solver_exp[isim],-(1.0-theta),simu_advec->dt);
      simu_advec->tnow=lbsimu->tnow+lbsimu->dt;
      for(int ie=0; ie < simu_advec->macromesh.nbelems; ++ie){
        simu_advec->fd[ie].tnow=simu_advec->tnow;
      }
      AssemblyImplicitLinearSolver(simu_advec, &solver_imp[isim],theta,simu_advec->dt);
      // compute residual
      MatVect(&solver_exp[isim], simu_advec->w, res_advec);
      //
      for(int i=0;i<solver_imp[isim].neq;i++){
        solver_imp[isim].rhs[i]=-solver_exp[isim].rhs[i]+solver_imp[isim].rhs[i]+res_advec[i];
      }
      //
      solver_imp[isim].solver_type=LU;
      Advanced_SolveLinearSolver(&solver_imp[isim],simu_advec);
      //
      for(int i=0;i<solver_imp[isim].neq;i++){
        simu_advec->w[i]=solver_imp[isim].sol[i];
      }
      // collect nodes ws to main w
      for(int ie = 0; ie < micsimu->macromesh.nbelems; ie++){
        f_glob=micsimu->fd+ie;
        for(int ipg=0; ipg < nb_ipg;ipg++){
          f_node= simu_advec->fd+ie;
          ig_glob=f_glob->varindex(f_glob->deg,f_glob->raf,f_glob->model.m,ipg,isim);
          ig_node=f_node->varindex(f_node->deg,f_node->raf,1,ipg,0);
          f_glob->wn[ig_glob]= f_node->wn[ig_node]; 
        }; // ipg end loop glops
      }; // ie end loop macroelements
    }; // isim end loop on velocity node 
    // post advec ops
    if (lbsimu->post_advec !=NULL){
      lbsimu->post_advec(lbsimu);
    }
    if (lbsimu->post_tstep !=NULL){
      lbsimu->post_tstep(lbsimu,macsimu->w);
    }
    if (iter==0){
      t_end= time(NULL);
      printf("First step duration %ld\n", t_end-t_start);
    }
    lbsimu->tnow += lbsimu->dt;
    simu_advec->tnow=lbsimu->tnow;
    micsimu->tnow=lbsimu->tnow;
    macsimu->tnow=lbsimu->tnow;
    for(int ie=0; ie < simu_advec->macromesh.nbelems; ++ie){
      micsimu->fd[ie].tnow=lbsimu->tnow;
      macsimu->fd[ie].tnow=lbsimu->tnow;
      simu_advec->fd[ie].tnow=lbsimu->tnow;
    }
    //
  }; // iter end of time loop 
  //
  if (res_advec !=NULL){
    free(res_advec);
  }
  if (simu_advec!= NULL){
    free(simu_advec);
  }
}
/***************************************************************************************************************/

/***********************************************************************************/
//********************************************************************************//
//********************************************************************************//
//********************************************************************************//
// per model specific routines and functions 
//******************************************************************************//
//********************************************************************************//
//********************************************************************************//
schnaps_real LBM_dummy_zeros_feq(int inode,int nb_macro,schnaps_real *w){
  return 0.0;
}
// D2Q9 isothermal
//********************************************************************************//
void LBM_Set_D2Q9_ISOTH_model(LBModelDescriptor *lb){
  assert(lb->d==2);
  assert(lb->q==9);
  assert(lb->nb_macro==3);
  NewMPolyDescriptorFromModel(&(lb->Moments[0]),&model_rho2D,NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[1]),&model_jx2D,NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[2]),&model_jy2D,NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[3]),&model_trMab2D,NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[4]),&model_DifMab2D,NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[5]),&model_Mxy2D,NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[6]),&model_Mabcx2D,NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[7]),&model_Mabcy2D,NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[8]),&model_V42D,NULL);
  //lb->cref= sqrt(1.0/3.0);
  sprintf(lb->macro_names[0],"rho");
  sprintf(lb->macro_names[1],"ux");
  sprintf(lb->macro_names[2],"uy");
  lb->cref= 1.0;
  for (int i=0;i< lb->q;i++){
    for (int j=0;j<lb->d;j++){
      lb->vi[i][j]=lb->cref * LBM_D2Q9_nodes[i][j];
    }
  }
  for (int i=0; i< lb->q;i++){
    lb->inode_min[i]=0;
    lb->inode_max[i]=lb->q-1;
  }
  ComputeLBModelDescriptorMomentMatrix(lb);
  //
  lb->feq=&LBM_feq_D2Q9_ISOTH;
  lb->f_to_macro=&LBM_f_to_macro_D2Q9_ISOTH;
  //
}
//
void LBM_f_to_macro_D2Q9_ISOTH(schnaps_real *f,schnaps_real *w){
    LatticeBoltzmannSimData *lsd=&schnaps_lbm_simdata;
    for (int i=0;i< lsd->lb_model->nb_macro;i++){
      w[i]=0.0;
      for (int inode=0;inode< lsd->lb_model->q;inode++){
        w[i]+=f[inode] * lsd->lb_model->M[i][inode];
      }
    }
    // normalize
    if (w[0] !=0.0){
      w[1]=w[1]/w[0];
      w[2]=w[2]/w[0];
    }
}
//
schnaps_real LBM_feq_D2Q9_ISOTH(int i_node,int nb_macro,schnaps_real *w){
    //
    LatticeBoltzmannSimData *lsd=&schnaps_lbm_simdata;
    schnaps_real rho=w[0];
    schnaps_real ux=w[1];
    schnaps_real uy=w[2];
    //schnaps_real temp=w[4];
    schnaps_real temp=1.0/3.0;
    //
    schnaps_real u2= (ux * ux + uy * uy)/temp;
    schnaps_real uv= (ux * lsd->lb_model->vi[i_node][0] + uy * lsd->lb_model->vi[i_node][1])/temp;
    schnaps_real feq= LBM_WEIGHTS_D2Q9_ISOTH[i_node] * rho * (1.0+uv+0.5 * (uv * uv- u2));
    return feq;
};
// D2Q9 isothermal linearized (2D wave equation)
void LBM_Set_D2Q9_ISOTH_LINEARIZED_model(LBModelDescriptor *lb){
  assert(lb->d==2);
  assert(lb->q==9);
  assert(lb->nb_macro==3);
  NewMPolyDescriptorFromModel(&(lb->Moments[0]),&model_rho2D,NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[1]),&model_jx2D,NULL);
  NewMPolyDescriptorFromModel(&(lb->Moments[2]),&model_jy2D,NULL);
  //
  sprintf(lb->macro_names[0],"rho");
  sprintf(lb->macro_names[1],"jx");
  sprintf(lb->macro_names[2],"jy");
  //
  lb->cref= 1.0;
  for (int i=0;i< lb->q;i++){
    for (int j=0;j<lb->d;j++){
      lb->vi[i][j]=lb->cref * LBM_D2Q9_nodes[i][j];
    }
  }
  for (int i=0; i< lb->q;i++){
    lb->inode_min[i]=0;
    lb->inode_max[i]=lb->q-1;
  }
  ComputeLBModelDescriptorMomentMatrix(lb);
  //
  lb->feq=&LBM_feq_D2Q9_ISOTH_LINEARIZED;
  lb->f_to_macro=&LBM_f_to_macro_D2Q9_ISOTH_LINEARIZED;
  //
}
schnaps_real LBM_feq_D2Q9_ISOTH_LINEARIZED(int i_node,int nb_macro,schnaps_real *w){
    //
    LatticeBoltzmannSimData *lsd=&schnaps_lbm_simdata;
    schnaps_real rho=w[0];
    schnaps_real jx=w[1];
    schnaps_real jy=w[2];
    //schnaps_real temp=w[4];
    schnaps_real temp=1.0/3.0;
    //
    schnaps_real jv= (jx * lsd->lb_model->vi[i_node][0] + jy * lsd->lb_model->vi[i_node][1])/temp;
    schnaps_real feq= LBM_WEIGHTS_D2Q9_ISOTH[i_node] * (rho+jv);
    return feq;
};
void LBM_f_to_macro_D2Q9_ISOTH_LINEARIZED(schnaps_real *f,schnaps_real *w){
    LatticeBoltzmannSimData *lsd=&schnaps_lbm_simdata;
    for (int i=0;i< lsd->lb_model->nb_macro;i++){
      w[i]=0.0;
      for (int inode=0;inode< lsd->lb_model->q;inode++){
        w[i]+=f[inode] * lsd->lb_model->M[i][inode];
      }
    }
}
//
