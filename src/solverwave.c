#include "solvercontinuous.h"
#include "geometry.h"
#include "quantities_vp.h"
#include "linear_solver.h"
#include "macromesh.h"
#include "solverwave.h"
#include "waterwave2d.h"

#include <math.h>

void Assembly_SW(Simulation* simu, void* solver, real* w, real theta, real dt){

  ContinuousSolver * ps=solver;
  real g=_GRAVITY;

  field* f0 = &ps->simu->fd[0];

  if (!ps->lsol.mat_is_assembly){
  for(int ie = 0; ie < ps->nbel; ie++){  


    int iemacro = ie / (f0->raf[0] * f0->raf[1] * f0->raf[2]);
    int isubcell = ie % (f0->raf[0] * f0->raf[1] * f0->raf[2]);

    // Inside the element ie

    //////////////////////////////////////////////
    ////// Init : Local vector of variables //////
    //////////////////////////////////////////////

    // value = 0, dx = 1, dy = 2, dz = 3
    real var[ps->nb_phy_vars][4][ps->nnodes]; 
    for(int iv = 0; iv<ps->nb_phy_vars; iv++){
      for(int i_der = 0; i_der<4; i_der++){
        for(int iloc = 0; iloc < ps->nnodes; iloc++){
          var[iv][i_der][iloc] = 0.0;
        }
      }
    }
    
    // Computation of the variable and its derivatives for each Gauss point
    for(int iloc = 0; iloc < ps->nnodes; iloc++){
      real wpg;
      real xref[3];
      real dtau[3][3],codtau[3][3];
      real dphiref_j[3];
      real dphi_j[3];
      real basisPhi_j[4];

      int ipgmacro = iloc + isubcell * ps->nnodes;
      ref_pg_vol(f0->deg,f0->raf,ipgmacro,xref,&wpg,NULL);
      
      // Loop on variables
      for(int iv=0; iv<ps->nb_phy_vars; iv++){
        int ino_dg = iloc + ie * ps->nnodes;
        int ino_fe = ps->dg_to_fe_index[ino_dg];
        var[iv][0][iloc] = w[iv+ino_fe*ps->nb_phy_vars];

        for(int jGauss = 0; jGauss < ps->nnodes; jGauss++){
          int jGaussmacro = jGauss + isubcell * ps->nnodes;
          int jGauss_dg = jGauss + ie * ps->nnodes;
          int jGauss_fe = ps->dg_to_fe_index[jGauss_dg];
          grad_psi_pg(f0->deg,f0->raf,jGaussmacro,ipgmacro,dphiref_j);
          Ref2Phy(ps->simu->fd[iemacro].physnode,
                  xref,dphiref_j,0,NULL,
                  dtau,codtau,dphi_j,NULL);
          real det = dot_product(dtau[0], codtau[0]);
          basisPhi_j[1]=dphi_j[0]/det;
          basisPhi_j[2]=dphi_j[1]/det;
          basisPhi_j[3]=dphi_j[2]/det;
          for (int i=1;i<4;i++){
            var[iv][i][iloc] += basisPhi_j[i]*w[iv+jGauss_fe*ps->nb_phy_vars];
          } // end for i (derivative)
        } // end for jGauss
      } // end for iv
    } // end for iloc


    // Init: local matrix and rhs
    real aloc[ps->nnodes*ps->nb_phy_vars][ps->nnodes*ps->nb_phy_vars];
    real rhsloc[ps->nnodes*ps->nb_phy_vars];
    for(int iloc = 0; iloc < ps->nnodes*ps->nb_phy_vars; iloc++){
      rhsloc[iloc] = 0.0;
      for(int jloc = 0; jloc < ps->nnodes*ps->nb_phy_vars; jloc++){
        aloc[iloc][jloc] = 0.0;
      }
    }

    /////////////////////////////////////////////////
    /////////// Now compute matrix and RHS //////////
    /////////////////////////////////////////////////
    for(int ipg = 0;ipg < ps->nnodes; ipg++){
      real wpg;
      real xref[3];
      int ipgmacro= ipg + isubcell * ps->nnodes;
      ref_pg_vol(f0->deg,f0->raf,ipgmacro,xref,&wpg,NULL);

      for(int iloc = 0; iloc < ps->nnodes; iloc++){


        real dtau[3][3],codtau[3][3];
        real dphiref_i[3],dphiref_j[3];
        real dphi_i[3],dphi_j[3];
        real basisPhi_i[4], basisPhi_j[4];
        int ilocmacro = iloc + isubcell * ps->nnodes;
        grad_psi_pg(f0->deg,f0->raf,ilocmacro,ipgmacro,dphiref_i);
        Ref2Phy(ps->simu->fd[iemacro].physnode,
                xref,dphiref_i,0,NULL,
                dtau,codtau,dphi_i,NULL);
      
        real det = dot_product(dtau[0], codtau[0]);
        if (ilocmacro==ipgmacro){
          basisPhi_i[0]=1;
        }
        else
        {
          basisPhi_i[0]=0;
        }
        basisPhi_i[1]=dphi_i[0]/det;
        basisPhi_i[2]=dphi_i[1]/det;
        basisPhi_i[3]=dphi_i[2]/det;

        real h  = var[0][0][iloc];
        real hx = var[0][1][iloc];
        real hy = var[0][2][iloc];
        real hz = var[0][3][iloc];
        real u  = var[1][0][iloc];
        real ux = var[1][1][iloc];
        real uy = var[1][2][iloc];
        real uz = var[1][3][iloc];
        real v  = var[2][0][iloc];
        real vx = var[2][1][iloc];
        real vy = var[2][2][iloc];
        real vz = var[2][3][iloc];

        // Building RHS
        rhsloc[iloc*ps->nb_phy_vars  ] += - basisPhi_i[0] * dt * (h*(ux+vy)+u*hx+v*hy+0.0*hz) * wpg * det;
        rhsloc[iloc*ps->nb_phy_vars+1] += - basisPhi_i[0] * dt * h * ( u*ux + v*uy + g*hx ) * wpg * det;
        rhsloc[iloc*ps->nb_phy_vars+2] += - basisPhi_i[0] * dt * h * ( u*vx + v*vy + g*hy ) * wpg * det;

        //printf("iloc=%d, h=%.12e\n",iloc,h);
        // Building Jacobian of the problem
        real Jac[9][4][4] = {{{1+dt*theta*(ux+uy),dt*theta*u,dt*theta*v,0},
                             {0,0,0,0},
                             {0,0,0,0},
                             {0,0,0,0}},
                            {{dt*theta*hx,dt*theta*h,0,0},
                             {0,0,0,0},
                             {0,0,0,0},
                             {0,0,0,0}},
                            {{dt*theta*hy,0,dt*theta*h,0},
                             {0,0,0,0},
                             {0,0,0,0},
                             {0,0,0,0}},
                            
                            {{dt*theta*(u*ux+v*uy+g*hx),dt*theta*g*h,0,0},
                             {0,0,0,0},
                             {0,0,0,0},
                             {0,0,0,0}},
                            {{h+dt*theta*h*ux,dt*theta*h*u,dt*theta*h*v,0},
                             {0,0,0,0},
                             {0,0,0,0},
                             {0,0,0,0}},
                            {{dt*theta*h*uy,0,0,0},
                             {0,0,0,0},
                             {0,0,0,0},
                             {0,0,0,0}},
                            
                            {{dt*theta*(u*vx+v*vy+g*hy),0,dt*theta*g*h,0},
                             {0,0,0,0},
                             {0,0,0,0},
                             {0,0,0,0}},
                            {{dt*theta*h*vx,0,0,0},
                             {0,0,0,0},
                             {0,0,0,0},
                             {0,0,0,0}},
                            {{h+dt*theta*h*vy,dt*theta*h*u,dt*theta*h*v,0},
                            {0,0,0,0},
                            {0,0,0,0},
                            {0,0,0,0}}};

        for(int jloc = 0; jloc < ps->nnodes; jloc++){
          int jlocmacro = jloc + isubcell * ps->nnodes;
          grad_psi_pg(f0->deg,f0->raf,jlocmacro,ipgmacro,dphiref_j);
          Ref2Phy(ps->simu->fd[iemacro].physnode,
                  xref,dphiref_j,0,NULL,
                  dtau,codtau,dphi_j,NULL);
          if (jlocmacro==ipgmacro){
            basisPhi_j[0]=1;
          }
          else
          {
            basisPhi_j[0]=0;
          }
          basisPhi_j[1]=dphi_j[0]/det;
          basisPhi_j[2]=dphi_j[1]/det;
          basisPhi_j[3]=dphi_j[2]/det;

          for (int iv1=0; iv1<ps->nb_phy_vars; iv1++){
            for (int iv2=0; iv2<ps->nb_phy_vars; iv2++){
              real res[4] = {0, 0, 0, 0};
              for (int i=0; i<4; i++){
                for (int j=0; j<4; j++){
                  res[i]+=basisPhi_j[j]*Jac[ps->nb_phy_vars*iv1+iv2][i][j];
                }
              }
              aloc[iv1+iloc*ps->nb_phy_vars][iv2+jloc*ps->nb_phy_vars] += dot_product(basisPhi_i, res) * wpg * det  ;
            }
          }
        } // end for jloc
      } // end for iloc
    } // end for ipg    

    // Building matrix
    for(int iloc = 0; iloc < ps->nnodes; iloc++){
      int ino_dg = iloc + ie * ps->nnodes;
      int ino_fe = ps->dg_to_fe_index[ino_dg];
      for (int iv=0; iv<ps->nb_phy_vars; iv++){
        int ipot_fe = ino_fe*ps->nb_phy_vars+iv;
        ps->lsol.rhs[ino_fe*ps->nb_phy_vars+iv] += rhsloc[iloc*ps->nb_phy_vars+iv];
      }
      for(int jloc = 0; jloc < ps->nnodes; jloc++){
        int jno_dg = jloc + ie * ps->nnodes;
        int jno_fe = ps->dg_to_fe_index[jno_dg];
        for (int iv1=0;iv1<ps->nb_phy_vars;iv1++){
          for (int iv2=0;iv2<ps->nb_phy_vars;iv2++){
            real val = aloc[iloc*ps->nb_phy_vars+iv1][jloc*ps->nb_phy_vars+iv2];
            AddLinearSolver(&ps->lsol,ino_fe*ps->nb_phy_vars+iv1,jno_fe*ps->nb_phy_vars+iv2,val);
          }
        }
      }
    }
  }
  }
}



void SW_test(ContinuousSolver* cs ,real* w, real theta, real dt){
  Assembly_SW(cs->simu, cs, w, theta, dt);

  real h=cs->simu->vmax*dt*theta;
  cs->FluxMatrix = calloc(cs->nb_phy_vars,sizeof(real));
  for (int i=0; i<cs->nb_phy_vars; i++){
    cs->FluxMatrix[i] = calloc(cs->nb_phy_vars,sizeof(real));
  }
  for (int i=0; i<cs->nb_phy_vars; i++){
    for (int j=0; j<cs->nb_phy_vars; j++){
      cs->FluxMatrix[i][j] = 0.0;
    }
  }

  cs->FluxMatrix[0][1]=h;
  cs->FluxMatrix[0][2]=h;
  cs->FluxMatrix[1][0]=h;
  cs->FluxMatrix[2][0]=h;
}



void Wave_test(ContinuousSolver* cs,real theta, real dt){
  real h=cs->simu->vmax*dt*theta;
  real tab[4][4];
  
  real waveMat[9][4][4] ={{{1.0,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0}},
                          {{0.0,0,0,0},
                           {-h,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0}},
                          {{0,0,0,0},
                           {0,0,0,0},
                           {-h,0,0,0},
                           {0,0,0,0}},
                          {{0,0,0,0},
                           {-h,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0}},
                          {{1.0,0,0,0}, 
                           {0,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0}},
                          {{0.0,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0}},
                          {{0,0,0,0},
                           {0,0,0,0},
                           {-h,0,0,0},
                           {0,0,0,0}},
                          {{0.0,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0}},
                          {{1.0,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0},
                           {0,0,0,0}}};
  
  cs->diff_op = calloc(cs->nb_phy_vars*cs->nb_phy_vars,sizeof(SDO));
  
  for (int i=0; i<cs->nb_phy_vars*cs->nb_phy_vars; i++){
    for (int j=0; j<4; j++){
      for (int k=0; k<4; k++){
        cs->diff_op[i].DO[j][k]=waveMat[i][j][k];
      }
    }
  }

  cs->FluxMatrix = calloc(cs->nb_phy_vars,sizeof(real));
  for (int i=0; i<cs->nb_phy_vars; i++){
    cs->FluxMatrix[i] = calloc(cs->nb_phy_vars,sizeof(real));
  }
  cs->FluxMatrix[0][1]=h;
  cs->FluxMatrix[0][2]=h;
  cs->FluxMatrix[1][0]=h;
  cs->FluxMatrix[2][0]=h;
  cs->FluxMatrix[0][0]=0.0;
  cs->FluxMatrix[1][1]=0.0;
  cs->FluxMatrix[1][2]=0.0;
  cs->FluxMatrix[2][1]=0.0;
  cs->FluxMatrix[2][2]=0.0; 
}


void BoundaryConditionFriedrichsAssembly(void * cs){
  ContinuousSolver * ps=cs;
  m=ps->nb_phy_vars;
  
  if (!ps->lsol.mat_is_assembly){

    for (int i=0; i<ps->simu->macromesh.nboundaryfaces;i++){
      int ifa = ps->simu->macromesh.boundaryface[i];
      int locfaL = ps->simu->macromesh.face2elem[4 * ifa + 1];
      int ie = ps->simu->macromesh.face2elem[4 * ifa ];
      int ieR = ps->simu->macromesh.face2elem[4 * ifa + 2];
      assert(ieR < 0);
      if (ieR<0){
	field *f = &ps->simu->fd[ie];
	
	for(int ipglf = 0;ipglf < NPGF(f->deg,f->raf,locfaL); ipglf++){
	  real xpgref[3], xpgref_in[3], wpg;


	  // Get the coordinates of the Gauss point and coordinates of a
	  // point slightly inside the opposite element in xref_in
	  int ipg = ref_pg_face(f->deg, f->raf, locfaL, ipglf, xpgref, &wpg, xpgref_in);
	  int ino_dg = ipg + ie * ps->npgmacrocell;
	  int ino_fe = ps->dg_to_fe_index[ino_dg];
	  // Normal vector at gauss point ipgL
	  real vnds[3], xpg[3];
	  {
	    real dtau[3][3], codtau[3][3];
	    Ref2Phy(f->physnode,
		    xpgref,
		    NULL, locfaL, // dpsiref, ifa
		    xpg, dtau,
		    codtau, NULL, vnds); // codtau, dpsi, vnds
	  }
	  // the boundary flux is an affine function
	  real flux0[m], wL[m];
	  for(int iv = 0; iv < m; iv++) {
	    wL[iv] = 0;
	  }
	  
	  ps->bc_flux(ps,xpg,wL,vnds,flux0);


	  for(int iv1 = 0; iv1 < m; iv1++) {
	    int ipot_fe1 = ino_fe*ps->nb_phy_vars + iv1;
	    for(int iv = 0; iv < m; iv++) {
	      wL[iv] = (iv == iv1);
	    }
	    
	    real flux[m];
	    ps->bc_flux(ps,xpg,wL,vnds,flux);
	    
	    for(int iv2 = 0; iv2 < m; iv2++) {
	      // The basis functions is also the gauss point index
	      int ipot_fe2 = ino_fe*ps->nb_phy_vars + iv2;
	      real val =  (flux[iv2]-flux0[iv2]) * wpg;
	      AddLinearSolver(&ps->lsol, ipot_fe2, ipot_fe1, val);
	    }
	  }	    
	}
      }
    }
 
    ps->lsol.mat_is_assembly=true;
  }
  
  for (int i=0; i<ps->simu->macromesh.nboundaryfaces;i++){
    int ifa = ps->simu->macromesh.boundaryface[i];
    int locfaL = ps->simu->macromesh.face2elem[4 * ifa + 1];
    int ie = ps->simu->macromesh.face2elem[4 * ifa ];
    int ieR = ps->simu->macromesh.face2elem[4 * ifa + 2];
    assert(ieR < 0);
    if (ieR<0){
      field *f = &ps->simu->fd[ie];

      for(int ipglf = 0;ipglf < NPGF(f->deg,f->raf,locfaL); ipglf++){
	real xpgref[3], xpgref_in[3], wpg;
          
	// Get the coordinates of the Gauss point and coordinates of a
	// point slightly inside the opposite element in xref_in
	int ipg = ref_pg_face(f->deg, f->raf, locfaL, ipglf, xpgref, &wpg, xpgref_in);
	int ino_dg = ipg + ie * ps->npgmacrocell;
	int ino_fe = ps->dg_to_fe_index[ino_dg];
	// Normal vector at gauss point ipgL
	real vnds[3], xpg[3];
	{
	  real dtau[3][3], codtau[3][3];
	  Ref2Phy(f->physnode,
		  xpgref,
		  NULL, locfaL, // dpsiref, ifa
		  xpg, dtau,
		  codtau, NULL, vnds); // codtau, dpsi, vnds
	}
          
	// the boundary flux is an affine function
	real w0[m],flux0[m];
	
	for(int ivv=0; ivv < ps->nb_phy_vars; ivv++){
	  w0[ivv]=0;
	}
	
	ps->bc_flux(ps,xpg,w0,vnds,flux0);
	
	for(int var=0; var < ps->nb_phy_vars; var++){
	  int ipot_fe = ino_fe*ps->nb_phy_vars + var;
	  //printf("clll %d %d  %.8e %.8e \n",ipot_fe,var,flux0[var],w0[var]);
	  real val = flux0[var] * wpg;
	  ps->lsol.rhs[ipot_fe] -= val;
	}
      }	    
    }
  }
}

void Wave_BC_normalvelocity_null(void * cs, real * xpg, real * w, real *vnorm, real * flux){

  ContinuousSolver * ps=cs;
  real M[3][3];
  real BC[3][3];
  real lambda=0;
  real mu=-1.e21;
  real p0=0,u0_1=0,u0_2=0;

  real h=ps->FluxMatrix[0][1]/ps->simu->vmax;
  
  M[0][0]=0;//lambda;
  M[0][1]=0;//-mu*vnorm[0];
  M[0][2]=0;//-mu*vnorm[1];
  M[1][0]=0;
  M[1][1]=mu*vnorm[0]*vnorm[0];
  M[1][2]=mu*vnorm[0]*vnorm[1];
  M[2][0]=0;
  M[2][1]=mu*vnorm[1]*vnorm[0];
  M[2][2]=mu*vnorm[1]*vnorm[1]; 


  BC[0][0]=h*M[0][0];
  BC[0][1]=(ps->FluxMatrix[0][1]*vnorm[0]+h*M[0][1]);
  BC[0][2]=(ps->FluxMatrix[0][2]*vnorm[1]+h*M[0][2]);
  BC[1][0]=(ps->FluxMatrix[1][0]*vnorm[0]+h*M[1][0]);
  BC[1][1]=h*M[1][1];
  BC[1][2]=h*M[1][2];
  BC[2][0]=(ps->FluxMatrix[2][0]*vnorm[1]+h*M[2][0]);
  BC[2][1]=h*M[2][1];
  BC[2][2]=h*M[2][2];

  flux[0]=(BC[0][0]*w[0]+BC[0][1]*w[1]+BC[0][2]*w[2])-h*(M[0][0]*p0+M[0][1]*u0_1+M[0][2]*u0_2);
  flux[1]=(BC[1][0]*w[0]+BC[1][1]*w[1]+BC[1][2]*w[2])-h*(M[1][0]*p0+M[1][1]*u0_1+M[1][2]*u0_2);
  flux[2]=(BC[2][0]*w[0]+BC[2][1]*w[1]+BC[2][2]*w[2])-h*(M[2][0]*p0+M[2][1]*u0_1+M[2][2]*u0_2);

  int angle=0;

  if(xpg[0]== 0.0 && xpg[1]== 0.0) { angle=1; }
  if(xpg[0]== 0.0 && xpg[1]== 1.0) { angle=1; }
  if(xpg[0]== 1.0 && xpg[1]== 0.0) { angle=1; }
  if(xpg[0]== 1.0 && xpg[1]== 1.0) { angle=1; }

  if(angle ==0){
    flux[0]= 0.0;
    flux[1]= mu * (w[1]*vnorm[0]*vnorm[0] + w[2]*vnorm[0]*vnorm[1])
      - mu * (u0_1*vnorm[0]*vnorm[0] + u0_2*vnorm[0]*vnorm[1]);
    flux[2]=  mu * (w[1]*vnorm[1]*vnorm[0] + w[2]*vnorm[1]*vnorm[1])
      -  mu * (u0_1*vnorm[1]*vnorm[0] + u0_2*vnorm[1]*vnorm[1]);
   }
  else {
    flux[0]= 0.0;
    flux[1]= mu * w[1] - mu * u0_1;
    flux[2]= mu * w[2] - mu * u0_2;
  }
}


void Wave_BC_pressure_null(void * cs, real * xpg, real * w, real *vnorm, real * flux){

  ContinuousSolver * ps=cs;
  real M[3][3];
  real BC[3][3];
  real lambda=+1e14;
  real mu=0;
  real p0=0,u0_1=0,u0_2=0;
  real x,y;
  x=xpg[0];
  y=xpg[1];

  real Win[ps->simu->fd[0].model.m];
  
  ps->simu->fd[0].model.ImposedData(xpg,ps->simu->dt,Win);

  real h=ps->FluxMatrix[0][1]/ps->simu->vmax;


  M[0][0]=lambda;
  M[0][1]=-mu*vnorm[0];
  M[0][2]=-mu*vnorm[1];
  M[1][0]=0;
  M[1][1]=0;
  M[1][2]=0;
  M[2][0]=0;
  M[2][1]=0;
  M[2][2]=0;

  flux[0]=M[0][0]*w[0]+M[0][1]*w[1]+M[0][2]*w[2] - M[0][0]*p0+M[0][1]*u0_1+M[0][2]*u0_2;
  flux[1]=0.0;
  flux[2]=0.0;
}


void Wave_BC_pressure_imposed(void * cs, real * xpg, real * w, real *vnorm, real * flux){

  ContinuousSolver * ps=cs;
  real M[3][3];
  real BC[3][3];
  real lambda=-1e7;
  real mu=0;
  real p0=0,u0_1=0,u0_2=0;
  real x,y;
  x=xpg[0];
  y=xpg[1];

  real Win[ps->simu->fd[0].model.m];
  
  ps->simu->fd[0].model.ImposedData(xpg,ps->simu->dt,Win);
  p0=Win[0];
  u0_1=Win[1];
  u0_2=Win[2];

  real h=ps->FluxMatrix[0][1]/ps->simu->vmax;


  M[0][0]=lambda;
  M[0][1]=-mu*vnorm[0];
  M[0][2]=-mu*vnorm[1];
  M[1][0]=0;
  M[1][1]=0;
  M[1][2]=0;
  M[2][0]=0;
  M[2][1]=0;
  M[2][2]=0;

  BC[0][0]=h*M[0][0];
  BC[0][1]=(ps->FluxMatrix[0][1]*vnorm[0]+h*M[0][1]);
  BC[0][2]=(ps->FluxMatrix[0][2]*vnorm[1]+h*M[0][2]);
  BC[1][0]=(ps->FluxMatrix[1][0]*vnorm[0]+h*M[1][0]);
  BC[1][1]=h*M[1][1];
  BC[1][2]=h*M[1][2];
  BC[2][0]=(ps->FluxMatrix[2][0]*vnorm[1]+h*M[2][0]);
  BC[2][1]=h*M[2][1];
  BC[2][2]=h*M[2][2];

  flux[0]=(BC[0][0]*w[0]+BC[0][1]*w[1]+BC[0][2]*w[2])-h*(M[0][0]*p0+M[0][1]*u0_1+M[0][2]*u0_2);
  flux[1]=(BC[1][0]*w[0]+BC[1][1]*w[1]+BC[1][2]*w[2])-h*(M[1][0]*p0+M[1][1]*u0_1+M[1][2]*u0_2);
  flux[2]=(BC[2][0]*w[0]+BC[2][1]*w[1]+BC[2][2]*w[2])-h*(M[2][0]*p0+M[2][1]*u0_1+M[2][2]*u0_2);

  //flux[0]=M[0][0]*w[0]+M[0][1]*w[1]+M[0][2]*w[2] - M[0][0]*p0+M[0][1]*u0_1+M[0][2]*u0_2;
  //flux[1]=0.0;
  //flux[2]=0.0;
}


  
void Source_Assembly(void * cs){
  ContinuousSolver * ps=cs;
  field* f0 = &ps->simu->fd[0];
  real dt = ps->simu->dt;
  
  m=ps->nb_phy_vars;
  for(int ie = 0; ie < ps->nbel; ie++){  

    int iemacro = ie / (f0->raf[0] * f0->raf[1] * f0->raf[2]);
    int isubcell = ie % (f0->raf[0] * f0->raf[1] * f0->raf[2]);
    
     for(int ipg = 0;ipg < ps->nnodes; ipg++){
      real wpg;
      real xref[3];
      int ipgmacro= ipg + isubcell * ps->nnodes;
      ref_pg_vol(f0->deg,f0->raf,ipgmacro,xref,&wpg,NULL);
      real dtau[3][3],codtau[3][3];
      
      Ref2Phy(ps->simu->fd[iemacro].physnode,
	      xref,NULL,0,NULL,
	      dtau,codtau,NULL,NULL);
        
      real det = dot_product(dtau[0], codtau[0]);
      real Source[f0->model.m];
      real w[f0->model.m];
      for(int var =0; var < f0->model.m; var++){
	w[var] = f0->wn[var];
      }
      f0->model.Source(xref, f0->tnow,w, Source);
    
      int ino_dg = ipg + ie * ps->nnodes;
      int ino_fe = ps->dg_to_fe_index[ino_dg];
      for(int var =0; var < ps->nb_phy_vars; var++){
	int ipot_fe = ino_fe * ps->nb_phy_vars + var;
        real val = dt * Source[ps->list_of_var[var]] * wpg * det;
	ps->lsol.rhs[ipot_fe] += val;
      }
    }
  }
}
