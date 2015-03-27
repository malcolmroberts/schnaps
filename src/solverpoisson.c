#include "solverpoisson.h"
#include "geometry.h"
#include "quantities_collision.h"


void SolvePoisson(Field *f){

  // for the moment, works only for the 1d case
  assert(f->is1d);

  // assembly of the rigidity matrix

  Skyline sky;

  // number of equation of the Poisson solver
  // = number of nodes in the mesh
  int degx=f->interp.interp_param[1];
  int nelx=f->interp.interp_param[4];
  printf("degx=%d nelx=%d _MV=%d\n",degx,nelx,_MV);
  double xmin=0;
  double xmax=1;  // TO DO: compute the maximal x coordinate
  int neq=degx*nelx+1;
  
  // number of conservatives variables
  // = number of velocity glops + 1 (potential)
  int m=f->model.m;

  InitSkyline(&sky,neq);

  // compute the profile of the matrix
  for(int ie=0;ie<nelx;ie++){
    for(int iloc=0;iloc<degx+1;iloc++){
      for(int jloc=0;jloc<degx+1;jloc++){
	int ino=iloc + ie * degx;
	int jno=jloc + ie * degx;
	SwitchOn(&sky,ino,jno);
      }
    }
  }

  AllocateSkyline(&sky);

  // local matrix (assuming identical finite elements)
  double aloc[degx+1][degx+1];
  for(int iloc=0;iloc<degx+1;iloc++){
    for(int jloc=0;jloc<degx+1;jloc++){
      aloc[iloc][jloc]=0;
    }
  }
  for(int ipg=0;ipg<degx+1;ipg++){
    double omega=wglop(degx,ipg);
    for(int iloc=0;iloc<degx+1;iloc++){
      for(int jloc=0;jloc<degx+1;jloc++){
	double dxi=dlag(degx,iloc,ipg);
	double dxj=dlag(degx,jloc,ipg);
	aloc[iloc][jloc]+=dxi*dxj*omega*nelx;
      }
    }
  }

  // assembly of the matrix
  for(int ie=0;ie<nelx;ie++){
    for(int iloc=0;iloc<degx+1;iloc++){
      for(int jloc=0;jloc<degx+1;jloc++){
	int ino=iloc + ie * degx;
	int jno=jloc + ie * degx;
	double val = aloc[iloc][jloc];
	SetSkyline(&sky,ino,jno,val);
      }
    }
  }


  // dirichlet boundary condition at the first and last location
  SetSkyline(&sky,0,0,1e20);
  SetSkyline(&sky,neq-1,neq-1,1e20);


  //DisplaySkyline(&sky);

  FactoLU(&sky);


  // charge computation
  Computation_charge_density(f);

    // source assembly 
  double source[neq];
  for(int i=0;i<neq;i++){
    source[i]=0;
  }

  for(int ie=0;ie<nelx;ie++){
    for(int iloc=0;iloc<degx+1;iloc++){
      double omega=wglop(degx,iloc);
      int ino=iloc + ie * degx;  
      int imem=f->varindex(f->interp_param,0,iloc+ie*(degx+1),_MV+2);
      double charge=f->wn[imem];
      printf("charge=%f\n",charge);
      source[ino]+= charge*omega/nelx;
    }
  }


  // assert(1==2);

  double sol[neq];
  SolveSkyline(&sky,source,sol);

  for(int i=0;i<neq;i++){
    printf("sol %d = %f\n",i,sol[i]);
  }

  // now put the solution at the right place
  for(int ie=0;ie<nelx;ie++){
     for(int ipg=0;ipg<degx+1;ipg++){
       // position in the continuous vector
       int ino=ipg + ie * degx;
       // position in the DG vector
       int imem=f->varindex(f->interp_param,0,ipg+ie*(degx+1),_MV);
       f->wn[imem]=sol[ino];
     }
  }
	
  FreeSkyline(&sky);


  Compute_electric_field(f);
}
