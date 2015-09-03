#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "test.h"
#include "schnaps.h"
#include <stdlib.h>

int main(void) {
  
  // unit tests
    
  int resu=TestWaveDF();
	 

  if (resu) printf("Wave DF test OK !\n");
  else printf("Linear solver test failed !\n");

  return !resu;
} 




int TestWaveDF(void){

  int test=0;
  Simulation simu;

  LinearSolver sky;
  real * solution=NULL;
  

#ifdef PARALUTION 
  paralution_begin();
#endif 


  
  int Nwave=20;
  int neq=2*Nwave;
  real hw=1.0/Nwave;
  real time=0;
  real dt=0,c=0;
  real pi=4.0*atan(1.0);
  real Coef=2.0*pi;
  real itermax=5;

  simu.dt=0.1;
  simu.vmax=1;
  dt=simu.dt;
  c=simu.vmax;
  

  solution=malloc(neq*sizeof(real));
 
  
  InitLinearSolver(&sky,neq,NULL,NULL);

  sky.solver_type=LU;
  sky.pc_type=NONE;//PHDF;//;JACOBI;
  sky.iter_max=50;
  sky.tol=1.e-11;
  sky.restart_gmres=20;
  
  for(int i=0;i<neq;i++){
    if (i==0){ 
      IsNonZero(&sky,0,0);
      IsNonZero(&sky,0,1);
      IsNonZero(&sky,0,2);
      IsNonZero(&sky,0,3);
    }
    else if (i==1){
      IsNonZero(&sky,1,0);
      IsNonZero(&sky,1,1);
      IsNonZero(&sky,1,2);
      IsNonZero(&sky,1,3);
    } 
    else if (i==neq-1){ 
      IsNonZero(&sky,neq-1,neq-1);
      IsNonZero(&sky,neq-1,neq-2);
      IsNonZero(&sky,neq-1,neq-3);
      IsNonZero(&sky,neq-1,neq-4);
    } 
    else if (i==neq-2){ 
      IsNonZero(&sky,neq-2,neq-1);
      IsNonZero(&sky,neq-2,neq-2);
      IsNonZero(&sky,neq-2,neq-3);
      IsNonZero(&sky,neq-2,neq-4);
    } 
    else if (i % 2 ==0) { 
      IsNonZero(&sky,i,i);
      IsNonZero(&sky,i,i+1);
    
      IsNonZero(&sky,i,i-1);
      IsNonZero(&sky,i,i-2);
    
      IsNonZero(&sky,i,i+2);
      IsNonZero(&sky,i,i+3);
    }
    else { 
      IsNonZero(&sky,i,i);
      IsNonZero(&sky,i,i-1);
    
      IsNonZero(&sky,i,i-2);
      IsNonZero(&sky,i,i-3);
    
      IsNonZero(&sky,i,i+1);
      IsNonZero(&sky,i,i+2);
    } 
  }
  // once the nonzero positions are known allocate memory
  AllocateLinearSolver(&sky);

  for(int i=0;i<neq;i++){
    if (i==0){ 
      AddLinearSolver(&sky,0,0,1.0+(c*dt)/(hw)); //pj
      AddLinearSolver(&sky,0,2,-(c*dt)/(2.0*hw)); //pj+1
	    
      AddLinearSolver(&sky,0,1,0.0); //uj
      AddLinearSolver(&sky,0,3,+(c*dt)/(2.0*hw)); //uj+1
    }
    else if (i==1){ 
      AddLinearSolver(&sky,1,1,1.0+(c*dt)/(hw)); //uj
      AddLinearSolver(&sky,1,3,-(c*dt)/(2.0*hw)); //uj+1
	    
      AddLinearSolver(&sky,1,0,0.0); //pj
      AddLinearSolver(&sky,1,2,+(c*dt)/(2.0*hw)); //pj+1
    } 
    else if (i==neq-1){ 
      AddLinearSolver(&sky,neq-1,neq-1,1.0+(c*dt)/(hw)); //uj
      AddLinearSolver(&sky,neq-1,neq-3,-(c*dt)/(2.0*hw)); //uj-1
       
      AddLinearSolver(&sky,neq-1,neq-2,0.0); //pj
      AddLinearSolver(&sky,neq-1,neq-4,-(c*dt)/(2.0*hw));//pj-1
    }
    else if (i==neq-2){ 
      AddLinearSolver(&sky,neq-2,neq-2,1.0+(c*dt)/(hw)); //pj
      AddLinearSolver(&sky,neq-2,neq-4,-(c*dt)/(2.0*hw)); //pj-1
       
      AddLinearSolver(&sky,neq-2,neq-1,0.0); //uj
      AddLinearSolver(&sky,neq-2,neq-3,-(c*dt)/(2.0*hw));//uj-1
    } 
    else if (i % 2 ==0) {
      AddLinearSolver(&sky,i,i-2,-(c*dt)/(2.0*hw)); //pj-1
      AddLinearSolver(&sky,i,i,1.0+(c*dt)/(hw)); //pj
      AddLinearSolver(&sky,i,i+2,-(c*dt)/(2.0*hw)); //pj+1
      
      AddLinearSolver(&sky,i,i-1,-(c*dt)/(2.0*hw)); //uj-1
      AddLinearSolver(&sky,i,i+1,0.0); //uj
      AddLinearSolver(&sky,i,i+3,+(c*dt)/(2.0*hw)); //uj+1
      
    }
    else { 
      AddLinearSolver(&sky,i,i-2,-(c*dt)/(2.0*hw)); //uj-1
      AddLinearSolver(&sky,i,i,1.0+(c*dt)/(hw)); //uj
      AddLinearSolver(&sky,i,i+2,-(c*dt)/(2.0*hw)); //uj+1
      
      AddLinearSolver(&sky,i,i-3,-(c*dt)/(2.0*hw)); //pj-1
      AddLinearSolver(&sky,i,i-1,0.0); //pj
      AddLinearSolver(&sky,i,i+1,+(c*dt)/(2.0*hw)); //pj+1
    }
    sky.sol[i]=0.0;
    sky.rhs[i]=0.0;
    
  }


  for(int i=0;i<neq;i++){
    SetLinearSolver(&sky,neq-1,i,0.0);
    SetLinearSolver(&sky,neq-2,i,0.0);
    SetLinearSolver(&sky,0,i,0.0);
    SetLinearSolver(&sky,1,i,0.0);
  }
  SetLinearSolver(&sky,neq-1,neq-1,1.0);
  SetLinearSolver(&sky,neq-2,neq-2,1.0);
  SetLinearSolver(&sky,0,0,1.0);
  SetLinearSolver(&sky,1,1,1.0);


   //DisplayLinearSolver(&sky);
  time=0; 
  for(int tstep=0;tstep<itermax;tstep++){

    for(int i=0;i<Nwave;i++){
      solution[2*i]=10;//-c*Coef*sqrt(2.0)*sin(c*Coef*sqrt(2.0)*time)*cos(Coef*(i+0.5*hw)); // presssure
      solution[2*i+1]=2;//c*Coef*cos(Coef*c*sqrt(2.0)*time)*sin(Coef*(i+0.5*hw)); //velocity
    }

    if(tstep==0){
      for(int i=0;i<neq;i++){
	sky.rhs[i]=solution[i];
      }
    }
    else
      {
	for(int i=0;i<neq;i++){
	  sky.rhs[i]=sky.sol[i];
	}
      }
    
    
    sky.rhs[0]=1.0*solution[0];
    sky.rhs[1]=1.0*solution[1];
    sky.rhs[neq-1]=1.0*solution[neq-1];
    sky.rhs[neq-2]=1.0*solution[neq-2];
     
    SolveLinearSolver(&sky,&simu);

     //PhyBasedPC_waveDF(&sky,&simu,sky.sol,sky.rhs);

    real error=0,norm=0;
    time=time+dt;
    for(int i=0;i<neq;i++){
      error=error+(sky.sol[i]-solution[i])*(sky.sol[i]-solution[i]);
      norm=norm+(solution[i])*(solution[i]);
    }
    printf("error time %.13e %.13e \n",time,sqrt(error/norm));
  }

  // deallocate memory
  FreeLinearSolver(&sky);



#ifdef PARALUTION 
  paralution_end();
#endif 

  return test;

}



