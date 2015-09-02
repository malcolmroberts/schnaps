#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "test.h"
#include "schnaps.h"

#define _NN 5

int main(void) {
  
  // unit tests
    
  int resu=TestWaveDF();
	 

  if (resu) printf("Wave DF test OK !\n");
  else printf("Linear solver test failed !\n");

  return !resu;
} 




int TestWaveDF(void){

  int test=0,test1=1,test2=1,test3=1,test4=1,test5=1;
  Simulation simu;

  LinearSolver sky;

#ifdef PARALUTION 
  paralution_begin();
#endif 


  
  int Nwave=60;
  real hw=1.0/Nwave;
  real c=1,dt=1;
  for(int i=0;i<Nwave;i++){
    if(i % 2 ==0){
      solution[i]=10; // presssure
    }
    else {
      solution[i]=2; //velocity
    }
  }
  
  InitLinearSolver(&sky,2*Nwave,NULL,NULL);

  sky.solver_type = LU;
  sky.pc_type=NONE;
  sky.iter_max=20000;
  sky.tol=1.e-7;
  
  for(int i=0;i<Nwave;i++){
    if (i==0){ 
      IsNonZero(&sky,0,0);
      IsNonZero(&sky,0,1);
      IsNonZero(&sky,0,2);
      IsNonZero(&sky,0,3);
    } 
    else if (i==2*NPoisson-1){ 
      IsNonZero(&sky,2*Nwave-1,NPoisson-1);
      IsNonZero(&sky,2*Nwave-1,NPoisson-2);
      IsNonZero(&sky,2*Nwave-1,NPoisson-3);
      IsNonZero(&sky,2*Nwave-1,NPoisson-4);
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

  for(int i=0;i<Nwave;i++){
    if (i==0){ 
      AddLinearSolver(&sky,0,0,1);
      AddLinearSolver(&sky,0,1,1);
      AddLinearSolver(&sky,0,2,1);
      AddLinearSolver(&sky,0,3,1);
    } 
    else if (i==2*NPoisson-1){ 
      AddLinearSolver(&sky,2*Nwave-1,NPoisson-1,1);
      AddLinearSolver(&sky,2*Nwave-1,NPoisson-2,1);
      AddLinearSolver(&sky,2*Nwave-1,NPoisson-3,1);
      AddLinearSolver(&sky,2*Nwave-1,NPoisson-4,1);
    } 
    else if (i % 2 ==0) { 
      AddLinearSolver(&sky,i,i,1);
      AddLinearSolver(&sky,i,i+1,1);
    
      AddLinearSolver(&sky,i,i-1,1);
      AddLinearSolver(&sky,i,i-2,1);
    
      AddLinearSolver(&sky,i,i+2,1);
      AddLinearSolver(&sky,i,i+3,1);
    }
    else { 
      AddLinearSolver(&sky,i,i,1);
      AddLinearSolver(&sky,i,i-1,1);
    
      AddLinearSolver(&sky,i,i-2,1);
      AddLinearSolver(&sky,i,i-3,1);

      AddLinearSolver(&sky,i,i+1,1);
      AddLinearSolver(&sky,i,i+2,1);
    }
    sky.sol[i]=0.0;
    sky.rhs[i]=solution[i];
    
  }

  real bigvalw=1.e+15;
  SetLinearSolver(&sky,2*Nwave-1,2*Nwave-1,bigvalw);
  SetLinearSolver(&sky,2*Nwave-2,2*Nwave-2,bigvalw);
  SetLinearSolver(&sky,0,0,bigvalw);
  SetLinearSolver(&sky,1,1,bigvalw);
  sky.rhs[0]=bigvalw*sky.rhs[0];
  sky.rhs[1]=bigvalw*sky.rhs[1];
  sky.rhs[2*Nwave-1]=bigvalw*sky.rhs[2*Nwave-1];
  sky.rhs[2*Nwave-2]=bigvalw*sky.rhs[2*Nwave-2];

  for(int tstep=0;tstep<2*dt;tstep++){
    time=tstep*dt;
    SolveLinearSolver(&sky,&simu);
    
    for(int i=0;i<Nwave;i++){
      sky.rhs[i]=sky.sol[i];
    }

    real error=0;
    for(int i=0;i<2*Nwave;i++){
      error=error+(sky.sol[i]-solution[i])*(sky.sol[i]-solution[i]);
    }
    prtintf("error time \n",time,sqrt(error));
  }


  // deallocate memory
  FreeLinearSolver(&sky);

  test5 = test5 && (verr<5.e-2);
  printf("Error =%.12e\n",verr);

  if(test1==1 &&  test2==1 && test3==1 && test4==1) test=1;


#ifdef PARALUTION 
  paralution_end();
#endif 

  return test;

}
