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


void Solution_computation(int Nwave,real * centers, real * solution,real eps,real c,real sigma,real time,int test);
void Rhs_computation(int Nwave,real * centers, real * solution,real eps,real c,real sigma,real time,int test);

real SolutionBC(int Nwave,int var,real x,real eps,real c,real sigma,real time,int test);

void Profil_matrix(LinearSolver *sky, int n);
void Assembly_matrix(LinearSolver *sky, int n);

int TestWaveDF(void){

  int test=0;
  Simulation simu;

  LinearSolver sky;

  

#ifdef PARALUTION 
  paralution_begin();
#endif 


  
  int Nwave=1000;
  int neq=2*Nwave;
  real L=1.0;
  real hw=L/Nwave;
  real time=0;
  real timefinal=0;
  real dt=0,c=0,eps=0,sigma=0;
  int numtest=0;
  int itermax=0;

  
  real * solution=NULL;  
  solution=malloc(neq*sizeof(real));
  real * rhs=NULL;  
  rhs=malloc(neq*sizeof(real));

  //------------- mesh ---------------//
  real * nodes;
  real *centers;
  real * volumes;
  real *dualvolumes;
  real * nodesBoundary;
  real * centersBoundary;

  nodes=malloc((Nwave+1)*sizeof(real));
  centers=malloc((Nwave)*sizeof(real));
  volumes=malloc((Nwave)*sizeof(real));
  dualvolumes=malloc((Nwave+1)*sizeof(real));
  centersBoundary=malloc(2*sizeof(real));
  nodesBoundary=malloc(2*sizeof(real));

  for(int i=0;i<Nwave+1;i++){
    nodes[i]=i*hw;
  }
  nodesBoundary[0]=-hw;
  nodesBoundary[1]=L+hw;

  for(int i=0;i<Nwave;i++){
    centers[i]=0.5*(nodes[i]+nodes[i+1]);
    volumes[i]=nodes[i+1]-nodes[i];
  }
  centersBoundary[0]=-0.5*hw;
  centersBoundary[1]=L+0.5*hw;

  for(int i=0;i<Nwave+1;i++){
    if(i==0){
      dualvolumes[i]=centers[0]-centersBoundary[0];
    }
    else if (i==Nwave){
      dualvolumes[i]=centersBoundary[1]-centers[Nwave];
      }
    else {
      dualvolumes[i]=centers[i+1]-centers[i];
    }
  }
  //------------- End mesh ---------------//

  //------------- Parameters ------------//
  simu.dt=0.001;
  eps=1.0;
  c=1.0;
  sigma=1.0;
  simu.vmax=c/eps;
  timefinal=0.004;
  itermax=timefinal/simu.dt;
  numtest=3;

  //------------- Parameters ------------//

 
  
  InitLinearSolver(&sky,neq,NULL,NULL);

  sky.solver_type=LU;
  sky.pc_type=NONE;//PHDF;//;JACOBI;
  sky.iter_max=50;
  sky.tol=1.e-11;
  sky.restart_gmres=20;

  Profil_matrix(&sky,neq);
 
  // once the nonzero positions are known allocate memory
  AllocateLinearSolver(&sky);

  for(int i=0;i<neq;i++){
    if (i==0){ 
      AddLinearSolver(&sky,0,0,1.0+(c*simu.dt)/(hw)); //pj
      AddLinearSolver(&sky,0,2,-(c*simu.dt)/(2.0*hw)); //pj+1
	    
      AddLinearSolver(&sky,0,1,0.0); //uj
      AddLinearSolver(&sky,0,3,+(c*simu.dt)/(2.0*hw)); //uj+1
    }
    else if (i==1){ 
      AddLinearSolver(&sky,1,1,1.0+(c*simu.dt)/(hw)); //uj
      AddLinearSolver(&sky,1,3,-(c*simu.dt)/(2.0*hw)); //uj+1
	   
      AddLinearSolver(&sky,1,0,0.0); //pj
      AddLinearSolver(&sky,1,2,+(c*simu.dt)/(2.0*hw)); //pj+1
    } 
    else if (i==neq-1){ 
      AddLinearSolver(&sky,neq-1,neq-1,1.0+(c*simu.dt)/(hw)); //uj
      AddLinearSolver(&sky,neq-1,neq-3,-(c*simu.dt)/(2.0*hw)); //uj-1
       
      AddLinearSolver(&sky,neq-1,neq-2,0.0); //pj
      AddLinearSolver(&sky,neq-1,neq-4,-(c*simu.dt)/(2.0*hw));//pj-1
    }
    else if (i==neq-2){ 
      AddLinearSolver(&sky,neq-2,neq-2,1.0+(c*simu.dt)/(hw)); //pj
      AddLinearSolver(&sky,neq-2,neq-4,-(c*simu.dt)/(2.0*hw)); //pj-1
       
      AddLinearSolver(&sky,neq-2,neq-1,0.0); //uj
      AddLinearSolver(&sky,neq-2,neq-3,-(c*simu.dt)/(2.0*hw));//uj-1
    } 
    else if (i % 2 ==0) {
      AddLinearSolver(&sky,i,i-2,-(c*simu.dt)/(2.0*hw)); //pj-1
      AddLinearSolver(&sky,i,i,1.0+(c*simu.dt)/(hw)); //pj
      AddLinearSolver(&sky,i,i+2,-(c*simu.dt)/(2.0*hw)); //pj+1
      
      AddLinearSolver(&sky,i,i-1,-(c*simu.dt)/(2.0*hw)); //uj-1
      AddLinearSolver(&sky,i,i+1,0.0); //uj
      AddLinearSolver(&sky,i,i+3,+(c*simu.dt)/(2.0*hw)); //uj+1
      
    }
    else { 
      AddLinearSolver(&sky,i,i-2,-(c*simu.dt)/(2.0*hw)); //uj-1
      AddLinearSolver(&sky,i,i,1.0+(c*simu.dt)/(hw)); //uj
      AddLinearSolver(&sky,i,i+2,-(c*simu.dt)/(2.0*hw)); //uj+1
      
      AddLinearSolver(&sky,i,i-3,-(c*simu.dt)/(2.0*hw)); //pj-1
      AddLinearSolver(&sky,i,i-1,0.0); //pj
      AddLinearSolver(&sky,i,i+1,+(c*simu.dt)/(2.0*hw)); //pj+1
    }
    sky.sol[i]=0.0;
    
  }

  
   //DisplayLinearSolver(&sky);
  time=0; 
  for(int tstep=0;tstep<itermax;tstep++){

    Solution_computation(Nwave,centers,solution,eps,c,sigma,time,numtest);
    Rhs_computation(Nwave,centers,rhs,eps,c,sigma,time,numtest);


    if(tstep==0){
      for(int i=0;i<neq;i++){
	sky.rhs[i]=solution[i]+simu.dt*rhs[i];
      }
    }
    else
      {
	for(int i=0;i<neq;i++){
	  sky.rhs[i]=sky.sol[i]+simu.dt*rhs[i];
	}
      }

    real pl=SolutionBC(Nwave,0,centersBoundary[0],eps,c,sigma,time,numtest);
    real pr=SolutionBC(Nwave,0,centersBoundary[1],eps,c,sigma,time,numtest);
    real ul=SolutionBC(Nwave,1,centersBoundary[0],eps,c,sigma,time,numtest);
    real ur=SolutionBC(Nwave,1,centersBoundary[1],eps,c,sigma,time,numtest);
    
    sky.rhs[0]=sky.rhs[0]+(c*simu.dt)/(2.0*hw)*pl+(c*simu.dt)/(2.0*hw)*ul;
    sky.rhs[1]=sky.rhs[1]+(c*simu.dt)/(2.0*hw)*ul+(c*simu.dt)/(2.0*hw)*pl;
    sky.rhs[neq-2]=sky.rhs[neq-2]+(c*simu.dt)/(2.0*hw)*pr-(c*simu.dt)/(2.0*hw)*ur; // pressure right BC
    sky.rhs[neq-1]=sky.rhs[neq-1]+(c*simu.dt)/(2.0*hw)*ur-(c*simu.dt)/(2.0*hw)*pr; // velocity right BC

    SolveLinearSolver(&sky,&simu);
    /*for(int i=0;i<Nwave;i++){
      printf("mmmmm p %d %f\n",i,sky.sol[2*i]);
    }
      for(int i=0;i<Nwave;i++){
      printf("mmmmm u %d %f \n",i,sky.sol[2*i+1]);
      }*/

    real error=0,norm=0;
    time=time+simu.dt;
    for(int i=0;i<Nwave;i++){
      error=error+volumes[i]*(sky.sol[2*i]-solution[2*i])*(sky.sol[2*i]-solution[2*i]);
      error=error+volumes[i]*(sky.sol[2*i+1]-solution[2*i+1])*(sky.sol[2*i+1]-solution[2*i+1]);
      norm=norm+volumes[i]*(solution[2*i])*(solution[2*i]);
      norm=norm+volumes[i]*(solution[2*i+1])*(solution[2*i+1]);
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

real SolutionBC(int Nwave,int var,real x,real eps,real c,real sigma,real time,int test){
  real pi=4.0*atan(1.0);
  real Coef=2.0*pi;
  real p,u;
  if(test == 1){
      p=10;
      u=2;
  }
  if(test == 2){
      p=1.0+x*x;
      u=2.0;
  }
  if(test == 3){
      p=-c*Coef*sqrt(2.0)*sin(c*Coef*sqrt(2.0)*time)*cos(Coef*x); // presssure
      u=c*Coef*cos(Coef*c*sqrt(2.0)*time)*sin(Coef*x); //velocity  
  }
  if(var==0){
    return p;
  }
  else if(var==1){
    return u;
  }
  else {
    return 0;
  }
}


void Solution_computation(int Nwave,real * centers, real * solution,real eps,real c,real sigma,real time,int test){
  real pi=4.0*atan(1.0);
  real Coef=2.0*pi;
  if(test == 1){
    for(int i=0;i<Nwave;i++){
      solution[2*i]=10;//-c*Coef*sqrt(2.0)*sin(c*Coef*sqrt(2.0)*time)*cos(Coef*(i+0.5*hw)); // presssure
      solution[2*i+1]=2;//c*Coef*cos(Coef*c*sqrt(2.0)*time)*sin(Coef*(i+0.5*hw)); //velocity
    }
  }
   if(test == 2){
    for(int i=0;i<Nwave;i++){
      solution[2*i]=1.0+centers[i]*centers[i];
      solution[2*i+1]=2.0;
    }
  }
  if(test == 3){
    for(int i=0;i<Nwave;i++){
      solution[2*i]=-c*Coef*sqrt(2.0)*sin(c*Coef*sqrt(2.0)*time)*cos(Coef*centers[i]); // presssure
      solution[2*i+1]=c*Coef*cos(Coef*c*sqrt(2.0)*time)*sin(Coef*centers[i]); //velocity
    }
  } 
}

void Rhs_computation(int Nwave,real * centers, real * rhs,real eps,real c,real sigma,real time,int test){
  real pi=4.0*atan(1.0);
  real Coef=2.0*pi;
  if(test == 1){
    for(int i=0;i<Nwave;i++){
      rhs[2*i]=0;
      rhs[2*i+1]=0;
    }
  }
  if(test == 2){
    for(int i=0;i<Nwave;i++){
      rhs[2*i]=0.0;
      rhs[2*i+1]=2.0*centers[i];
    }
  }
  if(test == 3){
    for(int i=0;i<Nwave;i++){
      rhs[2*i]=0; // presssure
      rhs[2*i+1]=0; //velocity
    }
  }
}


void Profil_matrix(LinearSolver *sky, int neq){

  for(int i=0;i<neq;i++){
    if (i==0){ 
      IsNonZero(sky,0,0);
      IsNonZero(sky,0,1);
      IsNonZero(sky,0,2);
      IsNonZero(sky,0,3);
    }
    else if (i==1){
      IsNonZero(sky,1,0);
      IsNonZero(sky,1,1);
      IsNonZero(sky,1,2);
      IsNonZero(sky,1,3);
    } 
    else if (i==neq-1){ 
      IsNonZero(sky,neq-1,neq-1);
      IsNonZero(sky,neq-1,neq-2);
      IsNonZero(sky,neq-1,neq-3);
      IsNonZero(sky,neq-1,neq-4);
    } 
    else if (i==neq-2){ 
      IsNonZero(sky,neq-2,neq-1);
      IsNonZero(sky,neq-2,neq-2);
      IsNonZero(sky,neq-2,neq-3);
      IsNonZero(sky,neq-2,neq-4);
    } 
    else if (i % 2 ==0) { 
      IsNonZero(sky,i,i);
      IsNonZero(sky,i,i+1);
    
      IsNonZero(sky,i,i-1);
      IsNonZero(sky,i,i-2);
    
      IsNonZero(sky,i,i+2);
      IsNonZero(sky,i,i+3);
    }
    else { 
      IsNonZero(sky,i,i);
      IsNonZero(sky,i,i-1);
    
      IsNonZero(sky,i,i-2);
      IsNonZero(sky,i,i-3);
    
      IsNonZero(sky,i,i+1);
      IsNonZero(sky,i,i+2);
    } 
  }
  
}


