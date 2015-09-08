#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "test.h"
#include "schnaps.h"
#include <stdlib.h>

typedef struct MeshDF{
  int N;
  real h;
  real L;
  real * nodes;
  real *centers;
  real * volumes;
  real *dualvolumes;
  real * nodesBoundary;
  real * centersBoundary;

} MeshDF;

typedef struct Data{
  real time;
  real timefinal;
  real dt;
  real c;
  real eps;
  real sigma;
  int test;
  int itermax;
  

} Data;

void Solution_computation(int Nwave,real * centers, real * solution,Data * d);
void Rhs_computation(int Nwave,real * centers, real * solution,Data * d);

real SolutionBC(int Nwave,int var,real x,Data * d);

void Profil_matrix(LinearSolver *sky, int n);
void Assembly_matrix(LinearSolver *sky, int n);

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

#ifdef PARALUTION 
  paralution_begin();
#endif 

  //------------- mesh ---------------//

  MeshDF m; 
  m.N=1000;
  m.L=1.0;
  m.h=m.L/m.N;
  int neq=2*m.N;

  
  real * solution=NULL;  
  solution=malloc(neq*sizeof(real));
  real * rhs=NULL;  
  rhs=malloc(neq*sizeof(real));

  m.nodes=malloc((m.N+1)*sizeof(real));
  m.centers=malloc((m.N)*sizeof(real));
  m.volumes=malloc((m.N)*sizeof(real));
  m.dualvolumes=malloc((m.N+1)*sizeof(real));
  m.centersBoundary=malloc(2*sizeof(real));
  m.nodesBoundary=malloc(2*sizeof(real));

  for(int i=0;i<m.N+1;i++){
    m.nodes[i]=i*m.h;
  }
  m.nodesBoundary[0]=-m.h;
  m.nodesBoundary[1]=m.L+m.h;

  for(int i=0;i<m.N;i++){
    m.centers[i]=0.5*(m.nodes[i]+m.nodes[i+1]);
    m.volumes[i]=m.nodes[i+1]-m.nodes[i];
  }
  m.centersBoundary[0]=-0.5*m.h;
  m.centersBoundary[1]=m.L+0.5*m.h;

  for(int i=0;i<m.N+1;i++){
    if(i==0){
      m.dualvolumes[i]=m.centers[0]-m.centersBoundary[0];
    }
    else if (i==m.N){
      m.dualvolumes[i]=m.centersBoundary[1]-m.centers[m.N];
      }
    else {
      m.dualvolumes[i]=m.centers[i+1]-m.centers[i];
    }
  }
  //------------- End mesh ---------------//

  //------------- Parameters ------------//
  Data d;
  
  d.dt=0.0001;
  d.eps=1.0;
  d.c=1.0;
  d.sigma=1.0;
  d.timefinal=0.1;
  d.itermax=d.timefinal/d.dt;
  d.test=3;
  d.time=0;

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
      AddLinearSolver(&sky,0,0,1.0+(d.c*d.dt)/(m.volumes[0])); //pj
      AddLinearSolver(&sky,0,2,-(d.c*d.dt)/(2.0*m.volumes[0])); //pj+1
	    
      AddLinearSolver(&sky,0,1,0.0); //uj
      AddLinearSolver(&sky,0,3,+(d.c*d.dt)/(2.0*m.volumes[0])); //uj+1
    }
    else if (i==1){ 
      AddLinearSolver(&sky,1,1,1.0+(d.c*d.dt)/(m.volumes[i])); //uj
      AddLinearSolver(&sky,1,3,-(d.c*d.dt)/(2.0*m.volumes[i])); //uj+1
	   
      AddLinearSolver(&sky,1,0,0.0); //pj
      AddLinearSolver(&sky,1,2,+(d.c*d.dt)/(2.0*m.volumes[i])); //pj+1
    } 
    else if (i==neq-1){ 
      AddLinearSolver(&sky,neq-1,neq-1,1.0+(d.c*d.dt)/(m.volumes[neq-1-m.N])); //uj
      AddLinearSolver(&sky,neq-1,neq-3,-(d.c*d.dt)/(2.0*m.volumes[neq-1-m.N])); //uj-1
       
      AddLinearSolver(&sky,neq-1,neq-2,0.0); //pj
      AddLinearSolver(&sky,neq-1,neq-4,-(d.c*d.dt)/(2.0*m.volumes[neq-1-m.N]));//pj-1
    }
    else if (i==neq-2){ 
      AddLinearSolver(&sky,neq-2,neq-2,1.0+(d.c*d.dt)/(m.volumes[neq-2-m.N])); //pj
      AddLinearSolver(&sky,neq-2,neq-4,-(d.c*d.dt)/(2.0*m.volumes[neq-2-m.N])); //pj-1
       
      AddLinearSolver(&sky,neq-2,neq-1,0.0); //uj
      AddLinearSolver(&sky,neq-2,neq-3,-(d.c*d.dt)/(2.0*m.volumes[neq-2-m.N]));//uj-1
    } 
    else if (i % 2 ==0) {
      AddLinearSolver(&sky,i,i-2,-(d.c*d.dt)/(2.0*m.h)); //pj-1
      AddLinearSolver(&sky,i,i,1.0+(d.c*d.dt)/(m.h)); //pj
      AddLinearSolver(&sky,i,i+2,-(d.c*d.dt)/(2.0*m.h)); //pj+1
      
      AddLinearSolver(&sky,i,i-1,-(d.c*d.dt)/(2.0*m.h)); //uj-1
      AddLinearSolver(&sky,i,i+1,0.0); //uj
      AddLinearSolver(&sky,i,i+3,+(d.c*d.dt)/(2.0*m.h)); //uj+1
      
    }
    else { 
      AddLinearSolver(&sky,i,i-2,-(d.c*d.dt)/(2.0*m.h)); //uj-1
      AddLinearSolver(&sky,i,i,1.0+(d.c*d.dt)/(m.h)); //uj
      AddLinearSolver(&sky,i,i+2,-(d.c*d.dt)/(2.0*m.h)); //uj+1
      
      AddLinearSolver(&sky,i,i-3,-(d.c*d.dt)/(2.0*m.h)); //pj-1
      AddLinearSolver(&sky,i,i-1,0.0); //pj
      AddLinearSolver(&sky,i,i+1,+(d.c*d.dt)/(2.0*m.h)); //pj+1
    }
    sky.sol[i]=0.0;
    
  }

  
   //DisplayLinearSolver(&sky);
  d.time=0; 
  for(int tstep=0;tstep<d.itermax;tstep++){

    Solution_computation(m.N,m.centers,solution,&d);
    Rhs_computation(m.N,m.centers,rhs,&d);


    if(tstep==0){
      for(int i=0;i<neq;i++){
	sky.rhs[i]=solution[i]+d.dt*rhs[i];
      }
    }
    else
      {
	for(int i=0;i<neq;i++){
	  sky.rhs[i]=sky.sol[i]+d.dt*rhs[i];
	}
      }

    real pl=SolutionBC(m.N,0,m.centersBoundary[0],&d);
    real pr=SolutionBC(m.N,0,m.centersBoundary[1],&d);
    real ul=SolutionBC(m.N,1,m.centersBoundary[0],&d);
    real ur=SolutionBC(m.N,1,m.centersBoundary[1],&d);
    //printf("mmp %d %f %f %f %f\n ",numtest,pl,pr,ul,ur);
    
    sky.rhs[0]=sky.rhs[0]+(d.c*d.dt)/(2.0*m.h)*pl+(d.c*d.dt)/(2.0*m.h)*ul;
    sky.rhs[1]=sky.rhs[1]+(d.c*d.dt)/(2.0*m.h)*ul+(d.c*d.dt)/(2.0*m.h)*pl;
    sky.rhs[neq-2]=sky.rhs[neq-2]+(d.c*d.dt)/(2.0*m.h)*pr-(d.c*d.dt)/(2.0*m.h)*ur; // pressure right BC
    sky.rhs[neq-1]=sky.rhs[neq-1]+(d.c*d.dt)/(2.0*m.h)*ur-(d.c*d.dt)/(2.0*m.h)*pr; // velocity right BC

    Advanced_SolveLinearSolver(&sky,&simu);

    real error=0,norm=0;
    d.time=d.time+d.dt;
    for(int i=0;i<m.N;i++){
      error=error+m.volumes[i]*(sky.sol[2*i]-solution[2*i])*(sky.sol[2*i]-solution[2*i]);
      error=error+m.volumes[i]*(sky.sol[2*i+1]-solution[2*i+1])*(sky.sol[2*i+1]-solution[2*i+1]);
      norm=norm+m.volumes[i]*(solution[2*i])*(solution[2*i]);
      norm=norm+m.volumes[i]*(solution[2*i+1])*(solution[2*i+1]);
    }
    printf("error time %.13e %.13e \n",d.time,sqrt(error/norm));
  }

  // deallocate memory
  FreeLinearSolver(&sky);



#ifdef PARALUTION 
  paralution_end();
#endif 

  return test;

}

real SolutionBC(int Nwave,int var,real x,Data * d){
  real pi=4.0*atan(1.0);
  real Coef=2.0*pi;
  real p,u;
  if(d->test == 1){
      p=10;
      u=2;
  }
  if(d->test == 2){
      p=1.0+x*x;
      u=2.0-x;
  }
  if(d->test == 3){
      p=-d->c*Coef*sin(d->c*Coef*d->time)*cos(Coef*x); // presssure
      u=d->c*Coef*cos(Coef*d->c*d->time)*sin(Coef*x); //velocity  
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


void Solution_computation(int Nwave,real * centers, real * solution, Data * d){
  real pi=4.0*atan(1.0);
  real Coef=2.0*pi;
  if(d->test == 1){
    for(int i=0;i<Nwave;i++){
      solution[2*i]=10;//-c*Coef*sqrt(2.0)*sin(c*Coef*sqrt(2.0)*time)*cos(Coef*(i+0.5*volumes[i])); // presssure
      solution[2*i+1]=2;//c*Coef*cos(Coef*c*sqrt(2.0)*time)*sin(Coef*(i+0.5*volumes[i])); //velocity
    }
  }
   if(d->test == 2){
    for(int i=0;i<Nwave;i++){
      solution[2*i]=1.0+centers[i]*centers[i];
      solution[2*i+1]=2.0-centers[i];
    }
  }
  if(d->test == 3){
    for(int i=0;i<Nwave;i++){
      solution[2*i]=-d->c*Coef*sin(d->c*Coef*d->time)*cos(Coef*centers[i]); // presssure
      solution[2*i+1]=d->c*Coef*cos(Coef*d->c*d->time)*sin(Coef*centers[i]); //velocity
    }
  } 
}

void Rhs_computation(int Nwave,real * centers, real * rhs,Data *d){
  real pi=4.0*atan(1.0);
  real Coef=2.0*pi;
  if(d->test == 1){
    for(int i=0;i<Nwave;i++){
      rhs[2*i]=0;
      rhs[2*i+1]=0;
    }
  }
  if(d->test == 2){
    for(int i=0;i<Nwave;i++){
      rhs[2*i]=-1.0;
      rhs[2*i+1]=2.0*centers[i];
    }
  }
  if(d->test == 3){
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


