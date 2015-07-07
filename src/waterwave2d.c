#include "waterwave2d.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>


void Wave_Upwind_NumFlux(real wL[],real wR[],real* vnorm,real* flux){
  real flux_temp=0;
  
  flux[0]=0.5*((wL[1]+wR[1])*vnorm[0] + (wL[2]+wR[2])*vnorm[1])+0.5*(wL[0]-wR[0]);
  
  flux_temp=0.5*(wL[0]+wR[0]) + 0.5*((wL[1]-wR[1])*vnorm[0] + (wL[2]-wR[2])*vnorm[1]);
  flux[1]=flux_temp*vnorm[0];
  flux[2]=flux_temp*vnorm[1];
 

  flux[0]=_SPEED_WAVE*flux[0];
  flux[1]=_SPEED_WAVE*flux[1];
  flux[2]=_SPEED_WAVE*flux[2];
  
};




void ShallowWater_Roe_NumFlux(real wL[],real wR[],real* vnorm,real* flux){
  real centered_flux[3];
  real e1[3], e2[3], e3[3];
  real viscosity[3];
  real alpha =0;
  real g=_GRAVITY;
  real alpha1=0, alpha2=0, alpha3=0,c=0;
  real uti=0,vti=0;
  real dh=0,dhu=0,dhv=0;
  real a1=0,a2=0,a3=0;
  real uL=0, uR=0, vL=0, vR=0;

  uL=wL[1]/wL[0];
  uR=wR[1]/wR[0];

  vL=wL[2]/wL[0];
  vR=wR[2]/wR[0];

  dh=(wR[0]-wL[0]);
  dhu=(wR[1]-wL[1]);
  dhv=(wR[2]-wL[2]);
  
  c=sqrt((g/2.0)*(wR[0]+wL[0]));
  uti=(uR*sqrt(wR[0])+uL*sqrt(wL[0]))/(sqrt(wR[0])+sqrt(wL[0]));
  vti=(vR*sqrt(wR[0])+vL*sqrt(wL[0]))/(sqrt(wR[0])+sqrt(wL[0]));

  a1=uti*vnorm[0]+vti*vnorm[1]+c;
  a2=uti*vnorm[0]+vti*vnorm[1];
  a3=uti*vnorm[0]+vti*vnorm[1]-c;

  e1[0]=1;
  e1[1]=uti+c*vnorm[0];
  e1[2]=vti+c*vnorm[1];

  e2[0]=0;
  e2[1]=-c*vnorm[1];
  e2[2]=c*vnorm[0];

  e3[0]=1;
  e3[1]=uti-c*vnorm[0];
  e3[2]=vti-c*vnorm[1];
  
  
  alpha1=0.5*dh+(1.0/(2.0*c))*(dhu*vnorm[0]+dhv*vnorm[1]-(uti*vnorm[0]-vti*vnorm[1])*dh);
  alpha3=alpha1;
  alpha2=(1.0/c)*((dhv-vti*dh)*vnorm[0]-(dhu-uti*dh)*vnorm[1]);

  viscosity[0]=alpha1*fabs(a1)*e1[0]+alpha2*fabs(a2)*e2[0]+alpha3*fabs(a3)*e3[0];
  viscosity[1]=alpha1*fabs(a1)*e1[1]+alpha2*fabs(a2)*e2[1]+alpha3*fabs(a3)*e3[1];
  viscosity[2]=alpha1*fabs(a1)*e1[2]+alpha2*fabs(a2)*e2[2]+alpha3*fabs(a3)*e3[2];
  
  centered_flux[0]= (wL[1]+wR[1])*vnorm[0] +(wL[2]+wR[2])*vnorm[1];
  
  centered_flux[1]= (wL[1]*uL+0.5*g*pow(wL[0],2.0)+wR[1]*uR+0.5*g*pow(wR[0],2.0))*vnorm[0] + (wL[1]*vL+wR[1]*vR)*vnorm[1];
  centered_flux[2]= (wL[1]*vL+wR[1]*vR)*vnorm[0] + (wL[2]*vL+0.5*g*pow(wL[0],2.0)+wR[2]*vR+0.5*g*pow(wR[0],2.0))*vnorm[1];
 
  flux[0]= 0.5*centered_flux[0]-viscosity[0];
  
  flux[1]= 0.5*centered_flux[1]-viscosity[1];
  flux[2]= 0.5*centered_flux[2]-viscosity[2];
  flux[3]= 0.0;
  flux[4]= 0.0;
  flux[5]= 0.0;

  //A comparative study of finite volume methods on unstructured meshes for simulation of 2D shallow water wave problems
  //Ji-Wen Wang a,b,∗, Ru-Xun Liu a
  
};


void ShallowWater_HLL_NumFlux(real wL[],real wR[],real* vnorm,real* flux){
  real centered_flux[3];
  real e1[3], e2[3], e3[3];
  real viscosity[3];
  real alpha =0;
  real g=_GRAVITY;
  real cL=0, cR=0, cstar=0,qn=0;
  real SL=0, SR=0, uL=0, uR=0, vL=0, vR=0;

  uL=wL[1]/wL[0];
  uR=wR[1]/wR[0];

  vL=wL[2]/wL[0];
  vR=wR[2]/wR[0];
 
  cL = sqrt(g*wL[0]);
  cR = sqrt(g*wR[0]);
  qn = 0.5 * ((uL+uR)*vnorm[0]+(vL+vR)*vnorm[1])+cL-cR;
  cstar = 0.5 * (cL + cR) + 0.25 * ((uL-uR)*vnorm[0]+(vL-vR)*vnorm[1]);

  SL=fmin((uL*vnorm[0]+vL*vnorm[1])-cL,qn-cstar);
  SR=fmax((uR*vnorm[0]+vR*vnorm[1])+cR,qn+cstar);
  
  viscosity[0]=SL*SR*(wR[0]-wL[0]);
  viscosity[1]=SL*SR*(wR[1]-wL[1]);
  viscosity[2]=SL*SR*(wR[2]-wL[2]);
  
  centered_flux[0]= (SR*wL[1]-SL*wR[1])*vnorm[0] +(SR*wL[2]-SL*wR[2])*vnorm[1];
  
  centered_flux[1]= (SR*wL[1]*uL+SR*0.5*g*pow(wL[0],2.0)-SL*wR[1]*uR-SL*0.5*g*pow(wR[0],2.0))*vnorm[0] + (SR*wL[1]*vL-SL*wR[1]*vR)*vnorm[1];
  centered_flux[2]= (SR*wL[1]*vL-SL*wR[1]*vR)*vnorm[0] + (SR*vL*wL[2]+SR*0.5*g*pow(wL[0],2.0)-SL*vR*wR[2]-SL*0.5*g*pow(wR[0],2.0))*vnorm[1];
 
  flux[0]= (1.0/(SR-SL))*(centered_flux[0]+viscosity[0]);
  flux[1]= (1.0/(SR-SL))*(centered_flux[1]+viscosity[1]);
  flux[2]= (1.0/(SR-SL))*(centered_flux[2]+viscosity[2]);
  flux[3]= 0.0;
  flux[4]= 0.0;
  flux[5]= 0.0;

  //A comparative study of finite volume methods on unstructured meshes for simulation of 2D shallow water wave problems
  //Ji-Wen Wang a,b,∗, Ru-Xun Liu a
  
};


void ShallowWater_Rusanov_NumFlux(real wL[],real wR[],real* vnorm,real* flux){
  real centered_flux[3];
  real viscosity[3];
  real g=_GRAVITY;
  real cL=0, cR=0, uL=0, uR=0, vL=0, vR=0;
  real hL=0, hR=0, pL=0, pR=0;
  real S=0;

  cL = sqrt(g*wL[0]);
  cR = sqrt(g*wR[0]);

  hL=wL[0];
  hR=wR[0];

  uL=wL[1]/wL[0];
  uR=wR[1]/wR[0];

  vL=wL[2]/wL[0];
  vR=wR[2]/wR[0];

  pL=0.5*g*pow(wL[0],2.0);
  pR=0.5*g*pow(wR[0],2.0);

  S=fmax(fabs(uR*vnorm[0]+vR*vnorm[1])+cR,fabs(uL*vnorm[0]+vL*vnorm[1])+cL);
  
  viscosity[0]=S*(hR-hL);
  viscosity[1]=S*(hR*uR-hL*uL);
  viscosity[2]=S*(hR*vR-hL*vL);

  
  centered_flux[0]= (hR*uR+hL*uL)*vnorm[0] +(hR*vR-hL*vL)*vnorm[1];
  
  centered_flux[1]= (hL*uL*uL+pL+hR*uR*uR+pR)*vnorm[0] + (hL*uL*vL+hR*uR*vR)*vnorm[1];
  centered_flux[2]= (hL*uL*vL+hR*uR*vR)*vnorm[0] + (hL*vL*vL+pL+hR*vR*vR+pR)*vnorm[1];
 
  flux[0]= 0.5*(centered_flux[0]-viscosity[0]);
  flux[1]= 0.5*(centered_flux[1]-viscosity[1]);
  flux[2]= 0.5*(centered_flux[2]-viscosity[2]);
  flux[3]= 0.0;
  flux[4]= 0.0;
  flux[5]= 0.0;
  
  //A comparative study of finite volume methods on unstructured meshes for simulation of 2D shallow water wave problems
  //Ji-Wen Wang a,b,∗, Ru-Xun Liu a
  
};


void ShallowWater_classical_SourceTerm(const real *x, const real t, const real *w, real *source){
  real g=_GRAVITY;
  real hL=0, hR=0, uL=0, uR=0, vL=0, vR=0;
  real S=0;

  hL = w[0];


  uL=w[1]/w[0];
  vL=w[2]/w[0];

 
  source[0]= 0.0;
  source[1]= -g*hL*w[4];
  source[2]= -g*hL*w[5];
  source[3]= 0.0;
  source[4]= 0.0;
  source[5]= 0.0;
 

};


void ShallowWater_HLLWB_NumFlux(real wL[],real wR[],real* vnorm,real* flux){
  real centered_flux[3];
  real viscosity[3];
  real g=_GRAVITY;
  real cL=0, cR=0, cstar=0,qn=0;
  real SL=0, SR=0, uL=0, uR=0, vL=0, vR=0;
  real bL=0,bR=0,hL=0,hR=0,pR=0,pL=0;

  hL=wL[0];
  hR=wR[0];

  uL=wL[1]/wL[0];
  uR=wR[1]/wR[0];

  vL=wL[2]/wL[0];
  vR=wR[2]/wR[0];

  pL=0.5*g*pow(wL[0],2.0);
  pR=0.5*g*pow(wR[0],2.0);

  cL = sqrt(g*wL[0]);
  cR = sqrt(g*wR[0]);
  
  bL=wL[3];
  bR=wR[3];

  printf("bl br %.15e %.15e \n",pL,pR);

 
  qn = 0.5 * ((uL+uR)*vnorm[0]+(vL+vR)*vnorm[1])+cL-cR;
  cstar = 0.5 * (cL + cR) + 0.25 * ((uL-uR)*vnorm[0]+(vL-vR)*vnorm[1]);

  SL=fmin((uL*vnorm[0]+vL*vnorm[1])-cL,qn-cstar);
  SR=fmax((uR*vnorm[0]+vR*vnorm[1])+cR,qn+cstar);
  
  viscosity[0]=SL*SR*(hR-hL)+SL*SR*(bR-bL);
  centered_flux[0]= (SR*wL[1]-SL*wR[1])*vnorm[0] +(SR*wL[2]-SL*wR[2])*vnorm[1];
  
   viscosity[1]=SL*SR*(hR*uR-hL*uL);//-SL*0.5*g*(hL+hR)*(bR-bL);//
   centered_flux[1]= (SR*hL*uL*uL+SR*pL-SL*uR*uR-SL*pR-SL*0.5*g*(hL+hR)*(bR-bL))*vnorm[0] + (SR*hL*uL*vL-SL*hR*uR*vR)*vnorm[1];

    viscosity[2]=SL*SR*(hR*vR-hL*vL);//-SL*0.5*g*(hL+hR)*(bR-bL);
    centered_flux[2]= (SR*hL*uL*vL-SL*hR*uR*vR)*vnorm[0] + (SR*hL*vL*vL+SR*pL-SL*vR*vR-SL*pR-SL*0.5*g*(hL+hR)*(bR-bL))*vnorm[1];
 


  
  flux[0]= (1.0/(SR-SL))*(centered_flux[0]+viscosity[0]); 
  flux[1]= (1.0/(SR-SL))*(centered_flux[1]+viscosity[1]); 
  flux[2]= (1.0/(SR-SL))*(centered_flux[2]+viscosity[2]); 
   flux[3]= 0.0; 
   flux[4]= 0.0; 
   flux[5]= 0.0; 

   //printf(" spleed %f %f\n",SL,SR); 
   //printf(" fluw %f %f %f\n",flux[0],flux[1],flux[2]);
   

  //A simple well-balanced and positive numerical scheme for the shallow-water system
//Emmanuel Audusse, Christophe Chalons, Philippe Ung
  
};


void ShallowWater_HLLWB_SourceTerm(const real *x, const real t, const real *w, real *source){
  real g=_GRAVITY;
  real hL=0, hR=0, uL=0, uR=0, vL=0, vR=0;
  real S=0;

  hL = w[0];


  uL=w[1]/w[0];
  vL=w[2]/w[0];

 
  source[0]= 0.0;
  source[1]= 0.0;
  source[2]= 0.0;
  source[3]= 0.0;
  source[4]= 0.0;
  source[5]= 0.0;
 

};
