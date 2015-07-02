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

  dh=(wR[0]-wL[0]);
  dhu=(wR[0]*wR[1]-wL[0]*wL[1]);
  dhv=(wR[0]*wR[2]-wL[0]*wL[2]);
  
  c=sqrt((g/2.0)*(wR[0]+wL[0]));
  uti=(wR[1]*sqrt(wR[0])+wL[1]*sqrt(wL[0]))/(sqrt(wR[0])+sqrt(wL[0]));
  vti=(wR[2]*sqrt(wR[0])+wL[2]*sqrt(wL[0]))/(sqrt(wR[0])+sqrt(wL[0]));

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
  
  centered_flux[0]= 0.5*(wL[0]*wL[1]+wR[0]*wR[1])*vnorm[0] + 0.5*(wL[0]*wL[2]+wR[0]*wR[2])*vnorm[1];
  
  centered_flux[1]= 0.5*(wL[0]*pow(wL[1],2.0)+0.5*g*pow(wL[0],2.0)+wR[0]*pow(wR[1],2.0)+0.5*g*pow(wR[0],2.0))*vnorm[0] + 0.5*(wL[0]*wL[1]*wL[2]+wR[0]*wR[1]*wR[2])*vnorm[1];
  centered_flux[2]= 0.5*(wL[0]*wL[1]*wL[2]+wR[0]*wR[1]*wR[2])*vnorm[0] + 0.5*(wL[0]*pow(wL[2],2.0)+0.5*g*pow(wL[0],2.0)+wR[0]*pow(wR[2],2.0)+0.5*g*pow(wR[0],2.0))*vnorm[1];
 
  flux[0]= centered_flux[0]-viscosity[0];
  
  flux[1]= centered_flux[1]-viscosity[1];
  flux[2]= centered_flux[2]-viscosity[2];

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
  real SL=0, SR=0;

  cL = sqrt(g*wL[0]);
  cR = sqrt(g*wR[0]);
  qn = 0.5 * ((wL[1]+wR[1])*vnorm[0]+(wL[2]+wR[2])*vnorm[1])+cL-cR;
  cstar = 0.5 * (cL + cR) + 0.25 * ((wL[1]-wR[1])*vnorm[0]+(wL[2]-wR[2])*vnorm[1]);

  SL=fmin((wL[1]*vnorm[0]+wL[2]*vnorm[1])-cL,qn-cstar);
  SR=fmax((wR[1]*vnorm[0]+wR[2]*vnorm[1])+cR,qn+cstar);
  
  viscosity[0]=SL*SR*(wR[0]-wL[0]);
  viscosity[1]=SL*SR*(wR[0]*wR[1]-wL[0]*wL[1]);
  viscosity[2]=SL*SR*(wR[0]*wR[2]-wL[0]*wL[2]);
  
  centered_flux[0]= (SR*wL[0]*wL[1]-SL*wR[0]*wR[1])*vnorm[0] +(SR*wL[0]*wL[2]-SL*wR[0]*wR[2])*vnorm[1];
  
  centered_flux[1]= (SR*wL[0]*pow(wL[1],2.0)+SR*0.5*g*pow(wL[0],2.0)-SL*wR[0]*pow(wR[1],2.0)-SL*0.5*g*pow(wR[0],2.0))*vnorm[0] + (SR*wL[0]*wL[1]*wL[2]-SL*wR[0]*wR[1]*wR[2])*vnorm[1];
  centered_flux[2]= (SR*wL[0]*wL[1]*wL[2]-SL*wR[0]*wR[1]*wR[2])*vnorm[0] + (SR*wL[0]*pow(wL[2],2.0)+SR*0.5*g*pow(wL[0],2.0)-SL*wR[0]*pow(wR[2],2.0)-SL*0.5*g*pow(wR[0],2.0))*vnorm[1];
 
  flux[0]= (1.0/(SR-SL))*(centered_flux[0]+viscosity[0]);
  flux[1]= (1.0/(SR-SL))*(centered_flux[1]+viscosity[1]);
  flux[2]= (1.0/(SR-SL))*(centered_flux[2]+viscosity[2]);

  //A comparative study of finite volume methods on unstructured meshes for simulation of 2D shallow water wave problems
  //Ji-Wen Wang a,b,∗, Ru-Xun Liu a
  
};


void ShallowWater_Rusanov_NumFlux(real wL[],real wR[],real* vnorm,real* flux){
  real centered_flux[3];
  real viscosity[3];
  real g=_GRAVITY;
  real cL=0, cR=0;
  real S=0;

  cL = sqrt(g*wL[0]);
  cR = sqrt(g*wR[0]);

  S=fmax(fabs(wR[1]*vnorm[0]+wR[2]*vnorm[1])+cR,fabs(wL[1]*vnorm[0]+wL[2]*vnorm[1])+cL);
  
  viscosity[0]=S*(wR[0]-wL[0]);
  viscosity[1]=S*(wR[0]*wR[1]-wL[0]*wL[1]);
  viscosity[2]=S*(wR[0]*wR[2]-wL[0]*wL[2]);
  
  centered_flux[0]= (wL[0]*wL[1]+wR[0]*wR[1])*vnorm[0] +(wL[0]*wL[2]+wR[0]*wR[2])*vnorm[1];
  
  centered_flux[1]= (wL[0]*pow(wL[1],2.0)+0.5*g*pow(wL[0],2.0)+wR[0]*pow(wR[1],2.0)+0.5*g*pow(wR[0],2.0))*vnorm[0] + (wL[0]*wL[1]*wL[2]+wR[0]*wR[1]*wR[2])*vnorm[1];
  centered_flux[2]= (wL[0]*wL[1]*wL[2]+wR[0]*wR[1]*wR[2])*vnorm[0] + (wL[0]*pow(wL[2],2.0)+0.5*g*pow(wL[0],2.0)+wR[0]*pow(wR[2],2.0)+0.5*g*pow(wR[0],2.0))*vnorm[1];
 
  flux[0]= 0.5*(centered_flux[0]-viscosity[0]);
  flux[1]= 0.5*(centered_flux[1]-viscosity[1]);
  flux[2]= 0.5*(centered_flux[2]-viscosity[2]);

  //A comparative study of finite volume methods on unstructured meshes for simulation of 2D shallow water wave problems
  //Ji-Wen Wang a,b,∗, Ru-Xun Liu a
  
};
