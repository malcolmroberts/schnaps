// vim: set sw=4 ts=4 sts=4 et tw=78 foldmarker={{{,}}} foldlevel=0 foldmethod=marker:

#include "model.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>

#define ONE_OVER_SQRT_3 (0.57735026918962584)
#define ONE_OVER_SQRT_2 (0.707106781186547524400844362105)

#define _CH (5)
#define _GAM (1.666666666666)

const double transport_v[] = {
  ONE_OVER_SQRT_3,
  ONE_OVER_SQRT_3,
  ONE_OVER_SQRT_3};

//const double transport_v[] = {1,0,0};

const double transport_v2d[] = {
  ONE_OVER_SQRT_2,
  ONE_OVER_SQRT_2,
  0};


// Transport {{{
void TransportNumFlux(double wL[],double wR[],double* vnorm,double* flux){

  double vn =
    transport_v[0] * vnorm[0] +
    transport_v[1] * vnorm[1] +
    transport_v[2] * vnorm[2];

   double vnp = vn>0 ? vn : 0;
   double vnm = vn-vnp;

   flux[0] = vnp * wL[0] + vnm * wR[0];

};

void TransportNumFlux2d(double wL[],double wR[],double* vnorm,double* flux){

  double vn =
    transport_v2d[0] * vnorm[0] +
    transport_v2d[1] * vnorm[1] +
    transport_v2d[2] * vnorm[2];

   double vnp = vn>0 ? vn : 0;
   double vnm = vn-vnp;

   flux[0] = vnp * wL[0] + vnm * wR[0];
   /* if (fabs(vnorm[2])>1e-6){ */
   /*   printf("vnds %lf %lf %lf \n",vnorm[0],vnorm[1],vnorm[2]); */
   /* } */
   // verify that 2d computations are actually
   // activated
   assert(fabs(vnorm[2])<1e-8);


};

void TransportBoundaryFlux(double x[3],double t,double wL[],double* vnorm,
			   double* flux){
  double wR[1];
  TransportImposedData(x,t,wR);
  TransportNumFlux(wL,wR,vnorm,flux);
};

void TransportBoundaryFlux2d(double x[3],double t,double wL[],double* vnorm,
			   double* flux){
  double wR[1];
  TransportImposedData2d(x,t,wR);
  TransportNumFlux2d(wL,wR,vnorm,flux);
};

void TransportInitData(double x[3],double w[]){

  double t=0;
  TransportImposedData(x,t,w);

};

void TransportInitData2d(double x[3],double w[]){

  double t=0;
  TransportImposedData2d(x,t,w);

};


void TransportImposedData(double x[3],double t,double w[]){

  double vx =
    transport_v[0] * x[0] +
    transport_v[1] * x[1] +
    transport_v[2] * x[2];

  double xx = vx - t;

  w[0]=cos(xx);
};

void TransportImposedData2d(double x[3],double t,double w[]){

  double vx =
    transport_v2d[0] * x[0] +
    transport_v2d[1] * x[1] +
    transport_v2d[2] * x[2];

  double xx = vx - t;

  w[0]=cos(xx);
};
// }}}

// Test Transport {{{
void TestTransportBoundaryFlux(double x[3],double t,double wL[],double* vnorm,
			   double* flux){
  double wR[1];
  TestTransportImposedData(x,t,wR);
  TransportNumFlux(wL,wR,vnorm,flux);
};

void TestTransportBoundaryFlux2d(double x[3],double t,double wL[],double* vnorm,
			   double* flux){
  double wR[1];
  TestTransportImposedData2d(x,t,wR);
  TransportNumFlux2d(wL,wR,vnorm,flux);
};

void TestTransportInitData(double x[3],double w[]){

  double t=0;
  TestTransportImposedData(x,t,w);

};


void TestTransportInitData2d(double x[3],double w[]){

  double t=0;
  TestTransportImposedData2d(x,t,w);

};

void TestTransportImposedData(double x[3],double t,double w[]){

  double vx =
    transport_v[0] * x[0] +
    transport_v[1] * x[1] +
    transport_v[2] * x[2];

  double xx = vx - t;

  w[0]=xx*xx;
  //w[0]=xx;
};

void TestTransportImposedData2d(double x[3],double t,double w[]){

  double vx =
    transport_v2d[0] * x[0] +
    transport_v2d[1] * x[1] +
    transport_v2d[2] * x[2];

  double xx = vx - t;

  w[0]=xx*xx;
};
// }}}

// MHD {{{

// conservative {{{
void conservatives(double* y, double* w){
    double gam = _GAM;

    w[0] = y[0];
    w[1] = y[0]*y[1];
    w[2] = y[2]/(gam-1) +  y[0]*(y[1]*y[1]+y[3]*y[3]+y[4]*y[4])/2
        + (y[7]*y[7]+y[5]*y[5]+y[6]*y[6])/2;
    w[3] = y[0]*y[3];
    w[4] = y[0]*y[4];
    w[5] = y[5];
    w[6] = y[6];
    w[7] = y[7];        // Bx
    w[8] = y[8];        // psi
}
// }}}

// fluxnum {{{
void fluxnum(double* W,double* vn, double* flux){

    double gam = _GAM;

    double un = W[1]/W[0]*vn[0]+W[3]/W[0]*vn[1]+W[4]/W[0]*vn[2];
    double bn = W[7]*vn[0]+W[5]*vn[1]+W[6]*vn[2];

    double p = (gam-1)*(W[2] - W[0]*(W[1]/W[0]*W[1]/W[0] + W[3]/W[0]*W[3]/W[0]\
                + W[4]/W[0]*W[4]/W[0])/2 - (W[7]*W[7]+W[5]*W[5]+W[6]*W[6])/2);

    flux[0] = W[0]*un;
    flux[1] = W[0]*un*W[1]/W[0] + (p + (W[7]*W[7] + W[5]*W[5] + W[6]*W[6])/2)\
              *vn[0] - bn*W[7];
    flux[2] = (W[2] + p + (W[7]*W[7] + W[5]*W[5] + W[6]*W[6])/2)*un\
              - (W[7]*W[1]/W[0] + W[5]*W[3]/W[0] + W[6]*W[4]/W[0])*bn;
    flux[3] = W[0]*un*W[3]/W[0] + (p + (W[7]*W[7] + W[5]*W[5]\
                + W[6]*W[6])/2)*vn[1] - bn*W[5];
    flux[4] = W[0]*un*W[4]/W[0] + (p + (W[7]*W[7] + W[5]*W[5]\
                + W[6]*W[6])/2)*vn[2] - bn*W[6];

    flux[5] = -bn*W[3]/W[0] + un*W[5] + W[8]*vn[1];
    flux[6] = -bn*W[4]/W[0] + un*W[6] + W[8]*vn[2];
    flux[7] = -bn*W[1]/W[0] + un*W[7] + W[8]*vn[0];

    flux[8] = _CH*_CH*bn;
}
// }}}

// MHDNumFlux {{{
void MHDNumFlux(double wL[],double wR[],double* vnorm,double* flux){
    double fluxL[9];
    double fluxR[9];
    fluxnum(wL,vnorm,fluxL);
    fluxnum(wR,vnorm,fluxR);

    for(int i=0; i<9; i++){
        flux[i] = (fluxL[i]+fluxR[i])/2 - _CH*(wR[i]-wL[i])/2;
    }
};
// }}}

// MHDBoundaryFlux {{{
void MHDBoundaryFlux(double x[3],double t,double wL[],double* vnorm,
			   double* flux){
  double wR[9];

  if(vnorm[0]!=0){
  MHDImposedData(x,t,wR);
  MHDNumFlux(wL,wR,vnorm,flux);
  }
  else if(vnorm[1]!=0){
      for(int i=0; i<9; i++){
          wR[i] = wL[i];
      }
  }
  else{
      printf("Error in MHDBoundaryFlux !\n");
      printf("vnorm = %f %f %f\n", vnorm[0], vnorm[1], vnorm[2]);
      assert(1==2);
  }
};
// }}}

// MHDInitData {{{
void MHDInitData(double x[3],double w[]){

  double t=0;
  MHDImposedData(x,t,w);

};
// }}}

// MHDImposedData {{{
void MHDImposedData(double x[3],double t,double w[]){
    double yL[9], yR[9];
    double wL[9], wR[9];

    yL[0] = 3.;
    yL[1] = 1.3;
    yL[3] = 0.;
    yL[4] = 0.;
    yL[2] = 3.;
    yL[5] = 1.;
    yL[6] = 1.;
    yL[7] = 1.5;
    yL[8] = 0.;

    yR[0] = 1.;
    yR[1] = 1.3;
    yR[3] = 0.;
    yR[4] = 0.;
    yR[2] = 1.;
    yR[5] = 0.0707372016677029;
    yR[6] = 0.9974949866040544;
    yR[7] = 1.5;
    yR[8] = 0.;

    conservatives(yL, wL);
    conservatives(yR, wR);

    if(x[0] < 0)
        for(int i=0; i<9; i++){
            w[i] = wL[i];
        }
    else
        for(int i=0; i<9; i++){
            w[i] = wR[i];
        }
};
// }}}

// }}}

