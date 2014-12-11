#include "model.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>

#define ONE_OVER_SQRT_3 (0.57735026918962584)
#define ONE_OVER_SQRT_2 (0.707106781186547524400844362105)

#define _CH (5)
#define _GAM (1.666666666666)

const double transport_v[] = {ONE_OVER_SQRT_3,
			      ONE_OVER_SQRT_3,
			      ONE_OVER_SQRT_3};

//const double transport_v[] = {1, 0, 0};

const double transport_v2d[] = {ONE_OVER_SQRT_2,
				ONE_OVER_SQRT_2,
				0};

// Transport {{{
void TransNumFlux(double wL[], double wR[], double* vnorm, double* flux) {
  double vn
    = transport_v[0] * vnorm[0]
    + transport_v[1] * vnorm[1]
    + transport_v[2] * vnorm[2];
  double vnp = vn > 0 ? vn : 0;
  double vnm = vn - vnp;
  flux[0] = vnp * wL[0] + vnm * wR[0];
};

void TransNumFlux2d(double wL[], double wR[], double* vnorm, double* flux) {
  double vn
    = transport_v2d[0] * vnorm[0]
    + transport_v2d[1] * vnorm[1]
    + transport_v2d[2] * vnorm[2];
  double vnp = vn > 0 ? vn : 0;
  double vnm = vn - vnp;
  flux[0] = vnp * wL[0] + vnm * wR[0];
  /* if (fabs(vnorm[2])>1e-6) { */
  /*   printf("vnds %lf %lf %lf \n", vnorm[0], vnorm[1], vnorm[2]); */
  /* } */
  // verify that 2d computations are actually
  // activated
  assert(fabs(vnorm[2]) < 1e-8);
};

void VecTransNumFlux2d(double wL[], double wR[], double* vnorm, double* flux) {
  double vn
    = transport_v2d[0] * vnorm[0]
    + transport_v2d[1] * vnorm[1]
    + transport_v2d[2] * vnorm[2];
  double vnp = vn > 0 ? vn : 0;
  double vnm = vn - vnp;
  flux[0] = vnp * wL[0] + vnm * wR[0];
  flux[1] = vnp * wL[1] + vnm * wR[1];
  /* if (fabs(vnorm[2])>1e-6) { */
  /*   printf("vnds %lf %lf %lf \n", vnorm[0], vnorm[1], vnorm[2]); */
  /* } */
  // Verify that 2d computations are actually activated
  assert(fabs(vnorm[2]) < 1e-8);
};

void TransBoundaryFlux(double x[3], double t, double wL[], double* vnorm,
		       double* flux) {
  double wR[1];
  TransImposedData(x, t, wR);
  TransNumFlux(wL, wR, vnorm, flux);
};

void TransBoundaryFlux2d(double x[3], double t, double wL[], double *vnorm,
			 double *flux) {
  double wR[1];
  TransImposedData2d(x, t, wR);
  TransNumFlux2d(wL, wR, vnorm, flux);
};

void VecTransBoundaryFlux2d(double x[3], double t,
			    double wL[], double *vnorm,
			    double* flux) {
  double wR[2];
  VecTransImposedData2d(x, t, wR);
  VecTransNumFlux2d(wL, wR, vnorm, flux);
};

void TransInitData(double x[3], double w[]) {
  double t = 0;
  TransImposedData(x, t, w);
};

void TransInitData2d(double x[3], double w[]) {
  double t = 0;
  TransImposedData2d(x, t, w);
};

void VecTransInitData2d(double x[3], double w[]) {
  double t = 0;
  VecTransImposedData2d(x, t, w);
};

void TransImposedData(double x[3], double t, double w[]) {
  double vx =
    transport_v[0] * x[0] +
    transport_v[1] * x[1] +
    transport_v[2] * x[2];
  double xx = vx - t;
  w[0]=cos(xx);
};

void TransImposedData2d(double x[3], double t, double w[]) {
  double vx
    = transport_v2d[0] * x[0]
    + transport_v2d[1] * x[1]
    + transport_v2d[2] * x[2];
  double xx = vx - t;
  w[0] = cos(xx);
};

// m = 2 test-case
void VecTransImposedData2d(double x[3], double t, double* w) {
  double vx
    = transport_v2d[0] * x[0]
    + transport_v2d[1] * x[1]
    + transport_v2d[2] * x[2];
  double xx = vx - t;
  w[0] = xx * xx;
  w[1] = 2 * xx * xx;
};
// }}}

// TestTransport {{{
void TestTransBoundaryFlux(double x[3], double t,
			       double wL[], double* vnorm,
			       double* flux) {
  double wR[1];
  TestTransImposedData(x, t, wR);
  TransNumFlux(wL, wR, vnorm, flux);
};

void TestTransBoundaryFlux2d(double x[3], double t, double wL[],
				 double* vnorm, double* flux) {
  double wR[1];
  TestTransImposedData2d(x, t, wR);
  TransNumFlux2d(wL, wR, vnorm, flux);
};

void TestTransInitData(double x[3], double w[]) {
  double t = 0;
  TestTransImposedData(x, t, w);
};

void TestTransInitData2d(double x[3], double w[]) {
  double t = 0;
  TestTransImposedData2d(x, t, w);
};

void TestTransImposedData(double x[3], double t, double w[]) {
  double vx
    = transport_v[0] * x[0]
    + transport_v[1] * x[1]
    + transport_v[2] * x[2];
  double xx = vx - t;
  w[0] = xx * xx;
};

void TestTransImposedData2d(double x[3], double t, double w[]) {
  double vx
    = transport_v2d[0] * x[0]
    + transport_v2d[1] * x[1]
    + transport_v2d[2] * x[2];
  double xx = vx - t;
  w[0] = xx * xx;
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
  MHDNumFlux(wL,wR,vnorm,flux);
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

// Vlasov 2D {{{
// Vlasov 2D transport equation functions

// Set the parameters for the Vlasov equation that are stored in the
// global space of model.h
void set_vlasov_params(int mx0, int my0, int mz0, double vmax0)
{
  mx = mx0;
  my = my0;
  mz = mz0;
  m = mx * my * mz;
  vmax = vmax0;
}

void vlaTransInitData2d(double x[3], double w[])
{
  double t = 0;
  vlaTransImposedData2d(x, t, w);
};

void vlaTransNumFlux2d(double wL[], double wR[], double* vnorm, double* flux)
{
  for(int ix = 0; ix < mx; ++ix) {
    double vx = vmax * (ix - (mx / 2));

    for(int iy = 0; iy < my; ++iy) {
      double vy = vmax * (iy - (my / 2));

      double vn = vx * vnorm[0]	+ vy * vnorm[1];
      double vnp = vn > 0 ? vn : 0;
      double vnm = vn - vnp;

      // NB: assumes a certain memory distribution for the velocity
      // components at each point.
      int im = ix * my + iy;
      flux[im] = vnp * wL[im] + vnm * wR[im];
    }
  }
  // Verify that 2d computations are actually activated
  assert(fabs(vnorm[2]) < 1e-8);
};

void vlaTransBoundaryFlux2d(double x[3], double t,
			    double wL[], double* vnorm,
			    double* flux)
{
  double wR[m];
  vlaTransImposedData2d(x, t, wR);
  vlaTransNumFlux2d(wL, wR, vnorm, flux);
};

// 6th-degree polynomial with compact support
double compact_poly6(double r)
{
  if (fabs(2 * r) > 1)
    return 0;
  double rrm1 = 2 * r - 1;
  double rrp1 = 2 * r + 1;
  return -35.0 / 16.0 * rrm1 * rrm1 * rrm1 * rrp1 * rrp1 * rrp1;
}

// Impose a compact-supporrt 6th degree polynomial in 2+2D
void vlaTransImposedData2d(double x[3], double t, double* w)
{
  double PI = 4.0 * atan(1.0);
  double s2pi = sqrt(2.0 * PI);
  double xval = 1.0;

  double r = sqrt(x[0] * x[0] + x[1] * x[1]);
  double pr = compact_poly6(r);

  for(int ix = 0; ix < mx; ++ix) {
    double vx = vmax * (ix - (mx / 2));

    for(int iy = 0; iy < my; ++iy) {
      double vy = vmax * (iy - (my / 2));

      double vr = sqrt(vx * vx + vy * vy);
      double pvr = compact_poly6(vr);

      // NB: assumes a certain memory distribution for the velocity
      // components at each point.
      int im = ix * my + iy;
      w[im] = pr * pvr;
    }
  }
};
// }}}
