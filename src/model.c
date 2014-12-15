#include "model.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>

#define ONE_OVER_SQRT_3 (0.57735026918962584)
#define ONE_OVER_SQRT_2 (0.707106781186547524400844362105)

// {{{ MHD
#define _CH (5)
#define _GAM (1.666666666666)
//#define PADE
//#define PADE2
#define P2
//#define P6
//#define DIFFPADE
//#define DIFFP2
//#define DIFFP6

const double transport_v[] = {ONE_OVER_SQRT_3,
			      ONE_OVER_SQRT_3,
			      ONE_OVER_SQRT_3};

//const double transport_v[] = {1, 0, 0};

const double transport_v2d[] = {ONE_OVER_SQRT_2,
				ONE_OVER_SQRT_2,
				0};

// {{{ Transport

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

// {{{ TestTransport
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

// {{{ MHD

// {{{   conservative
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
// {{{   primitives
void primitives(double* W, double* Y){

  double gam = _GAM;

  Y[0] = W[0];
  Y[1] = W[1]/W[0];
  Y[2] = (gam-1)*(W[2] - W[0]*(W[1]/W[0]*W[1]/W[0] + W[3]/W[0]*W[3]/W[0]
			       + W[4]/W[0]*W[4]/W[0])/2 - (W[7]*W[7]+W[5]*W[5]+W[6]*W[6])/2);
  Y[3] = W[3]/W[0];
  Y[4] = W[4]/W[0];
  Y[5] = W[5];
  Y[6] = W[6];
  Y[7] = W[7];        // Bx
  Y[8] = W[8];        // psi
}
// }}}
// {{{   jacobmhd

void jacobmhd(double* W,double* vn, double M[9][9]){

  double gam = _GAM;
  double Y[9];

  double rho, ux, uy, uz, by, bz, p, bx;

  int i,j;

  for(i=0; i<9; i++){
    for(j=0; j<9; j++){
      M[i][j] = 0;
    }
  }

  primitives(W,Y);

  rho = Y[0];
  ux = Y[1];
  p = Y[2];
  uy = Y[3];
  uz = Y[4];
  by = Y[5];
  bz = Y[6];
  bx = Y[7];

  M[0][0] = ux*vn[0]+uy*vn[1]+uz*vn[2];
  M[0][1] = rho*vn[0];
  M[0][3] = rho*vn[1];
  M[0][4] = rho*vn[2];

  M[1][1] = ux*vn[0]+uy*vn[1]+uz*vn[2];
  M[1][2] = 1/rho*vn[0];
  M[1][5] = -(vn[1]*bx-by*vn[0])/rho;
  M[1][6] = -(vn[2]*bx-bz*vn[0])/rho;
  M[1][7] = -(bx*vn[0]+by*vn[1]+bz*vn[2])/rho;

  M[2][1] = gam*p*vn[0];
  M[2][2] = ux*vn[0]+uy*vn[1]+uz*vn[2];
  M[2][3] = gam*p*vn[1];
  M[2][4] = gam*p*vn[2];
  M[2][5] = -ux*vn[1]*bx-uz*bz*vn[1]-uy*by*vn[1]+uy*by*vn[1]*gam+\
    ux*gam*vn[1]*bx+uz*gam*bz*vn[1];
  M[2][6] = -uz*vn[2]*bz+ux*gam*bx*vn[2]+uy*gam*by*vn[2]+uz*bz*vn[2]*gam-\
    ux*vn[2]*bx-uy*by*vn[2];
  M[2][7] = -ux*vn[0]*bx-uy*by*vn[0]-uz*bz*vn[0]+ux*bx*vn[0]*gam+\
    uy*gam*vn[0]*by+uz*gam*bz*vn[0];
  M[2][8] = -bx*vn[0]*gam+bx*vn[0]-by*vn[1]*gam+by*vn[1]\
    -bz*vn[2]*gam+bz*vn[2];

  M[3][2] = 1/rho*vn[1];
  M[3][3] = ux*vn[0]+uy*vn[1]+uz*vn[2];
  M[3][5] = -(bx*vn[0]+by*vn[1]+bz*vn[2])/rho;
  M[3][6] = -(vn[2]*by-bz*vn[1])/rho;
  M[3][7] = (vn[1]*bx-by*vn[0])/rho;

  M[4][2] = 1/rho*vn[2];
  M[4][4] = ux*vn[0]+uy*vn[1]+uz*vn[2];
  M[4][5] = (vn[2]*by-bz*vn[1])/rho;
  M[4][6] = -(bx*vn[0]+by*vn[1]+bz*vn[2])/rho;
  M[4][7] = (vn[2]*bx-bz*vn[0])/rho;

  M[5][1] = by*vn[0];
  M[5][3] = -bx*vn[0]-bz*vn[2];
  M[5][4] = vn[2]*by;
  M[5][5] = ux*vn[0]+uz*vn[2];
  M[5][6] = -vn[2]*uy;
  M[5][7] = -vn[0]*uy;
  M[5][8] = vn[1];

  M[6][1] = bz*vn[0];
  M[6][3] = bz*vn[1];
  M[6][4] = -vn[0]*bx-by*vn[1];
  M[6][5] = -vn[1]*uz;
  M[6][6] = ux*vn[0]+uy*vn[1];
  M[6][7] = -vn[0]*uz;
  M[6][8] = vn[2];

  M[7][1] = -vn[1]*by-bz*vn[2];
  M[7][3] = bx*vn[1];
  M[7][4] = bx*vn[2];
  M[7][5] = -vn[1]*ux;
  M[7][6] = -vn[2]*ux;
  M[7][7] = uy*vn[1]+uz*vn[2];
  M[7][8] = vn[0];

  M[8][5] = _CH*_CH*vn[1];
  M[8][6] = _CH*_CH*vn[2];
  M[8][7] = _CH*_CH*vn[0];

  for(int i=0; i<9; i++){
    for(int j=0; j<9; j++){
      M[i][j] /= _CH;
    }
  }
}

// }}}
// {{{   matmul
void matrix_vector(double A[9][9], double B[9], double* C){

  for(int i=0; i<9; i++){
    C[i] = 0;
  }

  for(int i=0; i<9; i++){
    for(int j=0; j<9; j++){
      C[i] += A[i][j]*B[j];
    }
  }
}

void matrix_matrix(double A[9][9],double B[9][9],double C[9][9]){
  for(int i=0; i<9; i++){
    for(int j=0; j<9; j++){
      C[i][j]=0;
    }
  }

  for(int i=0; i<9; i++){
    for(int j=0; j<9; j++){
      for(int k=0; k<9; k++){
	C[i][j] = C[i][j] + (A[i][k]*B[k][j]);
      }
    }
  }
}
// }}}
// {{{   gauss

void write_matrix(double A[9][9],double *second, double B[9][9+1]){
  
  for (int i = 0; i < 9; i++){
    for (int j = 0; j < 9 ; j++){
      B[i][j]=A[i][j];
    }
  }

  for (int i = 0; i < 9; i++){
    B[i][9]=second[i];
  }
}

void gauss(double A[9][9], double b[9], double *x){
  double B[9][9+1];
  write_matrix(A,b,B);
  // go down
  for (int i = 0; i < 9; i++){
    for (int j = 9; j >= i; --j){
      for (int k = i + 1; k < 9; k++){
	B[k][j] = B[k][j] - ((B[k][i]/B[i][i]) * B[i][j]);
      }
    }
  }

  // go up
  for (int i = 9 - 1; i >= 0; i--) {
    B[i][9] = B[i][9] / B[i][i];
    for (int j = i - 1; j >= 0; j--) {
      B[j][9] = B[j][9] - (B[j][i] * B[i][9]);
      B[j][i] = 0;
    }
  }

  for(int i=0; i<9; i++){
    x[i] = B[i][9];
  }
}

// }}}
// {{{   fluxnum
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
// {{{   MHDNumFlux
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
// {{{   MHDNumFlux_2
void MHDNumFlux_2(double wL[],double wR[],double* vn, double* flux){

  double wmil[9];
  double wRmwL[9];
  double Z[9];

  double M[9][9] = {{0}};
  double M2[9][9] = {{0}};

  for (int i=0; i< 9 ; i++){
    wmil[i] = (wL[i] + wR[i])/2;
  }

  // calcul de la matrice M
  jacobmhd(wmil,vn,M);

  // {{{ PADE
#ifdef PADE
  // calcul de la matrice M^2
  matrix_matrix(M,M,M2);

  // calcul de la matrice (I+3M^2) on la stock dans M
  for (int i = 0; i < 9; i++) {
    for (int j = 0; j < 9; j++) {
      if (i==j) {
	M[i][j] = 1 + 3*M2[i][j];
      } else {
	M[i][j] = 3*M2[i][j];
      }
    }
  }

  // calcul de (wR-wL)
  for (int i=0; i< 9 ; i++){
    wRmwL[i] = (wR[i] - wL[i]);
  }

  // calcul de (I+3M^2)*(wRmwL)
  matrix_vector(M,wRmwL,Z);

  // calcul de la matrice (3I+M^2) on le stock encore dans M
  for (int i = 0; i < 9; i++) {
    for (int j = 0; j < 9; j++) {
      if (i==j) {
	M[i][j] = 3 + M2[i][j];
      }
      else {
	M[i][j] = M2[i][j];
      }
    }
  }

  // resolution du systeme (3I + M^2)x = (I + 3M^2)(wR - wL)
  //                   <=> x = (3I + M^2)^-1(I + 3M^2)(wR - wL)
  //                   <=> x = |M|(wR - wL)
  double abs[9];
  gauss(M,Z,abs);
#endif
  // }}}
  // {{{ PADE2
#ifdef PADE2
  // calcul de la matrice M^2
  matrix_matrix(M,M,M2);

  // calcul de la matrice (I+12M^2) on la stock dans M
  for (int i = 0; i < 9; i++) {
    for (int j = 0; j < 9; j++) {
      if (i==j) {
	M[i][j] = 1 + 12*M2[i][j];
      }
      else {
	M[i][j] = 12*M2[i][j];
      }
    }
  }

  // calcul de (wR-wL)
  for (int i=0; i< 9 ; i++){
    wRmwL[i] = (wR[i] - wL[i]);
  }

  // calcul de (I+12M^2)*(wRmwL)
  matrix_vector(M,wRmwL,Z);

  // calcul de la matrice (6I+7M^2) on le stock encore dans M
  for (int i = 0; i < 9; i++) {
    for (int j = 0; j < 9; j++) {
      if (i==j) {
	M[i][j] = 6 + 7*M2[i][j];
      }
      else {
	M[i][j] = 7*M2[i][j];
      }
    }
  }

  // resolution du systeme (6I + 7M^2)x = (I + 12M^2)(wR - wL)
  //                   <=> x = (6I + 7M^2)^-1(I + 12M^2)(wR - wL)
  //                   <=> x = |M|(wR - wL)
  double abs[9];
  gauss(M,Z,abs);
#endif
  // }}}
  // {{{ P2
#ifdef P2
  double coef[3] = {1./2, 0., 1./2};

  // calcul de (wR-wL)
  for (int i=0; i< 9 ; i++){
    wRmwL[i] = (wR[i] - wL[i]);
  }

  double abs[9];

  for (int i=0; i< 9 ; i++){
    abs[i] = coef[2]*wRmwL[i];
  }

  double Mw[9];
  for(int i=1; i>=0; i--){
    matrix_vector(M,abs,Mw);
    for(int j=0; j<9; j++){
      abs[j] = coef[i]*wRmwL[j] + Mw[j];
    }
  }
#endif
  // }}}
  // {{{ P6
#ifdef P6
  double coef[7] = {5./16, 0., 15./16, 0., -5./16, 0., 1./16};

  // calcul de (wR-wL)
  for (int i=0; i< 9 ; i++){
    wRmwL[i] = (wR[i] - wL[i]);
  }

  double abs[9];

  for (int i=0; i< 9 ; i++){
    abs[i] = coef[6]*wRmwL[i];
  }

  double Mw[9];
  for(int i=5; i>=0; i--){
    matrix_vector(M,abs,Mw);
    for(int j=0; j<9; j++){
      abs[j] = coef[i]*wRmwL[j] + Mw[j];
    }
  }
#endif
  // }}}
  // {{{ DIFFPADE
#ifdef DIFFPADE
  // calcul de la matrice M^2
  matrix_matrix(M,M,M2);

  // calcul de la matrice (16M) on la stock dans M
  for (int i = 0; i < 9; i++) {
    for (int j = 0; j < 9; j++) {
      M[i][j] = 16*M[i][j];
    }
  }

  // calcul de (wR-wL)
  for (int i=0; i< 9 ; i++){
    wRmwL[i] = (wR[i] - wL[i]);
  }

  // calcul de (16M)*(wRmwL)
  matrix_vector(M,wRmwL,Z);

  // calcul de la matrice (3I+M^2) on le stock encore dans M
  for (int i = 0; i < 9; i++) {
    for (int j = 0; j < 9; j++) {
      if (i==j) {
	M[i][j] = 3 + M2[i][j];
      }
      else {
	M[i][j] = M2[i][j];
      }
    }
  }

  // calcul de la matrice (3I+M^2)^2 on le stock encore dans M2
  matrix_matrix(M,M,M2);

  // resolution du systeme (3I + M^2)^2x = (16M)(wR - wL)
  //                   <=> x = ((3I + M^2)^2)^-1(16M)(wR - wL)
  //                   <=> x = sgn(A)(wR - wL)
  double abs[9];
  gauss(M2,Z,abs);
#endif
  // }}}
  // {{{ DIFFP2
#ifdef DIFFP2
  // calcul de (wR-wL)
  for (int i=0; i< 9 ; i++){
    wRmwL[i] = (wR[i] - wL[i]);
  }

  double abs[9];

  matrix_vector(M,wRmwL,abs);
#endif
  // }}}
  // {{{ DIFFP6
#ifdef DIFFP6
  double coef[6] = {0., 15./8, 0., -10./8, 0., 3./8};

  // calcul de (wR-wL)
  for (int i=0; i< 9 ; i++){
    wRmwL[i] = (wR[i] - wL[i]);
  }

  double abs[9];

  for (int i=0; i< 9 ; i++){
    abs[i] = coef[5]*wRmwL[i];
  }

  double Mw[9];
  for(int i=4; i>=0; i--){
    matrix_vector(M,abs,Mw);
    for(int j=0; j<9; j++){
      abs[j] = coef[i]*wRmwL[j] + Mw[j];
    }
  }
#endif
  // }}}

  double fluxL[9];
  double fluxR[9];

  fluxnum(wL,vn,fluxL);
  fluxnum(wR,vn,fluxR);

  for(int i=0; i<9; i++){
    flux[i] = (fluxL[i]+fluxR[i])/2 - _CH*abs[i]/2;
  }
}
// }}}
// {{{   MHDBoundaryFlux
void MHDBoundaryFlux(double x[3],double t,double wL[],double* vnorm,
		     double* flux){
  double wR[9];

  if(vnorm[0] > 0.0001 || vnorm[0] < -0.0001){
    MHDImposedData(x,t,wR);
  }
  else if(vnorm[1] > 0.0001 || vnorm[1] < -0.0001){
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
// {{{   MHDInitData
void MHDInitData(double x[3],double w[]){

  double t=0;
  MHDImposedData(x,t,w);

};
// }}}
// {{{   MHDImposedData
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

// {{{ Vlasov 2D
// Vlasov 2D transport equation functions

// Set the parameters for the Vlasov equation that are stored in the
// global space of model.h
void set_vlasov_params(Model *mod) 
{
  m = mod->m;
  mx = mod->mx;
  my = mod->my;
  mz = mod->mz;
  assert(m == mx * my * mz);
  vmax = mod->vmax;
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
