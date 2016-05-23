#include "maxwell.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>

#pragma start_opencl
void Maxwell2DCleanNumFlux_centered(real *wL, real *wR, real *vnorm, real *flux)
{
  // w: (Ex, Ey, Hz, \lambda, rho, Jx, Jy)

  // source for E: -n \times (0, 0, Hz)
  // source for H: \hat{z} \cdot n \timex (Ex, Ey, 0)
  // source for lambda: n \cdot (Ex, Ey)
  // FIXME: check
  
  const real nx = vnorm[0];
  const real ny = vnorm[1];
  const real khi = 1.0;

  // mean of L and R
  const real MEx = 0.5 * ( wL[0] + wR[0] );
  const real MEy = 0.5 * ( wL[1] + wR[1] );
  const real MHz = 0.5 * ( wL[2] + wR[2] );
  const real Mlambda =  0.5 * ( wL[3] + wR[3] );

  // (Ex Ey) flux = (-ny , nx) * Mhz + khi * (nx , ny) * Mlambda 
  flux[0] = -ny * MHz + khi * nx * Mlambda;
  flux[1] =  nx * MHz + khi * ny * Mlambda;

  // Hz flux = - (Ex , Ey) x (nx , ny)
  flux[2] = -ny * MEx + nx * MEy;

  // Lambda flux
  flux[3] = khi * (nx * MEx + ny * MEy);

  // rho flux
  flux[4] = 0.0;

  // Jx, Jy flux
  flux[5] = 0.0;
  flux[6] = 0.0;
}
#pragma end_opencl

#pragma start_opencl
void Maxwell2DCleanNumFlux_upwind(real *wL, real *wR, real *vnorm, real *flux) 
{
  // w: (Ex, Ey, Hz, lambda, rho, Jx, Jy)

  // source for E: -n \times (0, 0, Hz)
  // source for H: \hat{z} \cdot n \timex (Ex, Ey, 0)
  // source for lambda: n \cdot (Ex, Ey)
  // FIXME: check
  
  // FIXME add documentation

  const real nx = vnorm[0];
  const real ny = vnorm[1];
  const real r = sqrt(nx * nx + ny * ny);
  const real overr = 1.0 / (r + 1e-16);
  const real khi = 1.0;

  // Mean of L and R
  const real MEx = 0.5 * ( wR[0] + wL[0] );
  const real MEy = 0.5 * ( wR[1] + wL[1] );
  const real MHz = 0.5 * ( wR[2] + wL[2] );
  const real Mlambda =  0.5 * ( wR[3] + wL[3] );

  // Half or R - L
  const real DEx = 0.5 * ( wR[0] - wL[0] );
  const real DEy = 0.5 * ( wR[1] - wL[1] );
  const real DHz = 0.5 * ( wR[2] - wL[2] );
  const real Dlambda =  0.5 * ( wR[3] - wL[3] );

  // (Ex , Ey) flux: (-ny, nx) Mhz + (nx , ny) Mlambda + FIXME
  flux[0] = 
    -ny * MHz + khi * nx * Mlambda
    -overr * ( DEx * (ny * ny + khi * nx * nx)
	       + DEy * nx * ny * (khi - 1) );
  flux[1] = 
    nx * MHz + khi * ny * Mlambda
    -overr * ((nx * ny * (khi - 1)) * DEx
	      +(nx * nx + khi * ny * ny) * DEy);

  // Hz fluxL 
  flux[2] = -ny * MEx + nx * MEy - r * DHz;

  // lambda flux: 
  flux[3] = khi * (nx * MEx + ny * MEy -r * Dlambda);

  // rho, Jx, Jy, flux:
  flux[4] = 0.0;
  flux[5] = 0.0;
  flux[6] = 0.0;
}
#pragma end_opencl

#pragma start_opencl
void Maxwell2DCleanNumFlux_unoptimised(real *wL, real *wR, real *vnorm,
				       real *flux) 
{
  // w: (Ex, Ey, Hz, \lambda, rho, Jx, Jy)

  // FIXME add documentation

  const real nx = vnorm[0];
  const real ny = vnorm[1];
  const real r = sqrt(nx * nx + ny * ny);
  const real overr = 1.0 / (r + 1e-16);
  // Centered flux if eps=0, upwind flux if eps=1
  const real eps = 1;
  const real khi = 1.0;

  flux[0] = 
    - ny * (wR[2] + wL[2]) + khi * nx * (wR[3] + wL[3])
    - eps * (ny * ny + khi * nx * nx) * overr * (wR[0] - wL[0])
    - eps * (ny * nx * (khi - 1)) * overr * (wR[1] - wL[1]);

  flux[1] =   
    nx * (wR[2] + wL[2]) + khi * ny * (wR[3] + wL[3])
    - eps * (ny * nx * (khi - 1)) * overr * (wR[0] - wL[0])
    - eps * (nx * nx + khi * ny * ny) * overr  * (wR[1] - wL[1]);

  flux[2] = - ny * (wR[0] + wL[0]) 
    + nx * (wR[1] + wL[1]) 
    - eps * r * (wR[2] - wL[2]);

  flux[3] = 
    khi * nx * (wR[0] + wL[0]) 
    + khi * ny * (wR[1] + wL[1]) 
    - eps * khi * r * (wR[3] - wL[3]);

  flux[0] *= 0.5;
  flux[1] *= 0.5;
  flux[2] *= 0.5;
  flux[3] *= 0.5;
  flux[4] = 0.0;
  flux[5] = 0.0;
  flux[6] = 0.0;
}
#pragma end_opencl

#pragma start_opencl
void Maxwell2DCleanImposedData(const real *x, const real t, real *w) 
{
  // w: (Ex, Ey, Hz, lambda, rho, Jx, Jy)

  // FIXME add documentation
  
  const real pi = 4.0 * atan(1.0);
  const real r = 1.0;
  const real theta = pi / 4.0;
  const real u = cos(theta);
  const real v = sin(theta); 
  const real k = 2.0 * pi / v;
  const real c = -cos(k * (u * x[0] + v * x[1] - t));
  
  w[0] = -v * c / r;
  w[1] = u * c / r;
  w[2] = c / r;
  w[3] = 0.0;
  w[4] = 0.0;
  w[5] = 0.0;
  w[6] = 0.0;
}
#pragma end_opencl

#pragma start_opencl
void Maxwell2DCleanBoundaryFlux_upwind(real *x, real t, 
				       real *wL, real *vnorm, real *flux)
{
  real wR[7];
  Maxwell2DCleanImposedData(x, t, wR);
  Maxwell2DCleanNumFlux_upwind(wL, wR, vnorm, flux);
}
#pragma end_opencl

void Maxwell2DCleanInitData(real *x, real *w) 
{
  real t = 0.0;
  Maxwell2DCleanImposedData(x, t, w);
}

#pragma start_opencl
void Maxwell2DCleanSource(const real *x, const real t, const real *w,
			  real *source, int m)
{
  // w: (Ex, Ey, Hz, \lambda, rho, Jx, Jy)
  
  // FIXME add documentation

  // FIXMEFIXMEFIXME: this looks totally wrong.
  
  const real khi = 1.0; // FIXME: what is khi ????
  source[0] = -w[4]; // Ex -= dt * lambda ????
  source[1] = -w[5]; // Ey -= dt * rho ????
  source[2] = 0.0;
  source[3] = khi * w[6]; // lambda += dt * khi * Jy ????
  source[4] = 0.0;
  source[5] = 0.0;
  source[6] = 0.0;
}
#pragma end_opencl

#pragma start_opencl
void Maxwell3DNumFlux_upwind(real *wL, real *wR, real *vnorm, real *flux) 
{
  // Upwind flux (upwind) for Maxwell's equations

  // Data layout: w = {Ex, Ey, Ez, Hx, Hy, Hz}

  // Let {{E}} = ( ER + EL ) / 2, [[E]] = ( ER - EL ) / 2 
  // The first three components of the flux are
  // - n x {{H}} + n x n x [[E]] / r 
  // and the last three are
  //   n x {{E}} + n x n x [[H]] / r 

  const real nx = vnorm[0];
  const real ny = vnorm[1];
  const real nz = vnorm[2];

  const real overr = 1.0 / ( sqrt( nx * nx + ny * ny + nz * nz ) + 1e-16 );
  const real nxy = overr * nx * ny;
  const real nxz = overr * nx * nz;
  const real nyz = overr * ny * nz;
  const real nxx = overr * nx * nx;
  const real nyy = overr * ny * ny;
  const real nzz = overr * nz * nz;
    
  const real Esx = 0.5 * ( wR[0] + wL[0] );
  const real Esy = 0.5 * ( wR[1] + wL[1] );
  const real Esz = 0.5 * ( wR[2] + wL[2] );

  const real Hsx = 0.5 * ( wR[3] + wL[3] );
  const real Hsy = 0.5 * ( wR[4] + wL[4] );
  const real Hsz = 0.5 * ( wR[5] + wL[5] );

  const real Edx = 0.5 * ( wR[0] - wL[0] );
  const real Edy = 0.5 * ( wR[1] - wL[1] );
  const real Edz = 0.5 * ( wR[2] - wL[2] );

  const real Hdx = 0.5 * ( wR[3] - wL[3] );
  const real Hdy = 0.5 * ( wR[4] - wL[4] );
  const real Hdz = 0.5 * ( wR[5] - wL[5] );

  // E flux
  flux[0] = nz * Hsy -ny * Hsz 
    -(nyy + nzz) * Edx + nxy * Edy + nxz * Edz;
  flux[1] = -nz * Hsx + nx * Hsz
    + nxy * Edx -(nxx + nzz) * Edy + nyz * Edz; 
  flux[2] = -nx * Hsy + ny * Hsx
    + nxz * Edx + nyz * Edy -(nxx + nyy) * Edz;

  // H flux
  flux[3] = -nz * Esy + ny * Esz
    -(nyy + nzz) * Hdx + nxy * Hdy + nxz * Hdz; 
  flux[4] = nz * Esx - nx * Esz
    + nxy * Hdx -(nxx + nzz) * Hdy + nyz * Hdz; 
  flux[5] = - ny * Esx + nx * Esy  
    + nxz * Hdx + nyz * Hdy -(nxx + nyy) * Hdz;
}
#pragma end_opencl

#pragma start_opencl
void Maxwell3DCleanNumFlux_upwind(real *wL, real *wR, real *vnorm, real *flux) 
{
  // Upwind flux for Maxwell's equations

  // Data layout: w = {Ex, Ey, Ez, Hx, Hy, Hz, lambda_E, lambda_H}

  // Let {{E}} = ( ER + EL ) / 2, [[E]] = ( ER - EL ) / 2 
  // The E flux is:
  // - n x {{H}} + n x n x [[E]] / r 
  //   + c1 n {{lambda_E}} + c1 n ( n . [[E]] ) / r 
  // The H flux is:
  //   n x {{E}} + n x n x [[H]] / r 
  //   + c2 n {{lambda_H}} + c2 n ( n . [[H]] ) / r 
  // The lambda_E flux is
  // c1 * ( n . {{E}} + r [[lambda_e]] )
  // The lambda_H flux is
  // c2 * ( n . {{H}} + r [[lambda_H]] )

  // We first compute the uncleaned flux, and then add the cleaning
  // (which is cleaner, though perhaps uses a few extra operations).
  Maxwell3DNumFlux_upwind(wL, wR, vnorm, flux);

  // FIXME: how do we set these?  What are good values?
  const real c1 = 0.1; // E-cleaning parameter
  const real c2 = 0.1; // H-cleaning parameter

  // Consts based on vnorm
  const real nx = vnorm[0];
  const real ny = vnorm[1];
  const real nz = vnorm[2];
  const real r = sqrt( nx * nx + ny * ny + nz * nz );
  
  // Consts based on the mean
  const real Esx = 0.5 * ( wR[0] + wL[0] );
  const real Esy = 0.5 * ( wR[1] + wL[1] );
  const real Esz = 0.5 * ( wR[2] + wL[2] );

  const real Hsx = 0.5 * ( wR[3] + wL[3] );
  const real Hsy = 0.5 * ( wR[4] + wL[4] );
  const real Hsz = 0.5 * ( wR[5] + wL[5] );

  const real lEs = 0.5 * ( wR[6] + wL[6] );
  const real lHs = 0.5 * ( wR[7] + wL[7] );

  // Consts based on the jump
  const real Edx = 0.5 * ( wR[0] - wL[0] );
  const real Edy = 0.5 * ( wR[1] - wL[1] );
  const real Edz = 0.5 * ( wR[2] - wL[2] );

  const real Hdx = 0.5 * ( wR[3] - wL[3] );
  const real Hdy = 0.5 * ( wR[4] - wL[4] );
  const real Hdz = 0.5 * ( wR[5] - wL[5] );

  const real lEd = 0.5 * ( wR[6] - wL[6] );
  const real lHd = 0.5 * ( wR[7] - wL[7] );

  // Add correction term to E flux
  // c_1 * ( n \cdot \jump{E} + \mean{\lambda_E} )
  const real Ec = c1 * (nx * Edx + ny * Edy + nz * Edz + lEs );
  flux[0] += nx * Ec;
  flux[1] += ny * Ec;
  flux[2] += nz * Ec;

  // Add correction term to H flux
  // c_2 * ( n \cdot \jump{E} + \mean{\lambda_E} )
  const real Hc = c2 * (nx * Hdx + ny * Hdy + nz * Hdz + lHs);
  flux[3] += nx * Hc;
  flux[4] += ny * Hc;
  flux[5] += nz * Hc;

  // Flux for lambda_E
  // c_1 * ( n \cdot \mean{E} + r \jump{\lambda_E} )
  flux[6] = c1 * (nx * Esx + ny * Esy + nz * Esz + r * lEd );

  // Flux for lambda_H
  // c_2 * ( n \cdot \mean{H} + r \jump{\lambda_H} )
  flux[7] = c2 * (nx * Hsx + ny * Hsy + nz * Hsz + r * lHd );
}
#pragma end_opencl

#pragma start_opencl
void Maxwell3DCleanImposedData(const real *x, const real t, real *w) 
{
  // Data layout: w = {Ex, Ey, Ez, Hx, Hy, Hz, lambda_E, lambda_H}

#if 0
  
#ifndef M_PI
#define M_PI (3.14159)
#endif

#define MAXWELL_THETA (M_PI / 4)
#define MAXWELL_PHI   (M_PI / 4)

 real Vn[] = {
    sin(MAXWELL_THETA) * cos(MAXWELL_PHI),
    sin(MAXWELL_THETA) * sin(MAXWELL_PHI),
    cos(MAXWELL_THETA)
  };
  real Ws[] = {
    cos(MAXWELL_THETA) * cos(MAXWELL_PHI),
    cos(MAXWELL_THETA) * sin(MAXWELL_PHI),
    -sin(MAXWELL_THETA),
    -sin(MAXWELL_PHI),
    cos(MAXWELL_PHI),
    0
  };
  real test_cos_frequency = 1;

  real xdotVn = Vn[0] * x[0] + Vn[1] * x[1] + Vn[2] * x[2];
  real magnitude = cos(M_PI * 2.0 * test_cos_frequency * (xdotVn - t));

  for(int ii = 0; ii < 6; ii++) {
    w[ii] = Ws[ii] * magnitude;
  }
#else

  const real pi = 4.0 * atan(1.0);
  const real theta = pi / 4.0;
  const real r = 1.0;

  const real u = cos(theta);
  const real v = sin(theta); 
  const real k = 2.0 * pi / v;
  const real c = -cos(k * (u * x[0] + v * x[1] - t));

  // set E
  w[0] = -v * c / r;
  w[1] = u * c / r;
  w[2] = 0.0;
  // set H
  w[3] = 0.0;
  w[4] = 0.0;
  w[5] = c / r;

#endif
  // set cleaners
  w[6] = 0.0;
  w[7] = 0.0;
}
#pragma end_opencl

#pragma start_opencl
void Maxwell3DCleanInitData(real *x, real *w) 
{
  real t = 0.0;
  Maxwell3DCleanImposedData(x, t, w);
}
#pragma end_opencl

#pragma start_opencl
void Maxwell3DCleanBoundaryFlux_upwind(real *x, real t, 
				       real *wL, real *vnorm, real *flux)
{
  real wR[8];
  Maxwell3DCleanImposedData(x, t, wR);
  Maxwell3DCleanNumFlux_upwind(wL, wR, vnorm, flux);
}
#pragma end_opencl

// TODO: add 3D clean source.
