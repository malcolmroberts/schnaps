#ifndef _MHD_H
#define _MHD_H
#include "model.h"

#pragma start_opencl
//! \brief computes the conservatives states from the primitives
//! \param[in] y : primitives states
//! \param[out] w : conservatives states
void conservatives(real *y, real *w);
void primitives(real *W, real *Y);


//! \brief Numerical flux for the MHD model
//! \param[in] w : states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void fluxnum(real *w, real *vn, real *flux);

//! \brief particular flux for the MHD model
//! \param[in] wL,wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void MHDNumFluxRusanov(real *wL, real *wR, real *vn, real *flux);
void MHDNumFluxP2(real *wL,real *wR,real *vn, real *flux);
void MHDNumFlux1D(real *wL,real *wR,real *vn, real *flux);


//! \brief particular boundary flux for the MHD model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void MHDBoundaryFlux(real *x, real t, real *wL, real *vn, real *flux);

//! \brief particular init data for the MHD model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void MHDInitData(real *x, real *w);

//! \brief particular imposed data for the MHD model
//! \param[in] x,t : space and time position
//! \param[out] w : imposed state at point x and time t
void MHDImposedData(real *x,real t, real *w);

// FIXME: using "real var[]" instead of "real *var" breaks OpenCL on
// certain platforms.
void jacobmhd(real* W,real* vn, real *M);
void matrix_vector(real *A, real B[9], real* C);

#pragma end_opencl


#endif
