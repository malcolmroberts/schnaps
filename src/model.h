#ifndef _MODEL_H
#define _MODEL_H

#include "global.h"

// typedefs for function pointers
// Numerical flux function pointer
typedef void (*fluxptr)(real*, real*, real*, real*);
// Boundary flux function pointer
typedef void (*bfluxptr)(real*, real, real*, real*, real*);
// Init data function pointer
typedef void (*initdataptr)(real*, real*);
// Imposed data function pointer
typedef void (*imposeddataptr)(const real*, const real, real* w); 

// Return a pointer to a numflux function based on the string given.
fluxptr numflux(const char *numfluxname);

// Return a pointer to a boundary blux function based on the string given.
bfluxptr bflux(const char *bfluxname);

// Return a pointer to an init data function based on the string given.
initdataptr initdata(const char *name);

// Return a pointer to an imposed data function based on the string given.
imposeddataptr imposeddata(const char *name);

//! \brief a unified framework for all physical models
typedef struct Model {
  //! Model name
  char name[50];

  //! Number of conservative variables
  int m;

  //! CFL coefficient for time-stepping
  real cfl; // 0.05

  //! Number of conservative variables in each dimension (NB: their
  //! product must equal m).
  int vlasov_mx, vlasov_my, vlasov_mz;

  //! The conservative variables have velocity in [-vlasov_vmax, vlasov_vmax]
  real vlasov_vmax;

  //! \brief A pointer to the numflux function
  //! \param[in] wL, wR : left and right states
  //! \param[in] vn : normal vector
  //! \param[out] flux : the flux
  void (*NumFlux)(real *wL, real *wR, real *vn, real *flux);

  //! \brief A pointer to the boundary flux function
  //! \param[in] x : space position
  //! \param[in] t : time
  //! \param[in] wL : left state
  //! \param[in] vn : normal vector
  //! \param[out] flux : the flux
  void (*BoundaryFlux)(real *x, real t, real *wL, real *vn, real *flux);

   //! \brief a pointer to the source function
  //! \param[in] x : space position
  //! \param[in] t : time
  //! \param[in] w :  state
  //! \param[out] source : the source
  void (*Source)(const real *x, const real t, const real *w, real *source,
		 int m);

  //! \brief A pointer to the init data function
  // !\param[in] x : space position
  //! \param[out] w : init state at point x
  void (*InitData)(real *x, real *w);

  //! \brief A pointer to the imposed data function
  //!\param[in] x, t : space and time position
  //! \param[out] w : imposed state at point x and time t
  void (*ImposedData)(const real *x, const real t, real *w);

} Model;

//! \brief The particular flux for the transport model
//! \param[in] wL, wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void TransNumFlux(real *wL, real *wR, real *vn, real *flux);

//! \brief The particular flux for the 2d transport model
//! \param[in] wL, wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
#pragma start_opencl
void TransNumFlux2d(real *wL, real *wR, real *vn, real *flux);
#pragma end_opencl

#pragma start_opencl
void VecTransNumFlux2d(__private real *wL, real *wR, real *vn, real *flux);
#pragma end_opencl

//! \brief The particular boundary flux for the transport model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void TransBoundaryFlux(real* x, real t, real* wL, real* vn, real* flux);

//! \brief The particular boundary flux for the 2d transport model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
#pragma start_opencl
void TransBoundaryFlux2d(real* x, real t, real* wL, real* vn, real* flux);
#pragma end_opencl

//! \brief The particular init data for the transport model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void TransInitData(real* x, real* w);

//! \brief The particular init data for the 2d transport model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void TransInitData2d(real *x, real *w);

void VecTransInitData2d(real *x, real *w);

//void vTransImposedData2d(real *x, real t, real *w);

//! \brief The particular imposed data for the transport model
//! \param[in] x, t : space and time position
//! \param[out] w : imposed state at point x and time t
void TransImposedData(const real* x, const real t, real* w);

//! \brief The particular imposed data for the 2d transport model
//! \param[in] x, t : space and time position
//! \param[out] w : imposed state at point x and time t
#pragma start_opencl
void TransImposedData2d(const real *x, const real t, real *w);
#pragma end_opencl

#pragma start_opencl
void VecTransImposedData2d(const real* x, const real t, real *w);
#pragma end_opencl

//! \brief The particular flux for testing the transport model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void TestTransBoundaryFlux(real* x, real t, real* wL, real* vn,
			   real* flux);

//! \brief The particular flux for testing the 2d transport model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void TestTransBoundaryFlux2d(real* x, real t, real* wL, real* vn,
			     real* flux);

#pragma start_opencl
void VecTransBoundaryFlux2d(real* x, real t, real* wL, real* vn,
			    real* flux);
#pragma end_opencl

//! \brief The particular init data for the transport model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void TestTransInitData(real* x, real* w);

//! \brief The particular init data for the 2d transport model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void TestTransInitData2d(real* x, real* w);

//! \brief The particular imposed data for the transport model
//! \param[in] x, t : space and time position
//! \param[out] w : imposed state at point x and time t
void TestTransImposedData(const real* x, const real t, real* w);

//! \brief The particular imposed data for the 2d transport model
//! \param[in] x, t : space and time position
//! \param[out] w : imposed state at point x and time t
void TestTransImposedData2d(const real* x, const real t, real* w);

//! Parameters used for computing the velocity for the Vlasov equation.
// FIXME: pass all of this instead of making it global.
static int m;
static int vlasov_mx;
static int vlasov_my;
static int vlasov_mz;
static real vlasov_vmax;

void set_global_m(int m0);
void set_vlasov_params(Model *mod);
real vlasov_vel(const int id, const int md, real vlasov_vmax);
void vlaTransInitData2d(real *x, real *w);
void vlaTransNumFlux2d(real *wL, real *wR, real* vnorm, real* flux);
void vlaTransBoundaryFlux2d(real *x, real t, 
			    real *wL, real* vnorm,
			    real* flux);
void vlaTransImposedData2d(const real *x, const real t, real* w);
real compact_bump(real r);
real icgaussian(real r, real c);

void cemracs2014_TransBoundaryFlux(real *x, real t, 
				   real *wL, real *vnorm,
				   real *flux);
void cemracs2014_TransInitData(real *x, real *w);
void cemcracs2014_imposed_data(const real *x, const real t, real *w); 
void cemracs2014a_TransInitData(real *x, real *w);
void cemcracs2014a_imposed_data(const real *x, const real t, real *w); 

void OneSource(const real *x, const real t, const real *w, real *source, int m);

#endif
