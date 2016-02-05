#ifndef _GYRO_H
#define _GYRO_H

// FIXME:

// The use of defines to set varialbes is going to cause massive
// problems with maintainability.  Please stop doing this.

#define _NB_ELEM_V 1
#define _DEG_V 1

#define _MV (_NB_ELEM_V *  _DEG_V + 1)
#define _INDEX_MAX_KIN (_MV-1)
#define _INDEX_PHI (_MV)
#define _INDEX_EX (_MV+1)
#define _INDEX_EY (_MV+2)
#define _INDEX_EZ (_MV+3)
#define _INDEX_MAX (_MV+4)
#define _VMAX 6.
#define _DV (2*_VMAX / _NB_ELEM_V)


#include "model.h"
#include "field.h"
// gyrokinetic models

//! \brief particular flux for the gyro model
//! \param[in] wL,wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void Gyro_Upwind_NumFlux(schnaps_real wL[],schnaps_real wR[],
			       schnaps_real vn[3],schnaps_real* flux);

//! \brief particular flux for the gyro model
//! \param[in] wL,wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void Gyro_Lagrangian_NumFlux(schnaps_real wL[],schnaps_real wR[],
			     schnaps_real vn[3],schnaps_real* flux);

//! \brief particular boundary flux for the gyro model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void Gyro_Lagrangian_BoundaryFlux(schnaps_real* x,schnaps_real t,schnaps_real* wL,schnaps_real* vn,
			   schnaps_real* flux);

//! \brief particular init data for the gyro model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void GyroInitData(schnaps_real* x,schnaps_real* w);

//! \brief particular imposed data for the  gyro model
//! \param[in] x  space position
//! \param[in] t time
//! \param[out] w  init state at point x
void GyroImposedData(const schnaps_real* x, const schnaps_real t,schnaps_real* w);

//! \brief particular imposed data for the  gyro model
//! \param[in] x  space
//! \param[in] t time
//! \param[in] v  velocity
//! \returns value of the distribution function
schnaps_real Gyro_ImposedKinetic_Data(const schnaps_real* x, const schnaps_real t,schnaps_real v);

//! \brief compute gyro L2 error in x and v
//! \param[in] f : a field
schnaps_real GyroL2_Kinetic_error(field* f);

//! \brief compute square of velocity L2 error
//! \param[in] x,t : space and time position
//! \param[in] w : values of f at glops
schnaps_real GyroL2VelError(schnaps_real* x,schnaps_real t,schnaps_real *w);

//! \brief compute compute the source term of the collision
//! model: electric force + true collisions
void GyroSource(const schnaps_real *x, const schnaps_real t, const schnaps_real *w, schnaps_real *source);

//! \brief gnuplot file for the distribution function
//! \param[in] w : values of f at glops
void Velocity_distribution_plot(schnaps_real *w);

#endif
