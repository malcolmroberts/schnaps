#ifndef _COLLISION_H
#define _COLLISION_H


#include "model.h"

// collision models

//! \brief particular flux for the collision model
//! \param[in] wL,wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void CollisionNumFlux(double wL[],double wR[],double vn[3],double* flux);

//! \brief particular boundary flux for the collision model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void CollisionBoundaryFlux(double* x,double t,double* wL,double* vn,
			   double* flux);

//! \brief particular init data for the collision model
//! \param[in] x : space position
//! \param[out] w : init state at point x

void CollisionInitData(double* x,double* w);
//! \brief particular init data for the  collision model
//! \param[in] x : space position
//! \param[out] w : init state at point x

void CollisionImposedData(double* x,double t,double* w);
//! \brief particular imposed data for the  collision model
//! \param[in] x,t : space and time position
//! \param[out] w : imposed state at point x and time t


#endif