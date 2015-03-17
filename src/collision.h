#ifndef _COLLISION_H
#define _COLLISION_H

#define _NB_ELEM_V 4
#define _DEG_V 1

#define _MV (_NB_ELEM_V *  _DEG_V + 1)
#define _VMAX 6
#define _DV (2*_VMAX / _NB_ELEM_V)


#include "model.h"
#include "field.h"
// collision models

//! \brief particular flux for the collision model
//! \param[in] wL,wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void Collision_Lagrangian_NumFlux(double wL[],double wR[],double vn[3],double* flux);

//! \brief particular boundary flux for the collision model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void Collision_Lagrangian_BoundaryFlux(double* x,double t,double* wL,double* vn,
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

double Collision_ImposedKinetic_Data(double* x,double t,double v);
//! \brief particular imposed data for the  collision model
//! \param[in] x,v,t : space, velocity and time position
//! \param[out] f : imposed state at point x, the velocity v and time t

double L2_Kinetic_error(field* f);

//! \brief compute square of velocity L2 error
//! \param[in] x,t : space and time position
//! \param[in] w : values of f at glops
double L2VelError(double* x,double t,double *w);

//! \brief compute compute the source term of the collision
//! model: electric force + true collisions
void CollisionSource(double* force, double* w, double* source);

#endif
