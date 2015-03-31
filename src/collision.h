#ifndef _COLLISION_H
#define _COLLISION_H

#define _NB_ELEM_V 16
#define _DEG_V 3

//! \brief number of conservative variables
//! values of  the distribution function at the velocity glops
//! and value of the potential 
#define _MV (_NB_ELEM_V *  _DEG_V + 1 )

#define _INDEX_MAX_KIN (_MV-1)
#define _INDEX_PHI (_MV)
#define _INDEX_EX (_MV+1)
#define _INDEX_RHO (_MV+2)
#define _INDEX_VELOCITY (_MV+3)
#define _INDEX_PRESSURE (_MV+4)
#define _INDEX_TEMP (_MV+5)
#define _INDEX_MAX (_MV+6)

#define _VMAX 6.
#define _DV (2*_VMAX / _NB_ELEM_V)
#define _Dirichlet_Poisson_BC (1)
#define _Periodic_Poisson_BC (2)



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

//! \brief particular imposed data for the  collision model
//! \param[in] x,t : space and time position
//! \param[out] w : imposed state at point x and time t
void CollisionImposedData(double* x,double t,double* w);



//! \brief compute compute the source term of the collision
//! model: electric force + true collisions
//! \param[in] w the distribution function
//! \param[in] f the force 
//! \param[out] source the source
void CollisionSource(double* x,double t, double* w, double* source);


#endif