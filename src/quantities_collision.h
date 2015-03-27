#ifndef _QUANTITIES_COLLISION_H
#define _QUANTITIES_COLLISION_H

#include "model.h"
#include "field.h"


//! \brief particular imposed data for the  collision model
//! \param[in] x,t : space and time position
//! \param[out] w : imposed state at point x and time t
double Collision_ImposedKinetic_Data(double* x,double t,double v);


//! \brief compute l2 error between imposed and actual solution of the collsion model
//! \param[in] f : a field
//! \return the error
double L2_Kinetic_error(Field* f);

//! \brief compute square of velocity L2 error
//! \param[in] x,t : space and time position
//! \param[in] w : values of f at glops
//! \return the velocity error
double L2VelError(double* x,double t,double *w);


//! \brief compute the charge density for t and x
//! \param[in] x,t : space and time position
//! \param[inout] w : values of f at glops
//! \return in w[_MV+2] the charge density at x and t 
void Computation_charge_density(double* x,double t,double *w);



//! \brief compute the electric field on the mesh
//! \param[inout] : field
void Compute_electric_field(Field* f);


#endif
