#ifndef _QUANTITIES_VP_H
#define _QUANTITIES_VP_H

#include "model.h"
#include "field.h"


//! \brief compute the charge density for t and x
//! \param[in] x,t : space and time position
//! \param[inout] w : values of f at glops
//! \return in w[_MV+2] the charge density at x and t 
void Computation_charge_density(Field* f);


//! \brief compute the electric field on the mesh
//! \param[inout] : field
void Compute_electric_field(Field* f);


#endif
