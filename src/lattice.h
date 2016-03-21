#ifndef _LATTICE_H
#define _LATTICE_H

//! \brief number of conservative variables
//! values of  the distribution function at the velocity glops
//! and value of the potential 

#include "model.h"
#include "field.h"
// collision models

//! \brief particular flux for the collision model
//! \param[in] wL,wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void Lattice_NumFlux(schnaps_real *wL, schnaps_real *wR, schnaps_real *vn, schnaps_real *flux);
		    
#endif
