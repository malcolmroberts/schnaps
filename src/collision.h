#ifndef _COLLISION_H
#define _COLLISION_H

#define _NB_ELEM_V 24
#define _DEG_V 2

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

#define _ENTROPY (1)
#define _EPS (1)
#define _LAMBDA (1)

#include "model.h"
#include "field.h"
// collision models

//! \brief particular flux for the collision model
//! \param[in] wL,wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void VlasovP_Lagrangian_NumFlux(schnaps_real *wL, schnaps_real *wR, schnaps_real *vn, schnaps_real *flux);

//! \brief  compute the source term of the collision
//! model: electric force + true collisions
//! \param[in] x space position
//! \param[in] t time
//! \param[in] w the distribution function
//! \param[out] source the source
void VlasovP_Lagrangian_Source(const schnaps_real *x, const schnaps_real t, const schnaps_real *w, schnaps_real *source);


//! \brief compute M^{-1 * M_f(v) * w for collision step
//! \param[in] f the field
//! \param[in] w the distribution function in entropic variable
//! \param[in] function the function which compute second order derivate of the adjoint entropic transformation
//! \param[out] product contains the result
void VlasovP_Mass_modified(field *f,schnaps_real * w,void (*function)(field *f,schnaps_real w,schnaps_real *tw),schnaps_real* product);

			    
#endif
