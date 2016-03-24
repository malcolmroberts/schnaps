#ifndef _LATTICE_H
#define _LATTICE_H

//! \brief number of conservative variables
//! values of  the distribution function at the velocity glops
//! and value of the potential 

#include "model.h"
#include "field.h"
#include "simulation.h"
#include "global.h"
// collision models

//! \brief particular flux for the collision model
//! \param[in] wL,wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux

void Lattice_NumFlux(schnaps_real *wL, schnaps_real *wR, schnaps_real *vn, schnaps_real *flux);


void Compute_distribution_eq(Simulation * simu,schnaps_real * w_eq);

void Compute_relaxation(Simulation * simu,schnaps_real * w_eq);

void Compute_moments(Simulation * simu);

// equilibrium functions
schnaps_real feq_isothermal_D2Q9(int i_node,void *lattice,schnaps_real rho,schnaps_real ux,schnaps_real uy,schnaps_real uz,schnaps_real temp,schnaps_real p);
schnaps_real feq_isothermal_linearwave_D2Q9(int i_node,void *lattice,schnaps_real rho,schnaps_real ux,schnaps_real uy,schnaps_real uz,schnaps_real temp,schnaps_real p);
#endif
