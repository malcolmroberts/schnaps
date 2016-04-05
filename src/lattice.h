#ifndef _LATTICE_H
#define _LATTICE_H

#include "model.h"
#include "field.h"
#include "simulation.h"
#include "global.h"

//! \brief flux for the lattice model
//! \param[in] wL,wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void Lattice_NumFlux(schnaps_real *wL, schnaps_real *wR, schnaps_real *vn, schnaps_real *flux);
void Lattice_OneNodeNumFlux(schnaps_real *wL, schnaps_real *wR, schnaps_real *vn, schnaps_real *flux);

//! \brief computation of equilibrium distribution (global wrapper)
//! \para[in] simu a simulation object
//! \para[inout] w_eq real array  
void Compute_distribution_eq(Simulation * simu,schnaps_real * w_eq);

void Compute_relaxation(Simulation * simu,schnaps_real * w_eq);

void Compute_moments(Simulation * simu);

// equilibrium functions

//! \brief equiilibrium distribution for the euler/navier stokes isothermal model
//! \paral[in] i_node index of velocity node 
//! \para[in]  lattice object 
//! \para[in] rho, ux, uy, iz, temp, p macroscopic variabbles
//! \returns
schnaps_real feq_isothermal_D2Q9(int i_node,void *lattice,schnaps_real rho,schnaps_real ux,schnaps_real uy,schnaps_real uz,schnaps_real temp,schnaps_real p);
schnaps_real feq_isothermal_linearwave_D2Q9(int i_node,void *lattice,schnaps_real rho,schnaps_real ux,schnaps_real uy,schnaps_real uz,schnaps_real temp,schnaps_real p);
#endif
