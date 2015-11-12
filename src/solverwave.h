#ifndef _SOLVERWAVE_H
#define _SOLVERWAVE_H

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "geometry.h"
#include "interpolation.h"
#include "advanced_linear_solver.h"
#include "solvercontinuous.h"

//! \brief compute the local operator fo the wave weak form
//! \param[inout] cs: a ContinuousSolver object
//! \param[in] theta: a coeffcient
//! \param[in] dt: a time step
void Wave_test(ContinuousSolver* cs, real theta, real dt);

//! \brief compute the boundary condition for Friedrics systems
//! \param[inout] cs: a ContinuousSolver object
void BoundaryConditionFriedrichsAssembly(void * cs);

//! \brief pointer on the function which compute the BC flux 
//! \param[in] cs a continuous solver
//! \param[in] xpg a point of the mesh
//! \param[in] w a vector of unknowns
//! \param[in] vnorm a vector of normal
//! \param[inout] flux a vector of flux
void Wave_BC_pressure_imposed(void * cs, real * xpg, real * w, real *vnorm, real * flux);

//! \brief pointer on the function which compute the BC flux 
//! \param[in] cs a continuous solver
//! \param[in] xpg a point of the mesh
//! \param[in] w a vector of unknowns
//! \param[in] vnorm a vector of normal
//! \param[inout] flux a vector of flux
void Wave_BC_normalvelocity_null(void * cs, real * xpg, real * w, real *vnorm, real * flux);

//! \brief construct the source associated to the source of the model 
//! \param[inout] cs a continuous solver
void Source_Assembly(void * cs);

#endif
