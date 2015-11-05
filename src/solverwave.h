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
//! \param[in] cs: a ContinuousSolver object
//! \param[in] theta: a coeffcient
//! \param[in] dt: a time step
void Wave_test(ContinuousSolver* cs, real theta, real dt);

//! \brief compute the source for Friedrics systems
//! \param[in] cs: a ContinuousSolver object
//! \param[in] lsol : a linear solver
void SourceFriedrichsAssembly(void * cs,LinearSolver* lsol);

//! \brief compute the boundary condition for Friedrics systems
//! \param[in] cs: a ContinuousSolver object
//! \param[in] lsol : a linear solver
void BoundaryConditionFriedrichsAssembly(void * cs,LinearSolver* lsol);

//! \brief pointer on the function which compute the BC flux 
//! \param[in] lsol a linear solver allocate
//! \param[in] cs a continuous solver
//! \param[in] xpg a point of the mesh
//! \param[in] w a vector of unknowns
//! \param[in] vnorm a vector of normal
//! \param[inout] flux a vector of flux
void Wave_BC_pressure_imposed(void * cs,LinearSolver* lsol, real * xpg, real * w, real *vnorm, real * flux);

//! \brief pointer on the function which compute the BC flux 
//! \param[in] lsol a linear solver allocate
//! \param[in] cs a continuous solver
//! \param[in] xpg a point of the mesh
//! \param[in] w a vector of unknowns
//! \param[in] vnorm a vector of normal
//! \param[inout] flux a vector of flux
void Wave_BC_normalvelocity_null(void * cs,LinearSolver* lsol, real * xpg, real * w, real *vnorm, real * flux);

#endif
