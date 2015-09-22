#ifndef _SOLVERPOISSON_H
#define _SOLVERPOISSON_H


#include "collision.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "geometry.h"
#include "interpolation.h"
#include "linear_solver.h"
#include "simulation.h"
#include "solvercontinuous.h"



//! \brief compute the electirc field for poisson
//! \param[inout] lsol a linear solver allocate
//! \param[in] a continuous solver
void Computation_ElectricField_Poisson(void * cs,LinearSolver* lsol);


//! \brief init the rhs for poisson solver
//! \param[inout] lsol a linear solver allocate
//! \param[in] a continuous solver
void RHSPoisson_Continuous(void * cs,LinearSolver* lsol);


void RobinBoundaryConditionAssembly(void * cs,LinearSolver* lsol);





#endif
