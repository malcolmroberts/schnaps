#ifndef _IMPLICIT_H
#define _IMPLICIT_H


#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "simulation.h"
#include "linear_solver.h"


//! \brief Assembly of the sparse matrix
//! for the generic implicit linear solver
//! \param[inout] simu a simulation
//! \param[inout] solver a linear solver
void AssemblyImplicitLinearSolver(Simulation *simu, LinearSolver *solver);




//! \brief Construct the profile of the linear solver
//! for the generic implicit linear solver
//! \param[inout] simu a simulation
//! \param[inout] solver a linear solver
void InitImplicitLinearSolver(Simulation *simu, LinearSolver *solver);



#endif
