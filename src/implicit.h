#ifndef _IMPLICIT_H
#define _IMPLICIT_H


#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "simulation.h"
#include "linear_solver.h"



//! \brief Construct the profile of the linear solver
//! for the generic implicit linear solver
//! \param[inout] simu a simulation
void InitImplicitLinearSolver(Simulation *simu, LinearSolver *solver);



#endif
