#ifndef _IMPLICIT_LS_H
#define _IMPLICIT_LS_H


#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "simulation.h"
#include "linear_solver.h"


//! \brief Assembly of the sparse matrix
//! for the generic implicit linear solver
//! \param[inout] simu a simulation
//! \param[inout] solver a linear solver
//void AssemblyImplicitLinearSolver(Simulation *simu, LinearSolver *solver);

void AssemblyImplicitLinearSolver(Simulation *simu, LinearSolver *solver,real theta, real dt);



//! \brief Construct the profile of the linear solver
//! for the generic implicit linear solver
//! \param[inout] simu a simulation
//! \param[inout] solver a linear solver
void InitImplicitLinearSolver(Simulation *simu, LinearSolver *solver);


void InternalCoupling(Simulation *simu,  LinearSolver *solver, int itest);
void FluxCoupling(Simulation *simu,  LinearSolver *solver,int itest);




#endif
