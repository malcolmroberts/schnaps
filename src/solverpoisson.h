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


//! \brief solve a 1D poisson problem
//! \param[in] simu a simulation
//! \param[in] w the field values (for computing the charge
//! , returning the potential and the electric field)
//! \param[in] type_bc the boundary condition type
//!  (1->dirichlet ; 2-> periodic)
//! \param[in] bc_l left boundary value (dirichlet case)
//! \param[in] bc_r right boundary value (dirichlet case)
//! \param[in] solver_sys linear solver
//! \param[in] precon preconditionner
void SolvePoisson1D(Simulation *simu,real * w,
		    int type_bc, real bc_l, real bc_r,Solver solver_sys, PC precon);


#endif
