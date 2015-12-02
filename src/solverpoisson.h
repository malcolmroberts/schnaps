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
//! \param[inout] a continuous solver
void Computation_ElectricField_Poisson(void * cs);

//! \brief init the rhs for poisson solver
//! \param[inout] a continuous solver
void RHSPoisson_Continuous(void * cs);

//! \brief init the matrix for 1D poisson solver
//! \param[inout] a continuous solver
void ContinuousOperator_Poisson1D(void * cs);

//! \brief init the matrix for 2D poisson solver
//! \param[inout] a continuous solver
void ContinuousOperator_Poisson2D(void * cs);

//! \brief init the rhs for poisson solver
//! \param[inout] a continuous solver
void Periodic_BoundaryCondition_Poisson1D(void * cs);

//! \brief compute the boundary condition for poisson solver
//! \param[inout] a continuous solver
void RobinBoundaryConditionAssembly(void * cs);

//! \brief pointer on the function which compute the BC flux for Robin condition (nabla p,n)+alpha p=beta g
//! \param[in] cs a continuous solver
//! \param[in] xpg a point of the mesh
//! \param[in] w a vector of unknowns
//! \param[in] vnorm a vector of normal
//! \param[inout] flux a vector of flux
void RobinFlux(void * cs, real * xpg, real * w, real *vnorm, real * flux);


//! \brief Lower order preconditioner
//! \param[in] lsol contains the matrices rhs and sol
void LowerOrderPC_Poisson(LinearSolver * lsol, Simulation *simu, real* globalSol, real*globalRHS);

void Interpolation1D_P1_Pq(ContinuousSolver *cs_LowOrder,int HighOrder, real * VecIn, real * VecOut);
void Restriction1D_Pq_P1(ContinuousSolver *cs_LowOrder,int HighOrder, real * VecIn, real * VecOut);






#endif
