#ifndef _IMPLICIT_LS_H
#define _IMPLICIT_LS_H


#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "linear_solver.h"

//! \brief Construct the profile of the linear solver
//! for the generic implicit linear solver
//! \param[inout] simu a simulation
//! \param[inout] solver a linear solver
void InitImplicitLinearSolver(field *f, LinearSolver *solver);

//! \brief Assembly of the DG operator into a sparse matrix
//! computations of all the terms
//! \param[inout] simu a simulation
//! \param[inout] solver a linear solver
//! \param[in] theta the crank nicholson parameter
//!  \param[in] dt time step
void AssemblyImplicitLinearSolver(field *f, LinearSolver *solver,real theta,
				  real tnow, real dt);


//! \brief Assembly of the DG operator into a sparse matrix
//! prepare the matrix structure of the differential terms inside the fields
//! \param[inout] simu a simulation
//! \param[inout] solver a linear solver
//! \param[in] itest should be 0 (1 is used for debugging purposes)
void InternalCoupling(field *f, LinearSolver *solver, int itest);

//! \brief Assembly of the DG operator into a sparse matrix
//! prepare the matrix structure of the fluxes inside the fields
//! \param[inout] simu a simulation
//! \param[inout] solver a linear solver
//! \param[in] theta the crank nicholson parameter
//! \param[in] itest should be 0 (1 for debugging purposes)
void FluxCoupling(field *f, LinearSolver *solver,int itest);

//! \brief Assembly of the DG operator into a sparse matrix
//! prepare the matrix structure of the interface fluxes between fields
//! \param[inout] simu a simulation
//! \param[inout] solver a linear solver
//! \param[in] theta the crank nicholson parameter
//! \param[in] itest should be 0 (1 for debugging purposes)
void InterfaceCoupling(field *f, LinearSolver *solver,int itest);

//! \brief Assembly of the DG operator into a sparse matrix
//! assembly of the differential terms
//! \param[inout] simu a simulation
//! \param[inout] solver a linear solver
//! \param[in] theta the crank nicholson parameter
//!  \param[in] dt time step
void InternalAssembly(field *f, LinearSolver *solver,real theta, real dt);

//! \brief Assembly of the DG operator into a sparse matrix
//! assembly of the internal flxes of the fields
//! \param[inout] simu a simulation
//! \param[inout] solver a linear solver
//! \param[in] theta the crank nicholson parameter
//!  \param[in] dt time step
void FluxAssembly(field *f, LinearSolver *solver,real theta, real dt);

//! ADD DESCRIPTION
void ThetaTimeScheme(field *f, real tmax, real dt);

//! \brief Assembly of the DG operator into a sparse matrix
//! assembly of the interface fluxes between the neighboring fields
//! \param[inout] simu a simulation
//! \param[inout] solver a linear solver
//! \param[in] theta the crank nicholson parameter
//!  \param[in] dt time step
void InterfaceAssembly(field *f, LinearSolver *solver,real theta, real dt);

//! \brief Assembly of the DG operator into a sparse matrix
//! assembly of the right hand side of the linear system:
//! volume terms and boundary terms
//! \param[inout] simu a simulation
//! \param[inout] solver a linear solver
//! \param[in] theta the crank nicholson parameter
//! \param[in] tnow
void SourceAssembly(field *f, LinearSolver *solver, real theta, real tnow,
		    real dt);

//! \brief Assembly of the DG operator into a sparse matrix
//! assembly of the mass terms
//! \param[inout] simu a simulation
//! \param[inout] solver a linear solver
void MassAssembly(field *f, LinearSolver *solver);




#endif
