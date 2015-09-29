#ifndef _IMPLICIT_LS_H
#define _IMPLICIT_LS_H


#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "simulation.h"
#include "linear_solver.h"


//! \brief Construct the profile of the linear solver
//! for the generic implicit linear solver
//! \param[inout] simu a simulation
//! \param[inout] solver a linear solver
void InitImplicitLinearSolver(Simulation *simu, LinearSolver *solver);

//! \brief Construct the profile of the linear solver
//! for the locally implicit scheme
//! \param[inout] field a field
void InitFieldImplicitSolver(field *fd);



//! \brief Assembly of the DG operator into a sparse matrix
//! computations of all the terms
//! \param[inout] simu a simulation
//! \param[inout] solver a linear solver
//! \param[in] theta the crank nicholson parameter
//!  \param[in] dt time step
void AssemblyImplicitLinearSolver(Simulation *simu, LinearSolver *solver,real theta, real dt);

//! \brief Assembly of the DG operator into a sparse matrix
//! computations of all the terms
//! case of the locally implicit scheme
//! \param[inout] field a field
//! \param[in] theta the crank nicholson parameter
//!  \param[in] dt time step
void AssemblyFieldImplicitSolver(field *fd,real theta, real dt);

//! \brief Assembly of the DG operator into a sparse matrix
//! prepare the matrix structure of the differential terms inside the fields
//! \param[inout] simu a simulation
//! \param[inout] solver a linear solver
//! \param[in] itest should be 0 (1 is used for debugging purposes)
void InternalCoupling(Simulation *simu,  LinearSolver *solver, int itest);

//! \brief Assembly of the DG operator into a sparse matrix
//! prepare the matrix structure of the differential terms inside the fields
//! case of the locally implicit scheme
//! \param[inout] field a field
//! \param[in] itest should be 0 (1 is used for debugging purposes)
void InternalLocalCoupling(field *fd, int itest);




//! \brief Assembly of the DG operator into a sparse matrix
//! prepare the matrix structure of the fluxes inside the fields
//! \param[inout] simu a simulation
//! \param[inout] solver a linear solver
//! \param[in] theta the crank nicholson parameter
//! \param[in] itest should be 0 (1 for debugging purposes)
void FluxCoupling(Simulation *simu,  LinearSolver *solver,int itest);

//! \brief Assembly of the DG operator into a sparse matrix
//! prepare the matrix structure of the fluxes inside the fields
//! case of the locally implicit scheme
//! \param[inout] field a field
//! \param[in] itest should be 0 (1 for debugging purposes)
void FluxLocalCoupling(field *fd,int itest);

//! \brief Assembly of the DG operator into a sparse matrix
//! prepare the matrix structure of the interface fluxes between fields
//! \param[inout] simu a simulation
//! \param[inout] solver a linear solver
//! \param[in] theta the crank nicholson parameter
//! \param[in] itest should be 0 (1 for debugging purposes)
void InterfaceCoupling(Simulation *simu,  LinearSolver *solver,int itest);

//! \brief Assembly of the DG operator into a sparse matrix
//! assembly of the differential terms
//! \param[inout] simu a simulation
//! \param[inout] solver a linear solver
//! \param[in] theta the crank nicholson parameter
//!  \param[in] dt time step
void InternalAssembly(Simulation *simu,  LinearSolver *solver,real theta, real dt);

//! \brief Assembly of DG differential terms
//! case of the locally implicit scheme
//! \param[inout] field a field
//! \param[in] theta the crank nicholson parameter
//!  \param[in] dt time step
void InternalLocalAssembly(field *fd, real theta, real dt);

//! \brief Assembly of the DG operator into a sparse matrix
//! assembly of the internal flxes of the fields
//! \param[inout] simu a simulation
//! \param[inout] solver a linear solver
//! \param[in] theta the crank nicholson parameter
//!  \param[in] dt time step
void FluxAssembly(Simulation *simu,  LinearSolver *solver,real theta, real dt);

//! \brief Assembly of the DG operator into a sparse matrix
//! assembly of the internal flxes of the fields
//! case of the local implicit solver
//! \param[inout] field a field
//! \param[in] theta the crank nicholson parameter
//!  \param[in] dt time step
void FluxLocalAssembly(field* fd,real theta, real dt);

//! \brief time-stepping by the crank nicholson scheme
//! \param[inout] simu a simulation
//! \param[in] tmax final time
//! \param[in] dt time step
void ThetaTimeScheme(Simulation *simu, real tmax, real dt);

//! \brief time-stepping by the crank nicholson scheme
//! macrocell local version
//! \param[inout] simu a simulation
//! \param[in] tmax final time
//! \param[in] dt time step
void LocalThetaTimeScheme(Simulation *simu, real tmax, real dt);

//! \brief time-stepping by the crank nicholson scheme StarPU version
//! macrocell local version 
//! \param[inout] simu a simulation
//! \param[in] tmax final time
//! \param[in] dt time step
void LocalThetaTimeScheme_SPU(Simulation *simu, real tmax, real dt);

///! \brief Assembly of the DG operator into a sparse matrix
//! assembly of the interface fluxes between the neighboring fields
//! \param[inout] simu a simulation
//! \param[inout] solver a linear solver
//! \param[in] theta the crank nicholson parameter
//!  \param[in] dt time step
void InterfaceAssembly(Simulation *simu,  LinearSolver *solver,real theta, real dt);

//! \brief Assembly of the DG operator into a sparse matrix
//! assembly of the right hand side of the linear system:
//! volume terms and boundary terms
//! \param[inout] simu a simulation
//! \param[inout] solver a linear solver
//! \param[in] theta the crank nicholson parameter
//!  \param[in] dt time step
void SourceAssembly(Simulation *simu,  LinearSolver *solver,real theta, real dt);

//! \brief Assembly of the DG operator into a sparse matrix
//! assembly of the right hand side of the linear system:
//! volume terms
//! case of the locally implicit scheme
//! \param[inout] field a field
//! \param[inout] solver a linear solver
//! \param[in] theta the crank nicholson parameter
//!  \param[in] dt time step
void SourceLocalAssembly(field *fd,real theta, real dt);

//! \brief Assembly of the DG operator into a sparse matrix
//! assembly of the mass terms
//! \param[inout] simu a simulation
//! \param[inout] solver a linear solver
void MassAssembly(Simulation *simu,  LinearSolver *solver);

//! \brief Assembly of the DG operator into a sparse matrix
//! assembly of the mass terms
//! case of the locally implicit scheme
//! \param[inout] field a field
void MassLocalAssembly(field *fd);



#endif
