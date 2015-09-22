#ifndef _PHYBASED_PC_H
#define _PHYBASED_PC_H


#include "collision.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "geometry.h"
#include "interpolation.h"
#include "linear_solver.h"
#include "simulation.h"
#include "solvercontinuous.h"
#include "solverwave.h"
#include "solverpoisson.h"

//! \brief Struct managing Physics Based Preconditioners
typedef struct PB_PC{

  //! \brief list containing integers (0 or 1).
  // Describes which matrices should be assembled (ordered as in the structure)
  int *list_mat2assemble;

  //! \brief D ContinuousSolver (matrix from Schur decomposition)
  ContinuousSolver D;
  //! \brief L1 ContinuousSolver (matrix from Schur decomposition)
  ContinuousSolver L1;
  //! \brief L2 ContinuousSolver (matrix from Schur decomposition)
  ContinuousSolver L2;
  //! \brief U1 ContinuousSolver (matrix from Schur decomposition)
  ContinuousSolver U1;
  //! \brief U2 ContinuousSolver (matrix from Schur decomposition)
  ContinuousSolver U2;
  //! \brief Schur ContinuousSolver (matrix from Schur decomposition)
  ContinuousSolver Schur;

  //! \brief Right-hand side for the prediciton step of the preconditioner
  real *rhs_prediction;
  //! \brief Right-hand side for the propagation (middle) step of the preconditioner
  real *rhs_propagation;
  //! \brief Right-hand side for the correction step of the preconditioner
  real *rhs_correction;

  Solver solver_prediction;
  Solver solver_propagation;
  Solver solver_correction;

  PC pc_prediction;
  PC pc_propagation;
  PC pc_correction;

  real tol_prediction;
  real tol_propagation;
  real tol_correction;

  int itermax_prediction;
  int itermax_propagation;
  int itermax_correction;
  
  int restart_prediction;
  int restart_propagation;
  int restart_correction;

  // \brief 0 if the system is linear and 1 if not
  int nonlinear;

  //! \brief provides the implementation of each operator of the corresponding system's Schur decomposition.
  //! \param[inout] pb_pc: a PB_PC object.
  //! \param[in] offset: An integer pointing to the "j" component of the variables.
  void (*mat_assembly)(void* pb_pc, int offset);

  LinearSolver Schur2;

} PB_PC;

// \brief Takes a vector in Discontinuous Galerkin, and returns the 
// same vector projected on the Continuous Galerkin discrete space.
// \param[in] cs: ContinuousSolver object.
// \param[in] rhsIn: Vector in DG.
// \param[out] rhsOut: Vector in CG.
void VectorDgToCg (ContinuousSolver * ps,real * rhsIn, real * rhsOut);
void PiDgToCg(ContinuousSolver * cs,real * rhsIn, real * rhsOut);

// \brief Takes a vector in Continuous Galerkin, and returns the 
// same vector projected on the Discontinuous Galerkin discrete space.
// \param[in] cs: ContinuousSolver object.
// \param[in] rhsIn: Vector in CG.
// \param[out] rhsOut: Vector in DG.
void VectorCgToDg(ContinuousSolver * cs, real * rhsIn, real * rhsOut);
void PiInvertCgToDg(ContinuousSolver * cs,real * rhsIn, real * rhsOut);


// \brief Initialize the data of the solvers.
// \param[inout] pb_pc: Physics-based Preconditioner object.
void Init_Parameters_PhyBasedPC(PB_PC* pb_pc);

// \brief Initialize the physics-based preconditioner with schur on the velocity
// \param[in] simu: Simulation object containing some run-related variables
// \param[inout] pb_pc: Physics-based Preconditioner object.
// \param[in] list_mat2assembly: Integer array. Tells which matrices shall be assembled.
void Init_PhyBasedPC_SchurVelocity_Wave(Simulation *simu, PB_PC* pb_pc, int* list_mat2assemble);

// \brief Initialize the physics-based preconditioner with schur on the pressure
// \param[in] simu: Simulation object containing some run-related variables
// \param[inout] pb_pc: Physics-based Preconditioner object.
// \param[in] list_mat2assembly: Integer array. Tells which matrices shall be assembled.
void Init_PhyBasedPC_SchurPressure_Wave(Simulation *simu, PB_PC* pb_pc, int* list_mat2assemble);

// \brief Initialize Generic matrices for differential opertors -> Has to be tuned to the considered problem.
// Schur decomposition is given by the following definition of the matrices:
// |  D   |  U1    U2 |
// |      |           |
// |------|-----------|
// |  L1  |           |
// |      |   Schur   |
// |      |           |
// |  L2  |           |
// \param[inout] pb_pc: Physics-based preconditioner
// \param[in] Dmat: DMatrix
// \param[in] L1mat: L1Matrix
// \param[in] L2mat: L2Matrix
// \param[in] U1mat: U1Matrix
// \param[in] U2mat: U2Matrix
// \param[in] Schurmat: SchurMatrix
void InitMat_ContinuousSolver(PB_PC* pb_pc, real Dmat[4][4], real L1Mat[4][4], real L2Mat[4][4], real U1Mat[4][4], real U2Mat[4][4], real Schurmat[4][4][4]);

// \brief Physics-based CG preconditioner for DG problem
// \param[in] pb_pc: Physics-based preconditioner (contains all the Schur decomposition)
// \param[out] globalSol: Stores the solution of the preconditioner.
// \param[in] globalRHS: Right-hand-side containing all explicit and source terms.
void PhyBased_PC_DG(PB_PC* pb_pc, Simulation *simu, real* globalSol, real*globalRHS);

// \brief Physics-based CG preconditioner for CG problem
// \param[in] pb_pc: Physics-based preconditioner (contains all the Schur decomposition)
// \param[out] globalSol: Stores the solution of the preconditioner.
// \param[in] globalRHS: Right-hand-side containing all explicit and source terms.
void PhyBased_PC_CG(PB_PC* pb_pc, Simulation *simu, real* globalSol, real*globalRHS);

// \brief Physics-based CG preconditioner for CG problem
// \param[in] pb_pc: Physics-based preconditioner (contains all the Schur decomposition)
// \param[out] globalSol: Stores the solution of the preconditioner.
// \param[in] globalRHS: Right-hand-side containing all explicit and source terms.
void PhyBased_PC_InvertSchur_CG(PB_PC* pb_pc, Simulation *simu, real* globalSol, real*globalRHS);

// \brief Frees any PB_PC object
// \param[inout] pb_pc a PhysicsBased_PreConditioner 
void freePB_PC(PB_PC* pb_pc);

// \brief Function assembling all differential operator matrices (schur of the velocity)
// \param[in] pb_pc: The working preconditioner.
void GenericOperator_PBPC_Velocity(PB_PC* pb_pc);

// \brief Function assembling all differential operator matrices (schur of the pressure)
// \param[in] pb_pc: The working preconditioner.
void GenericOperator_PBPC_Pressure(PB_PC* pb_pc);

// \brief Function resetting all but problem matrices.
// \param[in] pb_pc: The working preconditioner.
void reset(PB_PC* pb_pc);

void RobinFlux_pressure(void * cs,LinearSolver* lsol, real * xpg, real * w, real *vnorm, real * flux);


#endif
