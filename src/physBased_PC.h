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
  schnaps_real *rhs_prediction;
  //! \brief Right-hand side for the propagation (middle) step of the preconditioner
  schnaps_real *rhs_propagation;
  //! \brief Right-hand side for the correction step of the preconditioner
  schnaps_real *rhs_correction;

  Solver solver_prediction;
  Solver solver_propagation;
  Solver solver_correction;

  PC pc_prediction;
  PC pc_propagation;
  PC pc_correction;

  schnaps_real tol_prediction;
  schnaps_real tol_propagation;
  schnaps_real tol_correction;

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
void PiDgToCg(ContinuousSolver * cs,schnaps_real * rhsIn, schnaps_real * rhsOut);

// \brief Takes a vector in Continuous Galerkin, and returns the 
// same vector projected on the Discontinuous Galerkin discrete space.
// \param[in] cs: ContinuousSolver object.
// \param[in] rhsIn: Vector in CG.
// \param[out] rhsOut: Vector in DG.
void PiInvertCgToDg(ContinuousSolver * cs,schnaps_real * rhsIn, schnaps_real * rhsOut);


// \brief general initialisation of the solvers for eahc sub systems of the pc
// \param[inout] pb_pc: Physics-based Preconditioner object.
void Init_Parameters_PhyBasedPC(PB_PC* pb_pc);


// \brief Initialize the physics-based preconditioner with schur on the velocity boundary condition (u,n)=0
// \param[in] simu: Simulation object containing some run-related variables
// \param[inout] pb_pc: Physics-based Preconditioner object.
// \param[in] list_mat2assembly: Integer array. Tells which matrices shall be assembled.
void Init_PBPC_Wave_SchurVelocity_BCVelocity(Simulation *simu, PB_PC* pb_pc, int* list_mat2assemble);

// \brief Initialize the physics-based preconditioner with schur on the velocity boundary condition p=g
// \param[in] simu: Simulation object containing some run-related variables
// \param[inout] pb_pc: Physics-based Preconditioner object.
// \param[in] list_mat2assembly: Integer array. Tells which matrices shall be assembled.
void Init_PBPC_Wave_SchurVelocity_BCPressure(Simulation *simu, PB_PC* pb_pc, int* list_mat2assemble);

// \brief Initialize the physics-based preconditioner with schur on the pressure with boundary condition (u,n)=0
// \param[in] simu: Simulation object containing some run-related variables
// \param[inout] pb_pc: Physics-based Preconditioner object.
// \param[in] list_mat2assembly: Integer array. Tells which matrices shall be assembled.
void Init_PBPC_Wave_SchurPressure_BCVelocity(Simulation *simu, PB_PC* pb_pc, int* list_mat2assemble);

// \brief Initialize the physics-based preconditioner with schur on the pressure with boundary condition p=g
// \param[in] simu: Simulation object containing some run-related variables
// \param[inout] pb_pc: Physics-based Preconditioner object.
// \param[in] list_mat2assembly: Integer array. Tells which matrices shall be assembled.
void Init_PBPC_Wave_SchurPressure_BCPressure(Simulation *simu, PB_PC* pb_pc, int* list_mat2assemble);

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
void InitMat_ContinuousSolver(PB_PC* pb_pc, schnaps_real Dmat[4][4], schnaps_real L1Mat[4][4], schnaps_real L2Mat[4][4], schnaps_real U1Mat[4][4], schnaps_real U2Mat[4][4], schnaps_real Schurmat[4][4][4]);


// \brief Physics-based CG preconditioner for CG problem with schur on velocity
// \param[in] pb_pc a Physics-based preconditioner (Schur velocity)
// \param[inout] globalSol a solution of the preconditioner.
// \param[in] globalRHS a Right-hand-side containing all explicit and source terms.
// \param[in] simu a simulation
void PhyBased_PC_CG(PB_PC* pb_pc, Simulation *simu, schnaps_real* globalSol, schnaps_real*globalRHS);

// \brief Physics-based CG preconditioner for CG problem with schur on pressure
// \param[in] pb_pc a Physics-based preconditioner (Schur pressure)
// \param[inout] globalSol a solution of the preconditioner.
// \param[in] globalRHS a Right-hand-side containing all explicit and source terms.
// \param[in] simu a simulation
void PhyBased_PC_InvertSchur_CG(PB_PC* pb_pc, Simulation *simu, schnaps_real* globalSol, schnaps_real*globalRHS);


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

//! \brief pointer on the function which compute the BC flux for the condition (nabla p,n)+=0
//! \param[in] cs a continuous solver
//! \param[in] xpg a point of the mesh
//! \param[in] w a vector of unknowns
//! \param[in] vnorm a vector of normal
//! \param[inout] flux a vector of flux
void RobinFlux_SchurPressure(void * cs, schnaps_real * xpg, schnaps_real * w, schnaps_real *vnorm, schnaps_real * flux);

//! \brief pointer on the function which compute the BC flux for the condition (u,n)=0
//! \param[in] cs a continuous solver
//! \param[in] xpg a point of the mesh
//! \param[in] w a vector of unknowns
//! \param[in] vnorm a vector of normal
//! \param[inout] flux a vector of flux
void Dirichlet_Velocity(void * cs, schnaps_real * xpg, schnaps_real * w, schnaps_real *vnorm, schnaps_real * flux);

//! \brief pointer on the function which compute the BC flux for the operator dx 
//! \param[in] cs a continuous solver
//! \param[in] xpg a point of the mesh
//! \param[in] w a vector of unknowns
//! \param[in] vnorm a vector of normal
//! \param[inout] flux a vector of flux
void BoundaryTerm_Yderivative(void * cs, schnaps_real * xpg, schnaps_real * w, schnaps_real *vnorm, schnaps_real * flux);


//! \brief pointer on the function which compute the BC flux for the operator dy 
//! \param[in] cs a continuous solver
//! \param[in] xpg a point of the mesh
//! \param[in] w a vector of unknowns
//! \param[in] vnorm a vector of normal
//! \param[inout] flux a vector of flux
void BoundaryTerm_Xderivative(void * cs, schnaps_real * xpg, schnaps_real * w, schnaps_real *vnorm, schnaps_real * flux);

#endif
