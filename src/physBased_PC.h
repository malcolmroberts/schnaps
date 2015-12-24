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
  //! \param[in] cs: a ContinuousSolver object.
  void (*mat_assembly)(void* pb_pc, ContinuousSolver* cs);

  //! \brief provides the implementation of the local matrices of the corresponding system's Schur decomposition
  //! \param[inout] pb_pc: a PB_PC object.
  //! \param[in] var: an array of varialbes and their derivatives.
  void (*loc_mat_assembly)(void* pb_pc, real* var);

  //! \brief provides the implementation of the rhs of the corresponding system's Schur decomposition.
  //! \param[inout] pb_pc: a PB_PC object.
  //! \param[in] cs: a ContinuousSolver object.
  void (*rhs_assembly)(void* pb_pc, ContinuousSolver* cs);

  //! \brief provides the implementation of the local right-hand sides of the corresponding system's Schur decomposition
  //! \param[inout] pb_pc: a PB_PC object.
  //! \param[in] var: an array of varialbes and their derivatives.
  //! \param[in] coeff: a weighting coefficient
  //! \param[out] loc_rhs: a local array for right-hand side.
  void (*loc_rhs_assembly)(void* pb_pc, real* var, real coeff, real* loc_rhs);

  //! \brief provides the implementation of the boundary conditions for the rhs
  //! \param[inout] cs: a ContinuousSolver object.
  void (*bc_assembly)(void* cs);

  //! \brief provides the implementation of the local source terms for the rhs
  //! \param[in] cs: a ContinuousSolver object
  //! \param[in] xpg: the real coordinates of the Gauss point
  //! \param[in] w: the "imposedData" value at the Gauss point
  //! \param[in] vnorm: the normal to the cell
  //! \param[out] flux: the value of the flux on the boundary
  void (*loc_bc_assembly)(void* cs, real* xpg, real* w, real* vnorm, real* flux);

  //! \brief provides the implementation of the source terms for the rhs
  //! \param[inout] cs: a ContinuousSolver object.
  void (*source_assembly)(void* cs);

  //! \brief provides the implementation of the preconditioner main iteration process.
  //! \param[inout] pb_pc: a PB_PC object
  //! \param[in] simu: a simulation object
  //! \param[out] globalSol: the solution after one preconditioner solve
  //! \param[in] globalRHS: Right-hand side of the system at the current time.
  void (*solvePC)(void* pb_pc, Simulation* simu, real* globalSol, real* globalRHS);

  LinearSolver Schur2;

} PB_PC;

// \brief Takes a vector in Discontinuous Galerkin, and returns the 
// same vector projected on the Continuous Galerkin discrete space.
// \param[in] cs: ContinuousSolver object.
// \param[in] nbVarIn: Number of variable of rhsIn (should be >= cs->nb_phy_vars)
// \param[in] rhsIn: Vector in DG.
// \param[out] rhsOut: Vector in CG.
void PiDgToCg(ContinuousSolver * cs, int nbVarIn, real * rhsIn, real * rhsOut);

// \brief Takes a vector in Continuous Galerkin, and returns the 
// same vector projected on the Discontinuous Galerkin discrete space.
// \param[in] cs: ContinuousSolver object.
// \param[in] nbVarOut: Number of variable of rhsOut (should be >= cs->nb_phy_vars)
// \param[in] rhsIn: Vector in CG.
// \param[out] rhsOut: Vector in DG.
void PiInvertCgToDg(ContinuousSolver * cs, int nbVarOut, real * rhsIn, real * rhsOut);


// \brief general initialisation of the solvers for eahc sub systems of the pc
// \param[inout] pb_pc: Physics-based Preconditioner object.
void Init_Parameters_PhyBasedPC(PB_PC* pb_pc);

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


// \brief Physics-based CG preconditioner for CG problem with schur on velocity
// \param[in] pb_pc a Physics-based preconditioner (Schur velocity)
// \param[inout] globalSol a solution of the preconditioner.
// \param[in] globalRHS a Right-hand-side containing all explicit and source terms.
// \param[in] simu a simulation
void PhyBased_PC_CG(void* pb_pc, Simulation *simu, real* globalSol, real*globalRHS);

// \brief Physics-based CG preconditioner for CG problem with schur on pressure
// \param[in] pb_pc a Physics-based preconditioner (Schur pressure)
// \param[inout] globalSol a solution of the preconditioner.
// \param[in] globalRHS a Right-hand-side containing all explicit and source terms.
// \param[in] simu a simulation
void PhyBased_PC_InvertSchur_CG(void* pb_pc, Simulation *simu, real* globalSol, real*globalRHS);


// \brief Frees any PB_PC object
// \param[inout] pb_pc a PhysicsBased_PreConditioner 
void freePB_PC(PB_PC* pb_pc);

// \brief Function assembling all differential operator matrices for non linear cases (variable dependent)
// \param[inout] pb_pc: The working preconditioner.
// \param[in] cs: a ContinuousSolver object containing the solution at current time.
void GenericOperator_PBPC_NonLinear(void* pb_pc, ContinuousSolver* cs);

// \brief Function assembling RHS of the system for non linear cases (variable dependent)
// \param[inout] pb_pc: The working preconditioner.
// \param[in] cs: a ContinuousSolver object containing the solution at current time.
void GenericRHS_PBPC_NonLinear(void* pb_pc, ContinuousSolver* cs);

// \brief Function assembling all differential operator matrices
// \param[inout] pb_pc: The working preconditioner.
// \param[in] cs: a ContinuousSolver object.
void GenericOperator_PBPC(void* pb_pc, ContinuousSolver* cs);

// \brief Function resetting all but problem matrices.
// \param[in] pb_pc: The working preconditioner.
void reset(PB_PC* pb_pc);

//! \brief pointer on the function which compute the BC flux for the condition (nabla p,n)+=0
//! \param[in] cs a continuous solver
//! \param[in] xpg a point of the mesh
//! \param[in] w a vector of unknowns
//! \param[in] vnorm a vector of normal
//! \param[inout] flux a vector of flux
void RobinFlux_SchurPressure(void * cs, real * xpg, real * w, real *vnorm, real * flux);

//! \brief pointer on the function which compute the BC flux for the condition (u,n)=0
//! \param[in] cs a continuous solver
//! \param[in] xpg a point of the mesh
//! \param[in] w a vector of unknowns
//! \param[in] vnorm a vector of normal
//! \param[inout] flux a vector of flux
void Dirichlet_Velocity(void * cs, real * xpg, real * w, real *vnorm, real * flux);

//! \brief pointer on the function which compute the BC flux for the operator dx 
//! \param[in] cs a continuous solver
//! \param[in] xpg a point of the mesh
//! \param[in] w a vector of unknowns
//! \param[in] vnorm a vector of normal
//! \param[inout] flux a vector of flux
void BoundaryTerm_Yderivative(void * cs, real * xpg, real * w, real *vnorm, real * flux);


//! \brief pointer on the function which compute the BC flux for the operator dy 
//! \param[in] cs a continuous solver
//! \param[in] xpg a point of the mesh
//! \param[in] w a vector of unknowns
//! \param[in] vnorm a vector of normal
//! \param[inout] flux a vector of flux
void BoundaryTerm_Xderivative(void * cs, real * xpg, real * w, real *vnorm, real * flux);

#endif
