#ifndef _ADVANCEDLINEARSOLVER_H
#define _ADVANCEDLINEARSOLVER_H

#include "global.h"
#include <stdbool.h>
#include "field.h"
#include "paralution_c.h"
#include "simulation.h"
#include "linear_solver.h"

typedef struct JFLinearSolver{

  //! \brief number of equations
  int neq;

  //! brief eps for free Jacobian matrix
  real eps;

  //! solver type;
  Solver solver_type;

  //! name of the storage method;
  PC pc_type; 

  //! \brief solution of the linear system
  real* sol;
  //! \brief rhs of the linear system
  real* rhs;
  //! \brief sol at the time n
  real* soln;

    //! tolerance iterative solver
  real tol;

  //! restart for gmres
  int restart_gmres;
  //! number max of iteration
  int iter_max;

  //! \brief compute a matrix vector product
  //! \param[in] lsol the LinearSolver object containing matrix A
  //! \param[in] f the field
  //! \param[in] x a vector
  //! \param[out] prod Ax
  void (*MatVecProduct)(Simulation * simu,void* lsol,real x[],real prod[]);

  //! \brief compute the
  //! \param[in] simu the simulatio,n
  //! \param[in] lsol the LinearSolver object containing matrix A
  //! \param[in] solvector the solution at the time n
  //! \param[out] given the nonlinear vector for the free jacobian
  void (*NonlinearVector_computation)(Simulation * simu,void* lsol,real * solvector,real *nlvector);

} JFLinearSolver;



//! \brief init the LinearSolver structure with an empty matrix
//! \param[inout] lsol the LinearSolver object
//! \param[in] n number of equations
//! \param[in] solvtyp solver type (optional)
void InitJFLinearSolver(JFLinearSolver* lsol,int n,
		      Solver* solvtyp);

//! \brief free the allocated arrays
//! \param[inout] lsol the LinearSolver object
void FreeJFLinearSolver(JFLinearSolver* lsol);



//! \brief compute a matrix vector product
//! \param[in] system the LinearSolver object containing matrix A
//! \param[in] f a field
//! \param[in] x a vector
//! \param[out] prod Ax
void MatVecJacobianFree(Simulation * simu,void * system,real x[],real prod[]);

//! \brief solve the linear system
//! \param[inout] lsol the JFLinearSolver object
//! \param[in] simu contains the simuation
void SolveJFLinearSolver(JFLinearSolver* lsol,Simulation * simu);


//! \brief solve the linear system with the GMREs  and with a physic solver
//! \param[inout] lsol contains the matrices rhs and sol
//! \param[in] simu contains the simuation
void Advanced_GMRESSolver(LinearSolver* lsol, Simulation* simu);

//! \brief solve the linear system with a physic solver
//! \param[inout] lsol contains the matrices rhs and sol
//! \param[in] simu contains the simuation
void Advanced_SolveLinearSolver(LinearSolver* lsol, Simulation* simu);

//! \brief Lower order preconditioner
//! \param[in] lsol contains the matrices rhs and sol
void LowerOrderPC_Poisson(LinearSolver * lsol, Simulation *simu, real* globalSol, real*globalRHS);

#endif
