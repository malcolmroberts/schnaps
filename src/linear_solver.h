#ifndef _LINEARSOLVER_H
#define _LINEARSOLVER_H

#include "global.h"
#include <stdbool.h>
#include "field.h"
#include "paralution_c.h"
#include "simulation.h"



typedef enum MatrixStorage{SKYLINE,CSR} MatrixStorage;
typedef enum Solver{LU,GMRES,PAR_GMRES,PAR_FGMRES,PAR_CG,PAR_BICGSTAB,PAR_AMG,PAR_LU,PAR_QR} Solver;
typedef enum PC{NONE,JACOBI,PAR_JACOBI,PAR_ILU,PAR_MULTICOLOREDSGS,PAR_MULTICOLOREDGS,PAR_MULTICOLOREDILU,PAR_AMG_PC,PAR_ELIMI,PHY_BASED_P1,PHY_BASED_P2,PHY_BASED_EXACT,EXACT} PC;

//! class for managing linear solvers
typedef struct LinearSolver{

  //! \brief number of equations
  int neq;

  //! \brief storage struct for the matrix
  //! the actual type depends on the chosen format
  void* matrix;

  //! name of the storage method;
  MatrixStorage storage_type; 

  //! \brief true if the matrix is symmetric
  bool is_sym;

  //! \brief true if the struct is initialized
  bool is_init;

  //! \brief true if the arrays are allocated
  bool is_alloc;

   //! \brief true if the matrix is assembly
  bool mat_is_assembly;

   //! \brief true if the matrix is assembly
  bool rhs_is_assembly;

  //! solver type;
  Solver solver_type;

  //! solver type;
  bool is_CG;

  //! name of the storage method;
  PC pc_type; 

  //! \brief solution of the linear system
  real* sol;
  //! \brief rhs of the linear system
  real* rhs;

  //! tolerance iterative solver
  real tol;

  //! restart for gmres
  int restart_gmres;
  //! number max of iteration
  int iter_max;
  

  //! \brief compute a matrix vector product
  //! \param[in] lsol the LinearSolver object containing matrix A
  //! \param[in] x a vector
  //! \param[out] prod Ax
  void (*MatVecProduct)(void* lsol,real x[],real prod[]);

} LinearSolver;



//! \brief init the LinearSolver structure with an empty matrix
//! \param[inout] lsol the LinearSolver object
//! \param[in] n number of equations
//! \param[in] matstor storage type (optional)
//! \param[in] solvtyp solver type (optional)
void InitLinearSolver(LinearSolver* lsol,int n,
		      MatrixStorage* matstor,
		      Solver* solvtyp);

//! \brief free the allocated arrays
//! \param[inout] lsol the LinearSolver object
//! \param[in] freeAll: an integer whose purpose is to choose whether we free everything. (testing purposes)
void FreeLinearSolver(LinearSolver* lsol);


//! \brief indicates that elem (i,j) is nonzero
//! \param[inout] lsol the LinearSolver object
//! \param[in] i row index
//! \param[in] j column index
void IsNonZero(LinearSolver* lsol,int i,int j); 

//! \brief allocate the variable-size arrays
//! \brief the nonzero positions should first be "switched on"
//! \param[inout] lsol the LinearSolver object
void AllocateLinearSolver(LinearSolver* lsol);

//! \brief add  to elem (i,j)  value val
//! \param[inout] lsol the LinearSolver object
//! \param[in] i row index
//! \param[in] j column index
//! \param[in] val value
void AddLinearSolver(LinearSolver* lsol,int i,int j,real val);

//! \brief set  to elem (i,j)  value val
//! \param[inout] lsol the LinearSolver object
//! \param[in] i row index
//! \param[in] j column index
//! \param[in] val value
void SetLinearSolver(LinearSolver* lsol,int i,int j,real val); 

//! \brief get elem (i,j)
//! \param[inout] lsol the LinearSolver object
//! \param[in] i row index
//! \param[in] j column index
//! \return value at pos (i,j)
real GetLinearSolver(LinearSolver* lsol,int i,int j); 


//! \brief display the matrix
//! \param[inout] lsol the LinearSolver object
void DisplayLinearSolver(LinearSolver* lsol); 

//! \brief compute a matrix vector product
//! \param[in] system the LinearSolver object containing matrix A
//! \param[in] x a vector
//! \param[out] prod Ax
void MatVect(void * system,real x[],real prod[]);

//! \brief compute the inplace LU decomposition
//! \param[inout] lsol the LinearSolver object
void LUDecompLinearSolver(LinearSolver* lsol);

//! \brief solve the linear system
//! \param[inout] lsol the LinearSolver object
void SolveLinearSolver(LinearSolver* lsol);


//! \brief copy vector
//! \param[in] x vector
//! \param[inout] prod is a copy of x
//! \param[in] N size
void Vector_copy(real x[],real prod[],int N);

//! \brief return the dot product
//! \param[in] x vector
//! \param[in] y vector
//! \param[in] N size
real Vector_prodot(real x[],real y[],int N);

//! \brief return the l2 norm
//! \param[in] x vector
//! \param[in] N size
real Vector_norm2(real x[],int  N);

//! \brief solve the linear system with paralution
//! \param[inout] lsol contains the matrices rhs and sol
void Solver_Paralution(LinearSolver* lsol);

//! \brief solve the linear system with the GMREs of the cerfacs
//! \param[inout] lsol contains the matrices rhs and sol
void GMRESSolver(LinearSolver* lsol);

//! \brief Jacobi preconditioner
//! \param[in] lsol contains the matrices rhs and sol
void Jacobi_PC(LinearSolver* lsol, real* sol, real* rhs);

//! \brief Exact LU preconditioner
//! \param[in] lsol contains the matrices rhs and sol
void Exact_PC(LinearSolver* lsol, real* sol, real* rhs);

#endif
