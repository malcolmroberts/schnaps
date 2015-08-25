#ifndef _SOLVERCONTINUOUS_H
#define _SOLVERCONTINUOUS_H

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "geometry.h"
#include "interpolation.h"
#include "linear_solver.h"
#include "simulation.h"

#define _Dirichlet_Poisson_BC (1)
#define _Periodic_Poisson_BC (2)


//! \brief a struct for sorting and pasting the
//! nodes of the DG mesh for obtaining a FE mesh
typedef struct FatNode{

  //! \brief index in the dg mesh
  int dg_index;
  //! \brief index in the fe mesh
  int fe_index;

  //! \brief physical coordinates of the node
  real x[3];
  //! \brief int converted coordinates for sorting and searching
  int x_int[3];


} FatNode;

//! \brief a struct for the resolution of the poisson equation:
//! conversion between a DG and Finite Element FE mesh
//! FE assembly and resolution, etc. 
typedef struct ContinuousSolver{

  //! \brief a simulation (gives the mesh and the charge)
  Simulation* simu;

  //! linear solver
  LinearSolver lsol;
  
  //! \brief number of FE nodes
  int nb_fe_nodes;

  //! \brief number of DG nodes
  int nb_dg_nodes;

  //! \brief number of FE degrees of fredom
  int nb_fe_dof;

  //! \brief number of DG degrees of freedom
  int nb_dg_dof;

  //! \brief node list with coordinates converted to integer
  //! and DG/FE indices
  FatNode* fn_list;

  //! \brief connectivity DG node -> FE node
  int* dg_to_fe_index;

   //! \brief number of element
  int nbel;

  //! \brief number of local nodes
  int nnodes;

  //! \brief number of local nodes
  int npgmacrocell;

  //! \brief list that marks boundary nodes
  int* is_boundary_node;

  //! \brief number of FE nodes
  int nb_phy_vars;

  //! \brief list of index for the variables
  int * list_of_var;

  //! \brief for dirchilet homogeneous, 2 for periodic 
  int type_bc;

  //! \brief Differential operator for scalar problems
  real diff_op[4][4];

  //! \brief Differential operator for 2d vectorial problems
  real diff_op2vec[4][4][4];

  //! \brief pointer on the function which assembles the rhs
  //! \param[inout] lsol a linear solver allocate
  //! \param[in] a continuous solver
  void (*rhs_assembly)(void * cs,LinearSolver* lsol);

  //! \brief pointer on the function which assembles the matrix
  //! \param[inout] lsol a linear solver allocate
  //! \param[in] a continuous solver
  void (*matrix_assembly)(void * cs,LinearSolver* lsol);

  //! \brief pointer on the function which assembly the post computation
  //! \param[inout] lsol a linear solver allocate
  //! \param[in] a continuous solver
  void (*postcomputation_assembly)(void * cs,LinearSolver* lsol);

   //! \brief pointer on the function which assembly the post computation
  //! \param[inout] lsol a linear solver allocate
  //! \param[in] a continuous solver
  void (*bc_assembly)(void * cs,LinearSolver* lsol);
  

} ContinuousSolver;



//! \brief compare two nodes (used by quicksort).
//! Lexicographic order on the coordinates converted to integers.
//! \param[in] a first node
//! \param[in] b second node
//! \returns a value v, v<0 if a<b, v=0 if a==b, v>0 if a>b
int CompareFatNode(const void* a,const void* b);

//! \brief build the fat nodes list from a field
//! \param[in] simu an initialized simulation
//! \param[out] fn_list an allocated, prepared and sorted list of fat nodes
//! \returns the size of the list
int BuildFatNodeList(Simulation *simu,FatNode* fn_list);

//! \brief init a poisson solver
//! \param[inout] ps a PoissonSolver struct
//! \param[in] simu a simulation
//! \param[in] type_bc the number of bc type
//! \param[in] nb_phy_vars the number of variable for the solver
void InitContinuousSolver(void* cs, Simulation* simu,
			  int type_bc,int nb_phy_vars,int * listvar);


//! \brief compute the discontinuous unknown using the continuous one
//! \param[inout] lsol a linear solver allocate
//! \param[in] a continuous solver
void ContinuousToDiscontinuous_Copy(ContinuousSolver * cs,LinearSolver* lsol);


//! \brief allocate matrix for continuous solver
//! \param[inout] lsol a linear solver allocate
//! \param[in] a continuous solver
void AllocateContinuousMatrix(void * cs,LinearSolver* lsol);

//! \brief apply dirichlet homogeneous bc for continuous solver
//! \param[inout] lsol a linear solver allocate
//! \param[in] a continuous solver
void ExactDirichletContinuousMatrix(void * cs,LinearSolver* lsol);

//! \brief apply dirichlet homogeneous bc for continuous solver
//! \param[inout] lsol a linear solver allocate
//! \param[in] a continuous solver
void PenalizedDirichletContinuousMatrix(void * cs,LinearSolver* lsol);

//! \brief solve a 2D poisson problem
//! \param[inout] ps a Poisson solver (field + linear solver + other parameters)
//! \param[in] type_bc the boundary condition type
//!  (1->dirichlet ; 2-> periodic)
void SolveContinuous2D(void * cs);

//! TODOODODODODOOD
void GenericOperatorScalar_Continuous(void * cs,LinearSolver* lsol);

//! TODOODODODODOOD
void GenericOperator2Vec_Continuous(void * cs,LinearSolver* lsol);

//! TODOODODODODOOD
void cat2CGVectors(ContinuousSolver* L1Solver,ContinuousSolver* L2Solver, real *L1, real *L2, real *L);

//! TODOODODODODOOD
void catGradients(ContinuousSolver* L1Solver,ContinuousSolver* L2Solver, real *L1, real *L2, real *L);

//! TODOODODODODOOD
void extract2CGVectors(ContinuousSolver* L1Solver,ContinuousSolver* L2Solver, real *L, real *L1, real *L2);

//! TODOODODODODOOD
void MatrixPoisson_Continuous(void * cs,LinearSolver* lsol);

//! \brief frees a ContinuousSolver object
//! \param[in] cs: a ContinuousSolver object
//! \param[in] free_simu: an integer whose purpose is to decide whether we free the Simulation object.
void freeContinuousSolver(ContinuousSolver* cs, int free_simu);
#endif
