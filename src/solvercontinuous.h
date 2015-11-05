#ifndef _SOLVERCONTINUOUS_H
#define _SOLVERCONTINUOUS_H

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "geometry.h"
#include "interpolation.h"
#include "advanced_linear_solver.h"

#define _Dirichlet_Poisson_BC (1)
#define _Periodic_Poisson_BC (2)

//! \brief Scalar diffusion operator
typedef struct SDO{

  //! \brief physical coordinates of the node
  real DO[4][4];

} SDO;

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

  //! \brief for Neumann (or (u,n)=0 for wave), robin ( or p imposed for wave),
  int type_bc;

  //! \brief Differential operator for 2d vectorial (some variables) problems
  SDO * diff_op;

  //! \brief FluxMatrix is the matrix of the hyperbolic model
  real ** FluxMatrix;
  
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

  //! \brief pointer on the function which compute the BC flux 
  //! \param[in] lsol a linear solver allocate
  //! \param[in] cs a continuous solver
  //! \param[in] xpg a point of the mesh
  //! \param[in] w a vector of unknowns
  //! \param[in] vnorm a vector of normal
  //! \param[inout] flux a vector of flux
  void (*bc_flux)(void * cs,LinearSolver* lsol, real * xpg, real * w, real *vnorm, real * flux);

  

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

//! \brief apply dirichlet inhomogeneous bc for continuous solver
//! \param[inout] lsol a linear solver allocate
//! \param[in] a continuous solver
void ExactDirichletContinuousMatrix(void * cs,LinearSolver* lsol);


//! \brief apply dirichlet inhomogeneous bc for continuous solver
//! \param[inout] lsol a linear solver allocate
//! \param[in] a continuous solver
void PenalizedDirichletContinuousMatrix(void * cs,LinearSolver* lsol);

//! \brief construct the source 
//! \param[inout] lsol a linear solver allocate
//! \param[in] a continuous solver
void Source_Assembly(void * cs,LinearSolver* lsol);



//! \brief solve a 2D poisson problem
//! \param[inout] ps a Poisson solver (field + linear solver + other parameters)
//! \param[in] type_bc the boundary condition type
//!  (1->dirichlet ; 2-> periodic)
void SolveContinuous2D(void * cs);

//! TODOODODODODOOD
void GenericOperator_Continuous(void * cs,LinearSolver* lsol);

//! TODOODODODODOOD
void cat2CGVectors(ContinuousSolver* L1Solver,ContinuousSolver* L2Solver, real *L1, real *L2, real *L);

//! TODOODODODODOOD
void extract2CGVectors(ContinuousSolver* L1Solver,ContinuousSolver* L2Solver, real *L, real *L1, real *L2);


//! \brief frees a ContinuousSolver object
//! \param[in] cs: a ContinuousSolver object
void freeContinuousSolver(ContinuousSolver* cs);

#endif
