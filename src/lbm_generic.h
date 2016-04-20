#ifndef _LBMGENERIC_H
#define _LBMGENERIC_H
#include <string.h>
#include <stdbool.h>
#include <global.h>
#include "lbm_models_data.h"
#include "model.h"
#include "simulation.h"

//
static const int LBM_MAX_POLYSTRSIZE=1024;
static const int LBM_MAX_MACROVARNAMESIZE=32;
static const char LBM_POLYVARNAMES[3]={'X','Y','Z'};
//
//
//!\brief compact structure for encoding multivariate polynomials describing velocity moments. 
typedef struct MPolyDescriptor{
  int ndim ; //! space dimension i.e number of variables
  int nc ; //! number of multivariate monomials
  schnaps_real *c; //! [nc] array monomial coefficients
  int **e; //! [nc][nd] array exponents multiplets
  char *name; //! human readable name
  char *poly_string; //! 
} MPolyDescriptor;
//
static const MPolyDescriptor MPolyDescriptor_NULL; // zeroed out structure
//!\brief Constructor 
//!\param[inout] mom a MPolyDescriptor
//!\param[in] ndim dimension of velocity space
//!\param[in] nc number of monomials
//!\param[in] name a human readable name
void NewMPolyDescriptor(MPolyDescriptor *mpd, int ndim, int nc,char *name);
//!\brief MPolyDescriptor Destructor
//!\param[inout] mom a MPolyDescriptor
void DestroyMPolyDescriptor(MPolyDescriptor *mpd);
//!\brief Sets the actual values of coefficients and exponents. Creates the string polynomial representation.
//!\param[inout] mom a MPolyDescriptor
//!\param[in] c [nc] array of monomial coefficients
//!\param[in] e [nc][ndim] array of exponents 
void SetMPolyDescriptor(MPolyDescriptor *mpd, schnaps_real c_in[mpd->nc],int e_in[mpd->nc][mpd->ndim]);
//!\brief Output name and Polynomial form to standard output
//!\param[in] mom a MPolyDescriptor
void DisplayMPolyDescriptor(MPolyDescriptor *mpd);
//!\brief Computes the value of multivariate polynomial
//!\param[in] mom a MPolyDescriptor
//!\param[in] V [ndim] array of values
//!\returns value of the polynomial at V 
schnaps_real EvaluatePoly(MPolyDescriptor *mpd,schnaps_real *V);
//
//!\brief create MPolyDescriptor from commonly used hard-coded models
//!\param[inout] mpd MPolyDescriptor
//!\param[in] pointer to predefined model
//!\param[in] name a human readable name, if set to NULL the default name encoded in the model is used.
void NewMPolyDescriptorFromModel(MPolyDescriptor *mpd,const MPolyData *model,char *name);
//!\brief A structure describing a LBM model. 
typedef struct LBModelDescriptor{
  int d ; //!velocity space dimension
  int q ; //!number of velocity nodes
  int nb_macro;//!number of conserved macro quantities NB : for multi-lattice models nc>q may be true.
  char **macro_names; //! human readable names of macroscopic quantities
//  int nr;//!number of macro quantities to relax;
  schnaps_real cref; //! velocity scale to apply to reference nodes
  schnaps_real **vi; //! [q][d} array of velocity nodes;
  schnaps_real vmax; //! maximum velocity 
  MPolyDescriptor *Moments ;//![q] array of moments descriptors for multi-relaxation models;
  int *inode_min;//! [q] array ; for each moment starting index of node range for moment computation
  int *inode_max;//! [q] array ; for each moment end index of node range for moment computation
  schnaps_real **M; //! [q][q] moment matrix allowing to compute moments from f 
  //schnaps_real **Minv; //! [q][q] inverse of moment matrix
  schnaps_real *s; //! [q] array of relaxation parameters, will apply to either f or moments
//  bool  *is_relaxed; //! [q] switch relaxation on/off
//  bool  *is_conserved; //! [nb_macro] switch conservation of moments
//  bool  global_relax; //! special case, if true, relaxation does not need to be done in moment space
  //
  void (*f_to_macro)(schnaps_real *f, schnaps_real *w); // computation of macro-quantities from f
  void (*macro_to_f)(schnaps_real *w,schnaps_real *f); // computation of f from macro-quantities
  schnaps_real (*feq)(int inode,int nb_macro,schnaps_real *w); // computation of equilibrium function
  schnaps_real (*meq)(int imoment,int nb_macro,schnaps_real *w); // computation of equilibrium moment 
  //
  void *(model_spec_params);// pointer to model specific parameters each model needing it should have its specific struc to hold
  // such parameters and init them in the corresponding set routine; 
} LBModelDescriptor;
//
static const LBModelDescriptor LBModelDescriptor_NULL; //! zeroed out stucture
//
//!\brief LBModelDescriptor Constructor / Initializer
//!\param[inout] lb a LBModelDescriptor
//!\param[in] d dimension of velocity space
//!\param[in] q number of velocity nodes
void NewLBModelDescriptor(LBModelDescriptor *lb,int d, int nb_macro,int q);
//!\brief LBModelDescriptor Destructor
//!\param[inout] lb a LBModelDescriptor
//!\paral[in] d dimension of velocity space
//!\paral[in] q number of velocity nodes
void DestroyLBModelDescriptor(LBModelDescriptor *lb);
//!\brief Moment matric computation (assuming Moments,cref and nodes are set)
void ComputeLBModelDescriptorMomentMatrix(LBModelDescriptor *lb);
//!brief print moment basis associated polynomials
void DisplayLBModelDescriptorMomentPoly(LBModelDescriptor *lb);
//!\brief print moment matrix M to standard output
void DisplayLBModelDescriptorMomentMatrix(LBModelDescriptor *lb);
//!\brief compute max velocity of node grid
void ComputeLBModelDescriptorVmax(LBModelDescriptor *lb);
//!\brief check that macro quantities are conserved;
void CheckLBModelDescriptorMacroConservation(LBModelDescriptor *lb, bool verbose);
//******************************************************************************//
void LBM_Dummy_InitMacroData(schnaps_real x[3],schnaps_real w[]);
void LBM_Dummy_InitData(schnaps_real x[3],schnaps_real w[]);
void LBM_Dummy_InitData_OneNode(schnaps_real x[3],schnaps_real w[]);
//******************************************************************************//
//******************************************************************************//
typedef struct LBMSimulation{
  int d; //!velocity space dimension
  int nb_macro_fields; //! number of macro quantities
  int q; //! total number of velocity nodes 
  schnaps_real vmax;//!
  schnaps_real cfl;//!
  schnaps_real dt; //!
  schnaps_real tmax;//!
  schnaps_real itermax;//!
  schnaps_real tnow;
  schnaps_real iter_time;
  Model macro_model; // Model object to intiate macroscopic simulation
  Simulation macro_simu; //! a simulation object containing only macroscopic fields
  Model micro_model; //! array of nb_lattices Models object to initiate corresponding simulations.
  Simulation micro_simu; //! array of nb_lattices simulation object, one for each lattice.
  Model model_advec; //! model for advection of 1 velocity node, used in the implicit scheme;
  //
  void (*pre_advec)(void *lbsimu); //!
  void (*post_advec)(void *lbsimu); //! 
  void (*post_tstep) (void *lbsimu, schnaps_real *w);//! w to stay interface compatible with update_after_rk of simu
  //
  schnaps_real diag_2d_period;//! period for dumping of 2D diagnostics.
  void (*collect_diags) (void *lbsimu, schnaps_real *wmac,schnaps_real *wmic);//! routine for 1D time traces collection
  //
} LBMSimulation;
//*******************************************************************************//
void InitLBMSimulation( LBMSimulation *lbsimu,LatticeBoltzmannSimData *lsd, MacroMesh *mesh, int deg[3], int raf[3]);
void FreeBMSimulation(LBMSimulation *lbsimu);//
void LB_Relaxation_bgk_f( void *lbsimu);//
void LB_Relaxation_bgk_f_full( void *lbsimu);//
void LB_ComputeMacroFromMicro(void *lbsimu);//
// wrappers for compatibility with RK schemes operating on 1 simulation //
// in that case only the micsimu object is passed to rk. The global lbsim object, which is 
// necessary for pre/post advection steps is passed through a global pointer (in schnaps_lbm_simdata) to the current simu.
// the following wrappers get a pointer to the lbsimu object and call the corresponding routines (pre_advec,post_advec,post_tstep)
void LBM_pre_dtfields_wrapper( void * simu); // redirected to pre_advec;
void LBM_post_dtfields_wrapper( void * simu); // redirected to post_advec;
void LBM_update_after_rk_wrapper( void *simu,schnaps_real *w); // redirected to post_tstep;
//*******************************************************************************//
//***************** diagnostics routines ***************************************//

//******************************************************************************//
//******************************************************************************//
//******************************************************************************//
// flux functions for LB models
//******************************************************************************//
//******************************************************************************//
//! \brief flux for a lattice model; the lattice model index is specified by the global shared LatticeBoltzmannSimData structure
//! \param[in] wL,wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void LBM_OneLatticeNumFlux(schnaps_real *wL, schnaps_real *wR, schnaps_real *vnorm, schnaps_real *flux);
//! \brief flux for a lattice model and a specific velocity nodes model and node index specified by the global shared
//! LatticeBoltzmannSimData structure
//! \param[in] wL,wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void LBM_OneNodeNumFlux(schnaps_real *wL, schnaps_real *wR, schnaps_real *vnorm, schnaps_real *flux);
// time schemes
void LBMThetaTimeScheme(LBMSimulation *lbsimu,schnaps_real theta, schnaps_real tmax, schnaps_real dt);
// per model specific routines
//******************************************************************************//
//******************************************************************************//
//******************************************************************************//
//******************************************************************************//
//******************************************************************************//
schnaps_real LBM_dummy_zeros_feq(int inode,int nb_macro,schnaps_real *w);
//******* Hydrodynamic models (Euler, Navier-Stokes) ***************************//
// D2Q9 isothermal
void LBM_Set_D2Q9_ISOTH_model(LBModelDescriptor *lb,schnaps_real cref);
void LBM_f_to_macro_D2Q9_ISOTH(schnaps_real *f,schnaps_real *w);
schnaps_real LBM_feq_D2Q9_ISOTH(int inode, int nb_macro,schnaps_real *w);
// D2Q9 isothermal linearized (2D wave equation)
void LBM_Set_D2Q9_ISOTH_LINEARIZED_model(LBModelDescriptor *lb,schnaps_real cref);
void LBM_f_to_macro_D2Q9_ISOTH_LINEARIZED(schnaps_real *f,schnaps_real *w);
schnaps_real LBM_feq_D2Q9_ISOTH_LINEARIZED(int inode, int nb_macro,schnaps_real *w);
//
//***************************** MDH models **************************************//
#endif
