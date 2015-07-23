#ifndef _FIELD_H
#define _FIELD_H

#include "macromesh.h"
#include "interpolation.h"
#include "model.h"

#ifdef _WITH_OPENCL
#include "clinfo.h"
#endif


//! \brief Data structure for managing a  discrete vector field
//! solution of a DG approximation
typedef struct field {
  //! Physical nodes of the macrocell
  real physnode[20][3];
  
  //! Physical and numerical model
  Model model;
  //! Interpolation used for each component of the field
  Interpolation interp;

  //! Refinement of the macrocell in each direction
  int raf[3];

  //! Degrees of interpolation in each direction
  int deg[3];

  //! number og Gauss points in each direction
  int npg[3];
  
  //! Current time and time steps
  real tnow;
  real dt;

  //! ref length of the mesh subcells
  real hmin;

  //! period in each direction
  //! if negative: non-periodic computation (default)
  real period[3];


  //! PIC struct pointer (=NULL if not used)
  void *pic;

  //! Size of the field buffers
  int wsize;
  //! fields at current time step
  real *wn;
  //! Time derivative of the field
  real *dtwn;


  

  //! \brief Pointer to a generic function called before computing dtfield. 
  //! \param[inout] f a field (to be converted from void*)
  void (*pre_dtfield)(void *f, real *w);

  //! \brief Pointer to a generic function called after computing dtfield. 
  //! \param[inout] f a field (to be converted from void*)
  void (*post_dtfield)(void *f, real *w);

  //! \brief generic update function called 
  //! \brief called at each runge-kutta sustep
  //! \param[inout] f a field (to be converted from void*)
  //! \param[in] elem macro element index
  //! \param[in] ipg glop index
  //! \param[in] iv field component index
  void (*update_after_rk)(void *f, real *w);

  //! \brief Memory arrangement of field components
  //! \param[in] param interpolation parameters
  //! \param[in] elem macro element index
  //! \param[in] ipg glop index
  //! \param[in] iv field component index
  int (*varindex)(int* deg, int *ref, int m, int ipg, int iv);


#ifdef _WITH_OPENCL
  //! \brief opencl data
  CLInfo cli;
  //! \brief copy of the dtwn array
  cl_mem wn_cl;
  cl_mem dtwn_cl;
  //! \brief copy of the params
  cl_mem param_cl;
  //! \brief copy physnode
  cl_mem physnode_cl;
  cl_mem physnodes_cl; // The physnodes for all the macrocells
  real *physnode;

  cl_mem physnodeR_cl;
  real *physnodeR;

  bool use_source_cl;
  char *sourcename_cl;

  //! opencl kernels
  cl_kernel dgmass;
  cl_kernel dgflux;
  cl_kernel dgvolume;
  cl_kernel dgsource;
  cl_kernel dginterface;
  cl_kernel dgboundary;
  cl_kernel RK_out_CL;
  cl_kernel RK_in_CL;
  cl_kernel RK4_final_stage;
  cl_kernel zero_buf;

  // OpenCL events

  // set_buf_to_zero event
  cl_event clv_zbuf; 
  
  // Subcell mass events
  cl_event *clv_mass; 

  // Subcell flux events
  cl_event *clv_flux0, *clv_flux1, *clv_flux2;

  // Subcell volume events
  cl_event *clv_volume; 

  // Subcell volume events
  cl_event *clv_source; 

  // Macrocell interface events
  cl_event *clv_mci;
  // Boundary term events
  cl_event *clv_boundary;

  // OpenCL timing
  cl_ulong zbuf_time;
  cl_ulong mass_time;
  cl_ulong vol_time;
  cl_ulong flux_time;
  cl_ulong minter_time;
  cl_ulong boundary_time;
  cl_ulong source_time;
  cl_ulong rk_time;

  // OpenCL roofline measurements
  unsigned long int flops_vol, flops_flux, flops_mass; 
  unsigned long int reads_vol, reads_flux, reads_mass; 
#endif
} field;

//! \brief memory arrangement of field components.
//! Generic implementation.
//! \param[in] deg degrees parameters
//! \param[in] raf refinement parameters
//! \param[in] m number of conservative variables
//! \param[in] ipg glop index
//! \param[in] iv field component index
//! \returns the memory position in the arrays wn wnp1 or dtwn.
#pragma start_opencl
int GenericVarindex(__constant int *deg, __constant int *raf, int m,
		    int ipg, int iv);
#pragma end_opencl


//! \brief field initialization. Computation of the initial at each glop.
//! \param[inout] f a field
//! \param[in] m a model
//! \param[in] physnode list of geometrical nodes of the macroelement
//! \param[in] deg degrees parameters 
//! \param[in] raf refinements parameters 
//! \param[in] w a pointer to field value (if NULL memory will be allocated)
//! \param[inout] dtw a pointer to derivatives (if NULL memory will be allocated)
void Initfield(field *f, Model m, real physnode[][3], int *deg, int *raf, real *w, real* dtw);

void init_empty_field(field *f);

//! free the buffers created in Initfield.
//! \param[inout] f a field
void Freefield(field *f);

//! \brief  compute the Discontinuous Galerkin inter-macrocells boundary terms second implementation with a loop on the faces
//! \param[in] locfaL local index of the face in the left element
//! \param[inout] fL a left field
//! \param[in] offsetL left offset for accessing data in w and dtw
//! \param[inout] fR a right field
//! \param[in] offsetR right offset for accessing data in w and dtw
//! \param[in] w field data 
//! \param[out] dtw time derivative of the field data w
void DGMacroCellInterface(int locfaL,
			  field *fL, int offsetL, field *fR, int offsetR,
			  real *w, real *dtw);

//! \brief compute the Discontinuous Galerkin volume terms
//! \param[in] f a field
void DGVolume(field *f, real *w, real *dtw);

//! \brief compute the Discontinuous Galerkin inter-subcells terms
//! \param[in] f a field
void DGSubCellInterface(field *f, real *w, real *dtw);

//! \brief  apply the DG mass term
//! \param[in] f a field
void DGMass(field *f, real *w, real *dtw);

//! \brief Add the source term
//! \param[in] f a field
void DGSource(field *f, real *w, real *dtw);

//! \brief An out-of-place RK stage
//! \param[out] fwnp1 field at time n+1
//! \param[in] fwn field at time n
//! \param[in] fdtwn time derivative of the field
//! \param[in] dt time step
//! \param[in] sizew size of the field buffer
void RK_out(real *fwnp1, real *fwn, real *fdtwn, const real dt, 
	    const int sizew);

//! \brief An in-place RK stage
//! \param[inout] fwnp1 field at time n+1
//! \param[in] fdtwn time derivative of the field
//! \param[in] dt time step
//! \param[in] sizew size of the field buffer
void RK_in(real *fwnp1, real *fdtwn, const real dt, const int sizew);



#ifdef _WITH_OPENCL
//! \brief OpenCL version of RK2
//! time integration by a second order Runge-Kutta algorithm
//! \param[inout] f a field
//! \param[in] tmax physical duration of the simulation
void RK2_CL(field *f, real tmax, real dt,
	    cl_uint nwait, cl_event *wait, cl_event *done);
void RK4_CL(field *f, real tmax, real dt,
	    cl_uint nwait, cl_event *wait, cl_event *done);
#endif

/* //! \brief save the results in the gmsh format */
/* //! \param[in] typplot index of the field variable to plot. */
/* //! \param[in] compare if true, the numerical solution is compared */
/* //! with the analytical solution */
/* //! \param[in] f a field */
/* //! \param[in] fieldname name of the plotted data */
/* //! \param[in] filename the path to the gmsh visualization file. */
/* void Plotfield(int typplot, int compare, field *f, char *fieldname,  */
/* 	       char *filename); */

//! \brief interpolate field at a reference point a macrocell
//! \param[in] f a field
//! \param[in] ie the macrocell index
//! \param[in] xref reference coordinates
//! \param[out] w the m field values
void InterpField(field *f,int ie,real* xref,real* w);

//! \brief  display the field on screen
//! \param[in] f the field.
void Displayfield(field *f);




#endif
