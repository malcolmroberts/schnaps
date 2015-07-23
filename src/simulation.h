#ifndef _SIMULATION_H
#define _SIMULATION_H

#include "macromesh.h"
#include "field.h"

#ifdef _WITH_OPENCL
#include "clinfo.h"
#endif


//! \brief Data structure for managing a schnaps numerical simulation
typedef struct Simulation {
  //! Underlying mesh
  MacroMesh macromesh;
  //! List of fields for each macrocell
  field *fd;

  //! memory spaces for w and dtw
  real *w;
  real *dtw;

  //! sum of sizes of field data
  int wsize;

  //! Current time
  real tnow;
  //! CFL parameter min_i (vol_i / surf_i)
  real hmin;

  //! PIC struct pointer (=NULL if not used)
  void *pic;


  //! final time of simulation
  real tmax;

  //! time step and cfl
  real dt,cfl;

  //! current iteration of the RK algorithm
  int iter_time_rk;
  //! maximal number of iterations
  int itermax_rk;
  //! nb of diagnostics
  int nb_diags;
  //! table for diagnostics
  real *Diagnostics;

  //! vmax
  real vmax;


} Simulation;


//! \brief simulation initialization.
//! Computation of the initial data at each glop.
//! \param[inout] simu a simulation
//! \param[in] mesh a macromesh
//! \param[in] deg degrees parameters 
//! \param[in] raf refinements parameters 
//! \param[in] model a model
void InitSimulation(Simulation *simu, MacroMesh *mesh,
		    int *deg, int *raf, Model *model);


//! \brief apply the Discontinuous Galerkin approximation for computing
//! the time derivative of the fields. Works with several subcells.
//! \param[inout] simu A simulation
void DtFields(Simulation *simu, real *w, real *dtw);

//! \brief compute the time step of the RK scheme 
//! respecting a cfl condition
//! \param[inout] simu A simulation
real Get_Dt_RK(Simulation *simu);

//! \brief Time integration by a second order Runge-Kutta algorithm
//! \param[inout] simu a simulation
//! \param[in] tmax physical duration of the simulation
void RK2(Simulation *simu, real tmax);

//! \brief Time integration by a second order Runge-Kutta algorithm
//! \param[inout] simu a simulation
//! \param[in] tmax physical duration of the simulation
void RK4(Simulation *simu, real tmax);


// TODO: see how to manage opencl...
/* #ifdef _WITH_OPENCL */
/* //! \brief OpenCL version of RK2 */
/* //! time integration by a second order Runge-Kutta algorithm */
/* //! \param[inout] f a field */
/* //! \param[in] tmax physical duration of the simulation */
/* void RK2_CL(field *f, real tmax, real dt, */
/* 	    cl_uint nwait, cl_event *wait, cl_event *done); */
/* void RK4_CL(field *f, real tmax, real dt, */
/* 	    cl_uint nwait, cl_event *wait, cl_event *done); */
/* #endif */

//! \brief save the results in the gmsh format
//! \param[in] typplot index of the field variable to plot.
//! \param[in] compare if true, the numerical solution is compared
//! with the analytical solution
//! \param[in] simu a simulation
//! \param[in] fieldname name of the plotted data
//! \param[in] filename the path to the gmsh visualization file.
void PlotFields(int typplot, int compare, Simulation *simu, char *fieldname, 
	       char *filename);

//! \brief  list the valeus of the simulation
//! \param[in] simu a simulation.
void DisplaySimulation(Simulation *simu);

//! \brief Save 1D results in a text file
//! \param[in] simu a simulation.
//! \param[in] dir fixed direction to plot
//! \param[in] fixval fixed value to plot
//! \param[in] filename the path to the gmsh visualization file.
void Gnuplot(Simulation *simu,int dir, real fixval,char* filename);

//! \brief compute the normalized L2 distance with the imposed data
//! \param[in] simu a simulation.
//! \returns the error.
real L2error(Simulation *simu);


//! \brief compute the normalized L2 distance with the imposed data
//! \param[in] simu a simulation.
//! \param[in] nbvar index of one variable.
//! \returns the error.
real L2error_onefield(Simulation *simu, int nbvar);

#endif