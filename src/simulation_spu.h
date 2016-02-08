#ifndef _SIMULATION_SPU_H
#define _SIMULATION_SPU_H

#include "simulation.h"
#include <starpu.h>

//! \brief Display a starpu data handle on runtime
//! The function calls starpu_task_wait_for_all before display.
//! \param[in] handle a starpu handle
//! \param[in] name a starpu handle name (appears on every line: make it short)
void DisplayHandle_SPU(starpu_data_handle_t handle,
                       const char* name);


//! \brief Init and get ZeroBuffer codelet.
struct starpu_codelet* ZeroBuffer_codelet();

//! \brief Create a ZeroBuffer task and submit it.
//! \param[in,out] handle handle of the buffer to nullifiy
void ZeroBuffer_SPU(starpu_data_handle_t handle);


//! \brief Init and get AddBuffer codelet.
struct starpu_codelet* AddBuffer_codelet();

//! \brief Add a scaled buffer to another buffer.
//! handle_out = handle_out + alpha * handle_in
//! \param[in] alpha the scaling factor
//! \param[in] handle_in handle to the scaled buffer
//! \param[in,out] handle_out handle to the result buffer
void AddBuffer_SPU(schnaps_real alpha,
                   starpu_data_handle_t handle_in,
                   starpu_data_handle_t handle_out);


//! \brief Init and get DGVolume codelet.
struct starpu_codelet* DGVolume_codelet();

//! \brief Apply the "cross" derivative terms inside a macrocell
//! \param[in,out] fd field
void DGVolume_SPU(field* fd);


//! \brief apply the Discontinuous Galerkin approximation for computing
//! the time derivative of the fields. Works with several subcells.
//! starpu version
//! \param[inout] simu a simulation
//! \param[inout] w a starpu handle to the field value
//! \param[out] dtw a starpu handle to the time derivatives
void DtFields_SPU(Simulation *simu,
		  starpu_data_handle_t* w_handle,
		  starpu_data_handle_t* dtw_handle);

//! \brief RK2 integration of the DG approximation
//! starpu version
//! \param[inout] simu a simulation
//! \param[in] tmax tmax
void RK2_SPU(Simulation *simu, schnaps_real tmax);

//! \brief apply the flux terms inside a macrocell
//! StarPU version
//! \param[inout] fd a field
void DGSubCellInterface_SPU(field* fd);

//! \brief apply the source terms inside a macrocell
//! StarPU version
//! \param[inout] fd a field
void DGSource_SPU(field* fd);

//! \brief apply the inverse of the mass matrix in a macrocell
//! StarPU version
//! \param[inout] fd a field
void DGMass_SPU(field* fd);


//! \brief  apply the interface fluxes to a neighbouring field
//! \param[in] inter an Interface
//! \param[in] side the side: left if == 0 right if ==1
void InterfaceExplicitFlux_bis(Interface* inter, int side);

//! \brief  apply the interface fluxes to a neighbouring field
//! \param[in] inter an Interface
//! \param[in] side the side: left if == 0 right if ==1
void DGMacroCellInterface_SPU(Interface* inter, int side);

//! \brief  apply the interface fluxes to neighbouring fields
//! \param[in] inter an Interface
void DGMacroCellInterface_bis_SPU(Interface* inter);

//! \brief  apply the boundary flux to a neighbouring field
//! \param[in] inter an Interface
void DGMacroCellBoundaryFlux_SPU(Interface* inter);

//! \brief apply the Discontinuous Galerkin approximation for computing
//! the time derivative of the fields. Works with several subcells.
//! test version before starpu implementation
//! \param[in] simu a simulation
//! \param[inout] w an array to the field value
//! \param[out] dtw an array to the time derivatives
void DtFields_bis(Simulation *simu,
		  schnaps_real* w,
		  schnaps_real* dtw);


#endif