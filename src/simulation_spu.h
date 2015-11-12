#ifndef _SIMULATION_SPU_H
#define _SIMULATION_SPU_H

#include "simulation.h"
#include <starpu.h>

//! \brief apply the Discontinuous Galerkin approximation for computing
//! the time derivative of the fields. Works with several subcells.
//! starpu version
//! \param[in] simu a simulation
//! \param[inout] w a starpu handle to the field value
//! \param[out] dtw a starpu handle to the time derivatives
void DtFields_SPU(Simulation *simu,
		  starpu_data_handle_t* w_handle,
		  starpu_data_handle_t* dtw_handle);


void ZeroBuffer_SPU(starpu_data_handle_t w);


#endif
