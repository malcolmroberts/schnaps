#ifndef _WAVE_PC_H
#define _WAVE_PC_H

#include "linear_solver.h"
#include "simulation.h"
#include "solvercontinuous.h"
#include "physBased_PC.h"

//! \brief Initialize the physics-based preconditioner with schur on the velocity boundary condition (u,n)=0
//! \param[in] simu: Simulation object containing some run-related variables
//! \param[inout] pb_pc: Physics-based Preconditioner object.
//! \param[in] list_mat2assembly: Integer array. Tells which matrices shall be assembled.
void Init_PBPC_Wave_SchurVelocity_BCVelocity(Simulation *simu, PB_PC* pb_pc, int* list_mat2assemble);

//! \brief Initialize the physics-based preconditioner with schur on the velocity boundary condition p=g
//! \param[in] simu: Simulation object containing some run-related variables
//! \param[inout] pb_pc: Physics-based Preconditioner object.
//! \param[in] list_mat2assembly: Integer array. Tells which matrices shall be assembled.
void Init_PBPC_Wave_SchurVelocity_BCPressure(Simulation *simu, PB_PC* pb_pc, int* list_mat2assemble);

//! \brief Initialize the physics-based preconditioner with schur on the pressure with boundary condition (u,n)=0
//! \param[in] simu: Simulation object containing some run-related variables
//! \param[inout] pb_pc: Physics-based Preconditioner object.
//! \param[in] list_mat2assembly: Integer array. Tells which matrices shall be assembled.
void Init_PBPC_Wave_SchurPressure_BCVelocity(Simulation *simu, PB_PC* pb_pc, int* list_mat2assemble);

//! \brief Initialize the physics-based preconditioner with schur on the pressure with boundary condition p=g
//! \param[in] simu: Simulation object containing some run-related variables
//! \param[inout] pb_pc: Physics-based Preconditioner object.
//! \param[in] list_mat2assembly: Integer array. Tells which matrices shall be assembled.
void Init_PBPC_Wave_SchurPressure_BCPressure(Simulation *simu, PB_PC* pb_pc, int* list_mat2assemble);

#endif
