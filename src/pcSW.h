#ifndef _SW_PC_H
#define _SW_PC_H

#include "linear_solver.h"
#include "simulation.h"
#include "solvercontinuous.h"
#include "physBased_PC.h"

typedef enum SW_OP{ID, ASF, ESF, AAF, EAF} SW_OP;

//! \brief Initialize the physics-based preconditioner with schur on the velocity boundary condition (u,n)=0
//! \param[in] simu: Simulation object containing some run-related variables
//! \param[inout] pb_pc: Physics-based Preconditioner object.
//! \param[in] list_mat2assemble: Integer array. Tells which matrices shall be assembled.
void Init_PBPC_SW_SchurVelocity_BCVelocity(Simulation *simu, PB_PC* pb_pc, int* list_mat2assemble);

//! \brief Initialize the physics-based preconditioner with schur on the velocity boundary condition p=g
//! \param[in] simu: Simulation object containing some run-related variables
//! \param[inout] pb_pc: Physics-based Preconditioner object.
//! \param[in] list_mat2assembly: Integer array. Tells which matrices shall be assembled.
void Init_PBPC_SW_SchurVelocity_BCPressure(Simulation *simu, PB_PC* pb_pc, int* list_mat2assemble);

//! \brief Initialize the physics-based preconditioner with schur on the pressure with boundary condition (u,n)=0
//! \param[in] simu: Simulation object containing some run-related variables
//! \param[inout] pb_pc: Physics-based Preconditioner object.
//! \param[in] list_mat2assembly: Integer array. Tells which matrices shall be assembled.
void Init_PBPC_SW_SchurPressure_BCVelocity(Simulation *simu, PB_PC* pb_pc, int* list_mat2assemble);

//! \brief Initialize the physics-based preconditioner with schur on the pressure with boundary condition p=g
//! \param[in] simu: Simulation object containing some run-related variables
//! \param[inout] pb_pc: Physics-based Preconditioner object.
//! \param[in] list_mat2assembly: Integer array. Tells which matrices shall be assembled.
void Init_PBPC_SW_SchurPressure_BCPressure(Simulation *simu, PB_PC* pb_pc, int* list_mat2assemble);

//! \brief Implementation of local matrices in the Slow Flow Approximation paradigm
//! \param[inout] pb_pc: Physics-Based PreConditioner object.
//! \param[in] var: Local values of the variables and their derivatives.
void Schur_ASF(void* pb_pc, real* var);


//! \brief Implementation of the full RHS for SW
//! \param[inout] pb_pc: Physics-Based PreConditioner object.
//! \param[in] var: Local values of the variables and their derivatives.
//! \param[in] coeff = dt*det*wpg*basisPhi_i[0]
//! \param[out] loc_rhs: an array for the local values of the right-hand side.
void SW_RHS(void* pb_pc, real* var, real coeff, real* loc_rhs);
#endif

