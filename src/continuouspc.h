#ifndef _CONTINUOUSPC_H
#define _CONTINUOUSPC_H


#include "collision.h"
#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "geometry.h"
#include "interpolation.h"
#include "linear_solver.h"
#include "simulation.h"
#include "solvercontinuous.h"


// ! TOODODODODODOODODODOD
void VectorDgToCg (ContinuousSolver * ps,real * rhs, real * rhs_out);

// ! TOODODODODODOODODODOD
void VectorCgToDg(ContinuousSolver * cs, real * rhsIn, real * rhsOut);

// ! TOODODODODODOODODODOD
void physicPC_wave(Simulation *simu, real* globalSol, real* globalRHS);

// ! TOODODODODODOODODODOD
void test(void *cs);
#endif
