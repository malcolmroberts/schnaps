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
void VectorDgToCg (ContinuousSolver * ps,real * rhs);

#endif
