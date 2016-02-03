#ifndef _QUANTITIES_VP_H
#define _QUANTITIES_VP_H

#include "collision.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "interpolation.h"
#include "geometry.h"
#include "skyline.h"
#include "quantities_vp.h"
#include "simulation.h"

void Computation_charge_density(Simulation *simu);

//void Compute_electric_field(field * f, real * w);
void ComputeElectricField(field* f);
schnaps_real Computation_charge_average(Simulation *simu);

//void distribution_to_physic_entropy(field* f,real w,real *tw);

//void physic_entropy_to_distribution(field* f,real w,real *tw);

#endif
