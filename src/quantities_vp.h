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

void Computation_charge_density(field *f, real * w);

void Compute_electric_field(field* f, real * w);
void ComputeElectricField(field* f);
real Computation_charge_average(field *f,real * w);

void Entropy_transformation(field *f,real * w,void (*entropy_transform)(real f,real *ef));

void distribution_to_physic_entropy(real f,real *ef);

void physic_entropy_to_distribution(real f,real *ef);

#endif
