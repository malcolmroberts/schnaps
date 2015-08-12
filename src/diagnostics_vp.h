#ifndef _DIAGSNOTICS_VP_H
#define _DIAGSNOTICS_VP_H

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


//! \brief compute square of velocity L2 error
real L2VelError(field * f,real* x,real *w);

real L2_Kinetic_error(Simulation * simu);

real local_kinetic_energy(field * f,real* x,real *w);

void Energies(Simulation* simu,real * w,real k_energy, real e_energy, real t_energy,int first_diag);

void Charge_total(Simulation * simu,real * w, real t_charge,int first_diag);

void Plot_Energies(Simulation *simu, real dt);


#endif
