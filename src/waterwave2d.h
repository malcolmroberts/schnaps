#ifndef _WATERWAVE2D_H
#define _WATERWAVE2D_H

#include "model.h"


#define _SPEED_WAVE (1)

#define _LENGTH_DOMAIN (1.0)

#define _GRAVITY (1.0)


//! \brief boundardy flux based on the upwind scheme for wave
//! \param[in] t : current time
//! \param[in] x : current position
//! \param[in] wL : states
//! \param[in] vnorm : normal vector
//! \param[out] flux : the flux
void Wave_Upwind_BoundaryFlux(real *x, real t, real *wL, real *vnorm,real *flux);


//! \brief upwind flux for wave
//! \param[in] wL,wR : left and right states
//! \param[in] vnorm : normal vector
//! \param[out] flux : the flux
void Wave_Upwind_NumFlux(real wL[],real wR[],real* vnorm,real* flux);

//! \brief rusanov flux for wave
//! \param[in] wL,wR : left and right states
//! \param[in] vnorm : normal vector
//! \param[out] flux : the flux
void Wave_Rusanov_NumFlux(real wL[],real wR[],real* vnorm,real* flux);

//! \brief centered flux for wave
//! \param[in] wL,wR : left and right states
//! \param[in] vnorm : normal vector
//! \param[out] flux : the flux
void Wave_Centered_NumFlux(real wL[],real wR[],real* vnorm,real* flux);

//! \brief compute exact solution for x and t
//! \param[in] t : current time
//! \param[in] x : current position
//! \param[out] w : solution exact
void TestPeriodic_Wave_ImposedData(const real *x, const real t, real *w);


//! \brief init solution for x
//! \param[in] x : current position
//! \param[out] w : solution exact
void TestPeriodic_Wave_InitData(real *x, real *w);


//! \brief roe flux for shallows water
//! \param[in] wL,wR : left and right states
//! \param[in] vnorm : normal vector
//! \param[out] flux : the flux
void ShallowWater_Roe_NumFlux(real wL[],real wR[],real* vnorm,real* flux);


//! \brief HLL flux for shallow water
//! \param[in] wL,wR : left and right states
//! \param[in] vnorm : normal vector
//! \param[out] flux : the flux
void ShallowWater_HLL_NumFlux(real wL[],real wR[],real* vnorm,real* flux);

//! \brief Rusanov flux for Shallow water
//! \param[in] wL,wR : left and right states
//! \param[in] vnorm : normal vector
//! \param[out] flux : the flux
void ShallowWater_Rusanov_NumFlux(real wL[],real wR[],real* vnorm,real* flux);


 //! \brief a pointer to the source function
 //! \param[in] x : space position
 //! \param[in] t : time
  //! \param[in] w :  state
  //! \param[out] source : the source
void ShallowWater_classical_SourceTerm(const real *x, const real t, const real *w, real *source);


 //! \brief a pointer to the source function
 //! \param[in] x : space position
 //! \param[in] t : time
  //! \param[in] w :  state
  //! \param[out] source : the source
void ShallowWater_periodic_SourceTerm(const real *x, const real t, const real *w, real *source);

 //! \brief a pointer to the source function for the HLL scheme
 //! \param[in] x : space position
 //! \param[in] t : time
  //! \param[in] w :  state
  //! \param[out] source : the source
void ShallowWater_HLLWB_SourceTerm(const real *x, const real t, const real *w, real *source);


//! \brief init solution for x
//! \param[in] x : current position
//! \param[in] t : current time
//! \param[out] w : solution exact
void TestSH_equilibrium_ImposedData(const real *x, const real t, real *w);


//! \brief init solution for x
//! \param[in] x : current position
//! \param[out] w : solution exact
void TestSH_equilibrium_InitData(real *x, real *w);

//! \brief init solution for x
//! \param[in] x : current position
//! \param[in] t : current time
//! \param[out] w : solution exact
void TestSH_periodic_ImposedData(const real *x, const real t, real *w);


//! \brief init solution for x
//! \param[in] x : current position
//! \param[out] w : solution exact
void TestSH_periodic_InitData(real *x, real *w);

#endif
