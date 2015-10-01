#ifndef _INTERFACE_H
#define _INTERFACE_H

#include "field.h"


typedef struct Interface{

  //! \brief Left field
  field *fL;

  //! \brief Right field
  field *fR;

  //! \brief local interface index in Left field 
  int locfaL;

  //! \brief local interface index in Left field 
  int locfaR;
  
  //! \brief number of left glops
  int npgL;
  
  //! \brief number of right glops
  int npgR;

  //! \brief Left glops volume indices from the interface glop indices
  int* vol_indexL;
  starpu_data_handle_t vol_indexL_handle;
  
  //! \brief Right glops volume indices from the interface glop indices
  int * vol_indexR;
  starpu_data_handle_t vol_indexR_handle;

  //! \brief Size in memory of the stored conservative data at Left interface glops
  int wsizeL;

  //! \brief Size in memory of the stored conservative data at Left interface glops
  int wsizeR;

  //! \brief Left conservative glops data
  real *wL;
  starpu_data_handle_t wL_handle;
  
  //! \brief Right conservative glops data
  real *wR;
  starpu_data_handle_t wR_handle;
  
  //! \brief weighted normal vectors at interface glops 
  real *vnds;
  starpu_data_handle_t vnds_handle;
  
  //! \brief  glops coordinates 
  real *xpg;
  starpu_data_handle_t xpg_handle;


} Interface;

//! \brief  registration of starpu data for an interface
//! \param[inout] inter an Interface
void InitInterface_SPU(Interface* inter);


//! \brief  extract the values of the neighbouring fields to the interface
//! \param[inout] inter an Interface
//! \param[in] side the side: left if == 0 right if ==1
void ExtractInterface(Interface* inter, int side);

//! \brief  extract the values of the neighbouring fields to the interface
//! \brief StarPU version
//! \param[inout] inter an interface
//! \param[in] side the side: left if == 0 right if ==1
void ExtractInterface_SPU(Interface* inter, int side);


//! \brief  apply the interface fluxes to a neighbouring field
//! \param[in] inter an Interface
//! \param[in] side the side: left if == 0 right if ==1
void InterfaceExplicitFlux(Interface* inter, int side);


//! \brief  apply the interface fluxes to a neighbouring field
//! \brief StarPU version
//! \param[in] inter an Interface
//! \param[in] side the side: left if == 0 right if ==1
void InterfaceExplicitFlux_SPU(Interface* inter, int side);

//! \brief  apply the boundary fluxes to the Left field
//! \brief StarPU version
//! \param[in] inter an Interface
void InterfaceBoundaryFlux_SPU(Interface* inter);

//! \brief  assembly of the inter-fields matrix
//! \param[in] inter an Interface
//! \param[in] side the side: left if == 0 right if ==1
//! \param[in] theta Cranck-Nicolson parameter
//! \param[in] dt time step
void InterfaceLocalAssembly(Interface *inter,  real theta, real dt);

//! \brief  varindex function of the interface
//! \param[in] npg number of interface Gauss points
//! \param[in] m number of conservative variables
//! \param[in] ipgf Gauss point index
//! \param[in] iv conservative variable index
//! \returns the memory position of the variable
int VarindexFace(int npg, int m, int ipgf, int iv); 

#endif
