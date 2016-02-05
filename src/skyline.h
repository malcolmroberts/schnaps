#ifndef _SKYLINE_H
#define _SKYLINE_H

#include "global.h"



//! \brief a struct for managing skyline linear system
typedef struct Skyline{

  //! \brief number of equations
  int neq;

  //! \brief size in memory
  int nmem;

  //! \brief array for the upper part of  the matrix
  schnaps_real* vkgs;

  //! \brief array for the diagonal part of  the matrix (size neq)
  schnaps_real* vkgd;

  //! \brief array for the lower part of  the matrix
  schnaps_real* vkgi;

  //! \brief array for the upper part of  the matrix (copy used if factorized)
  schnaps_real* copy_vkgs;

  //! \brief array for the diagonal part of  the matrix (size neq) (copy used if factorized)
  schnaps_real* copy_vkgd;

  //! \brief array for the lower part of  the matrix (copy used if factorized)
  schnaps_real* copy_vkgi;

  //! \brief profile of the matrix (size neq)
  int* prof;

  //! \brief array of indices of column start (size neq+1)
  int* kld;

  //! \brief true if the matrix is symmetric
  bool is_sym;

  //! \brief true if the struct is initialized
  bool is_init;

  //! \brief true if the arrays are allocated
  bool is_alloc;

  //! \brief true if the arrays of the copied matrix are allocated
  bool copy_is_alloc;

  //! \brief true if the matrix is factorized
  bool is_lu;

} Skyline;

//! \brief init the skyline structure with an empty matrix
//! \param[inout] sky the skyline object
//! \param[in] n number of equations
void InitSkyline(Skyline* sky,int n);

//! \brief free the allocated arrays
//! \param[inout] sky the skyline object
void FreeSkyline(Skyline* sky);


//! \brief indicates that elem (i,j) is nonzero
//! \param[inout] sky the skyline object
//! \param[in] i row index
//! \param[in] j column index
void SwitchOn(Skyline* sky,int i,int j); 

//! \brief allocate the variable-size arrays
//! \brief the nonzero positions should first be "switched on"
//! \param[inout] sky the skyline object
void AllocateSkyline(Skyline* sky);

//! \brief allocate the variable-size arrays
//! \brief the nonzero positions should first be "switched on"
//! \param[inout] sky the skyline object
void AllocateCopySkyline(Skyline* sky);

//! \brief set elem (i,j) to value val
//! \param[inout] sky the skyline object
//! \param[in] i row index
//! \param[in] j column index
//! \param[in] val value
void SetSkyline(Skyline* sky,int i,int j,schnaps_real val);

//! \brief add elem (i,j) to value val
//! \param[inout] sky the skyline object
//! \param[in] i row index
//! \param[in] j column index
//! \param[in] val value
void AddSkyline(Skyline* sky,int i,int j,schnaps_real val); 

//! \brief get elem (i,j)
//! \param[inout] sky the skyline object
//! \param[in] i row index
//! \param[in] j column index
//! \return value at pos (i,j)
schnaps_real GetSkyline(Skyline* sky,int i,int j); 


//! \brief display the matrix
//! \param[inout] sky the skyline object
void DisplaySkyline(Skyline* sky);

//! \brief compute a matrix vector product
//! \param[in] sky a skyline matrix
//! \param[in] x a vector
//! \param[out] prod Ax
void MatVectSkyline(Skyline * sky, schnaps_real * x, schnaps_real * prod);


//! \brief compute the inplace LU decomposition
//! \param[inout] sky the skyline object
void FactoLU(Skyline* sky);

//! \brief solve the linear system
//! \param[in] sky the skyline object
//! \param[in] rhs the right hand side
//! \param[in] sol the solution
void SolveSkyline(Skyline* sky,schnaps_real* rhs,schnaps_real* sol);

#endif
