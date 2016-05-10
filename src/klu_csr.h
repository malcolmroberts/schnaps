#ifndef _KLU_CSR_H
#define _KLU_CSR_H

#include "global.h"
#include "cxsparse/cs.h"
#include "klu.h"


//! \brief a struct for managing KLU linear system
typedef struct KLU {

  //! \brief number of equations
  int neq;

  //! \brief true if the struct is initialized
  bool is_init;

  //! \brief true if the arrays are allocated
  bool is_alloc;

  //! \brief true if the arrays of the copied matrix are allocated
  bool copy_is_alloc;

  //! \brief true if the matrix is factorized
  bool is_lu;

  //! \brief symbolic and numeric objects
  klu_symbolic *symbolic;
  klu_numeric *numeric;
  klu_common common;

  //! csparse struct for triplet storage
  cs_di *T;

  //! csparse struct for csr storage
  cs_di *A;

  //! csparse struct for csr copy storage
  cs_di *Acopy;


} KLU;

//! \brief init the KLU structure with an empty matrix
//! \param[inout] sky the KLU object
//! \param[in] n number of equations
void InitKLU(KLU * sky, int n);

//! \brief free the allocated arrays
//! \param[inout] sky the KLU object
void FreeKLU(KLU * sky);


//! \brief indicates that elem (i,j) is nonzero
//! \param[inout] sky the KLU object
//! \param[in] i row index
//! \param[in] j column index
void SwitchOnKLU(KLU * sky, int i, int j);

//! \brief allocate the variable-size arrays
//! \brief the nonzero positions should first be "switched on"
//! \param[inout] sky the KLU object
void AllocateKLU(KLU * sky);

//! \brief allocate the variable-size arrays
//! \brief the nonzero positions should first be "switched on"
//! \param[inout] sky the KLU object
void AllocateCopyKLU(KLU * sky);

//! \brief set elem (i,j) to value val
//! \param[inout] sky the KLU object
//! \param[in] i row index
//! \param[in] j column index
//! \param[in] val value
void SetKLU(KLU * sky, int i, int j, schnaps_real val);

//! \brief add elem (i,j) to value val
//! \param[inout] sky the KLU object
//! \param[in] i row index
//! \param[in] j column index
//! \param[in] val value
void AddKLU(KLU * sky, int i, int j, schnaps_real val);

//! \brief get elem (i,j)
//! \param[inout] sky the KLU object
//! \param[in] i row index
//! \param[in] j column index
//! \return value at pos (i,j)
schnaps_real GetKLU(KLU * sky, int i, int j);


//! \brief display the matrix
//! \param[inout] sky the KLU object
void DisplayKLU(KLU * sky);

//! \brief compute a matrix vector product
//! \param[in] sky a KLU matrix
//! \param[in] x a vector
//! \param[out] prod Ax
void MatVectKLU(KLU * sky, schnaps_real * x, schnaps_real * prod);


//! \brief compute the inplace LU decomposition
//! \param[inout] sky the KLU object
void FactoKLU(KLU * sky);

//! \brief recompute the inplace LU decomposition
//! \param[inout] sky the KLU object
void ReFactoKLU(KLU * sky);

//! \brief solve the linear system
//! \param[in] sky the KLU object
//! \param[in] rhs the right hand side
//! \param[in] sol the solution
void SolveKLU(KLU * sky, schnaps_real * rhs, schnaps_real * sol);

#endif
