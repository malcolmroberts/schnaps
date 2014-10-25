#ifndef _INTERPOLATION_H
#define _INTERPOLATION_H

#include "geometry.h"

// utilitie functions for computing interpolation
// on a macrocell
//! \brief a struct for managing geometric mapping
typedef struct Interpolation{
  //! \brief interpolation parameters
  //! generally the convention is
  //! param[0] = M number of variables
  //! param[1] = deg x
  //! param[2] = deg y
  //! param[3] = deg z
  //! param[4] = raf x
  //! param[5] = raf y
  //! param[6] = raf z
  //! param[7..] = others param or return from interp
  //! functions
  int interp_param[8];

  //! \brief underlying geometry mapping
  Geom geo;
  //! \brief number of Gauss points in the volume
  int npgv;
  //! \brief current Gauss point index
  int ipgv;
  //! \brief current face id
  int ifa;
  //! \brief number of Gauss points on the current face
  int npgf;
  //! \brief current face Gauss point index
  int ipgf;
  //! \brief basis function index
  int ib;
  //! \brief current volume Gauss weight
  double wpgv;
  //! \brief  current face Gauss weight
  double wpgf;
  //! \brief Gauss point ref location
  double xpgref[3];
  //! \brief Gauss point physical location
  double xpg[3];
  //! \brief basis function values
  double phi;
  //! \brief basis function reference gradient
  double dphiref[3];
  //! \brief basis function physical gradient
  double dphi[3];

} Interpolation;

// in all these functions
// param is an input integer array
// param[0..2] : approximation degree in the 3 directions on 
// the reference cube
// param[3..5] : number of subcells in the 3 directions on 
// the reference cube
// param[6..nparam]: can be used for returning
// parameters (the user has to ensure that enough memory is reserved

//! \brief number of Gauss-LObatto Points (GLOPs) on the macro cell 
//! \param[in] param: the param list 
int NPG(int param[]);
// number of GLOPs on the face ifa of the macrocell
int NPGF(int* param,int ifa);

// return the reference coordinates xpg[3] and weight wpg of the GLOP ipg
void ref_pg_vol(int* param,int ipg,
		double* xpg,double* wpg,double* xpg_in);

// from a reference point find the nearest
// gauss point
int ref_ipg(int* param,double* xref);

// same function for the face 
// param[6] contains the volume GLOP index computed from face GLOP index.
//! \brief compute the position xpg of glop ipg in the local
//! numbering on face ifa. If xpgin is not NULL also compute
//! the position of point slightly inside the opposite subcell
void ref_pg_face(int* param,int ifa,int ipg,double* xpg,double* wpg,
		 double* xpgin);
// return the value psi  and the gradient dpsi[3] of the basis 
// function ib at point xref[3]. Warning: the value of the gradient is
// not reliable if xref is on the boundary 
//of a subcell (because the gradient is discontinuous)
void psi_ref(int* param, int ib, double* xref, double* psi, double* dpsiref);
// return the gradient dpsi[3] of the basis 
// function ib at GLOP ipg.
void grad_psi_pg(int* param,int ib,int ipg,double* dpsiref);

// same function but with specification of the subcell
// indices is[3] in the three directions
void psi_ref_subcell(int* param, int* is,int ib, double* xref, double* psi, double* dpsiref);


// return the 1d ith GLOP weight for degree deg
double wglop(int deg,int i);
// return the 1d derivative of lagrange polynomial ib at glop ipg
double dlag(int deg,int ib,int ipg);




#endif
