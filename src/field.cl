// Return the 1d derivative of lagrange polynomial ib at glop ipg

#ifndef _PERIODX
#define _PERIODX -1
#endif
#ifndef _PERIODY
#define _PERIODY -1
#endif
#ifndef _PERIODZ
#define _PERIODZ -1
#endif


// A kernel used solely for managing events (cf: SOCL)
__kernel
void empty_kernel()
{
}

real dlag(int deg, int ib, int ipg)
{
  return gauss_lob_dpsi[gauss_lob_dpsi_offset[deg] + ib * (deg + 1) + ipg];
}

#ifndef VARINDEX
#define VARINDEX GenericVarindex
#endif

int ref_ipg(__constant int *param, real *xref);

void compute_gradphi(const real xref[3], real gradphi[20][4]) 
{
  real x = xref[0];
  real y = xref[1];
  real z = xref[2];

  real t1 = -1 + z;
  real t2 = -1 + y;
  real t3 = t1 * t2;
  real t4 = 2 * y;
  real t5 = 2 * z;
  real t6 = 4 * x;
  real t9 = -1 + x;
  real t10 = t1 * t9;
  real t11 = 2 * x;
  real t12 = 4 * y;
  real t15 = t2 * t9;
  real t16 = 4 * z;
  real t24 = x * t1;
  real t27 = x * t2;
  real t33 = y * t1;
  real t38 = x * y;
  real t48 = y * t9;
  real t54 = z * t2;
  real t57 = z * t9;
  real t67 = x * z;
  real t75 = y * z;
  real t94 = t11 - 1;
  real t98 = 4 * t24 * t9;
  real t100 = 4 * t27 * t9;
  real t104 = 4 * t33 * t2;
  real t105 = t4 - 1;
  real t111 = 4 * y * t2 * t9;
  real t114 = z * t1;
  real t116 = 4 * t114 * t2;
  real t118 = 4 * t114 * t9;
  real t119 = t5 - 1;
  real t128 = 4 * t38 * t2;
  real t132 = 4 * t67 * t1;
  real t141 = 4 * t38 * t9;
  real t145 = 4 * t75 * t1;
  real t158 = 4 * t67 * t9;
  real t162 = 4 * t75 * t2;

  gradphi[0][0] = t3 * (t4 + t5 - 3 + t6);
  gradphi[0][1] = t10 * (t11 + t5 - 3 + t12);
  gradphi[0][2] = t15 * (t11 + t4 - 3 + t16);
  gradphi[0][3] = t3 * t9 * (t11 + t4 + t5 - 1);
  gradphi[1][0] = t3 * (-t4 - t5 - 1 + t6);
  gradphi[1][1] = t24 * (-t5 + 1 + t11 - t12);
  gradphi[1][2] = t27 * (-t4 + 1 + t11 - t16);
  gradphi[1][3] = t24 * t2 * (-t4 - t5 + t11 - 1);
  gradphi[2][0] = -t33 * (t4 - t5 - 3 + t6);
  gradphi[2][1] = -t24 * (-t5 - 3 + t11 + t12);
  gradphi[2][2] = -t38 * (t4 - 1 + t11 - t16);
  gradphi[2][3] = -t38 * t1 * (t4 - t5 - 3 + t11);
  gradphi[3][0] = -t33 * (-t4 + t5 - 1 + t6);
  gradphi[3][1] = -t10 * (t11 + t5 + 1 - t12);
  gradphi[3][2] = -t48 * (t11 - 1 - t4 + t16);
  gradphi[3][3] = -t33 * t9 * (t11 + t5 - t4 + 1);
  gradphi[4][0] = -t54 * (t4 - t5 - 1 + t6);
  gradphi[4][1] = -t57 * (t11 - t5 - 1 + t12);
  gradphi[4][2] = -t15 * (t11 + 1 + t4 - t16);
  gradphi[4][3] = -t54 * t9 * (t11 + t4 - t5 + 1);
  gradphi[5][0] = -t54 * (-t4 + t5 - 3 + t6);
  gradphi[5][1] = -t67 * (t5 - 1 + t11 - t12);
  gradphi[5][2] = -t27 * (-t4 - 3 + t11 + t16);
  gradphi[5][3] = -t67 * t2 * (-t4 + t5 - 3 + t11);
  gradphi[6][0] = t75 * (t4 + t5 - 5 + t6);
  gradphi[6][1] = t67 * (t5 - 5 + t11 + t12);
  gradphi[6][2] = t38 * (t4 - 5 + t11 + t16);
  gradphi[6][3] = t38 * z * (t4 + t5 - 5 + t11);
  gradphi[7][0] = t75 * (-t4 - t5 + 1 + t6);
  gradphi[7][1] = t57 * (t11 - t5 + 3 - t12);
  gradphi[7][2] = t48 * (t11 - t4 + 3 - t16);
  gradphi[7][3] = t75 * t9 * (t11 - t5 - t4 + 3);
  gradphi[8][0] = -4 * t3 * t94;
  gradphi[8][1] = -t98;
  gradphi[8][2] = -t100;
  gradphi[8][3] = -4 * t24 * t15;
  gradphi[9][0] = -t104;
  gradphi[9][1] = -4 * t1 * t105 * t9;
  gradphi[9][2] = -t111;
  gradphi[9][3] = -4 * t33 * t15;
  gradphi[10][0] = -t116;
  gradphi[10][1] = -t118;
  gradphi[10][2] = -4 * t119 * t2 * t9;
  gradphi[10][3] = -4 * t114 * t15;
  gradphi[11][0] = t104;
  gradphi[11][1] = 4 * t24 * t105;
  gradphi[11][2] = t128;
  gradphi[11][3] = 4 * t38 * t3;
  gradphi[12][0] = t116;
  gradphi[12][1] = t132;
  gradphi[12][2] = 4 * x * t119 * t2;
  gradphi[12][3] = 4 * t67 * t3;
  gradphi[13][0] = 4 * t33 * t94;
  gradphi[13][1] = t98;
  gradphi[13][2] = t141;
  gradphi[13][3] = 4 * t38 * t10;
  gradphi[14][0] = -t145;
  gradphi[14][1] = -t132;
  gradphi[14][2] = -4 * t38 * t119;
  gradphi[14][3] = -4 * t38 * t114;
  gradphi[15][0] = t145;
  gradphi[15][1] = t118;
  gradphi[15][2] = 4 * y * t119 * t9;
  gradphi[15][3] = 4 * t75 * t10;
  gradphi[16][0] = 4 * t54 * t94;
  gradphi[16][1] = t158;
  gradphi[16][2] = t100;
  gradphi[16][3] = 4 * t67 * t15;
  gradphi[17][0] = t162;
  gradphi[17][1] = 4 * z * t105 * t9;
  gradphi[17][2] = t111;
  gradphi[17][3] = 4 * t75 * t15;
  gradphi[18][0] = -t162;
  gradphi[18][1] = -4 * t67 * t105;
  gradphi[18][2] = -t128;
  gradphi[18][3] = -4 * t38 * t54;
  gradphi[19][0] = -4 * t75 * t94;
  gradphi[19][1] = -t158;
  gradphi[19][2] = -t141;
  gradphi[19][3] = -4 * t38 * t57;
}

void compute_xphy(__constant real *physnode,
		  real gradphi[20][4],
		  real xphy[3])
{
  for(int ii = 0; ii < 3; ++ii) {
    xphy[ii] = 0;
    for(int i = 0; i < 20; ++i) {
      xphy[ii] += physnode[3 * i + ii] * gradphi[i][3];
    }
  }
}

void compute_dtau(__constant real *physnode,
		  real gradphi[20][4],
		  real dtau[3][3])
{
  for(int ii = 0; ii < 3; ii++) {
    for(int jj = 0; jj < 3; jj++) {
      dtau[ii][jj] = 0;
    }
    for(int i = 0; i < 20; i++) {
      for(int jj = 0; jj < 3; jj++) {
	dtau[ii][jj] += physnode[3 * i + ii] * gradphi[i][jj];;
      }
    }
  }
}

void compute_codtau(real dtau[3][3], real codtau[3][3])
{
  codtau[0][0] =  dtau[1][1] * dtau[2][2] - dtau[1][2] * dtau[2][1];
  codtau[0][1] = -dtau[1][0] * dtau[2][2] + dtau[1][2] * dtau[2][0];
  codtau[0][2] =  dtau[1][0] * dtau[2][1] - dtau[1][1] * dtau[2][0];
  codtau[1][0] = -dtau[0][1] * dtau[2][2] + dtau[0][2] * dtau[2][1];
  codtau[1][1] =  dtau[0][0] * dtau[2][2] - dtau[0][2] * dtau[2][0];
  codtau[1][2] = -dtau[0][0] * dtau[2][1] + dtau[0][1] * dtau[2][0];
  codtau[2][0] =  dtau[0][1] * dtau[1][2] - dtau[0][2] * dtau[1][1];
  codtau[2][1] = -dtau[0][0] * dtau[1][2] + dtau[0][2] * dtau[1][0];
  codtau[2][2] =  dtau[0][0] * dtau[1][1] - dtau[0][1] * dtau[1][0];
}

void compute_dphi(real dphiref[3], real codtau[3][3], real dphi[3])
{
  for(int ii = 0; ii < 3; ii++) {
    dphi[ii]=0;
    for(int jj = 0; jj < 3; jj++) {
      dphi[ii] += codtau[ii][jj] * dphiref[jj];
    }
  }
}

void ComputeNormal(real codtau[3][3], int ifa, real vnds[3])
{
  int h20_refnormal[6][3]={ {0, -1,  0},
			    {1,  0,  0},
			    {0,  1,  0},
			    {-1, 0,  0},
			    {0,  0,  1},
			    {0,  0, -1}};
  for(int ii = 0; ii < 3; ii++) {
    vnds[ii]=0;
    for(int jj = 0; jj < 3; jj++) {
      vnds[ii] += codtau[ii][jj] * h20_refnormal[ifa][jj];
    }
  }
}

void Ref2Phy(__constant real *physnode,
             const real xref[3],
             real dphiref[3],       // can be NULL
             const int ifa,         // only needed for vnds calculation
             real xphy[3],          // can be NULL
             real dtau[3][3],       // can be NULL
             real codtau[3][3],     // can be NULL
             real dphi[3],          // can be NULL 
             real vnds[3])          // can be NULL
{
  // compute the mapping and its jacobian

  // gradient of the shape functions and value (4th component)
  // of the shape functions
  real gradphi[20][4];
  compute_gradphi(xref, gradphi);

  if (xphy != NULL)
    compute_xphy(physnode, gradphi, xphy);

  if (dtau != NULL)
    compute_dtau(physnode, gradphi, dtau);

  if (codtau != NULL)
    compute_codtau(dtau, codtau);
  
  if (dphi != NULL)
    compute_dphi(dphiref, codtau, dphi);

  if (vnds != NULL)
    ComputeNormal(codtau, ifa, vnds);
}

void Phy2Ref(__constant real *physnode, real *xphy, real *xref);

// Given parameters deg and nraf and input ipg, compute the reference
// coordinages (xpg) and the weght of the Gauss piont (wpg).
void ref_pg_vol(const int *deg, const int *nraf, 
		const int ipg, real *xpg, real *wpg)
{
  int ix[3], ic[3];
  ipg_to_xyz(nraf, deg, ic, ix, &ipg);

  real hx = 1 / (real) nraf[0];
  real hy = 1 / (real) nraf[1];
  real hz = 1 / (real) nraf[2];

  int offset[3];
  offset[0] = gauss_lob_offset[deg[0]] + ix[0];
  offset[1] = gauss_lob_offset[deg[1]] + ix[1];
  offset[2] = gauss_lob_offset[deg[2]] + ix[2];

  xpg[0] = hx * (ic[0] + gauss_lob_point[offset[0]]);
  xpg[1] = hy * (ic[1] + gauss_lob_point[offset[1]]);
  xpg[2] = hz * (ic[2] + gauss_lob_point[offset[2]]);
  
  *wpg = hx * hy * hz *
    gauss_lob_weight[offset[0]]*
    gauss_lob_weight[offset[1]]*
    gauss_lob_weight[offset[2]];
}

int ref_pg_face(const int *deg, const int *raf,
		const int ifa, // face index
		int ipgf,      // index of point in the face
                real *xpg, real *wpg, real *xpgin)
{
  // For each face, give the dimension index i
  int axis_permut[6][4] = { {0, 2, 1, 0},
			    {1, 2, 0, 1},
			    {2, 0, 1, 1},
			    {2, 1, 0, 0},
			    {0, 1, 2, 1},
			    {1, 0, 2, 0} };

  int paxis[4] = {axis_permut[ifa][0],
		  axis_permut[ifa][1],
		  axis_permut[ifa][2],
		  axis_permut[ifa][3]};
  
  // approximation degree in each permuted direction
  int pdeg[3] = {deg[paxis[0]],	deg[paxis[1]],	deg[paxis[2]]};

  // number of subcells in each permuted direction
  int praf[3] = {raf[paxis[0]], raf[paxis[1]], raf[paxis[2]]};

  // permuted point index in subcell
  int pix[3];
  pix[0] = ipgf % (pdeg[0] + 1);
  ipgf /= (pdeg[0] + 1);
  pix[1] = ipgf % (pdeg[1] + 1);
  ipgf /= (pdeg[1] + 1);
  pix[2] = paxis[3] * pdeg[2]; // Equals 0 or d depending on the face

  // Compute permuted subcell  indices of the subface
  int pic[3];
  pic[0] = ipgf % praf[0];
  pic[1] = ipgf / praf[0];
  pic[2] = paxis[3] * (praf[2] - 1); // Equals 0 or raf-1

  real h[3] = {1.0 / (real) praf[0],
	       1.0 / (real) praf[1],
	       1.0 / (real) praf[2] };
  
  // non-permuted subcell index
  int ic[3];
  ic[paxis[0]] = pic[0];
  ic[paxis[1]] = pic[1];
  ic[paxis[2]] = pic[2];

  // non-permuted point index
  int ix[3];
  ix[paxis[0]] = pix[0];
  ix[paxis[1]] = pix[1];
  ix[paxis[2]] = pix[2];
  
  // Compute the global index of the Gauss-Lobatto point in the volume
  int ipgv = xyz_to_ipg(raf, deg, ic, ix); 

  // Compute the reference coordinates of the Gauss-Lobatto point in
  // the volume
  int offset[2] = {gauss_lob_offset[pdeg[0]] + pix[0],
		   gauss_lob_offset[pdeg[1]] + pix[1]};

  xpg[paxis[0]] = h[0] * (pic[0] + gauss_lob_point[offset[0]]);
  xpg[paxis[1]] = h[1] * (pic[1] + gauss_lob_point[offset[1]]);
  xpg[paxis[2]] = paxis[3];

  *wpg = h[0] * h[1] *
    gauss_lob_weight[offset[0]] * gauss_lob_weight[offset[1]];

  // If xpgin exists, compute a point slightly INSIDE the opposite
  // subcell along the face.
  // NB: in OpenCL, we _always_ compute xpgin, so the test can be removed.
  //if(xpgin != NULL) {
  real small = 1e-3;  //0.001
  real vsmall = 1e-6; //0.000001;

  xpgin[paxis[0]] = h[0] * (pic[0] + gauss_lob_point[offset[0]]);
  xpgin[paxis[1]] = h[1] * (pic[1] + gauss_lob_point[offset[1]]);

  // TODO: can this be better handled with ifa?
  if(pix[0] == 0)
    xpgin[paxis[0]] = h[0] * (pic[0] + gauss_lob_point[offset[0]] + small);
  if(pix[0] == pdeg[0])
    xpgin[paxis[0]] = h[0] * (pic[0] + gauss_lob_point[offset[0]] - small);

  if(pix[1] == 0)
    xpgin[paxis[1]] = h[1] * (pic[1] + gauss_lob_point[offset[1]] + small);
  if(pix[1] == pdeg[1])
    xpgin[paxis[1]] = h[1] * (pic[1] + gauss_lob_point[offset[1]] - small);

  if(paxis[3] == 0)
    xpgin[paxis[2]] = -vsmall;
  if(paxis[3] == 1)
    xpgin[paxis[2]] = 1.0 + vsmall;
  //}

  return ipgv;
}

#ifndef NUMFLUX
#define NUMFLUX NumFlux
#endif

void NumFlux(real wL[], real wR[], real *vnorm, real *flux) {
  real vn = sqrt(0.5) * (vnorm[0] + vnorm[1]);

  real vnp = vn > 0 ? vn : 0;
  real vnm = vn - vnp;

  flux[0] = vnp * wL[0] + vnm * wR[0];
}

#ifndef vlasov_mx
#define vlasov_mx 1
#endif

#ifndef vlasov_my
#define vlasov_my 1
#endif

#ifndef vlasov_vmax
#define vlasov_vmax 0.5
#endif

// Return the component of the vlasov velocity with index id.
real vlasov_vel(const int id, const int md)
{
  int mid = md / 2;
  real dv = vlasov_vmax / mid;
  return (id - mid) * dv;
}

// Sample flux for 2D Vlasov equation
void vlaTransNumFlux2d(real wL[], real wR[], real *vnorm, real *flux)
{
  for(int ix = 0; ix < vlasov_mx; ++ix) {
    real vx = vlasov_vel(ix, vlasov_mx);

    for(int iy = 0; iy < vlasov_my; ++iy) {
      real vy = vlasov_vel(iy, vlasov_my);
      
      real vn = vx * vnorm[0]	+ vy * vnorm[1];
      real vnp = vn > 0 ? vn : 0;
      real vnm = vn - vnp;
      
      // NB: assumes a certain memory distribution for the velocity
      // components at each point.
      int im = ix * vlasov_my + iy;
      flux[im] = vnp * wL[im] + vnm * wR[im];
    }
  }
}

#ifndef _M
#define _M 1
#endif

void cemracs2014_TransBoundaryFlux(real x[3], real t, 
				   real wL[], real *vnorm,
				   real *flux)
{
  real wR[_M];
  int m = vlasov_mx * vlasov_my;
  for(int i = 0; i < m; ++i)
    wR[i] = 0;
  vlaTransNumFlux2d(wL, wR, vnorm, flux);
}

// Sample boundary flux
void BoundaryFlux(real x[3], real t, real *wL, real *vnorm,
                  real *flux) 
{
  real wR[_M];
  real vx = sqrt(0.5) * (x[0] + x[1]);
  wR[0] = cos(vx - t);

  NUMFLUX(wL, wR, vnorm, flux);
}

//! \brief 1d GLOP weights for a given degree
//! \param[in] deg degree
//! \param[in] i glop index
//! \returns the glop weight
real wglop(int deg, int i) 
{
  return gauss_lob_weight[gauss_lob_offset[deg] + i];
}

void get_dtau(real x, real y, real z,
	      __constant real *physnode, real dtau[][3]);

// Get the logical index of the gaussian point given the coordinate
// p[] of the point in the subcell and the index of the subcell icell.
int ipg(const int npg[], const int p[], const int icell) 
{
  return npg[0] * npg[1] * npg[2] * icell
    + p[0] + npg[0] * (p[1] + npg[1] * p[2]);
}

// Compute the surface terms inside one macrocell
__kernel
void DGFlux(__constant int *param,     // interp param
	    int dim0,                  // face direction
	    __constant real *physnode, // macrocell nodes
	    __global   real *wn,       // field values
	    __global   real *dtwn,     // time derivative
	    __local    real *wnloc     // wn and dtwn in local memory
	    )
{
  // Use __local memory in DGFlux kernel?
#define DGFLUX_LOCAL 1

  const int m = param[0];
  const int deg[3] = {param[1], param[2], param[3]};
  const int npg[3] = {deg[0] + 1, deg[1] + 1, deg[2] + 1};
  const int nraf[3] = {param[4], param[5], param[6]};

  int dim1 = (dim0 + 1) % 3;
  int dim2 = (dim1 + 1) % 3;

  // Subcell id
  int icL[3];
  int icR[3];
  int icell = get_group_id(0);

  icL[dim0] = icell % (nraf[dim0] - 1);
  icL[dim1] = (icell / (nraf[dim0] - 1)) % nraf[dim1];
  icL[dim2] = icell / (nraf[dim0]-1) / nraf[dim1];

  icR[dim0] = icL[dim0] + 1;
  icR[dim1] = icL[dim1];
  icR[dim2] = icL[dim2];

  __local real *wnlocL = wnloc;
  __local real *wnlocR = wnloc + get_local_size(0) * m;
  __local real *dtwnlocL = wnloc + 2 * get_local_size(0) * m;
  __local real *dtwnlocR = wnloc + 3 * get_local_size(0) * m;

  int pL[3];

  
#if DGFLUX_LOCAL
  for(int i = 0; i < m; i++) {
    int iread = get_local_id(0) + i * get_local_size(0);
    int iv = iread % m;
    int ipg = iread / m;
    
    int p[3];
    p[dim1] = ipg % npg[dim1];
    p[dim2] = ipg / npg[dim1];

    // Left point
    p[dim0] = deg[dim0];
    int ipgL = xyz_to_ipg(nraf, deg, icL, p);
    int imemL = VARINDEX(param, ipgL, iv);
    // wnlocL[iread] = wn[imemL];
    wnlocL[ipg * m + iv] = wn[imemL];

    // Right point
    p[dim0] = 0;
    int ipgR =  xyz_to_ipg(nraf, deg, icR, p);
    int imemR = VARINDEX(param, ipgR, iv);
    // wnlocR[iread] = wn[imemR];
    wnlocR[ipg * m + iv] = wn[imemR];
  }

  //barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
  barrier(CLK_LOCAL_MEM_FENCE);
#else
  // Gauss point id where we compute the jacobian

  int pR[3];
  
  //ipg_to_xyz(get_local_id(0), p, npg
  {
    int ipg = get_local_id(0);
    pL[dim0] = deg[dim0];
    pL[dim1] = ipg % npg[dim1];
    pL[dim2] = (ipg / npg[dim1]);

    pR[dim0] = 0;
    pR[dim1] = pL[dim1];
    pR[dim2] = pL[dim2];
  }
  
  real wL[_M], wR[_M];
  int ipgL = xyz_to_ipg(nraf, deg, icL, pL);
  int ipgR = xyz_to_ipg(nraf, deg, icR, pR);
  for(int iv = 0; iv < m; iv++) {
    int imemL = VARINDEX(param, ipgL, iv);
    wL[iv] = wn[imemL];
    int imemR = VARINDEX(param, ipgR, iv);
    wR[iv] = wn[imemR];
  }
#endif
 
  // Gauss point id where we compute the Jacobian
  //int pL[3], pR[3];
  //ipg_to_xyz(get_local_id(0), p, npg
  {
    int ipg = get_local_id(0);
    pL[dim0] = deg[dim0];
    pL[dim1] = ipg % npg[dim1];
    pL[dim2] = ipg / npg[dim1];
  }

  real hx = 1.0 / (real) nraf[0];
  real hy = 1.0 / (real) nraf[1];
  real hz = 1.0 / (real) nraf[2];

  int offset[3] = {gauss_lob_offset[deg[0]] + pL[0],
  		   gauss_lob_offset[deg[1]] + pL[1],
  		   gauss_lob_offset[deg[2]] + pL[2]};

  real x = hx * (icL[0] + gauss_lob_point[offset[0]]);
  real y = hy * (icL[1] + gauss_lob_point[offset[1]]);
  real z = hz * (icL[2] + gauss_lob_point[offset[2]]);

  /* real wpg = hx * hy * hz */
  /*   * gauss_lob_weight[offset[0]] */
  /*   * gauss_lob_weight[offset[1]] */
  /*   * gauss_lob_weight[offset[2]]; */

  real codtau[3][3];
  {
    real dtau[3][3];
    get_dtau(x, y, z, physnode, dtau);
    compute_codtau(dtau, codtau);
  }

  real h1h2 = 1.0 / nraf[dim1] / nraf[dim2];
  real vnds[3];
  vnds[0] = codtau[0][dim0] * h1h2;
  vnds[1] = codtau[1][dim0] * h1h2;
  vnds[2] = codtau[2][dim0] * h1h2;

#if DGFLUX_LOCAL
  real wL[_M], wR[_M]; // TODO: remove?
  __local real *wnL = wnlocL + get_local_id(0) * m;
  __local real *wnR = wnlocR + get_local_id(0) * m;
  for(int iv = 0; iv < m; iv++) {
    wL[iv] = wnL[iv];
    wR[iv] = wnR[iv];
  }
#else

  ipgL = xyz_to_ipg(nraf, deg, icL, pL);
  ipgR = xyz_to_ipg(nraf, deg, icR, pR);

  for(int iv = 0; iv < m; iv++) {
    int imemL = VARINDEX(param, ipgL, iv);
    wL[iv] = wn[imemL];
    int imemR = VARINDEX(param, ipgR, iv);
    wR[iv] = wn[imemR];
  }
#endif

  // TODO: wL and wR could be passed without a copy to __private.
  // (ie we can just pass *wnL and *wnR).
  real flux[_M];
  NUMFLUX(wL, wR, vnds, flux);

  real wpgs = wglop(deg[dim1], pL[dim1]) * wglop(deg[dim2], pL[dim2]);

#if DGFLUX_LOCAL
  __local real *dtwnL = dtwnlocL + get_local_id(0) * m;
  __local real *dtwnR = dtwnlocR + get_local_id(0) * m;
  for(int iv = 0; iv < m; ++iv) {
    // write flux to local memory
    dtwnL[iv] = -flux[iv] * wpgs;
    dtwnR[iv] =  flux[iv] * wpgs;
  }

  //barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
  barrier(CLK_LOCAL_MEM_FENCE);

  // Postfetch: 2m writes to global memory
  for(int i = 0; i < m; i++) {
    int p[3];
    int iread = get_local_id(0) + i * get_local_size(0);
    int iv = iread % m;
    int ipg = iread / m;
    
    p[dim1] = ipg % npg[dim1];
    p[dim2] = ipg / npg[dim1];

    // Left point
    p[dim0] = deg[dim0];
    int ipgL = xyz_to_ipg(nraf, deg, icL, p);
    int imemL = VARINDEX(param, ipgL, iv);
    // wnlocL[iread] = wn[imemL];
    dtwn[imemL] += dtwnlocL[ipg * m + iv];
    
    // Right point
    p[dim0] = 0;
    int ipgR =  xyz_to_ipg(nraf, deg, icR, p);
    int imemR = VARINDEX(param, ipgR, iv);
    // wnlocR[iread] = wn[imemR];
    dtwn[imemR] += dtwnlocR[ipg * m + iv];
  }
#else
  for(int iv = 0; iv < m; ++iv) {
    //int ipgL = ipg(npg, p, icell);
    //int imemL = VARINDEX(param, ie, ipgL, iv);

    int imemL = VARINDEX(param, ipgL, iv);
    dtwn[imemL] -= flux[iv] * wpgs;

    int imemR = VARINDEX(param, ipgR, iv);
    dtwn[imemR] += flux[iv] * wpgs;
  }
#endif
}

__kernel
void set_buffer_to_zero(__global real *w)
{
  w[get_global_id(0)] = 0.0;
}

#ifndef BOUNDARYFLUX
#define BOUNDARYFLUX BoundaryFlux
#endif

void prefetch_macrocell(const __global real *in,
			__local real *out,
			__constant int *param)
{
  const int m = param[0];
  int icell = get_group_id(0);
  for(int i = 0; i < m ; ++i){
    int iread = get_local_id(0) + i * get_local_size(0);
    int iv = iread % m;
    int ipgloc = iread / m;
    int ipgL = ipgloc + icell * get_local_size(0);
    int imem = VARINDEX(param, ipgL, iv);
    int imemloc = iv + ipgloc * m;
    
    out[imemloc] = in[imem];
  }
}

void zero_macrocell_buffer(__local real *out,
			   __constant int *param
			   )
{
  const int m = param[0];
  int icell = get_group_id(0);
  for(int i = 0; i < m ; ++i){
    int iread = get_local_id(0) + i * get_local_size(0);
    int iv = iread % m;
    int ipgloc = iread / m;
    int ipgL = ipgloc + icell * get_local_size(0);
    int imem = VARINDEX(param, ipgL, iv);
    int imemloc = iv + ipgloc * m;
    
    out[imemloc] = 0.0;
  }
}

void postfetch_macrocell(const __local real *in,
			 __global real *out,
			 __constant int *param
			 )
{
  const int m = param[0];
  int icell = get_group_id(0);
  for(int i = 0; i < m; ++i){
    int iread = get_local_id(0) + i * get_local_size(0);
    int iv = iread % m;
    int ipgloc = iread / m ;
    int ipg = ipgloc + icell * get_local_size(0);
    int imem = VARINDEX(param, ipg, iv);
    int imemloc = ipgloc * m + iv;
    out[imem] += in[imemloc];
  }
}

void postfetch_macrocell_overwrite(const __local real *in,
			 __global real *out,
			 __constant int *param
			 )
{
  const int m = param[0];
  int icell = get_group_id(0);
  for(int i = 0; i < m; ++i){
    int iread = get_local_id(0) + i * get_local_size(0);
    int iv = iread % m;
    int ipgloc = iread / m ;
    int ipg = ipgloc + icell * get_local_size(0);
    int imem = VARINDEX(param, ipg, iv);
    int imemloc = ipgloc * m + iv;
    out[imem] = in[imemloc];
  }
}

void compute_volume(__constant int *param,     // interp param
		    __constant real *physnode, // macrocell nodes
		    __local real *wnloc,       // cache for wn
		    __local real *dtwnloc      // cache for dtwn
		    )
{
  const int m = param[0];
  const int deg[3] = {param[1],param[2], param[3]};
  const int npg[3] = {deg[0] + 1, deg[1] + 1, deg[2] + 1};
  const int nraf[3] = {param[4], param[5], param[6]};
  
  // subcell id
  int icell = get_group_id(0);
  int icL[3];
  icL[0] = icell % nraf[0];
  icL[1] = (icell / nraf[0]) % nraf[1];
  icL[2] = icell / nraf[0] / nraf[1];

  // gauss point id where we compute the jacobian
  int p[3];
  //ipg_to_xyz(get_local_id(0), p, npg
  {
    int ipg = get_local_id(0);
    p[0] = ipg % npg[0];
    p[1] = (ipg / npg[0]) % npg[1];
    p[2] = ipg / npg[0] / npg[1];
  }

  // ref coordinates
  real hx = 1.0 / (real) nraf[0];
  real hy = 1.0 / (real) nraf[1];
  real hz = 1.0 / (real) nraf[2];

  int offset[3] = {gauss_lob_offset[deg[0]] + p[0],
		   gauss_lob_offset[deg[1]] + p[1],
		   gauss_lob_offset[deg[2]] + p[2]};

  real x = hx * (icL[0] + gauss_lob_point[offset[0]]);
  real y = hy * (icL[1] + gauss_lob_point[offset[1]]);
  real z = hz * (icL[2] + gauss_lob_point[offset[2]]);

  real wpg = hx * hy * hz
    * gauss_lob_weight[offset[0]]
    * gauss_lob_weight[offset[1]]
    * gauss_lob_weight[offset[2]];

  real codtau[3][3];
  {
    real dtau[3][3];
    get_dtau(x, y, z, physnode, dtau);
    compute_codtau(dtau, codtau);
  }

  real wL[_M];
  int ipgL = ipg(npg, p, 0);
  __local real *wnloc0 = wnloc + ipgL * m;

  for(int iv = 0; iv < m; iv++) {
    wL[iv] = wnloc0[iv];
  }

  real flux[_M];
  for(int dim0 = 0; dim0 < 3; dim0++) {
    int q[3] = {p[0], p[1], p[2]};

    // Loop on the "cross" points
    for(int iq = 0; iq < npg[dim0]; iq++) {
      q[dim0] = (p[dim0] + iq) % npg[dim0];
      real dphiref[3] = {0, 0, 0};
      dphiref[dim0] = dlag(deg[dim0], q[dim0], p[dim0]) * nraf[dim0];
      real dphi[3];
      for(int ii = 0; ii < 3; ii++) {
	/* dphi[ii] = 0; */
	/* for(int jj = 0; jj < 3; jj++) { */
	/*   dphi[ii] += codtau[ii][jj] * dphiref[jj]; */
	/* } */
	real *codtauii = codtau[ii];
	dphi[ii] 
	  = codtauii[0] * dphiref[0]
	  + codtauii[1] * dphiref[1]
	  + codtauii[2] * dphiref[2];
      }

      NUMFLUX(wL, wL, dphi, flux);

      int ipgR = ipg(npg, q, 0);

      int imemR0loc = ipgR * m;
      __local real *dtwnloc0 =  dtwnloc + imemR0loc;
      for(int iv = 0; iv < m; iv++) {
	dtwnloc0[iv] += flux[iv] * wpg;
      }
    }

  } // dim0 loop

}

void compute_volume_global(__constant int *param,     // interp param
			   __constant real *physnode, // macrocell nodes
			   __global real *wn,         // field values
			   __global real *dtwn       // time derivative
			   )
{
  
  const int m = param[0];
  const int deg[3] = {param[1],param[2], param[3]};
  const int npg[3] = {deg[0] + 1, deg[1] + 1, deg[2] + 1};
  const int nraf[3] = {param[4], param[5], param[6]};
  
  // subcell id
  int icell = get_group_id(0);
  int icL[3];
  icL[0] = icell % nraf[0];
  icL[1] = (icell / nraf[0]) % nraf[1];
  icL[2] = icell / nraf[0] / nraf[1];

  // gauss point id where we compute the jacobian
  int p[3];
  //ipg_to_xyz(get_local_id(0), p, npg
  {
    int ipg = get_local_id(0);
    p[0] = ipg % npg[0];
    p[1] = (ipg / npg[0]) % npg[1];
    p[2] = ipg / npg[0] / npg[1];
  }

  // ref coordinates
  real hx = 1.0 / (real) nraf[0];
  real hy = 1.0 / (real) nraf[1];
  real hz = 1.0 / (real) nraf[2];

  int offset[3] = {gauss_lob_offset[deg[0]] + p[0],
		   gauss_lob_offset[deg[1]] + p[1],
		   gauss_lob_offset[deg[2]] + p[2]};

  real x = hx * (icL[0] + gauss_lob_point[offset[0]]);
  real y = hy * (icL[1] + gauss_lob_point[offset[1]]);
  real z = hz * (icL[2] + gauss_lob_point[offset[2]]);

  real wpg = hx * hy * hz
    * gauss_lob_weight[offset[0]]
    * gauss_lob_weight[offset[1]]
    * gauss_lob_weight[offset[2]];

  real codtau[3][3];
  {
    real dtau[3][3];
    get_dtau(x, y, z, physnode, dtau);
    compute_codtau(dtau, codtau);
  }

  real wL[_M];
  int ipgL = ipg(npg, p, 0);
  for(int iv = 0; iv < m; iv++) {
    int ipgL = ipg(npg, p, icell);
    int imemL = VARINDEX(param, ipgL, iv);
    wL[iv] = wn[imemL];
  }

  real flux[_M];
  for(int dim0 = 0; dim0 < 3; dim0++) {
    int q[3] = {p[0], p[1], p[2]};

    // Loop on the "cross" points
    for(int iq = 0; iq < npg[dim0]; iq++) {
      q[dim0] = (p[dim0] + iq) % npg[dim0];
      real dphiref[3] = {0, 0, 0};
      dphiref[dim0] = dlag(deg[dim0], q[dim0], p[dim0]) * nraf[dim0];
      real dphi[3];
      for(int ii = 0; ii < 3; ii++) {
	real *codtauii = codtau[ii];
	dphi[ii] 
	  = codtauii[0] * dphiref[0]
	  + codtauii[1] * dphiref[1]
	  + codtauii[2] * dphiref[2];
      }

      NUMFLUX(wL, wL, dphi, flux);

      int ipgR = ipg(npg, q, icell);
      int imemR0 = VARINDEX(param, ipgR, 0);
      __global real *dtwn0 = dtwn + imemR0; 
      for(int iv = 0; iv < m; iv++) {
     	dtwn0[iv] += flux[iv] * wpg;
      }
    }

  } // dim0 loop

}
  
// Compute the volume  terms inside  one macrocell
__kernel
void DGVolume(__constant int *param,     // interp param
	      __constant real *physnode, // macrocell nodes
              __global real *wn,         // field values
	      __global real *dtwn,       // time derivative
	      __local real *wnloc,       // cache for wn
	      __local real *dtwnloc      // cache for dtwn
	      )
{
  // Use __local memory in DGVolume kernel?
#define DGVolume_LOCAL 1

#if DGVolume_LOCAL

  prefetch_macrocell(wn, wnloc, param);
  zero_macrocell_buffer(dtwnloc, param);
  
  barrier(CLK_LOCAL_MEM_FENCE);

  compute_volume(param, physnode, wnloc, dtwnloc);

  barrier(CLK_LOCAL_MEM_FENCE);

  postfetch_macrocell(dtwnloc, dtwn, param);

#else

  compute_volume_global(param, physnode, wn, dtwn);
  
#endif
}

void mass_division(__constant int *param,
		   const  __global real *mass,
		   __local real *dtwnloc)
{
  const int m = param[0];

  // The mass is __global
  int ipg = get_global_id(0);
  real overmassipg = 1.0 / mass[ipg];

  // dtwnloc is __local
  int ipgloc = get_local_id(0);
  __local real *dtwnloc0 =  dtwnloc + ipgloc * m;
  for(int iv = 0; iv < m; iv++) {
    dtwnloc0[iv] *= overmassipg;
  }
}

// Apply division by the mass matrix on one macrocell
__kernel
void DGMass(__constant int *param,      // interp param
            __constant real *physnode,  // macrocell nodes
	    const __global real *mass,  // macrocell masses
            __global real *dtwn,        // time derivative
	    __local real *dtwnloc       // cache for dtwn
	    )
{
  // TODO: if there's enough space, make mass __constant  

#if 0
  const int ipg = get_global_id(0);
  const int m = param[0];

  real overmass = 1.0 / mass[ipg];
  __global real *dtwn0 = dtwn + m * ipg;
  for(int iv = 0; iv < m; iv++) {
    dtwn0[iv] *= overmass;
  }

#else

    prefetch_macrocell(dtwn, dtwnloc, param);

    barrier(CLK_LOCAL_MEM_FENCE);
  
    mass_division(param, mass, dtwnloc);
  
    barrier(CLK_LOCAL_MEM_FENCE);

    postfetch_macrocell_overwrite(dtwnloc, dtwn, param);
#endif
}

__kernel
void ExtractInterface(const int d2,
		      const int stride0,   // stride for d0
		      const int stride1,   // stride for d1
		      const int stride2,   // stride for d2
		      const __global real *wn,   // volumic input
		      __global real *wface // output
		      )
{
  const int d0 = get_global_id(0); // first dimension
  const int d1 = get_global_id(1); // second dimension
  const int iv = get_global_id(2); // m index
  
  // We assume that the MacroCell's volumic points are given in the
  // standard C-ordering.

  // We must pre-compute the input strides because we don't know the
  // ordering for the face.
  int pin = d0 * stride0 + d1 * stride1 + d2 * stride2 + iv;

  const int m = get_global_size(2);
  const int n0 = get_global_size(0);
  // The output is in wface, with corrdinates (d0, d1, iv) in
  // [0, n0 -1] x [0, n1 -1] x [0, m - 1]
  int pout = d1 * n0 * m  + d0 * m + iv;
  
  wface[pout] = wn[pin];
}

__kernel
void DGInterfaceFlux(__constant int *param,        // interp param
		     __global real *wfaceL,
		     __global real *wfaceR
		     )
{
  const int d0 = get_global_id(0); // first dimension
  const int d1 = get_global_id(1); // second dimension
  const int iv = get_global_id(2); // m index

  const int n1 = get_global_size(1);

  int ipgf = d0 * n1 + d0;
  // FIXME: complete tis function.

  // find vnds

  // find wL and wR points

  // compute flux

  // put the computed flux into the correct __global_real buffer (or
  // just copy it back to the inputs????
}

// Compute the Discontinuous Galerkin inter-macrocells boundary terms.
// Second implementation with a loop on the faces.
__kernel
void DGMacroCellInterface(__constant int *param,        // interp param
                          int locfaL,                   // left face index
			  int locfaR,                   // right face index
                          __constant real *physnodeL,   // left physnode
			  __constant real *physnodeR,   // right physnode
                          __global real *wnL,           // field 
                          __global real *dtwnL,         // time derivative
                          __global real *wnR,           // field 
                          __global real *dtwnR,         // time derivative
			  __local real *cache           // local mem
			  )
{
  // TODO: use __local real *cache.

  // Index of the point on the face.
  int ipgfL = get_global_id(0);

  const int m = param[0];
  const int ndeg[3] = {param[1], param[2], param[3]};
  const int nraf[3] = {param[4], param[5], param[6]};

  real xpgref[3]; // reference point for L
  real xpgref_in[3]; // reference point slightly in R
  real wpg;
  int ipgL = ref_pg_face(ndeg, nraf, locfaL, ipgfL, xpgref, &wpg, xpgref_in);
  
  // Normal vector at gauss point based on xpgref
  real vnds[3];
  {
    real gradphi[20][4];
    compute_gradphi(xpgref, gradphi);

    real dtau[3][3];  
    compute_dtau(physnodeL, gradphi, dtau);

    real codtau[3][3];
    compute_codtau(dtau, codtau);

    ComputeNormal(codtau, locfaL, vnds);
  }

  // Find the volumic index for the opposite side:
  int ipgR;
  {
    real gradphi[20][4];
    compute_gradphi(xpgref_in, gradphi);

    real xpg_in[3];
    compute_xphy(physnodeL, gradphi, xpg_in);

    real period[3] = {_PERIODX,_PERIODY,_PERIODZ};
    PeriodicCorrection(xpg_in, period);

    real xrefL[3];
    Phy2Ref(physnodeR, xpg_in, xrefL);
    ipgR = ref_ipg(param + 1, xrefL);
  }
  
  // Test code
  /* { */
  /*   real xpgR[3], xrefR[3], wpgR; */
  /*   ref_pg_vol(param + 1, ipgR, xrefR, &wpgR, NULL); */
  /*   Ref2Phy(physnodeR, */
  /* 	  xrefR, */
  /* 	  NULL, -1, // dphiref, ifa */
  /* 	  xpgR, NULL,   */
  /* 	  NULL, NULL, NULL); // codtau, dphi,vnds */
  /*   assert(Dist(xpgR, xpg) < 1e-10); */
  /* }	 */
  
  real wL[_M];
  real wR[_M];

  int imemL0 = VARINDEX(param, ipgL, 0);
  int imemR0 = VARINDEX(param, ipgR, 0);

  __global real *wnL0 = wnL + imemL0;
  __global real *wnR0 = wnR + imemR0;
  for(int iv = 0; iv < m; iv++) {
    wL[iv] = wnL0[iv];
    wR[iv] = wnR0[iv];
  }

  real flux[_M];
  NUMFLUX(wL, wR, vnds, flux);

  __global real *dtwnL0 = dtwnL + imemL0;
  __global real *dtwnR0 = dtwnR + imemR0;
  for(int iv = 0; iv < m; ++iv) {
    real fluxwpg = flux[iv] * wpg;
    dtwnL0[iv] -= fluxwpg;
    dtwnR0[iv] += fluxwpg;
  }
}

// Compute the Discontinuous Galerkin inter-macrocells boundary terms.
// Second implementation with a loop on the faces.
__kernel
void DGBoundary(__constant int *param,      // interp param
		real tnow,                  // current time
		int locfa,                 // left face index
		__constant real *physnode, // geometry for all mcells
		__global real *wn,          // field 
		__global real *dtwn,        // time derivative
		__local real *cache         // local mem
		)
{
  // TODO: use __local real *cache.

  int ipgf = get_global_id(0);

  const int m = param[0];
  const int ndeg[3] = {param[1], param[2], param[3]};
  const int nraf[3] = {param[4], param[5], param[6]};

  real xref[3];    // reference coordinates
  real xref_in[3]; // unused
  real wpg;        // weighting for the GL point
  int ipgL = ref_pg_face(ndeg, nraf, locfa, ipgf, xref, &wpg, xref_in);

  // normal vector
  real vnds[3]; 
  // physical coordinates
  real xpg[3];  
  {
    real gradphi[20][4];
    compute_gradphi(xref, gradphi);
    
    compute_xphy(physnode, gradphi, xpg);

    real dtau[3][3];  
    compute_dtau(physnode, gradphi, dtau);

    real codtau[3][3];
    compute_codtau(dtau, codtau);

    ComputeNormal(codtau, locfa, vnds);
  }

  real wL[_M];
  
  int imemL0 = VARINDEX(param, ipgL, 0);
  __global real *wn0 = wn + imemL0;
  for(int iv = 0; iv < m; ++iv) {
    wL[iv] = wn0[iv];
  }

  real flux[_M];
  BOUNDARYFLUX(xpg, tnow, wL, vnds, flux);

  __global real *dtwn0 = dtwn + imemL0; 
  for(int iv = 0; iv < m; ++iv) {
    dtwn0[iv] -= flux[iv] * wpg;
  }
}

void get_dtau(real x, real y, real z,
	      __constant real *p, real dtau[][3]) 
{
  // Gradient of the shape functions and value (4th component) of the
  // shape functions

  /* real gradphi[20][3]; */
  /* //real x,y,z; */
  /* // this fills the values of gradphi */
  /* gradphi[0][0] = (-1 + z) * (-1 + y) * (2 * y + 2 * z + 4 * x - 3); */
  /* gradphi[0][1] = (-1 + z) * (-1 + x) * (2 * x + 2 * z - 3 + 4 * y); */
  /* gradphi[0][2] = (-1 + y) * (-1 + x) * (2 * x + 2 * y - 3 + 4 * z); */
  /* gradphi[1][0] = (-1 + z) * (-1 + y) * (-2 * y - 2 * z - 1 + 4 * x); */
  /* gradphi[1][1] = x * (-1 + z) * (-2 * z + 1 - 4 * y + 2 * x); */
  /* gradphi[1][2] = x * (-1 + y) * (-2 * y + 1 - 4 * z + 2 * x); */
  /* gradphi[2][0] = -y * (-1 + z) * (2 * y - 2 * z - 3 + 4 * x); */
  /* gradphi[2][1] = -x * (-1 + z) * (-2 * z - 3 + 4 * y + 2 * x); */
  /* gradphi[2][2] = -x * y * (2 * y - 4 * z - 1 + 2 * x); */
  /* gradphi[3][0] = -y * (-1 + z) * (-2 * y + 2 * z - 1 + 4 * x); */
  /* gradphi[3][1] = -(-1 + z) * (-1 + x) * (2 * x + 2 * z + 1 - 4 * y); */
  /* gradphi[3][2] = -y * (-1 + x) * (2 * x - 1 + 4 * z - 2 * y); */
  /* gradphi[4][0] = -z * (-1 + y) * (2 * y - 2 * z - 1 + 4 * x); */
  /* gradphi[4][1] = -z * (-1 + x) * (2 * x - 2 * z - 1 + 4 * y); */
  /* gradphi[4][2] = -(-1 + y) * (-1 + x) * (2 * x - 4 * z + 2 * y + 1); */
  /* gradphi[5][0] = -z * (-1 + y) * (-2 * y + 2 * z + 4 * x - 3); */
  /* gradphi[5][1] = -x * z * (2 * z - 4 * y - 1 + 2 * x); */
  /* gradphi[5][2] = -x * (-1 + y) * (-2 * y - 3 + 4 * z + 2 * x); */
  /* gradphi[6][0] = y * z * (2 * y + 2 * z + 4 * x - 5); */
  /* gradphi[6][1] = x * z * (2 * z + 4 * y - 5 + 2 * x); */
  /* gradphi[6][2] = x * y * (2 * y + 4 * z - 5 + 2 * x); */
  /* gradphi[7][0] = y * z * (-2 * y - 2 * z + 1 + 4 * x); */
  /* gradphi[7][1] = z * (-1 + x) * (2 * x - 2 * z + 3 - 4 * y); */
  /* gradphi[7][2] = y * (-1 + x) * (2 * x - 2 * y + 3 - 4 * z); */
  /* gradphi[8][0] = -4 * (-1 + z) * (-1 + y) * (2 * x - 1); */
  /* gradphi[8][1] = -4 * x * (-1 + z) * (-1 + x); */
  /* gradphi[8][2] = -4 * x * (-1 + y) * (-1 + x); */
  /* gradphi[9][0] = -4 * y * (-1 + z) * (-1 + y); */
  /* gradphi[9][1] = -4 * (-1 + z) * (2 * y - 1) * (-1 + x); */
  /* gradphi[9][2] = -4 * y * (-1 + y) * (-1 + x); */
  /* gradphi[10][0] = -4 * z * (-1 + z) * (-1 + y); */
  /* gradphi[10][1] = -4 * z * (-1 + z) * (-1 + x); */
  /* gradphi[10][2] = -4 * (2 * z - 1) * (-1 + y) * (-1 + x); */
  /* gradphi[11][0] = 4 * y * (-1 + z) * (-1 + y); */
  /* gradphi[11][1] = 4 * x * (-1 + z) * (2 * y - 1); */
  /* gradphi[11][2] = 4 * x * y * (-1 + y); */
  /* gradphi[12][0] = 4 * z * (-1 + z) * (-1 + y); */
  /* gradphi[12][1] = 4 * x * z * (-1 + z); */
  /* gradphi[12][2] = 4 * x * (2 * z - 1) * (-1 + y); */
  /* gradphi[13][0] = 4 * y * (-1 + z) * (2 * x - 1); */
  /* gradphi[13][1] = 4 * x * (-1 + z) * (-1 + x); */
  /* gradphi[13][2] = 4 * x * y * (-1 + x); */
  /* gradphi[14][0] = -4 * y * z * (-1 + z); */
  /* gradphi[14][1] = -4 * x * z * (-1 + z); */
  /* gradphi[14][2] = -4 * x * y * (2 * z - 1); */
  /* gradphi[15][0] = 4 * y * z * (-1 + z); */
  /* gradphi[15][1] = 4 * z * (-1 + z) * (-1 + x); */
  /* gradphi[15][2] = 4 * y * (2 * z - 1) * (-1 + x); */
  /* gradphi[16][0] = 4 * z * (-1 + y) * (2 * x - 1); */
  /* gradphi[16][1] = 4 * x * z * (-1 + x); */
  /* gradphi[16][2] = 4 * x * (-1 + y) * (-1 + x); */
  /* gradphi[17][0] = 4 * y * z * (-1 + y); */
  /* gradphi[17][1] = 4 * z * (2 * y - 1) * (-1 + x); */
  /* gradphi[17][2] = 4 * y * (-1 + y) * (-1 + x); */
  /* gradphi[18][0] = -4 * y * z * (-1 + y); */
  /* gradphi[18][1] = -4 * x * z * (2 * y - 1); */
  /* gradphi[18][2] = -4 * x * y * (-1 + y); */
  /* gradphi[19][0] = -4 * y * z * (2 * x - 1); */
  /* gradphi[19][1] = -4 * x * z * (-1 + x); */
  /* gradphi[19][2] = -4 * x * y * (-1 + x); */
  /* for(int ii=0;ii<3;ii++){ */
  /*   for(int jj=0;jj<3;jj++){ */
  /*     dtau[ii][jj]=0; */
  /*   } */
  /*   for(int i=0;i<20;i++){ */
  /*     //printf("xyzphy=%f %f %f \n",physnode[3*i+0],physnode[3*i+1],physnode[3*i+2]); */
  /*     for(int jj=0;jj<3;jj++){ */
  /*       dtau[ii][jj]+=physnode[3*i+ii]*gradphi[i][jj];; */
  /*     } */
  /*   } */
  /* } */

  dtau[0][0]=2*(-1+z)*(-1+y)*(-1+x)*p[0]+2*x*(-1+z)*(-1+y)*p[3]-2*x*y*(-1+z)*p[6]-2*y*(-1+z)*(-1+x)*p[9]-2*z*(-1+y)*(-1+x)*p[12]-2*x*z*(-1+y)*p[15]+2*x*y*z*p[18]+2*y*z*(-1+x)*p[21]+(-1+z)*(-1+y)*(2*x+2*y+2*z-1)*p[0]+(-1+z)*(-1+y)*(-2*y-2*z+2*x-1)*p[3]-y*(-1+z)*(2*y-2*z-3+2*x)*p[6]-y*(-1+z)*(2*x+2*z+1-2*y)*p[9]-z*(-1+y)*(2*x+2*y-2*z+1)*p[12]-z*(-1+y)*(-2*y+2*z-3+2*x)*p[15]+y*z*(2*y+2*z-5+2*x)*p[18]+y*z*(2*x-2*z+3-2*y)*p[21]-4*(-1+z)*(-1+y)*(-1+x)*p[24]-4*x*(-1+z)*(-1+y)*p[24]-4*y*(-1+z)*(-1+y)*p[27]-4*z*(-1+z)*(-1+y)*p[30]+4*y*(-1+z)*(-1+y)*p[33]+4*z*(-1+z)*(-1+y)*p[36]+4*y*(-1+z)*(-1+x)*p[39]+4*x*y*(-1+z)*p[39]-4*y*z*(-1+z)*p[42]+4*y*z*(-1+z)*p[45]+4*z*(-1+y)*(-1+x)*p[48]+4*x*z*(-1+y)*p[48]+4*y*z*(-1+y)*p[51]-4*y*z*(-1+y)*p[54]-4*y*z*(-1+x)*p[57]-4*x*y*z*p[57];

  dtau[0][1]=2*(-1+z)*(-1+y)*(-1+x)*p[0]-2*x*(-1+z)*(-1+y)*p[3]-2*x*y*(-1+z)*p[6]+2*y*(-1+z)*(-1+x)*p[9]-2*z*(-1+y)*(-1+x)*p[12]+2*x*z*(-1+y)*p[15]+2*x*y*z*p[18]-2*y*z*(-1+x)*p[21]+(-1+z)*(-1+x)*(2*x+2*y+2*z-1)*p[0]+x*(-1+z)*(-2*y-2*z+2*x-1)*p[3]-x*(-1+z)*(2*y-2*z-3+2*x)*p[6]-(-1+z)*(-1+x)*(2*x+2*z+1-2*y)*p[9]-z*(-1+x)*(2*x+2*y-2*z+1)*p[12]-x*z*(-2*y+2*z-3+2*x)*p[15]+x*z*(2*y+2*z-5+2*x)*p[18]+z*(-1+x)*(2*x-2*z+3-2*y)*p[21]-4*x*(-1+z)*(-1+x)*p[24]-4*(-1+z)*(-1+y)*(-1+x)*p[27]-4*y*(-1+z)*(-1+x)*p[27]-4*z*(-1+z)*(-1+x)*p[30]+4*x*(-1+z)*(-1+y)*p[33]+4*x*y*(-1+z)*p[33]+4*x*z*(-1+z)*p[36]+4*x*(-1+z)*(-1+x)*p[39]-4*x*z*(-1+z)*p[42]+4*z*(-1+z)*(-1+x)*p[45]+4*x*z*(-1+x)*p[48]+4*z*(-1+y)*(-1+x)*p[51]+4*y*z*(-1+x)*p[51]-4*x*z*(-1+y)*p[54]-4*x*y*z*p[54]-4*x*z*(-1+x)*p[57];

  dtau[0][2]=2*(-1+z)*(-1+y)*(-1+x)*p[0]-2*x*(-1+z)*(-1+y)*p[3]+2*x*y*(-1+z)*p[6]-2*y*(-1+z)*(-1+x)*p[9]+2*z*(-1+y)*(-1+x)*p[12]-2*x*z*(-1+y)*p[15]+2*x*y*z*p[18]-2*y*z*(-1+x)*p[21]+(-1+y)*(-1+x)*(2*x+2*y+2*z-1)*p[0]+x*(-1+y)*(-2*y-2*z+2*x-1)*p[3]-x*y*(2*y-2*z-3+2*x)*p[6]-y*(-1+x)*(2*x+2*z+1-2*y)*p[9]-(-1+y)*(-1+x)*(2*x+2*y-2*z+1)*p[12]-x*(-1+y)*(-2*y+2*z-3+2*x)*p[15]+x*y*(2*y+2*z-5+2*x)*p[18]+y*(-1+x)*(2*x-2*z+3-2*y)*p[21]-4*x*(-1+y)*(-1+x)*p[24]-4*y*(-1+y)*(-1+x)*p[27]-4*(-1+z)*(-1+y)*(-1+x)*p[30]-4*z*(-1+y)*(-1+x)*p[30]+4*x*y*(-1+y)*p[33]+4*x*(-1+z)*(-1+y)*p[36]+4*x*z*(-1+y)*p[36]+4*x*y*(-1+x)*p[39]-4*x*y*(-1+z)*p[42]-4*x*y*z*p[42]+4*y*(-1+z)*(-1+x)*p[45]+4*y*z*(-1+x)*p[45]+4*x*(-1+y)*(-1+x)*p[48]+4*y*(-1+y)*(-1+x)*p[51]-4*x*y*(-1+y)*p[54]-4*x*y*(-1+x)*p[57];

  dtau[1][0]=2*(-1+z)*(-1+y)*(-1+x)*p[1]+2*x*(-1+z)*(-1+y)*p[4]-2*x*y*(-1+z)*p[7]-2*y*(-1+z)*(-1+x)*p[10]-2*z*(-1+y)*(-1+x)*p[13]-2*x*z*(-1+y)*p[16]+2*x*y*z*p[19]+2*y*z*(-1+x)*p[22]+(-1+z)*(-1+y)*(2*x+2*y+2*z-1)*p[1]+(-1+z)*(-1+y)*(-2*y-2*z+2*x-1)*p[4]-y*(-1+z)*(2*y-2*z-3+2*x)*p[7]-y*(-1+z)*(2*x+2*z+1-2*y)*p[10]-z*(-1+y)*(2*x+2*y-2*z+1)*p[13]-z*(-1+y)*(-2*y+2*z-3+2*x)*p[16]+y*z*(2*y+2*z-5+2*x)*p[19]+y*z*(2*x-2*z+3-2*y)*p[22]-4*(-1+z)*(-1+y)*(-1+x)*p[25]-4*x*(-1+z)*(-1+y)*p[25]-4*y*(-1+z)*(-1+y)*p[28]-4*z*(-1+z)*(-1+y)*p[31]+4*y*(-1+z)*(-1+y)*p[34]+4*z*(-1+z)*(-1+y)*p[37]+4*y*(-1+z)*(-1+x)*p[40]+4*x*y*(-1+z)*p[40]-4*y*z*(-1+z)*p[43]+4*y*z*(-1+z)*p[46]+4*z*(-1+y)*(-1+x)*p[49]+4*x*z*(-1+y)*p[49]+4*y*z*(-1+y)*p[52]-4*y*z*(-1+y)*p[55]-4*y*z*(-1+x)*p[58]-4*x*y*z*p[58];

  dtau[1][1]=2*(-1+z)*(-1+y)*(-1+x)*p[1]-2*x*(-1+z)*(-1+y)*p[4]-2*x*y*(-1+z)*p[7]+2*y*(-1+z)*(-1+x)*p[10]-2*z*(-1+y)*(-1+x)*p[13]+2*x*z*(-1+y)*p[16]+2*x*y*z*p[19]-2*y*z*(-1+x)*p[22]+(-1+z)*(-1+x)*(2*x+2*y+2*z-1)*p[1]+x*(-1+z)*(-2*y-2*z+2*x-1)*p[4]-x*(-1+z)*(2*y-2*z-3+2*x)*p[7]-(-1+z)*(-1+x)*(2*x+2*z+1-2*y)*p[10]-z*(-1+x)*(2*x+2*y-2*z+1)*p[13]-x*z*(-2*y+2*z-3+2*x)*p[16]+x*z*(2*y+2*z-5+2*x)*p[19]+z*(-1+x)*(2*x-2*z+3-2*y)*p[22]-4*x*(-1+z)*(-1+x)*p[25]-4*(-1+z)*(-1+y)*(-1+x)*p[28]-4*y*(-1+z)*(-1+x)*p[28]-4*z*(-1+z)*(-1+x)*p[31]+4*x*(-1+z)*(-1+y)*p[34]+4*x*y*(-1+z)*p[34]+4*x*z*(-1+z)*p[37]+4*x*(-1+z)*(-1+x)*p[40]-4*x*z*(-1+z)*p[43]+4*z*(-1+z)*(-1+x)*p[46]+4*x*z*(-1+x)*p[49]+4*z*(-1+y)*(-1+x)*p[52]+4*y*z*(-1+x)*p[52]-4*x*z*(-1+y)*p[55]-4*x*y*z*p[55]-4*x*z*(-1+x)*p[58];

  dtau[1][2]=2*(-1+z)*(-1+y)*(-1+x)*p[1]-2*x*(-1+z)*(-1+y)*p[4]+2*x*y*(-1+z)*p[7]-2*y*(-1+z)*(-1+x)*p[10]+2*z*(-1+y)*(-1+x)*p[13]-2*x*z*(-1+y)*p[16]+2*x*y*z*p[19]-2*y*z*(-1+x)*p[22]+(-1+y)*(-1+x)*(2*x+2*y+2*z-1)*p[1]+x*(-1+y)*(-2*y-2*z+2*x-1)*p[4]-x*y*(2*y-2*z-3+2*x)*p[7]-y*(-1+x)*(2*x+2*z+1-2*y)*p[10]-(-1+y)*(-1+x)*(2*x+2*y-2*z+1)*p[13]-x*(-1+y)*(-2*y+2*z-3+2*x)*p[16]+x*y*(2*y+2*z-5+2*x)*p[19]+y*(-1+x)*(2*x-2*z+3-2*y)*p[22]-4*x*(-1+y)*(-1+x)*p[25]-4*y*(-1+y)*(-1+x)*p[28]-4*(-1+z)*(-1+y)*(-1+x)*p[31]-4*z*(-1+y)*(-1+x)*p[31]+4*x*y*(-1+y)*p[34]+4*x*(-1+z)*(-1+y)*p[37]+4*x*z*(-1+y)*p[37]+4*x*y*(-1+x)*p[40]-4*x*y*(-1+z)*p[43]-4*x*y*z*p[43]+4*y*(-1+z)*(-1+x)*p[46]+4*y*z*(-1+x)*p[46]+4*x*(-1+y)*(-1+x)*p[49]+4*y*(-1+y)*(-1+x)*p[52]-4*x*y*(-1+y)*p[55]-4*x*y*(-1+x)*p[58];

  dtau[2][0]=2*(-1+z)*(-1+y)*(-1+x)*p[2]+2*x*(-1+z)*(-1+y)*p[5]-2*x*y*(-1+z)*p[8]-2*y*(-1+z)*(-1+x)*p[11]-4*y*(-1+z)*(-1+y)*p[29]+(-1+z)*(-1+y)*(2*x+2*y+2*z-1)*p[2]+(-1+z)*(-1+y)*(-2*y-2*z+2*x-1)*p[5]-y*(-1+z)*(2*y-2*z-3+2*x)*p[8]-y*(-1+z)*(2*x+2*z+1-2*y)*p[11]-z*(-1+y)*(2*x+2*y-2*z+1)*p[14]-z*(-1+y)*(-2*y+2*z-3+2*x)*p[17]+y*z*(2*y+2*z-5+2*x)*p[20]+y*z*(2*x-2*z+3-2*y)*p[23]-4*(-1+z)*(-1+y)*(-1+x)*p[26]-4*x*(-1+z)*(-1+y)*p[26]-2*z*(-1+y)*(-1+x)*p[14]-2*x*z*(-1+y)*p[17]+2*x*y*z*p[20]+2*y*z*(-1+x)*p[23]-4*z*(-1+z)*(-1+y)*p[32]+4*y*(-1+z)*(-1+y)*p[35]+4*z*(-1+z)*(-1+y)*p[38]+4*y*(-1+z)*(-1+x)*p[41]+4*x*y*(-1+z)*p[41]-4*y*z*(-1+z)*p[44]+4*y*z*(-1+z)*p[47]+4*z*(-1+y)*(-1+x)*p[50]+4*x*z*(-1+y)*p[50]+4*y*z*(-1+y)*p[53]-4*y*z*(-1+y)*p[56]-4*y*z*(-1+x)*p[59]-4*x*y*z*p[59];

  dtau[2][1]=2*(-1+z)*(-1+y)*(-1+x)*p[2]-2*x*(-1+z)*(-1+y)*p[5]-2*x*y*(-1+z)*p[8]+2*y*(-1+z)*(-1+x)*p[11]-2*z*(-1+y)*(-1+x)*p[14]+2*x*z*(-1+y)*p[17]+2*x*y*z*p[20]-2*y*z*(-1+x)*p[23]+(-1+z)*(-1+x)*(2*x+2*y+2*z-1)*p[2]+x*(-1+z)*(-2*y-2*z+2*x-1)*p[5]-x*(-1+z)*(2*y-2*z-3+2*x)*p[8]-(-1+z)*(-1+x)*(2*x+2*z+1-2*y)*p[11]-z*(-1+x)*(2*x+2*y-2*z+1)*p[14]-x*z*(-2*y+2*z-3+2*x)*p[17]+x*z*(2*y+2*z-5+2*x)*p[20]+z*(-1+x)*(2*x-2*z+3-2*y)*p[23]-4*x*(-1+z)*(-1+x)*p[26]-4*(-1+z)*(-1+y)*(-1+x)*p[29]-4*y*(-1+z)*(-1+x)*p[29]-4*z*(-1+z)*(-1+x)*p[32]+4*x*(-1+z)*(-1+y)*p[35]+4*x*y*(-1+z)*p[35]+4*x*z*(-1+z)*p[38]+4*x*(-1+z)*(-1+x)*p[41]-4*x*z*(-1+z)*p[44]+4*z*(-1+z)*(-1+x)*p[47]+4*x*z*(-1+x)*p[50]+4*z*(-1+y)*(-1+x)*p[53]+4*y*z*(-1+x)*p[53]-4*x*z*(-1+y)*p[56]-4*x*y*z*p[56]-4*x*z*(-1+x)*p[59];

  dtau[2][2]=2*(-1+z)*(-1+y)*(-1+x)*p[2]-2*x*(-1+z)*(-1+y)*p[5]+2*x*y*(-1+z)*p[8]-2*y*(-1+z)*(-1+x)*p[11]+2*z*(-1+y)*(-1+x)*p[14]-2*x*z*(-1+y)*p[17]+2*x*y*z*p[20]-2*y*z*(-1+x)*p[23]+(-1+y)*(-1+x)*(2*x+2*y+2*z-1)*p[2]+x*(-1+y)*(-2*y-2*z+2*x-1)*p[5]-x*y*(2*y-2*z-3+2*x)*p[8]-y*(-1+x)*(2*x+2*z+1-2*y)*p[11]-(-1+y)*(-1+x)*(2*x+2*y-2*z+1)*p[14]-x*(-1+y)*(-2*y+2*z-3+2*x)*p[17]+x*y*(2*y+2*z-5+2*x)*p[20]+y*(-1+x)*(2*x-2*z+3-2*y)*p[23]-4*x*(-1+y)*(-1+x)*p[26]-4*y*(-1+y)*(-1+x)*p[29]-4*(-1+z)*(-1+y)*(-1+x)*p[32]-4*z*(-1+y)*(-1+x)*p[32]+4*x*y*(-1+y)*p[35]+4*x*(-1+z)*(-1+y)*p[38]+4*x*z*(-1+y)*p[38]+4*x*y*(-1+x)*p[41]-4*x*y*(-1+z)*p[44]-4*x*y*z*p[44]+4*y*(-1+z)*(-1+x)*p[47]+4*y*z*(-1+x)*p[47]+4*x*(-1+y)*(-1+x)*p[50]+4*y*(-1+y)*(-1+x)*p[53]-4*x*y*(-1+y)*p[56]-4*x*y*(-1+x)*p[59];

}

void Phy2Ref(__constant real *physnode, real *xphy, real *xref) 
{
#define ITERNEWTON 10
  real dxref[3], dxphy[3];
  xref[0] = 0.5;
  xref[1] = 0.5;
  xref[2] = 0.5;
  for(int iter = 0; iter < ITERNEWTON; ++iter ) {
    real dtau[3][3];
    real codtau[3][3];
#if 0
    int ifa =- 1;
    Ref2Phy(physnode, xref, 0, ifa, dxphy, dtau, codtau, 0, 0);
#else
    real gradphi[20][4];
    compute_gradphi(xref, gradphi);
    compute_xphy(physnode, gradphi, dxphy);
    compute_dtau(physnode, gradphi, dtau);
    compute_codtau(dtau, codtau);
#endif
    dxphy[0] -= xphy[0];
    dxphy[1] -= xphy[1];
    dxphy[2] -= xphy[2];
    real overdet = 1.0 / (  dtau[0][0] * codtau[0][0]
			    + dtau[0][1] * codtau[0][1]
			    + dtau[0][2] * codtau[0][2] );
    for(int ii = 0; ii < 3; ++ii) {
      dxref[ii] 
	= codtau[0][ii] * dxphy[0] 
	+ codtau[1][ii] * dxphy[1] 
	+ codtau[2][ii] * dxphy[2];
      xref[ii] -= dxref[ii] * overdet;
    }
  }
}

// From a reference point find the nearest gauss point
// Warning: works only  degree 1, 2, or 3 (FIXME: why?)
int ref_ipg(__constant int *param, real *xref) 
{
  // approximation degree in each direction
  int deg[3] = {param[0], param[1], param[2]};

  // number of subcells in each direction
  int nraf[3] = {param[3], param[4], param[5]};

  real hh[3] = {1.0 / nraf[0], 1.0 / nraf[1], 1.0 / nraf[2]};

  // get the subcell id
  int ncx = floor(xref[0] * nraf[0]);
  int ncy = floor(xref[1] * nraf[1]);
  int ncz = floor(xref[2] * nraf[2]);

  //printf("x=%f ncx=%d nrafx=%d\n",xref[0], ncx,nraf[0]);
  //printf("y=%f ncy=%d nrafy=%d\n",xref[1], ncy,nraf[1]);
  //printf("z=%f ncz=%d nrafz=%d\n",xref[2], ncz,nraf[2]);
  //assert(ncx >=0 && ncx<nraf[0]);
  //assert(ncy >=0 && ncy<nraf[1]);
  //assert(ncz >=0 && ncz<nraf[2]);

  // subcell index in the macrocell
  int nc = ncx + nraf[0] * (ncy + nraf[1] * ncz);
  int offset = (deg[0] + 1) * (deg[1] + 1) * (deg[2] + 1) * nc;

  // round to the nearest integer
  int ix = floor((xref[0] - ncx * hh[0]) / hh[0] * deg[0] + 0.5);
  int iy = floor((xref[1] - ncy * hh[1]) / hh[1] * deg[1] + 0.5);
  int iz = floor((xref[2] - ncz * hh[2]) / hh[2] * deg[2] + 0.5);
  //int iz=floor(xref[2]*deg[2]+0.5);

  //printf("xref %f %f %f ix=%d iy=%d iz=%d\n",
  //	 xref[0],xref[1],xref[2], ix, iy, iz);

  return ix + (deg[0] + 1) * (iy + (deg[1] + 1) * iz) + offset;
}

#ifndef _SOURCE_FUNC
#define _SOURCE_FUNC ZeroSource
#endif

// Sample source function
void ZeroSource(const real *x, const real t, const real *w, real *source,
		int m)
{
  for(int i = 0; i < m; ++i) 
    source[i] = 0.0;
}

// Sample source function
void OneSource(const real *x, const real t, const real *w, real *source,
	       int m)
{
  for(int i = 0; i < m; ++i) 
    source[i] = 1.0;
}

void compute_source(__constant int *param,     // interp param
		    __constant real *physnode, // macrocell nodes
		    const  __global real *mass,   // collocation point weights
		    const real tnow,           // the current time
		    __local real *wnloc,       // cache for wn
		    __local real *dtwnloc      // cache for dtwn
		    )
{
  // TODO: if there's enough space, make mass __constant

  const int m = param[0];
  const int deg[3] = {param[1], param[2], param[3]};
  const int nraf[3] = {param[4], param[5], param[6]};

  int ipg = get_global_id(0);
  int ipgloc = get_local_id(0);
  
  // Compute xref
  real xref[3];
  real wpg;
  ref_pg_vol(deg, nraf, ipg, xref, &wpg);

  // Compute xphy
  real xphy[3];
  real gradphi[20][4];
  compute_gradphi(xref, gradphi);
  compute_xphy(physnode, gradphi, xphy);
  
  // Copy w
  real w[_M];
  __local real *wnloc0 = wnloc + ipgloc * m;
  for(int iv = 0; iv < m; iv++) {
    w[iv] = wnloc0[iv];
  }

  // Compute source using w and xref
  real source[_M];
  _SOURCE_FUNC(xphy, tnow, w, source, _M);

  // The mass point is in __global memory, so get the global id.
  real massipg = mass[ipg];

  // Add the source buffer to dtw
  __local real *dtwnloc0 =  dtwnloc + ipgloc * m;;
  for(int iv = 0; iv < m; iv++) {
    dtwnloc0[iv] = source[iv] * massipg;
  }
}

// Compute the source terms inside  one macrocell
__kernel
void DGSource(__constant int *param,     // interp param
	      __constant real *physnode, // macrocell nodes
	      const __global real *mass, // collocation point weights
	      const real tnow,           // the current time
              __global real *wn,         // field values
	      __global real *dtwn,       // time derivative
	      __local real *wnloc,       // cache for wn
	      __local real *dtwnloc      // cache for dtwn
	      )
{
  prefetch_macrocell(wn, wnloc, param);
  prefetch_macrocell(dtwn, dtwnloc, param);
  
  barrier(CLK_LOCAL_MEM_FENCE);

  compute_source(param, physnode, mass, tnow, wnloc, dtwnloc);
  
  barrier(CLK_LOCAL_MEM_FENCE);

  postfetch_macrocell(dtwnloc, dtwn, param);
}

// Out-of-place RK stage
__kernel
void RK_out_CL(__global real *wnp1, 
	       __global const real *wn, 
	       __global const real *dtwn, 
	       const real dt)
{
  int ipg = get_global_id(0);
  wnp1[ipg] = wn[ipg] + dt * dtwn[ipg];
}

// In-place RK stage
__kernel
void RK_in_CL(__global real *wnp1, 
	      __global real *dtwn, 
	      const real dt)
{
  int ipg = get_global_id(0);
  wnp1[ipg] += dt * dtwn[ipg];
}


// Out-of-place RK stage
__kernel
void RK4_first_stages(__global real *wnp1, 
		     __global const real *wn, 
		     __global const real *dtwn, 
		      const real dt)
{
  int ipg = get_global_id(0);
  wnp1[ipg] = wn[ipg] + dt * dtwn[ipg];
}

// RK4 final stage
__kernel
void RK4_final_stage(__global real *w,
		     __global real *l1,
		     __global real *l2,
		     __global real *l3,
		     __global real *dtw, 
		     const real dt)
{
  const real b = -1.0 / 3.0;
  const real a[] = {1.0 / 3.0, 2.0 / 3.0, 1.0 / 3.0, dt / 6.0};
  int i = get_global_id(0);
  w[i] = 
    b * w[i] +
    a[0] * l1[i] +
    a[1] * l2[i] +
    a[2] * l3[i] +
    a[3] * dtw[i];
}
