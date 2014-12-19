#ifndef _MODEL_H
#define _MODEL_H

//! \brief a unified framework for all physical models
typedef struct Model {
  //! Number of conservative variables
  int m;

  //! CFL coefficient for time-stepping
  double cfl; // 0.05

  //! Number of conservative variables in each dimension (NB: their
  //! prouduct must equal m).
  int vlasov_mx, vlasov_my, vlasov_mz;

  //! The conservative variables have velocity in [-vlasov_vmax, vlasov_vmax]
  double vlasov_vmax;

  //! \brief A pointer to the numflux function
  //! \param[in] wL, wR : left and right states
  //! \param[in] vn : normal vector
  //! \param[out] flux : the flux
  void (*NumFlux)(double wL[], double wR[], double vn[3], double flux[]);

  //! \briefA pointer to the boundary flux function
  //! \param[in] x : space position
  //! \param[in] t : time
  //! \param[in] wL : left state
  //! \param[in] vn : normal vector
  //! \param[out] flux : the flux
  void (*BoundaryFlux)(double x[3], double t, double wL[], double vn[3],
		       double flux[]);

  //! \brief A pointer to the init data function
  // !\param[in] x : space position
  //! \param[out] w : init state at point x
  void (*InitData)(double x[3], double w[]);

  //! \brief A pointer to the imposed data function
  //!\param[in] x, t : space and time position
  //! \param[out] w : imposed state at point x and time t
  void (*ImposedData)(double x[3], double t, double w[]);

} Model;

//! \brief The particular flux for the transport model
//! \param[in] wL, wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void TransNumFlux(double wL[], double wR[], double vn[3], double* flux);

//! \brief The particular flux for the 2d transport model
//! \param[in] wL, wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void TransNumFlux2d(double wL[], double wR[], double vn[3], double* flux);

void VecTransNumFlux2d(double wL[], double wR[], double vn[3], double* flux);

//! \brief The particular boundary flux for the transport model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void TransBoundaryFlux(double* x, double t, double* wL, double* vn,
		       double* flux);

//! \brief The particular boundary flux for the 2d transport model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void TransBoundaryFlux2d(double* x, double t, double* wL, double* vn,
			 double* flux);

//! \brief The particular init data for the transport model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void TransInitData(double* x, double* w);

//! \brief The particular init data for the 2d transport model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void TransInitData2d(double *x, double *w);

void VecTransInitData2d(double *x, double *w);

void vTransImposedData2d(double *x, double t, double *w);

//! \brief The particular imposed data for the transport model
//! \param[in] x, t : space and time position
//! \param[out] w : imposed state at point x and time t
void TransImposedData(double* x, double t, double* w);

//! \brief The particular imposed data for the 2d transport model
//! \param[in] x, t : space and time position
//! \param[out] w : imposed state at point x and time t
void TransImposedData2d(double* x, double t, double* w);

void VecTransImposedData2d(double* x, double t, double* w);

//! \brief The particular flux for testing the transport model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void TestTransBoundaryFlux(double* x, double t, double* wL, double* vn,
			   double* flux);

//! \brief The particular flux for testing the 2d transport model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void TestTransBoundaryFlux2d(double* x, double t, double* wL, double* vn,
			     double* flux);

void VecTransBoundaryFlux2d(double* x, double t, double* wL, double* vn,
			    double* flux);

//! \brief The particular init data for the transport model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void TestTransInitData(double* x, double* w);

//! \brief The particular init data for the 2d transport model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void TestTransInitData2d(double* x, double* w);

//! \brief The particular imposed data for the transport model
//! \param[in] x, t : space and time position
//! \param[out] w : imposed state at point x and time t
void TestTransImposedData(double* x, double t, double* w);

//! \brief The particular imposed data for the 2d transport model
//! \param[in] x, t : space and time position
//! \param[out] w : imposed state at point x and time t
void TestTransImposedData2d(double* x, double t, double* w);

//! Parameters used for computing the velocity for the Vlasov equation.
static int m;
static int vlasov_mx;
static int vlasov_my;
static int vlasov_mz;
static double vlasov_vmax;

void set_vlasov_params(Model *mod);
double vlasov_vel(const int id, const int md, double vlasov_vmax);
void vlaTransInitData2d(double x[3], double w[]);
void vlaTransNumFlux2d(double wL[], double wR[], double* vnorm, double* flux);
void vlaTransBoundaryFlux2d(double x[3], double t, 
			    double wL[], double* vnorm,
			    double* flux);
void vlaTransImposedData2d(double x[3], double t, double* w);
#endif


// Format pour un model TRANSPORT
#pragma start_TRANSPORT

#pragma start_opencl
//! \brief The particular flux for the transport model
//! \param[in] wL, wR : left and right states
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void TRANSPORTNumFlux(double wL[], double wR[], double vn[3], double* flux);
#pragma end_opencl

#pragma start_opencl
//! \brief The particular flux for testing the transport model
//! \param[in] x : space position
//! \param[in] t : time
//! \param[in] wL : left state
//! \param[in] vn : normal vector
//! \param[out] flux : the flux
void TRANSPORTBoundaryFlux(double* x, double t, double* wL, double* vn,
			   double* flux);
#pragma end_opencl

//! \brief The particular init data for the transport model
//! \param[in] x : space position
//! \param[out] w : init state at point x
void TRANSPORTInitData(double* x, double* w);

#pragma start_opencl
//! \brief The particular imposed data for the transport model
//! \param[in] x, t : space and time position
//! \param[out] w : imposed state at point x and time t
void TRANSPORTImposedData(double* x, double t, double* w);
#pragma end_opencl

#pragma end_TRANSPORT
