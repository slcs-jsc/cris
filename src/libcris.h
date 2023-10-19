#include <netcdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_spline.h>
#include "jurassic.h"

/* ------------------------------------------------------------
   Dimensions...
   ------------------------------------------------------------ */

/*! Maximum number of data sets per granule. */
#define NDS 200000

/*! Maximum number of data points per granule. */
#define NPG 30

/*! Along-track size of CrIS radiance granule. */
#define L1_NTRACK 45

/*! Across-track size of CrIS radiance granule. */
#define L1_NXTRACK 30

/*! Number of field of views of CrIS radiance granule. */
#define L1_NFOV 9

/*! Number of CrIS longwave radiance channels. */
#define L1_NCHAN_LW 717

/*! Number of CrIS midwave radiance channels. */
#define L1_NCHAN_MW 869

/*! Number of CrIS shortwave radiance channels. */
#define L1_NCHAN_SW 637

/*! Along-track size of perturbation data. */
#define PERT_NTRACK 44000

/*! Across-track size of perturbation data. */
#define PERT_NXTRACK 120

/*! Number of field of views of perturbation data. */
#define PERT_NFOV 9

/*! Across-track size of wave analysis data. */
#define WX 300

/*! Along-track size of wave analysis data. */
#define WY 33000

/*! Maximum number of data points for spectral analysis. */
#define PMAX 512

/* ------------------------------------------------------------
   Macros...
   ------------------------------------------------------------ */

/*! Execute netCDF library command and check result. */
#define NC(cmd) {				     \
  int nc_result=(cmd);				     \
  if(nc_result!=NC_NOERR)			     \
    ERRMSG("%s", nc_strerror(nc_result));	     \
}

/* ------------------------------------------------------------
   Structs...
   ------------------------------------------------------------ */

/*! CrIS Level-1 data. */
typedef struct {

  /*! Time (seconds since 2000-01-01T00:00Z). */
  double time[L1_NTRACK][L1_NXTRACK];

  /*! Footprint longitude [deg]. */
  double lon[L1_NTRACK][L1_NXTRACK][L1_NFOV];

  /*! Footprint latitude [deg]. */
  double lat[L1_NTRACK][L1_NXTRACK][L1_NFOV];

  /*! Satellite altitude [km]. */
  double sat_z[L1_NTRACK];

  /*! Satellite longitude [deg]. */
  double sat_lon[L1_NTRACK];

  /*! Satellite latitude [deg]. */
  double sat_lat[L1_NTRACK];

  /*! Longwave channel frequencies [cm^-1]. */
  double nu_lw[L1_NCHAN_LW];

  /*! Longwave radiance [W/(m^2 sr cm^-1)]. */
  float rad_lw[L1_NTRACK][L1_NXTRACK][L1_NFOV][L1_NCHAN_LW];

  /*! Longwave radiance noise [W/(m^2 sr cm^-1)]. */
  float nedn_lw[L1_NFOV][L1_NCHAN_LW];

  /*! Longwave quality flag (0=best, 1=good, 2=don't use). */
  short qual_lw[L1_NTRACK][L1_NXTRACK][L1_NFOV];

  /*! Midwave channel frequencies [cm^-1]. */
  double nu_mw[L1_NCHAN_MW];

  /*! Midwave radiance [W/(m^2 sr cm^-1)]. */
  float rad_mw[L1_NTRACK][L1_NXTRACK][L1_NFOV][L1_NCHAN_MW];

  /*! Midwave radiance noise [W/(m^2 sr cm^-1)]. */
  float nedn_mw[L1_NFOV][L1_NCHAN_MW];

  /*! Midwave quality flag (0=best, 1=good, 2=don't use). */
  short qual_mw[L1_NTRACK][L1_NXTRACK][L1_NFOV];

  /*! Shortwave channel frequencies [cm^-1]. */
  double nu_sw[L1_NCHAN_SW];

  /*! Shortwave radiance [W/(m^2 sr cm^-1)]. */
  float rad_sw[L1_NTRACK][L1_NXTRACK][L1_NFOV][L1_NCHAN_SW];

  /*! Longwave radiance noise [W/(m^2 sr cm^-1)]. */
  float nedn_sw[L1_NFOV][L1_NCHAN_SW];

  /*! Shortwave quality flag (0=best, 1=good, 2=don't use). */
  short qual_sw[L1_NTRACK][L1_NXTRACK][L1_NFOV];

} cris_l1_t;

/*! Perturbation data. */
typedef struct {

  /*! Number of along-track values. */
  int ntrack;

  /*! Number of across-track values. */
  int nxtrack;

  /*! Number of field of views. */
  int nfov;

  /*! Time (seconds since 2000-01-01T00:00Z). */
  double time[PERT_NTRACK][PERT_NXTRACK][PERT_NFOV];

  /*! Longitude [deg]. */
  double lon[PERT_NTRACK][PERT_NXTRACK][PERT_NFOV];

  /*! Latitude [deg]. */
  double lat[PERT_NTRACK][PERT_NXTRACK][PERT_NFOV];

  /*! Brightness temperature (8.1 micron) [K]. */
  double dc[PERT_NTRACK][PERT_NXTRACK][PERT_NFOV];

  /*! Brightness temperature (4.3 or 15 micron) [K]. */
  double bt[PERT_NTRACK][PERT_NXTRACK][PERT_NFOV];

  /*! Brightness temperature perturbation (4 or 15 micron) [K]. */
  double pt[PERT_NTRACK][PERT_NXTRACK][PERT_NFOV];

  /*! Brightness temperature variance (4 or 15 micron) [K]. */
  double var[PERT_NTRACK][PERT_NXTRACK][PERT_NFOV];

} pert_t;

/*! Retrieval results. */
typedef struct {

  /*! Number of data sets. */
  int nds;

  /*! Number of data points. */
  int np;

  /*! Time (seconds since 2000-01-01T00:00Z). */
  double time[NDS][NPG];

  /*! Altitude [km]. */
  double z[NDS][NPG];

  /*! Longitude [deg]. */
  double lon[NDS][NPG];

  /*! Latitude [deg]. */
  double lat[NDS][NPG];

  /*! Pressure [hPa]. */
  double p[NDS][NPG];

  /*! Temperature [K]. */
  double t[NDS][NPG];

  /*! Temperature (a priori data) [K]. */
  double t_apr[NDS][NPG];

  /*! Temperature (total error) [K]. */
  double t_tot[NDS][NPG];

  /*! Temperature (noise error) [K]. */
  double t_noise[NDS][NPG];

  /*! Temperature (forward model error) [K]. */
  double t_fm[NDS][NPG];

  /*! Temperature (measurement content). */
  double t_cont[NDS][NPG];

  /*! Temperature (resolution). */
  double t_res[NDS][NPG];

  /*! Chi^2. */
  double chisq[NDS];

} ret_t;

/*! Wave analysis data. */
typedef struct {

  /*! Number of across-track values. */
  int nx;

  /*! Number of along-track values. */
  int ny;

  /*! Time (seconds since 2000-01-01T00:00Z). */
  double time;

  /*! Altitude [km]. */
  double z;

  /*! Longitude [deg]. */
  double lon[WX][WY];

  /*! Latitude [deg]. */
  double lat[WX][WY];

  /*! Across-track distance [km]. */
  double x[WX];

  /*! Along-track distance [km]. */
  double y[WY];

  /*! Temperature [K]. */
  double temp[WX][WY];

  /*! Background [K]. */
  double bg[WX][WY];

  /*! Perturbation [K]. */
  double pt[WX][WY];

  /*! Fit [K]. */
  double fit[WX][WY];

  /*! Variance [K]. */
  double var[WX][WY];

} wave_t;

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/*! Add variable attributes to netCDF file. */
void add_att(
  int ncid,
  int varid,
  const char *unit,
  const char *long_name);

/*! Add variable to netCDF file. */
void add_var(
  int ncid,
  const char *varname,
  const char *unit,
  const char *longname,
  int type,
  int dimid[],
  int *varid,
  int ndims);

/*! Get background based on polynomial fits. */
void background_poly(
  wave_t * wave,
  int dim_x,
  int dim_y);

/*! Get background based on polynomial fits. */
void background_poly_help(
  double *xx,
  double *yy,
  int n,
  int dim);

/*! Smooth background. */
void background_smooth(
  wave_t * wave,
  int npts_x,
  int npts_y);

/*! Set background... */
void create_background(
  wave_t * wave);

/*! Add noise to perturbations and temperatures... */
void create_noise(
  wave_t * wave,
  double nedt);

/*! Add linear wave pattern... */
void create_wave(
  wave_t * wave,
  double amp,
  double lx,
  double ly,
  double phi,
  double fwhm);

/*! Get day of year from date. */
void day2doy(
  int year,
  int mon,
  int day,
  int *doy);

/*! Get date from day of year. */
void doy2day(
  int year,
  int doy,
  int *mon,
  int *day);

/*! Evaluate wave fit... */
void fit_wave(
  wave_t * wave,
  double amp,
  double phi,
  double kx,
  double ky,
  double *chisq);

/*! Calculate 1-D FFT... */
void fft_help(
  double *fcReal,
  double *fcImag,
  int n);

/*! Calculate 2-D FFT... */
void fft(
  wave_t * wave,
  double *Amax,
  double *phimax,
  double *lhmax,
  double *kxmax,
  double *kymax,
  double *alphamax,
  double *betamax,
  char *filename);

/*! Apply Gaussian filter to perturbations... */
void gauss(
  wave_t * wave,
  double fwhm);

/*! Apply Hamming filter to perturbations... */
void hamming(
  wave_t * wave,
  int nit);

/*! Interpolate to regular grid in x-direction. */
void intpol_x(
  wave_t * wave,
  int n);

/*! Apply median filter to perturbations... */
void median(
  wave_t * wave,
  int dx);

/*! Merge wave structs in y-direction. */
void merge_y(
  wave_t * wave1,
  wave_t * wave2);

/*! Estimate noise. */
void noise(
  wave_t * wave,
  double *mu,
  double *sig);

/*! Estimate noise from perurbations. */
void noise_pert(
  pert_t * pert,
  int track0,
  int track1,
  double *mu,
  double *sig);

/*! Compute periodogram. */
void period(
  wave_t * wave,
  double lxymax,
  double dlxy,
  double *Amax,
  double *phimax,
  double *lhmax,
  double *kxmax,
  double *kymax,
  double *alphamax,
  double *betamax,
  char *filename);

/*! Convert radiance perturbation data to wave analysis struct. */
void pert2wave(
  pert_t * pert,
  wave_t * wave,
  int track0,
  int track1,
  int xtrack0,
  int xtrack1);

/*! Read CrIS Level-1 data. */
int read_cris_l1(
  char *filename,
  cris_l1_t * l1,
  int apo);

/*! Read radiance perturbation data. */
void read_pert(
  char *filename,
  char *pertname,
  pert_t * pert);

/*! Read CrIS retrieval data. */
void read_retr(
  char *filename,
  ret_t * ret);

/*! Convert array. */
void read_retr_help(
  double *help,
  int nds,
  int np,
  double mat[NDS][NPG]);

/*! Read wave analysis data. */
void read_wave(
  char *filename,
  wave_t * wave);

/*! Convert CrIS radiance data to wave analysis struct. */
void rad2wave(
  cris_l1_t * cris_l1,
  double *nu,
  int nd,
  wave_t * wave);

/*! Convert CrIS retrieval results to wave analysis struct. */
void ret2wave(
  ret_t * ret,
  wave_t * wave,
  int dataset,
  int ip);

/*! Compute local variance. */
void variance(
  wave_t * wave,
  double dh);

/*! Write wave analysis data. */
void write_wave(
  char *filename,
  wave_t * wave);
