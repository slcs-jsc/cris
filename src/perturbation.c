#include "libcris.h"

/* ------------------------------------------------------------
   Constants...
   ------------------------------------------------------------ */

/*! Number of 4.3 micron channels: */
#define N4 54

/*! Number of 15 micron low channels: */
#define N15_LOW 15

/*! Number of 15 micron high channels: */
#define N15_HIGH 1

/*! Number of 8.1 micron channels: */
#define N8 1

/*! Number of 10.4 micron channels: */
#define N10 1

/* ------------------------------------------------------------
   Macros...
   ------------------------------------------------------------ */

/*! Loop over all footprints. */
#define LOOP_ALL(ntrack)				\
  for (int track = 0; track < ntrack; track++)		\
    for (int xtrack = 0; xtrack < L1_NXTRACK; xtrack++)	\
      for (int ifov = 0; ifov < L1_NFOV; ifov++)

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  static cris_l1_t l1;

  static pert_t *pert_4mu, *pert_15mu_low, *pert_15mu_high;

  const double var_dh2 = 100. * 100.;

  static double help[PERT_NTRACK * PERT_NXTRACK * PERT_NFOV], rad, nu,
    x[L1_NXTRACK], y[L1_NXTRACK], x0[3], x1[3];

  const int list_4mu[N4]
    = { 270, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283,
    284, 285, 286, 288, 289, 291, 292, 294, 295, 296, 297, 298, 299, 300,
    302, 303, 304, 305, 306, 307, 318, 319, 321, 322, 323, 324, 325, 326,
    328, 330, 332, 333, 334, 335, 336, 337, 338, 339, 340, 341
  };

  const int list_15mu_low[N15_LOW]
  = { 5, 8, 10, 13, 15, 17, 23, 34, 36, 39, 41, 44, 46, 49, 51 };

  const int list_15mu_high[N15_HIGH] = { 31 };

  const int list_8mu[N8] = { 36 };

  const int list_10mu[N10] = { 502 };

  const int dtrack = 3, dxtrack = 3;

  static int dimid[3], n, ncid, track0,
    time_varid, lon_varid, lat_varid, bt_4mu_varid, init,
    bt_4mu_pt_varid, bt_4mu_var_varid, bt_8mu_varid, bt_10mu_varid,
    bt_15mu_low_varid, bt_15mu_low_pt_varid, bt_15mu_low_var_varid,
    bt_15mu_high_varid, bt_15mu_high_pt_varid, bt_15mu_high_var_varid;

  static double count_4mu[PERT_NXTRACK][PERT_NFOV],
    count_15mu_high[PERT_NXTRACK][PERT_NFOV],
    count_15mu_low[PERT_NXTRACK][PERT_NFOV];

  static double bias_4mu[PERT_NXTRACK][PERT_NFOV],
    bias_15mu_high[PERT_NXTRACK][PERT_NFOV],
    bias_15mu_low[PERT_NXTRACK][PERT_NFOV];

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <pert.nc> <l1b_file1> [<l1b_file2> ...]");

  /* Get control parameters... */
  const int apo = (int) scan_ctl(argc, argv, "APO", -1, "0", NULL);
  const int bias = (int) scan_ctl(argc, argv, "BIAS", -1, "0", NULL);

  /* Allocate... */
  ALLOC(pert_4mu, pert_t, 1);
  ALLOC(pert_15mu_low, pert_t, 1);
  ALLOC(pert_15mu_high, pert_t, 1);

  /* ------------------------------------------------------------
     Read CrIS Level-1B files...
     ------------------------------------------------------------ */

  /* Loop over files... */
  for (int iarg = 3; iarg < argc; iarg++) {

    /* Read CrIS data... */
    if (!read_cris_l1(argv[iarg], &l1, apo))
      continue;

    /* Write channel information... */
    if (!init) {
      init = 1;

      /* 4.3 micron channels... */
      LOG(2, "4.3 micron channels:");
      n = 0;
      nu = 0;
      for (int i = 0; i < N4; i++) {
	nu += l1.nu_sw[list_4mu[i]];
	n++;
	LOG(2, "  count= %2d | channel index= SW_%03d | nu= %.4f cm^-1", n,
	    list_4mu[i], l1.nu_sw[list_4mu[i]]);
      }
      LOG(2, "  nu_mean= %.4f cm^-1", nu / n);

      /* 15 micron low channels... */
      LOG(2, "15 micron low channels:");
      n = 0;
      nu = 0;
      for (int i = 0; i < N15_LOW; i++) {
	nu += l1.nu_lw[list_15mu_low[i]];
	n++;
	LOG(2, "  count= %2d | channel index= LW_%03d | nu= %.4f cm^-1", n,
	    list_15mu_low[i], l1.nu_lw[list_15mu_low[i]]);
      }
      LOG(2, "  nu_mean= %.4f cm^-1", nu / n);

      /* 15 micron high channels... */
      LOG(2, "15 micron high channels:");
      n = 0;
      nu = 0;
      for (int i = 0; i < N15_HIGH; i++) {
	nu += l1.nu_lw[list_15mu_high[i]];
	n++;
	LOG(2, "  count= %2d | channel index= LW_%03d | nu= %.4f cm^-1", n,
	    list_15mu_high[i], l1.nu_lw[list_15mu_high[i]]);
      }
      LOG(2, "  nu_mean= %.4f cm^-1", nu / n);

      /* 8.1 micron channels... */
      LOG(2, "8.1 micron channels:");
      n = 0;
      nu = 0;
      for (int i = 0; i < N8; i++) {
	nu += l1.nu_mw[list_8mu[i]];
	n++;
	LOG(2, "  count= %2d | channel index= MW_%03d | nu= %.4f cm^-1", n,
	    list_8mu[i], l1.nu_mw[list_8mu[i]]);
      }
      LOG(2, "  nu_mean= %.4f cm^-1", nu / n);

      /* 10.4 micron channels... */
      LOG(2, "10.4 micron channels:");
      n = 0;
      nu = 0;
      for (int i = 0; i < N10; i++) {
	nu += l1.nu_lw[list_10mu[i]];
	n++;
	LOG(2, "  count= %2d | channel index= LW_%03d | nu= %.4f cm^-1", n,
	    list_10mu[i], l1.nu_lw[list_10mu[i]]);
      }
      LOG(2, "  nu_mean= %.4f cm^-1", nu / n);
    }

    /* Save geolocation... */
    pert_4mu->ntrack += L1_NTRACK;
    if (pert_4mu->ntrack > PERT_NTRACK)
      ERRMSG("Too many granules!");
    pert_4mu->nxtrack = L1_NXTRACK;
    if (pert_4mu->nxtrack > PERT_NXTRACK)
      ERRMSG("Too many tracks!");
    pert_4mu->nfov = L1_NFOV;
    if (pert_4mu->nfov > PERT_NFOV)
      ERRMSG("Too many field of views!");
    LOOP_ALL(L1_NTRACK) {
      pert_4mu->time[track0 + track][xtrack][ifov]
	= l1.time[track][xtrack] - 220838400.;
      pert_4mu->lon[track0 + track][xtrack][ifov]
	= l1.lon[track][xtrack][ifov];
      pert_4mu->lat[track0 + track][xtrack][ifov]
	= l1.lat[track][xtrack][ifov];
    }

    /* Get 8.1 micron brightness temperature... */
    LOOP_ALL(L1_NTRACK)
      pert_4mu->dc[track0 + track][xtrack][ifov]
      = BRIGHT(l1.rad_mw[track][xtrack][ifov][list_8mu[0]] * 1e-3,
	       l1.nu_mw[list_8mu[0]]);

    /* Get 10.4 micron brightness temperature... */
    LOOP_ALL(L1_NTRACK)
      pert_15mu_high->dc[track0 + track][xtrack][ifov]
      = BRIGHT(l1.rad_lw[track][xtrack][ifov][list_10mu[0]] * 1e-3,
	       l1.nu_lw[list_10mu[0]]);

    /* Get 4.3 micron brightness temperature... */
    LOOP_ALL(L1_NTRACK) {
      n = 0;
      nu = rad = 0;
      for (int i = 0; i < N4; i++)
	if (gsl_finite(l1.rad_sw[track][xtrack][ifov][list_4mu[i]])) {
	  rad += l1.rad_sw[track][xtrack][ifov][list_4mu[i]] * 1e-3;
	  nu += l1.nu_sw[list_4mu[i]];
	  n++;
	}
      if (n > 0.9 * N4)
	pert_4mu->bt[track0 + track][xtrack][ifov] = BRIGHT(rad / n, nu / n);
      else
	pert_4mu->bt[track0 + track][xtrack][ifov] = GSL_NAN;
    }

    /* Get 15 micron low brightness temperature... */
    LOOP_ALL(L1_NTRACK) {
      n = 0;
      nu = rad = 0;
      for (int i = 0; i < N15_LOW; i++)
	if (gsl_finite(l1.rad_lw[track][xtrack][ifov][list_15mu_low[i]])) {
	  rad += l1.rad_lw[track][xtrack][ifov][list_15mu_low[i]] * 1e-3;
	  nu += l1.nu_lw[list_15mu_low[i]];
	  n++;
	}
      if (n > 0.9 * N15_LOW)
	pert_15mu_low->bt[track0 + track][xtrack][ifov] =
	  BRIGHT(rad / n, nu / n);
      else
	pert_15mu_low->bt[track0 + track][xtrack][ifov] = GSL_NAN;
    }

    /* Get 15 micron high brightness temperature... */
    LOOP_ALL(L1_NTRACK) {
      n = 0;
      nu = rad = 0;
      for (int i = 0; i < N15_HIGH; i++)
	if (gsl_finite(l1.rad_lw[track][xtrack][ifov][list_15mu_high[i]])) {
	  rad += l1.rad_lw[track][xtrack][ifov][list_15mu_high[i]] * 1e-3;
	  nu += l1.nu_lw[list_15mu_high[i]];
	  n++;
	}
      if (n > 0.9 * N15_HIGH)
	pert_15mu_high->bt[track0 + track][xtrack][ifov] =
	  BRIGHT(rad / n, nu / n);
      else
	pert_15mu_high->bt[track0 + track][xtrack][ifov] = GSL_NAN;
    }

    /* Increment track counter... */
    track0 += L1_NTRACK;
  }

  /* ------------------------------------------------------------
     Calculate perturbations...
     ------------------------------------------------------------ */

  /* Write info... */
  LOG(1, "Calculate perturbations...");

  /* Loop over scans and field of views... */
  for (int track = 0; track < pert_4mu->ntrack; track++)
    for (int ifov = 0; ifov < L1_NFOV; ifov++) {

      /* Skip scan edges... */
      pert_4mu->dc[track][0][ifov] = GSL_NAN;
      pert_4mu->dc[track][L1_NXTRACK - 1][ifov] = GSL_NAN;
      pert_15mu_high->dc[track][0][ifov] = GSL_NAN;
      pert_15mu_high->dc[track][L1_NXTRACK - 1][ifov] = GSL_NAN;
      pert_4mu->bt[track][0][ifov] = GSL_NAN;
      pert_4mu->bt[track][L1_NXTRACK - 1][ifov] = GSL_NAN;
      pert_15mu_low->bt[track][0][ifov] = GSL_NAN;
      pert_15mu_low->bt[track][L1_NXTRACK - 1][ifov] = GSL_NAN;
      pert_15mu_high->bt[track][0][ifov] = GSL_NAN;
      pert_15mu_high->bt[track][L1_NXTRACK - 1][ifov] = GSL_NAN;

      /* Get 4 micron perturbations... */
      for (int xtrack = 0; xtrack < L1_NXTRACK; xtrack++) {
	x[xtrack] = (double) xtrack;
	y[xtrack] = pert_4mu->bt[track][xtrack][ifov];
      }
      background_poly_help(x, y, L1_NXTRACK, 5);
      for (int xtrack = 0; xtrack < L1_NXTRACK; xtrack++)
	pert_4mu->pt[track][xtrack][ifov] =
	  pert_4mu->bt[track][xtrack][ifov] - y[xtrack];

      /* Get 15 micron low perturbations... */
      for (int xtrack = 0; xtrack < L1_NXTRACK; xtrack++) {
	x[xtrack] = (double) xtrack;
	y[xtrack] = pert_15mu_low->bt[track][xtrack][ifov];
      }
      background_poly_help(x, y, L1_NXTRACK, 7);
      for (int xtrack = 0; xtrack < L1_NXTRACK; xtrack++)
	pert_15mu_low->pt[track][xtrack][ifov] =
	  pert_15mu_low->bt[track][xtrack][ifov] - y[xtrack];

      /* Get 15 micron high perturbations... */
      for (int xtrack = 0; xtrack < L1_NXTRACK; xtrack++) {
	x[xtrack] = (double) xtrack;
	y[xtrack] = pert_15mu_high->bt[track][xtrack][ifov];
      }
      background_poly_help(x, y, L1_NXTRACK, 5);
      for (int xtrack = 0; xtrack < L1_NXTRACK; xtrack++)
	pert_15mu_high->pt[track][xtrack][ifov] =
	  pert_15mu_high->bt[track][xtrack][ifov] - y[xtrack];
    }

  /* ------------------------------------------------------------
     Bias correction...
     ------------------------------------------------------------ */

  /* Check flag... */
  if (bias) {

    /* Write info... */
    LOG(1, "Calculate bias correction...");

    /* Loop over scans and field of views... */
    for (int xtrack = 0; xtrack < L1_NXTRACK; xtrack++)
      for (int ifov = 0; ifov < L1_NFOV; ifov++) {

	for (int track = 0; track < pert_4mu->ntrack; track++)
	  if (isfinite(pert_4mu->pt[track][xtrack][ifov])) {
	    bias_4mu[xtrack][ifov] += pert_4mu->pt[track][xtrack][ifov];
	    count_4mu[xtrack][ifov]++;
	  }

	for (int track = 0; track < pert_15mu_low->ntrack; track++)
	  if (isfinite(pert_15mu_low->pt[track][xtrack][ifov])) {
	    bias_15mu_low[xtrack][ifov] +=
	      pert_15mu_low->pt[track][xtrack][ifov];
	    count_15mu_low[xtrack][ifov]++;
	  }

	for (int track = 0; track < pert_15mu_high->ntrack; track++)
	  if (isfinite(pert_15mu_high->pt[track][xtrack][ifov])) {
	    bias_15mu_high[xtrack][ifov] +=
	      pert_15mu_high->pt[track][xtrack][ifov];
	    count_15mu_high[xtrack][ifov]++;
	  }

	bias_4mu[xtrack][ifov] /= count_4mu[xtrack][ifov];
	bias_15mu_low[xtrack][ifov] /= count_15mu_low[xtrack][ifov];
	bias_15mu_high[xtrack][ifov] /= count_15mu_high[xtrack][ifov];
      }

    /* Write log messages... */
    LOG(2, "4.3 micron channels:");
    for (int ifov = 0; ifov < L1_NFOV; ifov++)
      for (int xtrack = 0; xtrack < L1_NXTRACK; xtrack++)
	LOG(2, "  xtrack= %d | ifov= %d | bias_4mu= %g K", xtrack, ifov,
	    bias_4mu[xtrack][ifov]);

    LOG(2, "15 micron low channels:");
    for (int ifov = 0; ifov < L1_NFOV; ifov++)
      for (int xtrack = 0; xtrack < L1_NXTRACK; xtrack++)
	LOG(2, "  xtrack= %d | ifov= %d | bias_15mu_low= %g K", xtrack, ifov,
	    bias_15mu_low[xtrack][ifov]);

    LOG(2, "15 micron high channels:");
    for (int ifov = 0; ifov < L1_NFOV; ifov++)
      for (int xtrack = 0; xtrack < L1_NXTRACK; xtrack++)
	LOG(2, "  xtrack= %d | ifov= %d | bias_15mu_high= %g K", xtrack, ifov,
	    bias_15mu_high[xtrack][ifov]);

    /* Subtract bias... */
    for (int track = 0; track < pert_4mu->ntrack; track++)
      for (int xtrack = 0; xtrack < L1_NXTRACK; xtrack++)
	for (int ifov = 0; ifov < L1_NFOV; ifov++) {
	  pert_4mu->pt[track][xtrack][ifov]
	    -= bias_4mu[xtrack][ifov];
	  pert_15mu_low->pt[track][xtrack][ifov]
	    -= bias_15mu_low[xtrack][ifov];
	  pert_15mu_high->pt[track][xtrack][ifov]
	    -= bias_15mu_high[xtrack][ifov];
	}
  }

  /* ------------------------------------------------------------
     Calculate variances...
     ------------------------------------------------------------ */

  /* Write info... */
  LOG(1, "Calculate variances...");

  /* Loop over footprints... */
  LOOP_ALL(pert_4mu->ntrack) {

    /* Get geolocation... */
    geo2cart(0.0, pert_4mu->lon[track][xtrack][ifov],
	     pert_4mu->lat[track][xtrack][ifov], x0);

    /* Init... */
    int n_4mu = 0, n_15mu_low = 0, n_15mu_high = 0;
    double mean_4mu = 0, mean_15mu_low = 0, mean_15mu_high = 0;
    double var_4mu = 0, var_15mu_low = 0, var_15mu_high = 0;

    /* Loop over neighbouring footprints... */
    for (int track2 = track - dtrack; track2 <= track + dtrack; track2++)
      for (int xtrack2 = xtrack - dxtrack; xtrack2 <= xtrack + dxtrack;
	   xtrack2++)
	for (int ifov2 = 0; ifov2 < L1_NFOV; ifov2++)
	  if (track2 >= 0 && track2 < pert_4mu->ntrack
	      && xtrack2 >= 0 && xtrack2 < L1_NXTRACK) {

	    /* Check horizontal distance... */
	    geo2cart(0.0, pert_4mu->lon[track2][xtrack2][ifov2],
		     pert_4mu->lat[track2][xtrack2][ifov2], x1);
	    if (DIST2(x0, x1) > var_dh2)
	      continue;

	    /* Calculate variance... */
	    if (gsl_finite(pert_4mu->pt[track2][xtrack2][ifov2])) {
	      mean_4mu += pert_4mu->pt[track2][xtrack2][ifov2];
	      var_4mu += gsl_pow_2(pert_4mu->pt[track2][xtrack2][ifov2]);
	      n_4mu++;
	    }

	    if (gsl_finite(pert_15mu_low->pt[track2][xtrack2][ifov2])) {
	      mean_15mu_low += pert_15mu_low->pt[track2][xtrack2][ifov2];
	      var_15mu_low +=
		gsl_pow_2(pert_15mu_low->pt[track2][xtrack2][ifov2]);
	      n_15mu_low++;
	    }

	    if (gsl_finite(pert_15mu_high->pt[track2][xtrack2][ifov2])) {
	      mean_15mu_high += pert_15mu_high->pt[track2][xtrack2][ifov2];
	      var_15mu_high +=
		gsl_pow_2(pert_15mu_high->pt[track2][xtrack2][ifov2]);
	      n_15mu_high++;
	    }
	  }

    /* Save variance data... */
    if (n_4mu > 0 && xtrack > 0 && xtrack < L1_NXTRACK - 1)
      pert_4mu->var[track][xtrack][ifov] =
	var_4mu / n_4mu - gsl_pow_2(mean_4mu / n_4mu);
    else
      pert_4mu->var[track][xtrack][ifov] = GSL_NAN;

    if (n_15mu_low > 0 && xtrack > 0 && xtrack < L1_NXTRACK - 1)
      pert_15mu_low->var[track][xtrack][ifov] =
	var_15mu_low / n_15mu_low - gsl_pow_2(mean_15mu_low / n_15mu_low);
    else
      pert_15mu_low->var[track][xtrack][ifov] = GSL_NAN;

    if (n_15mu_high > 0 && xtrack > 0 && xtrack < L1_NXTRACK - 1)
      pert_15mu_high->var[track][xtrack][ifov] =
	var_15mu_high / n_15mu_high - gsl_pow_2(mean_15mu_high / n_15mu_high);
    else
      pert_15mu_high->var[track][xtrack][ifov] = GSL_NAN;
  }

  /* ------------------------------------------------------------
     Write netCDF file...
     ------------------------------------------------------------ */

  /* Write info... */
  LOG(1, "Write perturbation data file: %s", argv[2]);

  /* Create netCDF file... */
  NC(nc_create(argv[2], NC_CLOBBER, &ncid));

  /* Set dimensions... */
  NC(nc_def_dim(ncid, "NTRACK", (size_t) pert_4mu->ntrack, &dimid[0]));
  NC(nc_def_dim(ncid, "NXTRACK", L1_NXTRACK, &dimid[1]));
  NC(nc_def_dim(ncid, "NFOV", L1_NFOV, &dimid[2]));

  /* Add variables... */
  NC(nc_def_var(ncid, "time", NC_DOUBLE, 3, dimid, &time_varid));
  add_att(ncid, time_varid, "s", "time (seconds since 2000-01-01T00:00Z)");
  NC(nc_def_var(ncid, "lon", NC_DOUBLE, 3, dimid, &lon_varid));
  add_att(ncid, lon_varid, "deg", "footprint longitude");
  NC(nc_def_var(ncid, "lat", NC_DOUBLE, 3, dimid, &lat_varid));
  add_att(ncid, lat_varid, "deg", "footprint latitude");

  NC(nc_def_var(ncid, "bt_8mu", NC_FLOAT, 3, dimid, &bt_8mu_varid));
  add_att(ncid, bt_8mu_varid, "K", "brightness temperature at 8.1 micron");

  NC(nc_def_var(ncid, "bt_10mu", NC_FLOAT, 3, dimid, &bt_10mu_varid));
  add_att(ncid, bt_10mu_varid, "K", "brightness temperature at 10.4 micron");

  NC(nc_def_var(ncid, "bt_4mu", NC_FLOAT, 3, dimid, &bt_4mu_varid));
  add_att(ncid, bt_4mu_varid, "K", "brightness temperature" " at 4.3 micron");
  NC(nc_def_var(ncid, "bt_4mu_pt", NC_FLOAT, 3, dimid, &bt_4mu_pt_varid));
  add_att(ncid, bt_4mu_pt_varid, "K", "brightness temperature perturbation"
	  " at 4.3 micron");
  NC(nc_def_var(ncid, "bt_4mu_var", NC_FLOAT, 3, dimid, &bt_4mu_var_varid));
  add_att(ncid, bt_4mu_var_varid, "K^2", "brightness temperature variance"
	  " at 4.3 micron");

  NC(nc_def_var(ncid, "bt_15mu_low", NC_FLOAT, 3, dimid, &bt_15mu_low_varid));
  add_att(ncid, bt_15mu_low_varid, "K", "brightness temperature"
	  " at 15 micron (low altitudes)");
  NC(nc_def_var(ncid, "bt_15mu_low_pt", NC_FLOAT, 3, dimid,
		&bt_15mu_low_pt_varid));
  add_att(ncid, bt_15mu_low_pt_varid, "K",
	  "brightness temperature perturbation"
	  " at 15 micron (low altitudes)");
  NC(nc_def_var
     (ncid, "bt_15mu_low_var", NC_FLOAT, 3, dimid, &bt_15mu_low_var_varid));
  add_att(ncid, bt_15mu_low_var_varid, "K^2",
	  "brightness temperature variance" " at 15 micron (low altitudes)");

  NC(nc_def_var(ncid, "bt_15mu_high", NC_FLOAT, 3, dimid,
		&bt_15mu_high_varid));
  add_att(ncid, bt_15mu_high_varid, "K", "brightness temperature"
	  " at 15 micron (high altitudes)");
  NC(nc_def_var(ncid, "bt_15mu_high_pt", NC_FLOAT, 3, dimid,
		&bt_15mu_high_pt_varid));
  add_att(ncid, bt_15mu_high_pt_varid, "K",
	  "brightness temperature perturbation"
	  " at 15 micron (high altitudes)");
  NC(nc_def_var
     (ncid, "bt_15mu_high_var", NC_FLOAT, 3, dimid, &bt_15mu_high_var_varid));
  add_att(ncid, bt_15mu_high_var_varid, "K^2",
	  "brightness temperature variance" " at 15 micron (high altitudes)");

  /* Leave define mode... */
  NC(nc_enddef(ncid));

  /* Write data... */
  n = 0;
  LOOP_ALL(pert_4mu->ntrack)
    help[n++] = pert_4mu->time[track][xtrack][ifov];
  NC(nc_put_var_double(ncid, time_varid, help));

  n = 0;
  LOOP_ALL(pert_4mu->ntrack)
    help[n++] = pert_4mu->lon[track][xtrack][ifov];
  NC(nc_put_var_double(ncid, lon_varid, help));

  n = 0;
  LOOP_ALL(pert_4mu->ntrack)
    help[n++] = pert_4mu->lat[track][xtrack][ifov];
  NC(nc_put_var_double(ncid, lat_varid, help));

  n = 0;
  LOOP_ALL(pert_4mu->ntrack)
    help[n++] = pert_4mu->dc[track][xtrack][ifov];
  NC(nc_put_var_double(ncid, bt_8mu_varid, help));

  n = 0;
  LOOP_ALL(pert_4mu->ntrack)
    help[n++] = pert_15mu_high->dc[track][xtrack][ifov];
  NC(nc_put_var_double(ncid, bt_10mu_varid, help));

  n = 0;
  LOOP_ALL(pert_4mu->ntrack)
    help[n++] = pert_4mu->bt[track][xtrack][ifov];
  NC(nc_put_var_double(ncid, bt_4mu_varid, help));

  n = 0;
  LOOP_ALL(pert_4mu->ntrack)
    help[n++] = pert_4mu->pt[track][xtrack][ifov];
  NC(nc_put_var_double(ncid, bt_4mu_pt_varid, help));

  n = 0;
  LOOP_ALL(pert_4mu->ntrack)
    help[n++] = pert_4mu->var[track][xtrack][ifov];
  NC(nc_put_var_double(ncid, bt_4mu_var_varid, help));

  n = 0;
  LOOP_ALL(pert_4mu->ntrack)
    help[n++] = pert_15mu_low->bt[track][xtrack][ifov];
  NC(nc_put_var_double(ncid, bt_15mu_low_varid, help));

  n = 0;
  LOOP_ALL(pert_4mu->ntrack)
    help[n++] = pert_15mu_low->pt[track][xtrack][ifov];
  NC(nc_put_var_double(ncid, bt_15mu_low_pt_varid, help));

  n = 0;
  LOOP_ALL(pert_4mu->ntrack)
    help[n++] = pert_15mu_low->var[track][xtrack][ifov];
  NC(nc_put_var_double(ncid, bt_15mu_low_var_varid, help));

  n = 0;
  LOOP_ALL(pert_4mu->ntrack)
    help[n++] = pert_15mu_high->bt[track][xtrack][ifov];
  NC(nc_put_var_double(ncid, bt_15mu_high_varid, help));

  n = 0;
  LOOP_ALL(pert_4mu->ntrack)
    help[n++] = pert_15mu_high->pt[track][xtrack][ifov];
  NC(nc_put_var_double(ncid, bt_15mu_high_pt_varid, help));

  n = 0;
  LOOP_ALL(pert_4mu->ntrack)
    help[n++] = pert_15mu_high->var[track][xtrack][ifov];
  NC(nc_put_var_double(ncid, bt_15mu_high_var_varid, help));

  /* Close file... */
  NC(nc_close(ncid));

  /* Free... */
  free(pert_4mu);
  free(pert_15mu_low);
  free(pert_15mu_high);

  return EXIT_SUCCESS;
}
