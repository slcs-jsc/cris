#include "libcris.h"

int main(
  int argc,
  char *argv[]) {

  static pert_t *pert;

  static double d, dmin, dmax, dmu, dx[L1_NXTRACK], x0[3], x1[3], x2[3];

  static int i, itrack, ixtrack, n, nx[L1_NXTRACK];

  /* Check arguments... */
  if (argc < 3)
    ERRMSG("Give parameters: <ctl> <pert.nc>");

  /* Allocate... */
  ALLOC(pert, pert_t, 1);

  /* Read perturbation data... */
  read_pert(argv[2], "4mu", 0, pert);

  /* Init... */
  dmin = 1e100;
  dmax = -1e100;
  dmu = 0;
  n = 0;

  /* Get swath width... */
  for (itrack = 0; itrack < pert->ntrack; itrack++) {
    geo2cart(0, pert->lon[itrack][0][4], pert->lat[itrack][0][4], x0);
    geo2cart(0, pert->lon[itrack][pert->nxtrack - 1][4],
	     pert->lat[itrack][pert->nxtrack - 1][4], x1);
    d = 2. * RE * asin(DIST(x0, x1) / (2. * RE));
    dmin = GSL_MIN(dmin, d);
    dmax = GSL_MAX(dmax, d);
    dmu += d;
    n++;
  }

  /* Write output... */
  printf("\nmean_swath_width=    %.1f km\n", dmu / n);
  printf("minimum_swath_width= %.1f km\n", dmin);
  printf("maximum_swath_width= %.1f km\n", dmax);

  /* Init... */
  dmin = 1e100;
  dmax = -1e100;
  dmu = 0;
  n = 0;

  /* Get across-track sampling distances... */
  for (itrack = 0; itrack < pert->ntrack; itrack++) {
    for (ixtrack = 0; ixtrack < pert->nxtrack - 1; ixtrack++) {
      geo2cart(0, pert->lon[itrack][ixtrack][4],
	       pert->lat[itrack][ixtrack][4], x0);
      geo2cart(0, pert->lon[itrack][ixtrack + 1][4],
	       pert->lat[itrack][ixtrack + 1][4], x1);
      d = 2. * RE * asin(DIST(x0, x1) / (2. * RE));
      dmin = GSL_MIN(dmin, d);
      dmax = GSL_MAX(dmax, d);
      dmu += d;
      n++;
    }
  }

  /* Write output... */
  printf("\nmean_across_track_sampling_distance=    %.1f km | %.1f km\n",
	 dmu / n, dmu / n / 3.);
  printf("minimum_across_track_sampling_distance= %.1f km | %.1f km\n", dmin,
	 dmin / 3.);
  printf("maximum_across_track_sampling_distance= %.1f km | %.1f km\n", dmax,
	 dmax / 3.);

  /* Init... */
  dmin = 1e100;
  dmax = -1e100;
  dmu = 0;
  n = 0;

  /* Get along-track sampling distances... */
  for (itrack = 0; itrack < pert->ntrack - 1; itrack++) {
    for (ixtrack = 0; ixtrack < pert->nxtrack; ixtrack++) {
      geo2cart(0, pert->lon[itrack][ixtrack][4],
	       pert->lat[itrack][ixtrack][4], x0);
      geo2cart(0, pert->lon[itrack + 1][ixtrack][4],
	       pert->lat[itrack + 1][ixtrack][4], x1);
      d = 2. * RE * asin(DIST(x0, x1) / (2. * RE));
      dmin = GSL_MIN(dmin, d);
      dmax = GSL_MAX(dmax, d);
      dmu += d;
      n++;
    }
  }

  /* Write output... */
  printf("\nmean_along_track_sampling_distance=    %.1f km | %.1f km\n",
	 dmu / n, dmu / n / 3.);
  printf("minimum_along_track_sampling_distance= %.1f km | %.1f km\n", dmin,
	 dmin / 3.);
  printf("maximum_along_track_sampling_distance= %.1f km | %.1f km\n", dmax,
	 dmax / 3.);

  /* Init... */
  dmin = 1e100;
  dmax = -1e100;
  dmu = 0;
  n = 0;

  /* Get angle between along-track and across-track direction... */
  for (itrack = 0; itrack < pert->ntrack - 1; itrack++) {
    geo2cart(0, pert->lon[itrack][pert->nxtrack / 2][4],
	     pert->lat[itrack][pert->nxtrack / 2][4], x0);
    geo2cart(0, pert->lon[itrack][pert->nxtrack / 2 + 1][4],
	     pert->lat[itrack][pert->nxtrack / 2 + 1][4], x1);
    geo2cart(0, pert->lon[itrack + 1][pert->nxtrack / 2][4],
	     pert->lat[itrack + 1][pert->nxtrack / 2][4], x2);
    for (i = 0; i < 3; i++) {
      x1[i] -= x0[i];
      x2[i] -= x0[i];
    }
    d = acos(DOTP(x1, x2) / (NORM(x1) * NORM(x2))) * 180. / M_PI;
    dmin = GSL_MIN(dmin, d);
    dmax = GSL_MAX(dmax, d);
    dmu += d;
    n++;
  }

  /* Write output... */
  printf("\nmean_across_track_angle=    %.1f deg\n", dmu / n);
  printf("minimum_across_track_angle= %.1f deg\n", dmin);
  printf("maximum_across_track_angle= %.1f deg\n", dmax);

  /* Get across-track distances... */
  for (itrack = 0; itrack < pert->ntrack; itrack++) {
    for (ixtrack = 0; ixtrack < pert->nxtrack; ixtrack++) {
      geo2cart(0, pert->lon[itrack][0][4], pert->lat[itrack][0][4], x0);
      geo2cart(0, pert->lon[itrack][ixtrack][4],
	       pert->lat[itrack][ixtrack][4], x1);
      d = 2. * RE * asin(DIST(x0, x1) / (2. * RE));
      dx[ixtrack] += d;
      nx[ixtrack]++;
    }
  }

  /* Write output... */
  printf("\n");
  for (ixtrack = 0; ixtrack < pert->nxtrack; ixtrack++)
    printf("ixtrack= %2d | dx= %.1f km | cx= %.1f km\n", ixtrack,
	   dx[ixtrack] / nx[ixtrack],
	   0.5 * dx[pert->nxtrack - 1] / nx[pert->nxtrack - 1] -
	   dx[ixtrack] / nx[ixtrack]);

  /* Free... */
  free(pert);

  return EXIT_SUCCESS;
}
