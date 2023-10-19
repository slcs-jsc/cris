#include "libcris.h"

int main(
  int argc,
  char *argv[]) {

  static pert_t *pert;

  char set[LEN], pertname[LEN];

  double t230 = 230.0, tbg, nesr, nedt = 0, sza2 = 0;

  int asc, itrack, ixtrack, ifov, orb = 0;

  FILE *out;

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <pert.nc> <map.tab>");

  /* Get control parameters... */
  scan_ctl(argc, argv, "PERTNAME", -1, "4mu", pertname);
  scan_ctl(argc, argv, "SET", -1, "full", set);
  int orbit = (int) scan_ctl(argc, argv, "ORBIT", -1, "-999", NULL);
  int track0 = (int) scan_ctl(argc, argv, "TRACK0", -1, "0", NULL);
  int track1 = (int) scan_ctl(argc, argv, "TRACK1", -1, "100000", NULL);
  int xtrack0 = (int) scan_ctl(argc, argv, "XTRACK0", -1, "0", NULL);
  int xtrack1 = (int) scan_ctl(argc, argv, "XTRACK1", -1, "29", NULL);
  int ifov0 = (int) scan_ctl(argc, argv, "IFOV0", -1, "0", NULL);
  int ifov1 = (int) scan_ctl(argc, argv, "IFOV1", -1, "8", NULL);
  double orblat = scan_ctl(argc, argv, "ORBLAT", -1, "0", NULL);
  double t0 = scan_ctl(argc, argv, "T0", -1, "-1e100", NULL);
  double t1 = scan_ctl(argc, argv, "T1", -1, "1e100", NULL);
  double sza0 = scan_ctl(argc, argv, "SZA0", -1, "-1e100", NULL);
  double sza1 = scan_ctl(argc, argv, "SZA1", -1, "1e100", NULL);
  double dt230 = scan_ctl(argc, argv, "DT230", -1, "-999", NULL);
  double nu = scan_ctl(argc, argv, "NU", -1, "-999", NULL);
  int dc = (int) scan_ctl(argc, argv, "DC", -1, "0", NULL);

  /* Allocate... */
  ALLOC(pert, pert_t, 1);

  /* Read perturbation data... */
  read_pert(argv[2], pertname, dc, pert);

  /* Check ranges... */
  track0 = GSL_MIN(GSL_MAX(track0, 0), pert->ntrack - 1);
  track1 = GSL_MIN(GSL_MAX(track1, 0), pert->ntrack - 1);
  xtrack0 = GSL_MIN(GSL_MAX(xtrack0, 0), pert->nxtrack - 1);
  xtrack1 = GSL_MIN(GSL_MAX(xtrack1, 0), pert->nxtrack - 1);
  ifov0 = GSL_MIN(GSL_MAX(ifov0, 0), pert->nfov - 1);
  ifov1 = GSL_MIN(GSL_MAX(ifov1, 0), pert->nfov - 1);

  /* Create output file... */
  LOG(1, "Write perturbation data: %s", argv[3]);
  if (!(out = fopen(argv[3], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = time (seconds since 01-JAN-2000, 00:00 UTC)\n"
	  "# $2 = along-track index\n"
	  "# $3 = longitude [deg]\n"
	  "# $4 = latitude [deg]\n"
	  "# $5 = cloud channel brightness temperature [K]\n"
	  "# $6 = %s brightness temperature [K]\n"
	  "# $7 = %s brightness temperature perturbation [K]\n"
	  "# $8 = %s brightness temperature variance [K^2]\n",
	  pertname, pertname, pertname);

  /* Write data... */
  for (itrack = track0; itrack <= track1; itrack++) {

    /* Count orbits... */
    if (itrack > 0)
      if (pert->lat[itrack - 1][pert->nxtrack / 2][0] <= orblat
	  && pert->lat[itrack][pert->nxtrack / 2][0] >= orblat)
	orb++;

    /* Write output... */
    fprintf(out, "\n");

    /* Check for data gaps... */
    if (itrack > 0 && pert->time[itrack][pert->nxtrack / 2][0]
	- pert->time[itrack - 1][pert->nxtrack / 2][0] >= 10)
      fprintf(out, "\n");

    /* Loop over scans and field of views... */
    for (ixtrack = xtrack0; ixtrack <= xtrack1; ixtrack++)
      for (ifov = ifov0; ifov <= ifov1; ifov++) {

	/* Check data... */
	if (pert->lon[itrack][ixtrack][ifov] < -180
	    || pert->lon[itrack][ixtrack][ifov] > 180
	    || pert->lat[itrack][ixtrack][ifov] < -90
	    || pert->lat[itrack][ixtrack][ifov] > 90)
	  continue;

	/* Get ascending/descending flag... */
	asc = (pert->lat[itrack > 0 ? itrack : itrack + 1][pert->nxtrack / 2]
	       [ifov]
	       > pert->lat[itrack >
			   0 ? itrack - 1 : itrack][pert->nxtrack / 2][ifov]);

	/* Calculate solar zenith angle... */
	if (sza0 >= -1e10 && sza0 <= 1e10 && sza1 >= -1e10 && sza1 <= 1e10)
	  sza2 = sza(pert->time[itrack][ixtrack][ifov],
		     pert->lon[itrack][ixtrack][ifov],
		     pert->lat[itrack][ixtrack][ifov]);

	/* Estimate noise... */
	if (dt230 > 0 && nu > 0) {
	  nesr = planck(t230 + dt230, nu) - planck(t230, nu);
	  tbg =
	    pert->bt[itrack][ixtrack][ifov] - pert->pt[itrack][ixtrack][ifov];
	  nedt = brightness(planck(tbg, nu) + nesr, nu) - tbg;
	}

	/* Write data... */
	if (orbit < 0 || orb == orbit)
	  if (set[0] == 'f' || (set[0] == 'a' && asc)
	      || (set[0] == 'd' && !asc))
	    if (pert->time[itrack][ixtrack][ifov] >= t0
		&& pert->time[itrack][ixtrack][ifov] <= t1
		&& sza2 >= sza0 && sza2 <= sza1)
	      fprintf(out, "%.2f %d %g %g %g %g %g %g\n",
		      pert->time[itrack][ixtrack][ifov], itrack,
		      pert->lon[itrack][ixtrack][ifov],
		      pert->lat[itrack][ixtrack][ifov],
		      pert->dc[itrack][ixtrack][ifov],
		      pert->bt[itrack][ixtrack][ifov],
		      pert->pt[itrack][ixtrack][ifov],
		      pert->var[itrack][ixtrack][ifov] - gsl_pow_2(nedt));
      }
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(pert);

  return EXIT_SUCCESS;
}
