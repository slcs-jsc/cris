#include "libcris.h"

int main(
  int argc,
  char *argv[]) {

  static pert_t *pert;

  char set[LEN], pertname[LEN];

  double orblat, nu, t230 = 230.0, dt230, tbg, nesr, nedt = 0,
    t0, t1, sza0, sza1, sza2 = 0;

  int asc, itrack, ixtrack, ifov, orb = 0, orbit;

  FILE *out;

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <pert.nc> <map.tab>");

  /* Get control parameters... */
  scan_ctl(argc, argv, "PERTNAME", -1, "4mu", pertname);
  scan_ctl(argc, argv, "SET", -1, "full", set);
  orbit = (int) scan_ctl(argc, argv, "ORBIT", -1, "-999", NULL);
  orblat = scan_ctl(argc, argv, "ORBLAT", -1, "0", NULL);
  t0 = scan_ctl(argc, argv, "T0", -1, "-1e100", NULL);
  t1 = scan_ctl(argc, argv, "T1", -1, "1e100", NULL);
  sza0 = scan_ctl(argc, argv, "SZA0", -1, "-1e100", NULL);
  sza1 = scan_ctl(argc, argv, "SZA1", -1, "1e100", NULL);
  dt230 = scan_ctl(argc, argv, "DT230", -1, "0.16", NULL);
  nu = scan_ctl(argc, argv, "NU", -1, "2345.0", NULL);

  /* Allocate... */
  ALLOC(pert, pert_t, 1);

  /* Read perturbation data... */
  read_pert(argv[2], pertname, pert);

  /* Create output file... */
  printf("Write perturbation data: %s\n", argv[3]);
  if (!(out = fopen(argv[3], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = time (seconds since 01-JAN-2000, 00:00 UTC)\n"
	  "# $2 = along-track index\n"
	  "# $3 = longitude [deg]\n"
	  "# $4 = latitude [deg]\n"
	  "# $5 = 8mu brightness temperature [K]\n"
	  "# $6 = %s brightness temperature [K]\n"
	  "# $7 = %s brightness temperature perturbation [K]\n"
	  "# $8 = %s brightness temperature variance [K^2]\n",
	  pertname, pertname, pertname);

  /* Write data... */
  for (itrack = 0; itrack < pert->ntrack; itrack++) {

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
    for (ixtrack = 0; ixtrack < pert->nxtrack; ixtrack++)
      for (ifov = 0; ifov < pert->nfov; ifov++) {

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
	if (dt230 > 0) {
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
