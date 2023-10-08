#include "libcris.h"

int main(
  int argc,
  char *argv[]) {

  static pert_t *pert;

  FILE *out;

  char pertname[LEN];

  double mu, nedt = -1e99, nedt_old;

  int itrack;

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <pert.nc> <noise.tab>");

  /* Read control parameters... */
  scan_ctl(argc, argv, "PERTNAME", -1, "4mu", pertname);
  int bsize = (int) scan_ctl(argc, argv, "BSIZE", -1, "-999", NULL);
  int maxvar = (int) scan_ctl(argc, argv, "MAXVAR", -1, "-999", NULL);

  /* Allocate... */
  ALLOC(pert, pert_t, 1);

  /* Read perturbation data... */
  read_pert(argv[2], pertname, pert);

  /* Set block size... */
  if (bsize < 0)
    bsize = pert->nxtrack;

  /* Create file... */
  LOG(1, "Write noise data: %s", argv[3]);
  if (!(out = fopen(argv[3], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = longitude [deg]\n"
	  "# $2 = latitude [deg]\n"
	  "# $3 = mean brightness temperature [K]\n"
	  "# $4 = noise estimate [K]\n\n");

  /* Loop over granules... */
  for (itrack = 0; itrack < pert->ntrack - bsize; itrack += bsize) {

    /* Estimate noise... */
    nedt_old = nedt;
    noise_pert(pert, itrack, itrack + bsize, &mu, &nedt);

    /* Write output... */
    if (maxvar <= 0
	|| fabs(200 * (nedt - nedt_old) / (nedt + nedt_old)) < maxvar)
      fprintf(out, "%g %g %g %g\n",
	      pert->lon[itrack + bsize / 2][pert->nxtrack / 2][4],
	      pert->lat[itrack + bsize / 2][pert->nxtrack / 2][4], mu, nedt);
  }

  /* Close file... */
  fclose(out);

  /* Free... */
  free(pert);

  return EXIT_SUCCESS;
}
