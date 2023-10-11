#include "libcris.h"

int main(
  int argc,
  char *argv[]) {

  static cris_l1_t l1;

  FILE *out;

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <map.tab> <l1b_file1> {<l1b_file2> ...]");

  /* Get control parameters... */
  int apo = (int) scan_ctl(argc, argv, "APO", -1, "0", NULL);
  double nu = scan_ctl(argc, argv, "NU", -1, "1231.2500", NULL);

  /* Create file... */
  LOG(1, "Write map: %s", argv[2]);
  if (!(out = fopen(argv[2], "w")))
    ERRMSG("Cannot create file!");

  /* Write header... */
  fprintf(out,
	  "# $1 = time (seconds since 01-JAN-2000, 00:00 UTC)\n"
	  "# $2 = satellite longitude [deg]\n"
	  "# $3 = satellite latitude [deg]\n"
	  "# $4 = footprint longitude [deg]\n"
	  "# $5 = footprint latitude [deg]\n"
	  "# $6 = wavenumber [cm^-1]\n"
	  "# $7 = brightness temperature [K]\n"
	  "# $8 = radiance [W/(m^2 sr cm^-1)]\n"
	  "# $9 = noise [W/(m^2 sr cm^-1)]\n" "# $10 = channel index\n\n");

  /* Loop over files... */
  for (int iarg = 3; iarg < argc; iarg++) {

    /* Read CrIS data... */
    if (!read_cris_l1(argv[iarg], &l1, apo))
      continue;

    /* Write data... */
    for (int ichan = 0; ichan < L1_NCHAN_LW; ichan++)
      if (fabs(l1.nu_lw[ichan] - nu) < 0.1)
	for (int track = 0; track < L1_NTRACK; track++)
	  for (int xtrack = 0; xtrack < L1_NXTRACK; xtrack++)
	    for (int ifov = 0; ifov < L1_NFOV; ifov++)
	      fprintf(out, "%.2f %g %g %g %g %.4f %g %g %g LW_%03d\n",
		      l1.time[track][xtrack] - 220838400,
		      l1.sat_lon[track], l1.sat_lat[track],
		      l1.lon[track][xtrack][ifov],
		      l1.lat[track][xtrack][ifov], l1.nu_lw[ichan],
		      brightness(l1.rad_lw[track][xtrack][ifov][ichan] * 1e-3,
				 l1.nu_lw[ichan]),
		      l1.rad_lw[track][xtrack][ifov][ichan] * 1e-3,
		      l1.nedn_lw[ifov][ichan] * 1e-3, ichan);

    for (int ichan = 0; ichan < L1_NCHAN_MW; ichan++)
      if (fabs(l1.nu_mw[ichan] - nu) < 0.1)
	for (int track = 0; track < L1_NTRACK; track++)
	  for (int xtrack = 0; xtrack < L1_NXTRACK; xtrack++)
	    for (int ifov = 0; ifov < L1_NFOV; ifov++)
	      fprintf(out, "%.2f %g %g %g %g %.4f %g %g %g MW_%03d\n",
		      l1.time[track][xtrack] - 220838400,
		      l1.sat_lon[track], l1.sat_lat[track],
		      l1.lon[track][xtrack][ifov],
		      l1.lat[track][xtrack][ifov], l1.nu_mw[ichan],
		      brightness(l1.rad_mw[track][xtrack][ifov][ichan] * 1e-3,
				 l1.nu_mw[ichan]),
		      l1.rad_mw[track][xtrack][ifov][ichan] * 1e-3,
		      l1.nedn_mw[ifov][ichan] * 1e-3, ichan);

    for (int ichan = 0; ichan < L1_NCHAN_SW; ichan++)
      if (fabs(l1.nu_sw[ichan] - nu) < 0.1)
	for (int track = 0; track < L1_NTRACK; track++)
	  for (int xtrack = 0; xtrack < L1_NXTRACK; xtrack++)
	    for (int ifov = 0; ifov < L1_NFOV; ifov++)
	      fprintf(out, "%.2f %g %g %g %g %.4f %g %g %g SW_%03d\n",
		      l1.time[track][xtrack] - 220838400,
		      l1.sat_lon[track], l1.sat_lat[track],
		      l1.lon[track][xtrack][ifov],
		      l1.lat[track][xtrack][ifov], l1.nu_sw[ichan],
		      brightness(l1.rad_sw[track][xtrack][ifov][ichan] * 1e-3,
				 l1.nu_sw[ichan]),
		      l1.rad_sw[track][xtrack][ifov][ichan] * 1e-3,
		      l1.nedn_sw[ifov][ichan] * 1e-3, ichan);
  }

  /* Close file... */
  fclose(out);

  return EXIT_SUCCESS;
}
