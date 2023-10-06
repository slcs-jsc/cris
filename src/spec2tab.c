#include "libcris.h"

int main(
  int argc,
  char *argv[]) {

  static cris_l1_t l1;

  FILE *out;

  double dmin = 1e100, x0[3], x1[3];

  int ichan, track = -1, track2, xtrack = -1, xtrack2, ifov = -1, ifov2;

  /* Check arguments... */
  if (argc != 7)
    ERRMSG("Give parameters: <cris_l1b_file> "
	   "[index <track> <xtrack> <ifov> | geo <lon> <lat> <dummy>] <spec.tab>");

  /* Read CrIS data... */
  printf("Read CrIS Level-1B data file: %s\n", argv[1]);
  read_cris_l1(argv[1], &l1);

  /* Get indices... */
  if (argv[2][0] == 'i') {
    track = atoi(argv[3]);
    xtrack = atoi(argv[4]);
    ifov = atoi(argv[5]);
  }

  /* Find nearest footprint... */
  else {
    geo2cart(0, atof(argv[3]), atof(argv[4]), x0);
    for (track2 = 0; track2 < L1_NTRACK; track2++)
      for (xtrack2 = 0; xtrack2 < L1_NXTRACK; xtrack2++)
	for (ifov2 = 0; ifov2 < L1_NFOV; ifov2++) {
	  geo2cart(0, l1.lon[track2][xtrack2][ifov2],
		   l1.lat[track2][xtrack2][ifov2], x1);
	  if (DIST2(x0, x1) < dmin) {
	    dmin = DIST2(x0, x1);
	    track = track2;
	    xtrack = xtrack2;
	    ifov = ifov2;
	  }
	}
    if (dmin > 2500)
      ERRMSG("Geolocation not covered by granule!");
    printf("nearest footprint: lon= %g, lat= %g, track= %d, xtrack=%d\n",
	   l1.lon[track][xtrack][ifov],
	   l1.lat[track][xtrack][ifov], track, xtrack);
  }

  /* Check indices... */
  if (track < 0 || track >= L1_NTRACK)
    ERRMSG("Along-track index out of range!");
  if (xtrack < 0 || xtrack >= L1_NXTRACK)
    ERRMSG("Across-track index out of range!");

#if 0
  /* Flag bad observations... */
  for (ichan = 0; ichan < AIRS_RAD_CHANNEL; ichan++)
    if ((airs_rad_gran.state[track][xtrack] != 0)
	|| (airs_rad_gran.ExcludedChans[ichan] > 2)
	|| (airs_rad_gran.CalChanSummary[ichan] & 8)
	|| (airs_rad_gran.CalChanSummary[ichan] & (32 + 64))
	|| (airs_rad_gran.CalFlag[track][ichan] & 16))
      airs_rad_gran.radiances[track][xtrack][ichan]
	= (float) sqrt(-1.0);
#endif

  /* Create file... */
  printf("Write spectrum: %s\n", argv[6]);
  if (!(out = fopen(argv[6], "w")))
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
	  "# $9 = waveband\n" "# $10 = channel number\n\n");

  /* Write data... */
  for (ichan = 0; ichan < L1_NCHAN_LW; ichan++)
    fprintf(out, "%.2f %g %g %g %g %g %g %g LW %d\n",
	    l1.time[track][xtrack] - 220838400,
	    l1.sat_lon[track], l1.sat_lat[track],
	    l1.lon[track][xtrack][ifov], l1.lat[track][xtrack][ifov],
	    l1.nu_lw[ichan],
	    brightness(l1.rad_lw[track][xtrack][ifov][ichan] * 1e-3,
		       l1.nu_lw[ichan]),
	    l1.rad_lw[track][xtrack][ifov][ichan] * 1e-3, ichan);

  fprintf(out, "\n");
  for (ichan = 0; ichan < L1_NCHAN_MW; ichan++)
    fprintf(out, "%.2f %g %g %g %g %g %g %g MW %d\n",
	    l1.time[track][xtrack] - 220838400,
	    l1.sat_lon[track], l1.sat_lat[track],
	    l1.lon[track][xtrack][ifov], l1.lat[track][xtrack][ifov],
	    l1.nu_mw[ichan],
	    brightness(l1.rad_mw[track][xtrack][ifov][ichan] * 1e-3,
		       l1.nu_mw[ichan]),
	    l1.rad_mw[track][xtrack][ifov][ichan] * 1e-3, ichan);

  fprintf(out, "\n");
  for (ichan = 0; ichan < L1_NCHAN_SW; ichan++)
    fprintf(out, "%.2f %g %g %g %g %g %g %g SW %d\n",
	    l1.time[track][xtrack] - 220838400,
	    l1.sat_lon[track], l1.sat_lat[track],
	    l1.lon[track][xtrack][ifov], l1.lat[track][xtrack][ifov],
	    l1.nu_sw[ichan],
	    brightness(l1.rad_sw[track][xtrack][ifov][ichan] * 1e-3,
		       l1.nu_sw[ichan]),
	    l1.rad_sw[track][xtrack][ifov][ichan] * 1e-3, ichan);

  /* Close file... */
  fclose(out);

  return EXIT_SUCCESS;
}
