#include "libcris.h"

int main(
  int argc,
  char *argv[]) {

  static cris_l1_t l1;

  FILE *out;

  double dmin = 1e100, x0[3], x1[3];

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <l1b_file> <spec.tab>");

  /* Get control parameters... */
  const int apo = (int) scan_ctl(argc, argv, "APO", -1, "0", NULL);
  int track = (int) scan_ctl(argc, argv, "TRACK", -1, "0", NULL);
  int xtrack = (int) scan_ctl(argc, argv, "XTRACK", -1, "0", NULL);
  int ifov = (int) scan_ctl(argc, argv, "IFOV", -1, "0", NULL);
  const double lon = (int) scan_ctl(argc, argv, "LON", -1, "-999", NULL);
  const double lat = (int) scan_ctl(argc, argv, "LAT", -1, "-999", NULL);

  /* Read CrIS data... */
  if (!read_cris_l1(argv[2], &l1, apo))
    ERRMSG("Cannot read CrIS Level-1B file!");

  /* Find nearest footprint... */
  if (lon >= -180 && lon <= 180 && lat >= -90 && lat <= 90) {
    geo2cart(0, lon, lat, x0);
    for (int track2 = 0; track2 < L1_NTRACK; track2++)
      for (int xtrack2 = 0; xtrack2 < L1_NXTRACK; xtrack2++)
	for (int ifov2 = 0; ifov2 < L1_NFOV; ifov2++) {
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
    LOG(1, "nearest footprint: lon= %g, lat= %g, track= %d, xtrack=%d",
	l1.lon[track][xtrack][ifov], l1.lat[track][xtrack][ifov],
	track, xtrack);
  }

  /* Check indices... */
  if (track < 0 || track >= L1_NTRACK)
    ERRMSG("Along-track index out of range!");
  if (xtrack < 0 || xtrack >= L1_NXTRACK)
    ERRMSG("Across-track index out of range!");
  if (ifov < 0 || ifov >= L1_NFOV)
    ERRMSG("Field of view index out of range!");

  /* Create file... */
  LOG(1, "Write spectrum: %s", argv[3]);
  if (!(out = fopen(argv[3], "w")))
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

  /* Write data... */
  for (int ichan = 0; ichan < L1_NCHAN_LW; ichan++)
    fprintf(out, "%.2f %g %g %g %g %.4f %g %g %g LW_%03d\n",
	    l1.time[track][xtrack] - 220838400,
	    l1.sat_lon[track], l1.sat_lat[track],
	    l1.lon[track][xtrack][ifov], l1.lat[track][xtrack][ifov],
	    l1.nu_lw[ichan],
	    BRIGHT(l1.rad_lw[track][xtrack][ifov][ichan] * 1e-3,
		   l1.nu_lw[ichan]),
	    l1.rad_lw[track][xtrack][ifov][ichan] * 1e-3,
	    l1.nedn_lw[ifov][ichan] * 1e-3, ichan);

  fprintf(out, "\n");
  for (int ichan = 0; ichan < L1_NCHAN_MW; ichan++)
    fprintf(out, "%.2f %g %g %g %g %.4f %g %g %g MW_%03d\n",
	    l1.time[track][xtrack] - 220838400,
	    l1.sat_lon[track], l1.sat_lat[track],
	    l1.lon[track][xtrack][ifov], l1.lat[track][xtrack][ifov],
	    l1.nu_mw[ichan],
	    BRIGHT(l1.rad_mw[track][xtrack][ifov][ichan] * 1e-3,
		   l1.nu_mw[ichan]),
	    l1.rad_mw[track][xtrack][ifov][ichan] * 1e-3,
	    l1.nedn_mw[ifov][ichan] * 1e-3, ichan);

  fprintf(out, "\n");
  for (int ichan = 0; ichan < L1_NCHAN_SW; ichan++)
    fprintf(out, "%.2f %g %g %g %g %.4f %g %g %g SW_%03d\n",
	    l1.time[track][xtrack] - 220838400,
	    l1.sat_lon[track], l1.sat_lat[track],
	    l1.lon[track][xtrack][ifov], l1.lat[track][xtrack][ifov],
	    l1.nu_sw[ichan],
	    BRIGHT(l1.rad_sw[track][xtrack][ifov][ichan] * 1e-3,
		   l1.nu_sw[ichan]),
	    l1.rad_sw[track][xtrack][ifov][ichan] * 1e-3,
	    l1.nedn_sw[ifov][ichan] * 1e-3, ichan);

  /* Close file... */
  fclose(out);

  return EXIT_SUCCESS;
}
