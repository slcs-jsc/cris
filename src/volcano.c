#include "libcris.h"

/* ------------------------------------------------------------
   Functions...
   ------------------------------------------------------------ */

/* Estimate noise. */
double get_noise(
  double rad,
  double nesr,
  double nu);

/* ------------------------------------------------------------
   Main...
   ------------------------------------------------------------ */

int main(
  int argc,
  char *argv[]) {

  FILE *out;

  static cris_l1_t l1;

  static double ci, ci_err,
    ai_low, ai_low_err, ai_low_bt1, ai_low_bt2, ai_high, ai_high_err,
    ai_high_bt1, ai_high_bt2, ai_old, ai_old_err, ai_old_bt1, ai_old_bt2,
    si_high, si_high_err, si_high_bt1, si_high_bt2, si_low, si_low_err,
    si_low_bt1, si_low_bt2, si_old, si_old_err, si_old_bt1, si_old_bt2,
    si_oper, si_oper_err, si_oper_bt1, si_oper_bt2;

  static int iarg, ai_low_nu1 = 341, ai_low_nu2 = 499, ai_high_nu1 =
    40, ai_high_nu2 = 698, ai_old_nu1 = 295, ai_old_nu2 = 499, ci_nu =
    36, si_low_nu1 = 326, si_low_nu2 = 260, si_high_nu1 = 328, si_high_nu2 =
    282, si_old_nu1 = 318, si_old_nu2 = 260, si_oper_nu1 = 359, si_oper_nu2 =
    244;

  /* Check arguments... */
  if (argc < 4)
    ERRMSG("Give parameters: <ctl> <out.tab> <l1b_file1> [<l1b_file2> ...]");

  /* Get control parameters... */
  int apo = (int) scan_ctl(argc, argv, "APO", -1, "0", NULL);

  /* Create file... */
  printf("Write volcanic emission data: %s\n", argv[2]);
  if (!(out = fopen(argv[2], "w")))
    ERRMSG("Cannot create file!");

  /* Loop over CrIS files... */
  for (iarg = 3; iarg < argc; iarg++) {

    /* Read CrIS data... */
    if (!read_cris_l1(argv[iarg], &l1, apo))
      continue;

    /* Write header... */
    if (iarg == 3) {
      fprintf(out,
	      "# $1  = time [s]\n"
	      "# $2  = footprint longitude [deg]\n"
	      "# $3  = footprint latitude [deg]\n"
	      "# $4  = satellite altitude [km]\n"
	      "# $5  = satellite longitude [deg]\n"
	      "# $6  = satellite latitude [deg]\n");
      fprintf(out,
	      "# $7  = cloud index, BT(%.4f/cm) [K]\n"
	      "# $8  = cloud index error [K]\n"
	      "# $9  = ash index (low wavenumbers),"
	      " BT(%.4f/cm) - BT(%.4f/cm) [K]\n"
	      "# $10 = ash index (low wavenumbers) error [K]\n"
	      "# $11 = ash index (high wavenumbers),"
	      " BT(%.4f/cm) - BT(%.4f/cm) [K]\n"
	      "# $12 = ash index (high wavenumbers) error [K]\n"
	      "# $13 = ash index (Hoffmann et al., 2014),"
	      " BT(%.4f/cm) - BT(%.4f/cm) [K]\n"
	      "# $14 = ash index (Hoffmann et al., 2014) error [K]\n",
	      l1.nu_mw[ci_nu],
	      l1.nu_lw[ai_low_nu1], l1.nu_lw[ai_low_nu2],
	      l1.nu_mw[ai_high_nu1], l1.nu_lw[ai_high_nu2],
	      l1.nu_lw[ai_old_nu1], l1.nu_lw[ai_old_nu2]);
      fprintf(out,
	      "# $15 = SO2 index (low concentrations),"
	      " BT(%.4f/cm) - BT(%.4f/cm) [K]\n"
	      "# $16 = SO2 index (low concentrations) error [K]\n"
	      "# $17 = SO2 index (high concentrations),"
	      " BT(%.4f/cm) - BT(%.4f/cm) [K]\n"
	      "# $18 = SO2 index (high concentrations) error [K]\n"
	      "# $19 = SO2 index (operational),"
	      " BT(%.4f/cm) - BT(%.4f/cm) [K]\n"
	      "# $20 = SO2 index (operational) error [K]\n"
	      "# $21 = SO2 index (Hoffmann et al., 2014),"
	      " BT(%.4f/cm) - BT(%.4f/cm) [K]\n"
	      "# $22 = SO2 index (Hoffmann et al., 2014) error [K]\n\n",
	      l1.nu_mw[si_low_nu1], l1.nu_mw[si_low_nu2],
	      l1.nu_mw[si_high_nu1], l1.nu_mw[si_high_nu2],
	      l1.nu_mw[si_oper_nu1], l1.nu_mw[si_oper_nu2],
	      l1.nu_mw[si_old_nu1], l1.nu_mw[si_old_nu2]);
    }

    /* Loop over footprints... */
    for (int track = 0; track < L1_NTRACK; track++)
      for (int xtrack = 0; xtrack < L1_NXTRACK; xtrack++)
	for (int ifov = 0; ifov < L1_NFOV; ifov++) {

	  /* cloud index... */
	  ci = BRIGHT(l1.rad_mw[track][xtrack][ifov][ci_nu] * 0.001,
		      l1.nu_mw[ci_nu]);
	  ci_err = get_noise(l1.rad_mw[track][xtrack][ifov][ci_nu] * 0.001,
			     l1.nedn_mw[ifov][ci_nu] * 0.001,
			     l1.nu_mw[ci_nu]);

	  /* ash index (low wavenumbers)... */
	  ai_low_bt1 =
	    BRIGHT(l1.rad_lw[track][xtrack][ifov][ai_low_nu1] *
		   0.001, l1.nu_lw[ai_low_nu1]);
	  ai_low_bt2 =
	    BRIGHT(l1.rad_lw[track][xtrack][ifov][ai_low_nu2] *
		   0.001, l1.nu_lw[ai_low_nu2]);
	  ai_low = ai_low_bt1 - ai_low_bt2;
	  ai_low_err =
	    sqrt(gsl_pow_2
		 (get_noise
		  (l1.rad_lw[track][xtrack][ifov][ai_low_nu1] * 0.001,
		   l1.nedn_lw[ifov][ai_low_nu1] * 0.001,
		   l1.nu_lw[ai_low_nu1]))
		 +
		 gsl_pow_2(get_noise
			   (l1.rad_lw[track][xtrack][ifov][ai_low_nu2] *
			    0.001, l1.nedn_lw[ifov][ai_low_nu2] * 0.001,
			    l1.nu_lw[ai_low_nu2])));

	  /* ash index (high wavenumbers)... */
	  ai_high_bt1 =
	    BRIGHT(l1.rad_mw[track][xtrack][ifov][ai_high_nu1] *
		   0.001, l1.nu_mw[ai_high_nu1]);
	  ai_high_bt2 =
	    BRIGHT(l1.rad_lw[track][xtrack][ifov][ai_high_nu2] *
		   0.001, l1.nu_lw[ai_high_nu2]);
	  ai_high = ai_high_bt1 - ai_high_bt2;
	  ai_high_err =
	    sqrt(gsl_pow_2
		 (get_noise
		  (l1.rad_mw[track][xtrack][ifov][ai_high_nu1] * 0.001,
		   l1.nedn_mw[ifov][ai_high_nu1] * 0.001,
		   l1.nu_mw[ai_high_nu1]))
		 +
		 gsl_pow_2(get_noise
			   (l1.rad_lw[track][xtrack][ifov][ai_high_nu2] *
			    0.001, l1.nedn_lw[ifov][ai_high_nu2] * 0.001,
			    l1.nu_lw[ai_high_nu1])));

	  /* ash index (old)... */
	  ai_old_bt1 =
	    BRIGHT(l1.rad_lw[track][xtrack][ifov][ai_old_nu1] *
		   0.001, l1.nu_lw[ai_old_nu1]);
	  ai_old_bt2 =
	    BRIGHT(l1.rad_lw[track][xtrack][ifov][ai_old_nu2] *
		   0.001, l1.nu_lw[ai_old_nu2]);
	  ai_old = ai_old_bt1 - ai_old_bt2;
	  ai_old_err =
	    sqrt(gsl_pow_2
		 (get_noise
		  (l1.rad_lw[track][xtrack][ifov][ai_old_nu1] * 0.001,
		   l1.nedn_lw[ifov][ai_old_nu1] * 0.001,
		   l1.nu_lw[ai_old_nu1]))
		 +
		 gsl_pow_2(get_noise
			   (l1.rad_lw[track][xtrack][ifov][ai_old_nu2] *
			    0.001, l1.nedn_lw[ifov][ai_old_nu2] * 0.001,
			    l1.nu_lw[ai_old_nu2])));

	  /* SO2 index (low concentrations)... */
	  si_low_bt1 =
	    BRIGHT(l1.rad_mw[track][xtrack][ifov][si_low_nu1] *
		   0.001, l1.nu_mw[si_low_nu1]);
	  si_low_bt2 =
	    BRIGHT(l1.rad_mw[track][xtrack][ifov][si_low_nu2] *
		   0.001, l1.nu_mw[si_low_nu2]);
	  si_low = si_low_bt1 - si_low_bt2;
	  si_low_err =
	    sqrt(gsl_pow_2
		 (get_noise
		  (l1.rad_mw[track][xtrack][ifov][si_low_nu1] * 0.001,
		   l1.nedn_mw[ifov][si_low_nu1] * 0.001,
		   l1.nu_mw[si_low_nu1]))
		 +
		 gsl_pow_2(get_noise
			   (l1.rad_mw[track][xtrack][ifov][si_low_nu2] *
			    0.001, l1.nedn_mw[ifov][si_low_nu2] * 0.001,
			    l1.nu_mw[si_low_nu2])));

	  /* SO2 index (high concentrations)... */
	  si_high_bt1 =
	    BRIGHT(l1.rad_mw[track][xtrack][ifov][si_high_nu1] *
		   0.001, l1.nu_mw[si_high_nu1]);
	  si_high_bt2 =
	    BRIGHT(l1.rad_mw[track][xtrack][ifov][si_high_nu2] *
		   0.001, l1.nu_mw[si_high_nu2]);
	  si_high = si_high_bt1 - si_high_bt2;
	  si_high_err =
	    sqrt(gsl_pow_2
		 (get_noise
		  (l1.rad_mw[track][xtrack][ifov][si_high_nu1] * 0.001,
		   l1.nedn_mw[ifov][si_high_nu1] * 0.001,
		   l1.nu_mw[si_high_nu1]))
		 +
		 gsl_pow_2(get_noise
			   (l1.rad_mw[track][xtrack][ifov][si_high_nu2] *
			    0.001, l1.nedn_mw[ifov][si_high_nu2] * 0.001,
			    l1.nu_mw[si_high_nu2])));

	  /* SO2 index (operational)... */
	  si_oper_bt1 =
	    BRIGHT(l1.rad_mw[track][xtrack][ifov][si_oper_nu1] *
		   0.001, l1.nu_mw[si_oper_nu1]);
	  si_oper_bt2 =
	    BRIGHT(l1.rad_mw[track][xtrack][ifov][si_oper_nu2] *
		   0.001, l1.nu_mw[si_oper_nu2]);
	  si_oper = si_oper_bt1 - si_oper_bt2;
	  si_oper_err =
	    sqrt(gsl_pow_2
		 (get_noise
		  (l1.rad_mw[track][xtrack][ifov][si_oper_nu1] * 0.001,
		   l1.nedn_mw[ifov][si_oper_nu1] * 0.001,
		   l1.nu_mw[si_oper_nu1]))
		 +
		 gsl_pow_2(get_noise
			   (l1.rad_mw[track][xtrack][ifov][si_oper_nu2] *
			    0.001, l1.nedn_mw[ifov][si_oper_nu2] * 0.001,
			    l1.nu_mw[si_oper_nu2])));

	  /* SO2 index (old)... */
	  si_old_bt1 =
	    BRIGHT(l1.rad_mw[track][xtrack][ifov][si_old_nu1] *
		   0.001, l1.nu_mw[si_old_nu1]);
	  si_old_bt2 =
	    BRIGHT(l1.rad_mw[track][xtrack][ifov][si_old_nu2] *
		   0.001, l1.nu_mw[si_old_nu2]);
	  si_old = si_old_bt1 - si_old_bt2;
	  si_old_err =
	    sqrt(gsl_pow_2
		 (get_noise
		  (l1.rad_mw[track][xtrack][ifov][si_old_nu1] * 0.001,
		   l1.nedn_mw[ifov][si_old_nu1] * 0.001,
		   l1.nu_mw[si_old_nu1]))
		 +
		 gsl_pow_2(get_noise
			   (l1.rad_mw[track][xtrack][ifov][si_old_nu2] *
			    0.001, l1.nedn_mw[ifov][si_old_nu2] * 0.001,
			    l1.nu_mw[si_old_nu2])));

	  /* Write output... */
	  fprintf(out,
		  "%.2f %.4f %.4f %.3f %.4f %.4f %.2f %.2f %.2f %.2f %.2f %.2f "
		  "%.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n",
		  l1.time[track][xtrack] - 220838400,
		  l1.lon[track][xtrack][ifov],
		  l1.lat[track][xtrack][ifov],
		  l1.sat_z[track],
		  l1.sat_lon[track],
		  l1.sat_lat[track],
		  ci, ci_err, ai_low, ai_low_err, ai_high, ai_high_err,
		  ai_old, ai_old_err, si_low, si_low_err, si_high,
		  si_high_err, si_oper, si_oper_err, si_old, si_old_err);
	}
  }

  /* Close file... */
  fclose(out);

  return EXIT_SUCCESS;
}

/************************************************************************/

double get_noise(
  double rad,
  double nesr,
  double nu) {

  return BRIGHT(rad + nesr, nu) - BRIGHT(rad, nu);
}
