#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include "piez.h"
#include "constants.h"

int	n_key = 27;
char *key[] = {
	"z_obs",
	"lambda",
	"mu",
	"beta",
	"exf_inc",
	"exf_dec",
	"mgz_int",
	"mgz_inc",
	"mgz_dec",
	"dcurier",
	"u1",
	"u2",
	"u3",
	"fstrike",
	"fdip",
	"flength1",
	"flength2",
	"fwidth1",
	"fwidth2",
	"fdepth",
	"x01",
	"x02",
	"dx",
	"y01",
	"y02",
	"dy",
	"output_comp"
};

char *key_string[] = {
	"z coordinates of observation point",
	"lame's constants (lambda)",
	"lame's constants (mu)",
	"stress sensitivity",
	"inclination of external geomagnetic field",
	"declination of external geomagnetic field",
	"intensity of initial crustal magnetization",
	"inclination of initial crustal magnetization",
	"declination of initial crustal magnetization",
	"depth of curier point isotherm",
	"displacement (strike slip)",
	"displacement (dip slip)",
	"displacement (tensile opening)",
	"strike angle of fault",
	"dip angle of fault plane",
	"fault length1",
	"fault length2",
	"fault width1",
	"fault width2",
	"fault depth",
	"x-range left",
	"x-range right",
	"x-grid spacing",
	"y-range bottom",
	"y-range upper",
	"y-grid spacing",
	"component of output result"
};

void
usage (char *toolname)
{
	char	*p = strrchr (toolname, '/');
	if (p) p++;
	else p = toolname;
	fprintf (stderr, "USAGE	 : %s -f <parameter file name>\n", p);
	fprintf (stderr, "optional: -v (verbos mode)\n");
	return;
}

extern char	*optarg;

bool
initialize (int argc, char **argv)
{
	char	in_fn[80];
	char	c;
	FILE	*fp;

	if (argc <= 1) {
		usage (argv[0]);
		return false;
	}

	verbos = false;
	while ((c = getopt (argc, argv, "f:v")) != -1) {

		switch (c) {
			case 'f':
			case 'F':
				strcpy (in_fn, optarg);
				break;

			case 'v':
			case 'V':
				verbos = true;
				break;

			default:
				break;
		}

	}

	if ((fp = fopen (in_fn, "r")) == NULL) {
		fprintf (stderr, "ERROR: cannot open input file %s\nprogram aborted.\n", in_fn);
		return false;
	}
	if (!fread_params (fp)) return false;
	fclose (fp);

	if (verbos) fwrite_params (stderr);
	if (!set_constants ()) return false;

	return true;
}

void
set_params (double *items)
{
	int i;

	for (i = 0; i < n_key; i++) {
		double val = items[i];
		switch (i) {
		case 0:
			z_obs = val;
			break;

		case 1:
			lambda = val;
			break;

		case 2:
			mu = val;
			break;

		case 3:
			beta = val;
			break;

		case 4:
			exf_inc = val;
			break;

		case 5:
			exf_dec = val;
			break;

		case 6:
			mgz_int = val;
			break;

		case 7:
			mgz_inc = val;
			break;

		case 8:
			mgz_dec = val;
			break;

		case 9:
			dcurier = val;
			break;

		case 10:
			u1 = val;
			break;

		case 11:
			u2 = val;
			break;

		case 12:
			u3 = val;
			break;

		case 13:
			fstrike = val;
			break;

		case 14:
			fdip = val;
			break;

		case 15:
			flength1 = val;
			break;

		case 16:
			flength2 = val;
			break;

		case 17:
			fwidth1 = val;
			break;

		case 18:
			fwidth2 = val;
			break;

		case 19:
			fdepth = val;
			break;

		case 20:
			x01 = val;
			break;

		case 21:
			x02 = val;
			break;

		case 22:
			dx = val;
			break;

		case 23:
			y01 = val;
			break;

		case 24:
			y02 = val;
			break;

		case 25:
			dy = val;
			break;

		case 26:
			output_comp = (int) val;
			break;

		default:
			break;
		}
	}
	return;
}

void
fwrite_params (FILE *stream)
{
	int i;
	for (i = 0; i < n_key; i++) {
		double val = 0.;
		switch (i) {
		case 0:
			val = z_obs;
			break;

		case 1:
			val = lambda;
			break;

		case 2:
			val = mu;
			break;

		case 3:
			val = beta;
			break;

		case 4:
			val = exf_inc;
			break;

		case 5:
			val = exf_dec;
			break;

		case 6:
			val = mgz_int;
			break;

		case 7:
			val = mgz_inc;
			break;

		case 8:
			val = mgz_dec;
			break;

		case 9:
			val = dcurier;
			break;

		case 10:
			val = u1;
			break;

		case 11:
			val = u2;
			break;

		case 12:
			val = u3;
			break;

		case 13:
			val = fstrike;
			break;

		case 14:
			val = fdip;
			break;

		case 15:
			val = flength1;
			break;

		case 16:
			val = flength2;
			break;

		case 17:
			val = fwidth1;
			break;

		case 18:
			val = fwidth2;
			break;

		case 19:
			val = fdepth;
			break;

		case 20:
			val = x01;
			break;

		case 21:
			val = x02;
			break;

		case 22:
			val = dx;
			break;

		case 23:
			val = y01;
			break;

		case 24:
			val = y02;
			break;

		case 25:
			val = dy;
			break;

		case 26:
			val = output_comp;
			break;

		default:
			break;
		}
		fprintf (stream, "%s\t : %s", key_string[i], key[i]);
		if (i == 26) {
			if (output_comp == X_COMP) fprintf (stream, "= X_COMP\n");
			else if (output_comp == Y_COMP) fprintf (stream, "= Y_COMP\n");
			else if (output_comp == Z_COMP) fprintf (stream, "= Z_COMP\n");
			else if (output_comp == TOTAL_FORCE) fprintf (stream, "= TOTAL_FORCE\n");
		} else fprintf (stream, "= %f\n", val);
	}
	return;
}

#define SPEC_COMP 26

bool
fread_params (FILE *fp)
{
	int		i;
	bool	is_set_item[n_key];
	double items[n_key];
	char	 buf[BUFSIZ];

	for (i = 0; i < n_key; i++) is_set_item[i] = false;

	while (fgets (buf, BUFSIZ, fp) != NULL) {

		if (buf[0] == '#' || buf[0] == '\n') continue;

		for (i = 0; i < n_key; i++) {
			if (strncmp (buf, key[i], (size_t) strlen (key[i])) == 0) {
				char	*ptr;
				if ((ptr = strrchr (buf, '=')) == NULL) continue;
				while (ptr[0] == '=' || ptr[0] == ' ') ptr++;

				items[i] = (double) atof (ptr);
				is_set_item[i] = true;

				break;
			}
		}
	}
	if ((int) items[SPEC_COMP] < 0 || (int) items[SPEC_COMP] >= 4) {
		fprintf (stderr, "ERROR: specified component invalid.\n");
		return false;
	}
	for (i = 0; i < n_key; i++) {
		if (!is_set_item[i]) {
			fprintf (stderr, "ERROR: following parameter is not specified.\n");
			fprintf (stderr, "			 %s : %s\n", key[i], key_string[i]);
			return false;
		}
	}
	set_params (items);
	return true;
}

void
coordinates_transform (double theta, double *x, double *y)
{
	double theta_rad = deg2rad (theta);
	double x1 = (*x) * cos (theta_rad) - (*y) * sin (theta_rad);
	double y1 = (*y) * cos (theta_rad) + (*x) * sin (theta_rad);
	*x = x1;
	*y = y1;
	return;
}

void
clear_singular_flag (int i)
{
	singular_R[i] = false;
	singular_RE[i] = false;
	return;
}

void
clear_all_singular_flag (void)
{
	int	i;
	for (i = 0; i < 4; i++) clear_singular_flag (i);
	return;
}

void
set_singular_flag (int i)
{
	if (i <= 0 || i > 3) return;
	singular_R[0] = singular_R[i];
	singular_RE[0] = singular_RE[i];
	return;
}

bool
is_singular_point (bool *flag)
{
	return flag[1] + flag[2] + flag[3];
}

void
check_singular_point (double x, double y, double eps)
{
	int	i;

	if (fabs (x + flength1) < eps || fabs (x - flength2) < eps) {
		for (i = 0; i < 3; i++) {
			double p = y * cd - d[i + 1] * sd;
			double q = y * sd + d[i + 1] * cd;
			if ((fabs (p + fwidth1) < eps || fabs (p - fwidth2) < eps) && fabs (q) < eps) singular_R[i + 1] = true;
			if ((p + fwidth1 < 0.0 || p < fwidth2) && fabs (q) < eps) singular_RE[i + 1] = true;
		}
	}
	return;
}
