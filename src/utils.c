#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "piezomag.h"
#include "private.h"

#define DUMMY 0

const int	n_key = 21;
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
	"component of output result"
};

static bool
set_params (double *items, fault_params *fault, magnetic_params *mag)
{
	int i;

	if (fault == NULL || mag == NULL) return false;

	for (i = 0; i < n_key; i++) {
		double val = items[i];
		switch (i) {
			case 0:
				z_obs = val;
				break;

			case 1:
				fault->lambda = val;
				break;

			case 2:
				fault->mu = val;
				break;

			case 3:
				mag->beta = val;
				break;

			case 4:
				mag->exf_inc = val;
				break;

			case 5:
				mag->exf_dec = val;
				break;

			case 6:
				mag->mgz_int = val;
				break;

			case 7:
				mag->mgz_inc = val;
				break;

			case 8:
				mag->mgz_dec = val;
				break;

			case 9:
				mag->dcurier = val;
				break;

			case 10:
				fault->u1 = val;
				break;

			case 11:
				fault->u2 = val;
				break;

			case 12:
				fault->u3 = val;
				break;

			case 13:
				fault->fstrike = val;
				break;

			case 14:
				fault->fdip = val;
				break;

			case 15:
				fault->flength1 = val;
				break;

			case 16:
				fault->flength2 = val;
				break;

			case 17:
				fault->fwidth1 = val;
				break;

			case 18:
				fault->fwidth2 = val;
				break;

			case 19:
				fault->fdepth = val;
				break;

			case 20:
				output_comp = (int) val;
				break;

			default:
				break;
		}
	}
	return true;
}

void
fwrite_params (FILE *stream, const fault_params *fault, const magnetic_params *mag)
{
	int i;
	for (i = 0; i < n_key; i++) {
		double val = 0.;
		switch (i) {
			case 0:
				val = z_obs;
				break;

			case 1:
				val = fault->lambda;
				break;

			case 2:
				val = fault->mu;
				break;

			case 3:
				val = mag->beta;
				break;

			case 4:
				val = mag->exf_inc;
				break;

			case 5:
				val = mag->exf_dec;
				break;

			case 6:
				val = mag->mgz_int;
				break;

			case 7:
				val = mag->mgz_inc;
				break;

			case 8:
				val = mag->mgz_dec;
				break;

			case 9:
				val = mag->dcurier;
				break;

			case 10:
				val = fault->u1;
				break;

			case 11:
				val = fault->u2;
				break;

			case 12:
				val = fault->u3;
				break;

			case 13:
				val = fault->fstrike;
				break;

			case 14:
				val = fault->fdip;
				break;

			case 15:
				val = fault->flength1;
				break;

			case 16:
				val = fault->flength2;
				break;

			case 17:
				val = fault->fwidth1;
				break;

			case 18:
				val = fault->fwidth2;
				break;

			case 19:
				val = fault->fdepth;
				break;

			case 20:
				val = output_comp;
				break;

			default:
				break;
		}
		fprintf (stream, "%s\t : %s", key_string[i], key[i]);
		if (i == 20) {
			if (output_comp == X_COMP) fprintf (stream, "= X_COMP\n");
			else if (output_comp == Y_COMP) fprintf (stream, "= Y_COMP\n");
			else if (output_comp == Z_COMP) fprintf (stream, "= Z_COMP\n");
			else if (output_comp == TOTAL_FORCE) fprintf (stream, "= TOTAL_FORCE\n");
		} else fprintf (stream, "= %f\n", val);
	}
	return;
}

#define SPEC_COMP 20

bool
fread_params (FILE *fp, fault_params *fault, magnetic_params *mag)
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
	set_params (items, fault, mag);
	return true;
}

bool
set_constants (fault_params *fault, magnetic_params *mag)
{
	double	jx, jy, jz;

	fault->fstrike = 90.0 - fault->fstrike;

	sd = sin (deg2rad (fault->fdip));
	cd = cos (deg2rad (fault->fdip));
	td = tan (deg2rad (fault->fdip));
	secd = 1.0 / cd;
	sd2 = pow (sd, 2.0);
	cd2 = pow (cd, 2.0);
	s2d = sin (deg2rad (2.0 * fault->fdip));
	c2d = cos (deg2rad (2.0 * fault->fdip));

	d[0] = DUMMY;
	d[1] = fault->fdepth - z_obs;
	d[2] = fault->fdepth - 2.0 * mag->dcurier + z_obs;
	d[3] = fault->fdepth + 2.0 * mag->dcurier - z_obs;

	fault->alpha = (fault->lambda + fault->mu) / (fault->lambda + 2.0 * fault->mu);

	alpha0 = 4.0 * fault->alpha - 1.0;
	alpha1 = 3.0 * fault->alpha / alpha0;

	alpha2 = 6.0 * pow (fault->alpha, 2.) / alpha0;
	alpha3 = 2.0 * fault->alpha * (1.0 - fault->alpha) / alpha0;
	alpha4 = fault->alpha * (2.0 * fault->alpha + 1.0) / alpha0;
	alpha5 = fault->alpha * (2.0 * fault->alpha - 5.0) / alpha0;
	alpha6 = 3.0 * fault->alpha * (1.0 - 2.0 * fault->alpha) / alpha0;

	jx = mag->mgz_int * cos (deg2rad (mag->mgz_inc)) * sin (deg2rad (mag->mgz_dec));
	jy = - mag->mgz_int * cos (deg2rad (mag->mgz_inc)) * cos (deg2rad (mag->mgz_dec));
	jz = mag->mgz_int * sin (deg2rad (mag->mgz_inc));
	coordinates_transform (fault->fstrike, &jx, &jy);

	mag->c0 = 0.25 * mag->beta * fault->mu * (3.0 * fault->lambda + 2.0 * fault->mu) / (fault->lambda + fault->mu);
	mag->cx = mag->c0 * jx;
	mag->cy = mag->c0 * jy;
	mag->cz = mag->c0 * jz;

	return true;
}

/***************************************
 * coordinate transform
 ***************************************/
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

/***************************************
 * treatment of singular point
 ***************************************/
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
check_singular_point (const fault_params *fault, double x, double y, double eps)
{
	int	i;

	if (fabs (x + fault->flength1) < eps || fabs (x - fault->flength2) < eps) {
		for (i = 0; i < 3; i++) {
			double p = y * cd - d[i + 1] * sd;
			double q = y * sd + d[i + 1] * cd;
			if ((fabs (p + fault->fwidth1) < eps || fabs (p - fault->fwidth2) < eps) && fabs (q) < eps) singular_R[i + 1] = true;
			if ((p + fault->fwidth1 < 0.0 || p < fault->fwidth2) && fabs (q) < eps) singular_RE[i + 1] = true;
		}
	}
	return;
}
