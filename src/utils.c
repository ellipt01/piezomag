#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "piezomag.h"
#include "private.h"

#define DUMMY 0

/* total force
 * double	hx, hy, hz: x(E+W-), y(N+S-) and z(Up-Down+) components
 * double	exf_inc, exf_dec:	inclination and declination of external field
 * double	exf_dec:	declination of external field */
double
total_force (double hx, double hy, double hz, double exf_inc, double exf_dec)
{
	double	f = hx * cos (deg2rad (exf_inc)) * sin (deg2rad (exf_dec))
		- hy * cos (deg2rad (exf_inc)) * cos (deg2rad (exf_dec))
		+ hz * sin (deg2rad (exf_inc));
	return f;
}

/* coordinate rotation */
void
rotate (double theta, double *x, double *y)
{
	double theta_rad = deg2rad (theta);
	double x1 = (*x) * cos (theta_rad) - (*y) * sin (theta_rad);
	double y1 = (*y) * cos (theta_rad) + (*x) * sin (theta_rad);
	*x = x1;
	*y = y1;
	return;
}

/* allocate structure */
fault_params *
fault_params_alloc (void)
{
	return (fault_params *) malloc (sizeof (fault_params));
}

/* allocate structure */
magnetic_params *
magnetic_params_alloc (void)
{
	return (magnetic_params *) malloc (sizeof (magnetic_params));
}

/* keywords for parameters */
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

/* descriptions of keyword */
char *key_string[] = {
	"z coordinates of observation point",
	"lame's constants (lambda)",
	"rigidity (mu)",
	"stress sensitivity",
	"inclination of external geomagnetic field",
	"declination of external geomagnetic field",
	"intensity of initial crustal magnetization",
	"inclination of initial crustal magnetization",
	"declination of initial crustal magnetization",
	"depth of Curier point isotherm",
	"dislocation (strike slip)",
	"dislocation (dip slip)",
	"dislocation (tensile opening)",
	"strike angle of fault",
	"dip angle of fault plane",
	"fault length1",
	"fault length2",
	"fault width1",
	"fault width2",
	"fault depth",
	"component of output result"
};

/* store parameters to global variables and structures */
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

/* calculate constants and store them to global variables or members of structure */
static bool
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

	d[0] = DUMMY;	// not referred
	d[1] = fault->fdepth - z_obs;	// source depth
	d[2] = fault->fdepth - 2.0 * mag->dcurier + z_obs;	// depth of mirror image
	d[3] = fault->fdepth + 2.0 * mag->dcurier - z_obs;	// depth of sub-mirror image

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
	rotate (fault->fstrike, &jx, &jy);

	mag->c0 = 0.25 * mag->beta * fault->mu * (3.0 * fault->lambda + 2.0 * fault->mu) / (fault->lambda + fault->mu);
	mag->cx = mag->c0 * jx;
	mag->cy = mag->c0 * jy;
	mag->cz = mag->c0 * jz;

	return true;
}

/* write parameters to stream */
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
			if (output_comp == X_COMP) fprintf (stream, " = X_COMP\n");
			else if (output_comp == Y_COMP) fprintf (stream, " = Y_COMP\n");
			else if (output_comp == Z_COMP) fprintf (stream, " = Z_COMP\n");
			else if (output_comp == TOTAL_FORCE) fprintf (stream, " = TOTAL_FORCE\n");
		} else fprintf (stream, " = %f\n", val);
	}
	return;
}

#define SPEC_COMP 20

/*c***************************************************
 * read parameters from file
 * and store them in beforehand allocated structures
 ** INPUT **
 * FILE *fp: input file descriptor
 ** OUTPUT **
 * fault_params *fault:	fault parameters
 * magnetic_params *mag:	magnetic parameters
 *c***************************************************/
bool
fread_params (FILE *fp, fault_params *fault, magnetic_params *mag)
{
	int		i;
	bool	is_set_item[n_key];
	double items[n_key];
	char	 buf[BUFSIZ];

	if (fault == NULL || mag == NULL) {
		fprintf (stderr, "ERROR: fread_params: structure is empty.");
		return false;
	}

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

	// todo: In here, check whether all parameters are valid

	// z_obs must be < 0
	if (items[0] >= 0) {
		fprintf (stderr, "ERROR: fread_params: z_obs must be < 0.\n");
		fprintf (stderr, "observation point must be located outside the medium.\n");
		return false;
	}
	// output_comp must be X_COMP(0), Y_COMP(1), Z_COMP(2) or TOTAL_FORCE(3)
	if ((int) items[SPEC_COMP] < 0 || (int) items[SPEC_COMP] >= 4) {
		fprintf (stderr, "ERROR: fread_params: output_comp is invalid.\n");
		return false;
	}
	for (i = 0; i < n_key; i++) {
		if (!is_set_item[i]) {
			fprintf (stderr, "ERROR: fread_params: following parameter is not specified.\n");
			fprintf (stderr, "       %s : %s\n", key[i], key_string[i]);
			return false;
		}
	}
	if (!set_params (items, fault, mag)) return false;
	if (!set_constants (fault, mag)) return false;

	return true;
}

/*c***************************************************
 * calculate some arithmetic constants
 * and store them in global variables
 ** INPUT **
 * sign: + or -
 * double xi, eta and qq: coordinates on fault plane
 *c***************************************************/
void
set_geometry_variables (double sign, double xi, double et, double qq)
{
	double	r3, r5;
	double	rx2, re2, rc2;
	double	rx3, re3, rc3;
	double	r3x2, r3e2, r3c2;
	double	r5x3, r5e3, r5c3;

	r2 = pow (xi, 2.0) + pow (et, 2.0) + pow (qq, 2.0);
	r = sqrt (r2);
	r3 = pow (r, 3.0);
	r5 = pow (r, 5.0);

	yy = et * cd + qq * sd;
	cc = sign * (qq * cd - et * sd);

	rx = r + xi;
	rx2 = pow (rx, 2.0);	// (r + xi)^2
	rx3 = pow (rx, 3.0);	// (r + xi)^2
	r3x2 = r3 * rx2;		// r^3 * (r + xi)^2
	r5x3 = r5 * rx3;		// r^5 * (r + xi)^2

	re = r + et;
	re2 = pow (re, 2.0);	// (r + et)^2
	re3 = pow (re, 3.0);	// (r + et)^2
	r3e2 = r3 * re2;		// r^3 * (r + et)^2
	r5e3 = r5 * re3;		// r^5 * (r + et)^2

	rc = r + cc;
	rc2 = pow (rc, 2.0); // (r + cc)^2
	rc3 = pow (rc, 3.0); // (r + cc)^2
	r3c2 = r3 * rc2;		// r^3 * (r + cc)^2
	r5c3 = r5 * rc3;		// r^5 * (r + cc)^2

	ir  = 0.0;
	ir3 = 0.0;
	ir5 = 0.0;
	if (!singular_R[0]) {
		if (fabs (r) > DBL_EPSILON)  ir  = 1.0 / r;
		if (fabs (r3) > DBL_EPSILON) ir3 = 1.0 / r3;
		if (fabs (r5) > DBL_EPSILON) ir5 = 1.0 / r5;
	}

	irx   = 0.0;
	irx2  = 0.0;
	irx3  = 0.0;
	ir3x2 = 0.0;
	ir5x3 = 0.0;
	if (fabs (rx) > DBL_EPSILON) {
		irx   = 1.0 / rx;
		irx2  = 1.0 / rx2;
		irx3  = 1.0 / rx3;
		ir3x2 = 1.0 / r3x2;
		ir5x3 = 1.0 / r5x3;
	}

	ire   = 0.0;
	ire2  = 0.0;
	ire3  = 0.0;
	ir3e2 = 0.0;
	ir5e3 = 0.0;
	if (!singular_RE[0]) {
		if (fabs (re) > DBL_EPSILON)   ire   = 1.0 / re;
		if (fabs (re2) > DBL_EPSILON)  ire2  = 1.0 / re2;
		if (fabs (r3e2) > DBL_EPSILON) ir3e2 = 1.0 / r3e2;
		if (fabs (re3) > DBL_EPSILON)  ire3  = 1.0 / re3;
		if (fabs (r5e3) > DBL_EPSILON) ir5e3 = 1.0 / r5e3;
	}

	irc   = 0.0;
	irc2  = 0.0;
	irc3  = 0.0;
	ir3c2 = 0.0;
	ir5c3 = 0.0;
	if (fabs (rc) > DBL_EPSILON) {
		irc   = 1.0 / rc;
		irc2  = 1.0 / rc2;
		irc3  = 1.0 / rc3;
		ir3c2 = 1.0 / r3c2;
		ir5c3 = 1.0 / r5c3;
	}

	return;
}

/* treatment of singular points */
static void
clear_singular_flag (int i)
{
	singular_R[i] = false;
	singular_RE[i] = false;
	return;
}

void
clear_all_singular_flags (void)
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
	return (flag[1] || flag[2] || flag[3]);
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
