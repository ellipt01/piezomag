#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "piezomag.h"
#include "private.h"

#define DUMMY 0

/* allocate fault_params structure */
fault_params *
fault_params_alloc (void)
{
	return (fault_params *) malloc (sizeof (fault_params));
}

/* allocate magnetic_params structure */
magnetic_params *
magnetic_params_alloc (void)
{
	return (magnetic_params *) malloc (sizeof (magnetic_params));
}

/* keywords for parameters */
const int	n_key = 19;
const char *key[] = {
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
};

/* descriptions of keyword */
const char *key_string[] = {
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
};

/* store parameters to global variables and structures */
static bool
set_params (double *items, fault_params *fault, magnetic_params *mag)
{
	int 	i;
	bool	status = true;

	if (fault == NULL || mag == NULL) return false;

	for (i = 0; i < n_key; i++) {
		double val = items[i];
		switch (i) {
			case 0:
				fault->lambda = val;
				break;

			case 1:
				fault->mu = val;
				break;

			case 2:
				mag->beta = val;
				break;

			case 3:
				mag->exf_inc = val;
				break;

			case 4:
				mag->exf_dec = val;
				break;

			case 5:
				mag->mgz_int = val;
				break;

			case 6:
				mag->mgz_inc = val;
				break;

			case 7:
				mag->mgz_dec = val;
				break;

			case 8:
				mag->dcurier = val;
				break;

			case 9:
				fault->u1 = val;
				break;

			case 10:
				fault->u2 = val;
				break;

			case 11:
				fault->u3 = val;
				break;

			case 12:
				fault->fstrike = val;
				break;

			case 13:
				fault->fdip = val;
				break;

			case 14:
				fault->flength1 = val;
				break;

			case 15:
				fault->flength2 = val;
				break;

			case 16:
				fault->fwidth1 = val;
				break;

			case 17:
				fault->fwidth2 = val;
				break;

			case 18:
				fault->fdepth = val;
				break;

			default:
				break;
		}
	}
	return status;
}

/* calculate constants and store them to global variables or members of structure */
static bool
set_constants (fault_params *fault, magnetic_params *mag)
{
	double	jx, jy, jz;

	sd = sin (deg2rad_ (fault->fdip));
	cd = cos (deg2rad_ (fault->fdip));
	td = tan (deg2rad_ (fault->fdip));
	secd = 1.0 / cd;
	sd2 = pow (sd, 2.0);
	cd2 = pow (cd, 2.0);
	s2d = sin (deg2rad_ (2.0 * fault->fdip));
	c2d = cos (deg2rad_ (2.0 * fault->fdip));

	d[0] = DUMMY;	// not referred
	d[1] = fault->fdepth;	// source depth
	d[2] = fault->fdepth - 2.0 * mag->dcurier;	// depth of mirror image
	d[3] = fault->fdepth + 2.0 * mag->dcurier;	// depth of sub-mirror image

	fault->alpha = (fault->lambda + fault->mu) / (fault->lambda + 2.0 * fault->mu);

	alpha0 = 4.0 * fault->alpha - 1.0;
	alpha1 = 3.0 * fault->alpha / alpha0;

	alpha2 = 6.0 * pow (fault->alpha, 2.) / alpha0;
	alpha3 = 2.0 * fault->alpha * (1.0 - fault->alpha) / alpha0;
	alpha4 = fault->alpha * (2.0 * fault->alpha + 1.0) / alpha0;
	alpha5 = fault->alpha * (2.0 * fault->alpha - 5.0) / alpha0;
	alpha6 = 3.0 * fault->alpha * (1.0 - 2.0 * fault->alpha) / alpha0;

	jx = mag->mgz_int * cos (deg2rad_ (mag->mgz_inc)) * cos (deg2rad_ (mag->mgz_dec));
	jy = - mag->mgz_int * cos (deg2rad_ (mag->mgz_inc)) * sin (deg2rad_ (mag->mgz_dec));
	jz = mag->mgz_int * sin (deg2rad_ (mag->mgz_inc));
	rotate (fault->fstrike, &jx, &jy);

	// seismomagnetic moment
	// c0 : intensity
	{
		double	_c1 = fault->mu * (3.0 * fault->lambda + 2.0 * fault->mu);
		double	_c2 = fault->lambda + fault->mu;
		mag->c0 = 0.25 * mag->beta * _c1 / _c2;
	}
	// cx, cy, cz : x(N+S-), y(E+W-) and z(Down+Up-) components
	mag->cx = mag->c0 * jx;
	mag->cy = mag->c0 * jy;
	mag->cz = mag->c0 * jz;

	return true;
}

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
		fprintf (stderr, "ERROR: fread_params: structure is empty.\n");
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

	for (i = 0; i < n_key; i++) {
		if (!is_set_item[i]) {
			fprintf (stderr, "ERROR: fread_params: following parameter is not specified.\n");
			fprintf (stderr, "       %s : %s\n", key[i], key_string[i]);
			return false;
		}
	}
	if (!set_params (items, fault, mag)) return false;
	if (!set_constants (fault, mag)) return false;

	// todo: In here, check whether all parameters are valid

	// check fault dip
	if (fault->fdip < 0. || 90. < fault->fdip) {
		fprintf (stderr, "ERROR: fread_params: fdip must be in [0, 90] (deg.).\n");
		return false;
	}
	fault_is_vertical = (fabs (fault->fdip - 90.) < 1.e-3) ? true : false;

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
				val = fault->lambda;
				break;

			case 1:
				val = fault->mu;
				break;

			case 2:
				val = mag->beta;
				break;

			case 3:
				val = mag->exf_inc;
				break;

			case 4:
				val = mag->exf_dec;
				break;

			case 5:
				val = mag->mgz_int;
				break;

			case 6:
				val = mag->mgz_inc;
				break;

			case 7:
				val = mag->mgz_dec;
				break;

			case 8:
				val = mag->dcurier;
				break;

			case 9:
				val = fault->u1;
				break;

			case 10:
				val = fault->u2;
				break;

			case 11:
				val = fault->u3;
				break;

			case 12:
				val = fault->fstrike;
				break;

			case 13:
				val = fault->fdip;
				break;

			case 14:
				val = fault->flength1;
				break;

			case 15:
				val = fault->flength2;
				break;

			case 16:
				val = fault->fwidth1;
				break;

			case 17:
				val = fault->fwidth2;
				break;

			case 18:
				val = fault->fdepth;
				break;

			default:
				break;
		}
		fprintf (stream, "%s\t : %s = %f\n", key_string[i], key[i], val);
	}
	return;
}

