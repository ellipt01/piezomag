#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "piezomag.h"
#include "private.h"

typedef struct {
	char	*keyword;
	char	*description;
} param_keyword;

/* keyword and description of parameters */
#define n_keys 19
const param_keyword	keys[n_keys] = {
		{"lambda",   "lame's constants (lambda)"},
		{"mu",       "rigidity (mu)"},
		{"u1",       "dislocation (strike-slip)"},
		{"u2",       "dislocation (dip-slip)"},
		{"u3",       "dislocation (tensile-opening)"},
		{"fstrike",  "strike angle of fault"},
		{"fdip",     "dip angle of fault"},
		{"flength1", "fault length1"},
		{"flength2", "fault length2"},
		{"fwidth1",  "fault width1"},
		{"fwidth2",  "fault width2"},
		{"fdepth",   "fault depth"},
		{"beta",     "stress sensitivity"},
		{"exf_inc",  "inclination of external geomagnetic field"},
		{"exf_dec",  "declination of external geomagnetic field"},
		{"mgz_int",  "intensity of initial crustal magnetization"},
		{"mgz_inc",  "inclination of initial crustal magnetization"},
		{"mgz_dec",  "declination of initial crustal magnetization"},
		{"dcurier",  "depth of Curier point isotherm"}
};

/* store parameters to structures */
static bool
set_params (double *items, fault_params *fault, magnetic_params *mag)
{
	int 	i;
	bool	status = true;

	if (fault == NULL || mag == NULL) return false;

	for (i = 0; i < n_keys; i++) {
		double val = items[i];
		switch (i) {
			// fault parameters
			case 0:
				fault->lambda = val;
				break;

			case 1:
				fault->mu = val;
				break;

			case 2:
				fault->u1 = val;
				break;

			case 3:
				fault->u2 = val;
				break;

			case 4:
				fault->u3 = val;
				break;

			case 5:
				fault->fstrike = val;
				break;

			case 6:
				fault->fdip = val;
				break;

			case 7:
				fault->flength1 = val;
				break;

			case 8:
				fault->flength2 = val;
				break;

			case 9:
				fault->fwidth1 = val;
				break;

			case 10:
				fault->fwidth2 = val;
				break;

			case 11:
				fault->fdepth = val;
				break;

			case 12:
				mag->beta = val;
				break;

			// magnetic properties
			case 13:
				mag->exf_inc = val;
				break;

			case 14:
				mag->exf_dec = val;
				break;

			case 15:
				mag->mgz_int = val;
				break;

			case 16:
				mag->mgz_inc = val;
				break;

			case 17:
				mag->mgz_dec = val;
				break;

			case 18:
				mag->dcurier = val;
				break;

			default:
				break;
		}
	}
	return status;
}

/* calculate constants and store them to global variables and structures */
static bool
set_constants (fault_params *fault, magnetic_params *mag)
{
	double	jx, jy, jz;

	/* trigonometric functions */
	sd = sin (deg2rad (fault->fdip));
	cd = cos (deg2rad (fault->fdip));
	td = tan (deg2rad (fault->fdip));
	secd = 1.0 / cd;
	sd2 = pow (sd, 2.0);
	cd2 = pow (cd, 2.0);
	s2d = sin (deg2rad (2.0 * fault->fdip));
	c2d = cos (deg2rad (2.0 * fault->fdip));

	/* depth of source and mirror images */
	d[0] = _PIEZOMAG_DUMMY_;	// not referred
	d[1] = fault->fdepth;	// d1 : source depth
	d[2] = fault->fdepth - 2.0 * mag->dcurier;	// d2: depth of mirror image
	d[3] = fault->fdepth + 2.0 * mag->dcurier;	// d3: depth of sub-mirror image

	/* alpha */
	fault->alpha = (fault->lambda + fault->mu) / (fault->lambda + 2.0 * fault->mu);

	alpha0 = 4.0 * fault->alpha - 1.0;
	alpha1 = 3.0 * fault->alpha / alpha0;

	alpha2 = 6.0 * pow (fault->alpha, 2.) / alpha0;
	alpha3 = 2.0 * fault->alpha * (1.0 - fault->alpha) / alpha0;
	alpha4 = fault->alpha * (2.0 * fault->alpha + 1.0) / alpha0;
	alpha5 = fault->alpha * (2.0 * fault->alpha - 5.0) / alpha0;
	alpha6 = 3.0 * fault->alpha * (1.0 - 2.0 * fault->alpha) / alpha0;

	/* seismomagnetic moment on fault coordinate system */
	// c0 : intensity
	{
		double	_c1 = fault->mu * (3.0 * fault->lambda + 2.0 * fault->mu);
		double	_c2 = fault->lambda + fault->mu;
		mag->c0 = 0.25 * mag->beta * _c1 / _c2;
	}
	// initial crustal magnetization on fault coordinate system
	jx = mag->mgz_int * cos (deg2rad (mag->mgz_inc)) * cos (deg2rad (mag->mgz_dec));
	jy = - mag->mgz_int * cos (deg2rad (mag->mgz_inc)) * sin (deg2rad (mag->mgz_dec));
	jz = mag->mgz_int * sin (deg2rad (mag->mgz_inc));
	rotate (fault->fstrike, &jx, &jy);
	// (cx, cy, cz) : seismomagnetic moment vector on fault coordinate system
	mag->cx = mag->c0 * jx;
	mag->cy = mag->c0 * jy;
	mag->cz = mag->c0 * jz;

	return true;
}

/*c***************************************************
 * read parameters from file
 * and store them in previously allocated structures
 ** INPUT **
 * FILE *fp: input file descriptor
 * fault_params *fault:	fault parameters
 * magnetic_params *mag:	magnetic parameters
 *c***************************************************/
bool
fread_params (FILE *fp, fault_params *fault, magnetic_params *mag)
{
	int		i;
	bool	is_set_item[n_keys];
	double	items[n_keys];
	char	buf[BUFSIZ];

	if (fault == NULL || mag == NULL) {
		fprintf (stderr, "ERROR: fread_params: structure is empty.\n");
		return false;
	}

	for (i = 0; i < n_keys; i++) is_set_item[i] = false;

	while (fgets (buf, BUFSIZ, fp) != NULL) {
		char	*ptr_buf = buf;

		// skip comment and blank lines
		if (ptr_buf[0] == '#' || ptr_buf[0] == '\n') continue;
		// skip blanks of head of line
		while (ptr_buf[0] == ' ' || ptr_buf[0] == '\t') ptr_buf++;

		// read buffer
		for (i = 0; i < n_keys; i++) {
			// if first word match the keyword, read values after '='
			if (strncmp (ptr_buf, keys[i].keyword, (size_t) strlen (keys[i].keyword)) == 0) {
				char	*ptr_val;
				if ((ptr_val = strrchr (buf, '=')) == NULL) continue;
				while (ptr_val[0] == '=' || ptr_val[0] == ' ' || ptr_val[0] == '\t') ptr_val++;

				items[i] = (double) atof (ptr_val);
				is_set_item[i] = true;

				break;
			}
		}
	}
	// check whether all parameters are specified
	{
		bool	status = true;
		for (i = 0; i < n_keys; i++) {
			if (!is_set_item[i]) {
				if (status) fprintf (stderr, "ERROR: fread_params: following parameter(s) is not specified.\n");
				fprintf (stderr, "       %s : %s\n", keys[i].keyword, keys[i].description);
				if (status) status = false;
			}
		}
		if (!status) return false;
	}
	// set parameters, calculate constants
	if (!set_params (items, fault, mag)) return false;
	if (!set_constants (fault, mag)) return false;

	// todo: In here, check whether all parameters are valid

	// check elastic properties of the crust
	if (fault->lambda < 0.) {
		fprintf (stderr, "ERROR: Lame's constant lambda must be >= 0.\n");
		return false;
	}
	if (fault->mu < 0.) {
		fprintf (stderr, "ERROR: rigidity mu must be >= 0.\n");
		return false;
	}

	// fault must be buried in the ground
	if (fault->fdepth + fault->fwidth1 * sd < 0.0) {
		fprintf (stderr, "ERROR: fault must be buried in the ground\n");
		fprintf (stderr, "       i.e. fdepth + fwidth1 * sin (delta) must be >= 0.\n");
		return false;
	}

	// check fault dip (0 <= fdip <= 90 deg.)
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
	for (i = 0; i < n_keys; i++) {
		double val = 0.;
		switch (i) {
			case 0:
				val = fault->lambda;
				break;

			case 1:
				val = fault->mu;
				break;

			case 2:
				val = fault->u1;
				break;

			case 3:
				val = fault->u2;
				break;

			case 4:
				val = fault->u3;
				break;

			case 5:
				val = fault->fstrike;
				break;

			case 6:
				val = fault->fdip;
				break;

			case 7:
				val = fault->flength1;
				break;

			case 8:
				val = fault->flength2;
				break;

			case 9:
				val = fault->fwidth1;
				break;

			case 10:
				val = fault->fwidth2;
				break;

			case 11:
				val = fault->fdepth;
				break;

			case 12:
				val = mag->beta;
				break;

			case 13:
				val = mag->exf_inc;
				break;

			case 14:
				val = mag->exf_dec;
				break;

			case 15:
				val = mag->mgz_int;
				break;

			case 16:
				val = mag->mgz_inc;
				break;

			case 17:
				val = mag->mgz_dec;
				break;

			case 18:
				val = mag->dcurier;
				break;

			default:
				break;
		}
		fprintf (stream, "%s\t : %s = %f\n", keys[i].description, keys[i].keyword, val);
	}
	return;
}

