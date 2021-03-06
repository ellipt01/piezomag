/*
 * main.c
 *
 *  Created on: 2014/08/14
 *      Author: utsugi
 *
 * An example program
 *
 * Description:
 * reads fault and magnetic parameters from file and calculates seismomagnetic field
 * on the grid in the specified range.
 *
 * Usage:
 * see the function usage()
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <float.h>

#include "piezomag.h"

/*** range of interest. By default, x=[-10:0.1:10], y=[-10:0.1:10] (km) ***/
double		xwest = -10.;
double		xeast = 10.;
double		dx = 0.1;

double		ysouth = -10.;
double		ynorth = 10.;
double		dy = 0.1;

/* z-coordinate of observation point. By default, z = -0.001 (km) ***/
double		zobs = -0.001;

/* output magnetic component:
 * MAG_COMP_F(default), MAG_COMP_X, MAG_COMP_Y or MAG_COMP_Z */
MagComp	output_comp = MAG_COMP_F;

/* verbos mide */
bool		verbos = false;

/* convert int -> MagComp */
static void
set_output_comp (int val)
{
	switch (val) {
	case 0:
		output_comp = MAG_COMP_F;
		break;
	case 1:
		output_comp = MAG_COMP_X;
		break;
	case 2:
		output_comp = MAG_COMP_Y;
		break;
	case 3:
		output_comp = MAG_COMP_Z;
		break;
	default:
		output_comp = MAG_COMP_NONE;
		break;
	}
	return;
}

static void
usage (char *toolname)
{
	char	*p = strrchr (toolname, '/');
	if (p) p++;
	else p = toolname;
	fprintf (stderr, "USAGE: %s -f <parameter file name> -r <x0(S)/x1(N)/y0(W)/y1(E)> -i <dx/dy>\n", p);
	fprintf (stderr, "          -z <zobs> -o <output component>\n");
	fprintf (stderr, "          [ -v -h ]\n");
	fprintf (stderr, "-f:  gives input parameter file name.\n");
	fprintf (stderr, "-r:  specifies the min/max coordinates of region of interest.\n");
	fprintf (stderr, "-i:  gives grid increments.\n");
	fprintf (stderr, "-z:  gives z-coordinates of observation point (< 0).\n");
	fprintf (stderr, "-o:  specifies output magnetic component (0=F,1=X,2=Y,3=Z).\n");
	fprintf (stderr, "=== optional ===\n");
	fprintf (stderr, "-v:  verbos mode.\n");
	fprintf (stderr, "-h:  show this message.\n\n");
	fprintf (stderr, "The seismomagnetic field is calculated on the grid x=[x0:dx:x1], y=[y0:dy:y1], z=zobs.\n\n");
	return;
}


/*** set parameters from command-line options ***/
static bool
initialize (int argc, char **argv, fault_params **fault, magnetic_params **mag)
{
	extern char	*optarg;
	char			in_fn[80] = "\0";
	char			c;
	FILE			*fp;

	fault_params		*_fault = fault_params_alloc ();
	magnetic_params	*_mag = magnetic_params_alloc ();

	if (argc <= 1) {
		usage (argv[0]);
		return false;
	}

	verbos = false;
	while ((c = getopt (argc, argv, "f:r:i:z:o:vh")) != -1) {

		switch (c) {
			case 'f':
				strcpy (in_fn, optarg);
				break;

			case 'r':
				sscanf (optarg, "%lf/%lf/%lf/%lf", &xwest, &xeast, &ysouth, &ynorth);
				break;

			case 'i':
				sscanf (optarg, "%lf/%lf", &dx, &dy);
				break;

			case 'z':
				zobs = (double) atof (optarg);
				break;

			case 'o':
				set_output_comp (atoi (optarg));
				break;

			case 'v':
				verbos = true;
				break;

			case 'h':
				usage (argv[0]);
				exit (1);

			default:
				break;
		}

	}

	if (xwest >= xeast || ysouth >= ynorth) {
		fprintf (stderr, "ERROR: -r : specified range is invalid.\n");
		return false;
	}

	if (fabs (dx) < DBL_EPSILON || fabs (dy) < DBL_EPSILON) {
		fprintf (stderr, "ERROR: -i : specified grid interval is invalid.\n");
		return false;
	}

	if ((fp = fopen (in_fn, "r")) == NULL) {
		fprintf (stderr, "ERROR: cannot open input file %s.\n", in_fn);
		return false;
	}
	if (!fread_params (fp, _fault, _mag)) return false;
	fclose (fp);

	if (verbos) fwrite_params (stderr, _fault, _mag);

	if (fault) *fault = _fault;
	if (mag) *mag = _mag;

	return true;
}

int
main (int argc, char **argv)
{
	fault_params		*fault;
	magnetic_params	*mag;

	if (!initialize (argc, argv, &fault, &mag)) {
		fprintf (stderr, "initialization of program failed, aborted.\n");
		exit (1);
	}

	fprintf_seismomagnetic_field (stdout, output_comp, fault, mag, xwest, xeast, dx, ysouth, ynorth, dy, zobs);

	free (fault);
	free (mag);

	return EXIT_SUCCESS;
}
