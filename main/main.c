/*
 * piezomag.c
 *
 *  Created on: 2014/08/14
 *      Author: utsugi
 *
 * Description:
 * reads fault and magnetic parameters from file and calculate seismomagnetic field
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

extern char	*optarg;

/*** range of calculation. By default, x=[-10:0.1:10], y=[-10:0.1:10](km) ***/
double		xwest = -10.;
double		xeast = 10.;
double		dx = 0.1;

double		ysouth = -10.;
double		ynorth = 10.;
double		dy = 0.1;

void
usage (char *toolname)
{
	char	*p = strrchr (toolname, '/');
	if (p) p++;
	else p = toolname;
	fprintf (stderr, "USAGE	 : %s -f <parameter file name> -r <x0/x1/y0/y1> -i <dx/dy>\n", p);
	fprintf (stderr, "optional: -v (verbos mode)\n");
	fprintf (stderr, "-r and -i specify the range of calculation.\n");
	fprintf (stderr, "The seismomagnetic field is calculated on the grid [x0:dx:x1][y0:dy:y1].\n");
	return;
}

/*** set parameters from command-line options ***/
bool
initialize (int argc, char **argv, fault_params **fault, magnetic_params **mag)
{
	char	in_fn[80] = "\0";
	char	c;
	FILE	*fp;

	fault_params		*_fault = (fault_params *) malloc (sizeof (fault_params));
	magnetic_params	*_mag = (magnetic_params *) malloc (sizeof (magnetic_params));

	if (argc <= 1) {
		usage (argv[0]);
		return false;
	}

	verbos = false;
	while ((c = getopt (argc, argv, "f:r:i:v")) != -1) {

		switch (c) {
			case 'f':
			case 'F':
				strcpy (in_fn, optarg);
				break;

			case 'r':
			case 'R':
				sscanf (optarg, "%lf/%lf/%lf/%lf", &xwest, &xeast, &ysouth, &ynorth);
				break;

			case 'i':
			case 'I':
				sscanf (optarg, "%lf/%lf", &dx, &dy);
				break;

			case 'v':
			case 'V':
				verbos = true;
				break;

			default:
				break;
		}

	}

	if (xwest >= xeast || ysouth >= ynorth) {
		fprintf (stderr, "ERROR: -r : specified range is invalid.\n");
		return false;
	}

	if (fabs (dx) < DBL_EPSILON || fabs (dy) < DBL_EPSILON) {
		fprintf (stderr, "ERROR: -r : specified range is invalid.\n");
		return false;
	}


	if ((fp = fopen (in_fn, "r")) == NULL) {
		fprintf (stderr, "ERROR: cannot open input file %s\nprogram aborted.\n", in_fn);
		return false;
	}
	if (!fread_params (fp, _fault, _mag)) return false;
	fclose (fp);

	if (verbos) fwrite_params (stderr, _fault, _mag);
	if (!set_constants (_fault, _mag)) return false;

	if (fault) *fault = _fault;
	if (mag) *mag = _mag;

	return true;
}

int
main (int argc, char **argv)
{
	fault_params		*fault;
	magnetic_params	*mag;

	if (!initialize (argc, argv, &fault, &mag)) exit (1);

  fprintf_piezomagnetic_effect (stdout, output_comp, fault, mag, xwest, xeast, dx, ysouth, ynorth, dy, z_obs);

  free (fault);
  free (mag);
  return EXIT_SUCCESS;
}
