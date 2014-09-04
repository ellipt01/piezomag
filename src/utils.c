/*
 * utils.c
 *
 *  Created on: 2014/09/04
 *      Author: utsugi
 */

#include <stdio.h>
#include <math.h>

#include "piezomag.h"

/* total force
 * double	hx, hy, hz: x(N+S-), y(E+W-) and z(Down+Up-) components
 * double	exf_inc, exf_dec:	inclination and declination of external field
 * double	exf_dec:	declination of external field (clockwise +). */
double
total_force (double hx, double hy, double hz, double exf_inc, double exf_dec)
{
	double	f = hx * cos (deg2rad_ (exf_inc)) * cos (deg2rad_ (exf_dec))
		- hy * cos (deg2rad_ (exf_inc)) * sin (deg2rad_ (exf_dec))
		+ hz * sin (deg2rad_ (exf_inc));
	return f;
}

/* coordinate rotation */
void
rotate (double theta, double *x, double *y)
{
	double theta_rad = deg2rad_ (theta);
	double x1 = (*x) * cos (theta_rad) + (*y) * sin (theta_rad);
	double y1 = (*y) * cos (theta_rad) - (*x) * sin (theta_rad);
	*x = x1;
	*y = y1;
	return;
}

/* check specified magnetic component is valid */
bool
check_mag_component (MagComp component)
{
	if (MAG_COMP_F <= component && component <= MAG_COMP_Z) return true;
	return false;
}

/* check specified seismomagnetic term is valid */
bool
check_seismo_mag_term (SeismoMagTerm term)
{
	if (term & SEISMO_MAG_MAIN) return true;
	if (term & SEISMO_MAG_MIRROR) return true;
	if (term & SEISMO_MAG_SUBMIRROR) return true;
	return false;
}
