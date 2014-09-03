/*
 * private.c
 *
 *  Created on: 2014/09/02
 *      Author: utsugi
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "piezomag.h"
#include "private.h"

/* total force
 * double	hx, hy, hz: x(N+S-), y(E+W-) and z(Down+Up-) components
 * double	exf_inc, exf_dec:	inclination and declination of external field
 * double	exf_dec:	declination of external field (clockwise +). */
double
total_force (double hx, double hy, double hz, double exf_inc, double exf_dec)
{
	double	f = hx * cos (_deg2rad_ (exf_inc)) * cos (_deg2rad_ (exf_dec))
		- hy * cos (_deg2rad_ (exf_inc)) * sin (_deg2rad_ (exf_dec))
		+ hz * sin (_deg2rad_ (exf_inc));
	return f;
}

/* coordinate rotation */
void
rotate (double theta, double *x, double *y)
{
	double theta_rad = _deg2rad_ (theta);
	double x1 = (*x) * cos (theta_rad) + (*y) * sin (theta_rad);
	double y1 = (*y) * cos (theta_rad) - (*x) * sin (theta_rad);
	*x = x1;
	*y = y1;
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
check_singular_point (const fault_params *fault, double x, double y, double z, double eps)
{
	if (fabs (x + fault->flength1) < eps || fabs (x - fault->flength2) < eps) {
		int		i;
		for (i = 1; i <= 3; i++) {
			double	di = d[i] + ((i == 2) ? z : -z);
			double p = y * cd - di * sd;
			double q = y * sd + di * cd;
			if ((fabs (p + fault->fwidth1) < eps || fabs (p - fault->fwidth2) < eps) && fabs (q) < eps) singular_R[i] = true;
			if ((p + fault->fwidth1 < 0.0 || p < fault->fwidth2) && fabs (q) < eps) singular_RE[i] = true;
		}
	}
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
		ir  = 1.0 / r;
		ir3 = 1.0 / r3;
		ir5 = 1.0 / r5;
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
		ire   = 1.0 / re;
		ire2  = 1.0 / re2;
		ir3e2 = 1.0 / r3e2;
		ire3  = 1.0 / re3;
		ir5e3 = 1.0 / r5e3;
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
