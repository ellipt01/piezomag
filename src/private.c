/*
 * private.c
 *
 *  Created on: 2014/09/02
 *      Author: utsugi
 */

#include <stdio.h>
#include <math.h>

#include "piezomag.h"
#include "private.h"

/* clear all singular flags */
void
clear_all_singular_flags (void)
{
	singular_iR  = false;
	singular_iRX = false;
	singular_iRE = false;
	singular_iRC = false;
	return;
}

/* check all singular flags */
bool
is_singular_point (void)
{
	return (singular_iR || singular_iRX || singular_iRE || singular_iRC);
}

const double	eps = 5.e-4;

/*c*************************************************************
 * calculate some arithmetic constants
 * and store them in global variables
 ** INPUT **
 * sign: + or -
 * double xi, eta and qq: coordinates obs. pont on fault plane
 *c*************************************************************/
void
calc_geometry_variables (double sign, double xi, double et, double qq)
{
	double	r3, r5;
	double	rx2, re2, rc2;
	double	rx3, re3, rc3;
	double	r3x2, r3e2, r3c2;

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

	re = r + et;
	re2 = pow (re, 2.0);	// (r + et)^2
	re3 = pow (re, 3.0);	// (r + et)^2
	r3e2 = r3 * re2;		// r^3 * (r + et)^2

	rc = r + cc;
	rc2 = pow (rc, 2.0); // (r + cc)^2
	rc3 = pow (rc, 3.0); // (r + cc)^2
	r3c2 = r3 * rc2;		// r^3 * (r + cc)^2

	if (fabs (r) > eps) {
		ir  = 1.0 / r;
		ir3 = 1.0 / r3;
		ir5 = 1.0 / r5;
	} else {
		ir  = 0.0;
		ir3 = 0.0;
		ir5 = 0.0;
		if (!singular_iR) singular_iR = true;
}

	if (fabs (rx) > eps) {
		irx   = 1.0 / rx;
		irx2  = 1.0 / rx2;
		irx3  = 1.0 / rx3;
		ir3x2 = 1.0 / r3x2;
	} else {
		irx   = 0.0;
		irx2  = 0.0;
		irx3  = 0.0;
		ir3x2 = 0.0;
		if (!singular_iRX) singular_iRX = true;
	}

	if (fabs (re) > eps) {
		ire   = 1.0 / re;
		ire2  = 1.0 / re2;
		ir3e2 = 1.0 / r3e2;
		ire3  = 1.0 / re3;
	} else {
		ire   = 0.0;
		ire2  = 0.0;
		ire3  = 0.0;
		ir3e2 = 0.0;
		if (!singular_iRE) singular_iRE = true;
	}

	if (fabs (rc) > eps) {
		irc   = 1.0 / rc;
		irc2  = 1.0 / rc2;
		irc3  = 1.0 / rc3;
		ir3c2 = 1.0 / r3c2;
	} else {
		irc   = 0.0;
		irc2  = 0.0;
		irc3  = 0.0;
		ir3c2 = 0.0;
		if (!singular_iRC) singular_iRC = true;
	}

	return;
}
