/*
 * piezomag.c
 *
 *  Created on: 2014/08/08
 *      Author: utsugi
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "piezomag.h"
#include "private.h"

/*
 * calculate total
 * double	hx:			x(EW) component
 * double	hy:			y(NS) component
 * double	hz:			z component
 * double	exf_inc:	inclination of external field
 * double	exf_dec:	declination of external field
 */
static double
total_force (double hx, double hy, double hz, double exf_inc, double exf_dec)
{
	double	f = hx * cos (deg2rad (exf_inc)) * sin (deg2rad (exf_dec))
		- hy * cos (deg2rad (exf_inc)) * cos (deg2rad (exf_dec))
		+ hz * sin (deg2rad (exf_inc));
	return f;
}

/*
 * calculate some constants which are used in piezomagnetic field estimation
 */
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

	if (singular_R[0]) {
		ir  = 0.0;
		ir3 = 0.0;
		ir5 = 0.0;
	} else {
		if (fabs (r) > DBL_EPSILON)  ir  = 1.0 / r;
		if (fabs (r3) > DBL_EPSILON) ir3 = 1.0 / r3;
		if (fabs (r5) > DBL_EPSILON) ir5 = 1.0 / r5;
	}

	if (fabs (rx) > DBL_EPSILON) {
		irx   = 1.0 / rx;
		irx2  = 1.0 / rx2;
		irx3  = 1.0 / rx3;
		ir3x2 = 1.0 / r3x2;
		ir5x3 = 1.0 / r5x3;
	} else {
		irx   = 0.0;
		irx2  = 0.0;
		irx3  = 0.0;
		ir3x2 = 0.0;
		ir5x3 = 0.0;
	}

	if (singular_RE[0]) {
		ire   = 0.0;
		ire2  = 0.0;
		ire3  = 0.0;
		ir3e2 = 0.0;
		ir5e3 = 0.0;
	} else {
		if (fabs (re) > DBL_EPSILON)   ire   = 1.0 / re;
		if (fabs (re2) > DBL_EPSILON)  ire2  = 1.0 / re2;
		if (fabs (r3e2) > DBL_EPSILON) ir3e2 = 1.0 / r3e2;
		if (fabs (re3) > DBL_EPSILON)  ire3  = 1.0 / re3;
		if (fabs (r5e3) > DBL_EPSILON) ir5e3 = 1.0 / r5e3;
	}

	if (fabs (rc) > DBL_EPSILON) {
		irc   = 1.0 / rc;
		irc2  = 1.0 / rc2;
		irc3  = 1.0 / rc3;
		ir3c2 = 1.0 / r3c2;
		ir5c3 = 1.0 / r5c3;
	} else {
		irc   = 0.0;
		irc2  = 0.0;
		irc3  = 0.0;
		ir3c2 = 0.0;
		ir5c3 = 0.0;
	}
	return;
}

/*
 * calculates specified component of seismo-magnetic effect on obs. point (xobs, yobs, zobs)
 * int	component:		X_COMP(0), Y_COMP(1), Z_COMP(2) or TOTAL_FORCE(3)
 */
static double
piezomagnetic_effect_component (int component, const fault_params *fault, const magnetic_params *mag, double xobs, double yobs, double zobs)
{
	double	val = 0.0;

	if (fabs (fault->u1) > DBL_EPSILON) val += fault->u1 * strike_slip (component, fault, mag, xobs, yobs, zobs);
	if (fabs (fault->u2) > DBL_EPSILON) val += fault->u2 * dip_slip (component, fault, mag, xobs, yobs, zobs);
	if (fabs (fault->u3) > DBL_EPSILON) val += fault->u3 * tensile_opening (component, fault, mag, xobs, yobs, zobs);

	return val;
}

/*
 * calculates specified X, Y, Z component or total force of seismo-magnetic effect on obs. point (xobs, yobs, zobs)
 * int	component:		X_COMP(0), Y_COMP(1), Z_COMP(2) or TOTAL_FORCE(3)
 */
double
piezomagnetic_effect (int component, const fault_params *fault, const magnetic_params *mag, double xobs, double yobs, double zobs)
{
	double	val;
	if (component < 0 || component >= 4) {
		fprintf (stderr, "ERROR: component of magnetic field must be X_COMP, Y_COMP, Z_COMP or TOTAL_FORCE\n");
		exit (1);
	}
	if (component == TOTAL_FORCE) {
		double	hx = piezomagnetic_effect_component (X_COMP, fault, mag, xobs, yobs, zobs);
		double	hy = piezomagnetic_effect_component (Y_COMP, fault, mag, xobs, yobs, zobs);
		double hz = piezomagnetic_effect_component (Z_COMP, fault, mag, xobs, yobs, zobs);
		if (fabs (fault->fstrike) > DBL_EPSILON) coordinates_transform (-fault->fstrike, &hx, &hy);
		val = total_force (hx, hy, hz, mag->exf_inc, mag->exf_dec);
	} else val = piezomagnetic_effect_component (component, fault, mag, xobs, yobs, zobs);

	return val;
}

/*
 * calculates and outputs seismo-magnetic effect in the range [xobs1 : dx : xobs2], [yobs1 : dy : yobs2] and z = zobs
 *
 */
void
fprintf_piezomagnetic_effect (FILE *stream, int component, const fault_params *fault, const magnetic_params *mag,
		double xobs1, double xobs2, double dx, double yobs1, double yobs2, double dy, double zobs)
{
	int		i, j;
	double x, y;
	int		n_grid_x = (int) floor ((xobs2 - xobs1) / dx);
	int		n_grid_y = (int) floor ((yobs2 - yobs1) / dy);
	double eps, grid = MIN (dx, dy);

	eps = 2.0 * grid;
	for (i = 0, x = xobs1; i <= n_grid_x; i++, x += dx) {
		for (j = 0, y = yobs1; j <= n_grid_y; j++, y += dy) {
			double tx, ty;
			double val;

			tx = x;
			ty = -y;
			if (fabs (fault->fstrike) > DBL_EPSILON) coordinates_transform (fault->fstrike, &tx, &ty);

			clear_all_singular_flag ();
			check_singular_point (fault, tx, ty, eps);
			if (is_singular_point (singular_R)) {
				if (verbos) {
					fprintf (stderr, "## SINGULAR: R: ");
					fprintf (stderr, "x = %f, y = %f", x, y);
					fprintf (stderr, " evaluation skiped ##\n");
				}
				continue;
			}
			if (is_singular_point (singular_RE)) {
				if (verbos) {
					fprintf (stderr, "## SINGULAR: R + ETA: ");
					fprintf (stderr, "x = %f, y = %f", x, y);
					fprintf (stderr, " evaluation skiped ##\n");
				}
				continue;
			}
			val = piezomagnetic_effect (component, fault, mag, tx, ty, z_obs);
			fprintf (stream, "%.4f\t%.4f\t%.8f\n", x, y, val);
		}
	}
	return;
}

