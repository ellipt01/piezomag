#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "constants.h"
#include "piez.h"

#define DUMMY 0
#define MIN(x, y) ((x) <= (y)) ? (x) : (y)

bool
set_constants (void)
{
	double	c0, jx, jy, jz;

	fstrike = 90.0 - fstrike;

	sd = sin (deg2rad (fdip));
	cd = cos (deg2rad (fdip));
	td = tan (deg2rad (fdip));
	secd = 1.0 / cd;
	sd2 = pow (sd, 2.0);
	cd2 = pow (cd, 2.0);
	s2d = sin (deg2rad (2.0 * fdip));
	c2d = cos (deg2rad (2.0 * fdip));

	d[0] = DUMMY;
	d[1] = fdepth - z_obs;
	d[2] = fdepth - 2.0 * dcurier + z_obs;
	d[3] = fdepth + 2.0 * dcurier - z_obs;

	alpha = (lambda + mu) / (lambda + 2.0 * mu);

	alpha0 = 4.0 * alpha - 1.0;
	alpha1 = 3.0 * alpha / alpha0;

	alpha2 = 6.0 * alpha * alpha / alpha0;
	alpha3 = 2.0 * alpha * (1.0 - alpha) / alpha0;
	alpha4 = alpha * (2.0 * alpha + 1.0) / alpha0;
	alpha5 = alpha * (2.0 * alpha - 5.0) / alpha0;
	alpha6 = 3.0 * alpha * (1.0 - 2.0 * alpha) / alpha0;

	jx = mgz_int * cos (deg2rad (mgz_inc)) * sin (deg2rad (mgz_dec));
	jy = - mgz_int * cos (deg2rad (mgz_inc)) * cos (deg2rad (mgz_dec));
	jz = mgz_int * sin (deg2rad (mgz_inc));
	coordinates_transform (fstrike, &jx, &jy);

	c0 = 0.25 * beta * mu * (3.0 * lambda + 2.0 * mu) / (lambda + mu);
	cx = c0 * jx;
	cy = c0 * jy;
	cz = c0 * jz;

	return true;
}

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

double
total_force (double hx, double hy, double hz, double exf_inc, double exf_dec)
{
	double	f = hx * cos (deg2rad (exf_inc)) * sin (deg2rad (exf_dec))
		- hy * cos (deg2rad (exf_inc)) * cos (deg2rad (exf_dec))
		+ hz * sin (deg2rad (exf_inc));
	return f;
}

double
piezomagnetic_effect_component (int component, double u1, double u2, double u3, double x, double y, double z)
{
	double	val = 0.0;

	if (fabs (u1) > DBL_EPSILON) val += u1 * strike_slip (component, x, y, z);
	if (fabs (u2) > DBL_EPSILON) val += u2 * dip_slip (component, x, y, z);
	if (fabs (u3) > DBL_EPSILON) val += u3 * tensile_opening (component, x, y, z);

	return val;
}

double
piezomagnetic_effect (int component, double u1, double u2, double u3, double x, double y, double z)
{
	double	val;
	if (component < 0 || component >= 4) {
		fprintf (stderr, "ERROR: component of magnetic field must be X_COMP, Y_COMP, Z_COMP or TOTAL_FORCE\n");
		exit (1);
	}
	if (component == TOTAL_FORCE) {
		double	hx = piezomagnetic_effect_component (X_COMP, u1, u2, u3, x, y, z);
		double	hy = piezomagnetic_effect_component (Y_COMP, u1, u2, u3, x, y, z);
		double hz = piezomagnetic_effect_component (Z_COMP, u1, u2, u3, x, y, z);
		if (fabs (fstrike) > DBL_EPSILON) coordinates_transform (-fstrike, &hx, &hy);
		val = total_force (hx, hy, hz, exf_inc, exf_dec);
	} else
		val = piezomagnetic_effect_component (component, u1, u2, u3, x, y, z);

	return val;
}

void
fprintf_piezomagnetic_effect (FILE *stream, int component, double x1, double x2, double dx, double y1, double y2, double dy, double z)
{
	int		i, j;
	double x, y;
	int		n_grid_x = (int) floor ((x02 - x01) / dx);
	int		n_grid_y = (int) floor ((y02 - y01) / dy);
	double eps, grid = MIN (dx, dy);

	eps = 2.0 * grid;
	for (i = 0, x = x01; i <= n_grid_x; i++, x += dx) {
		for (j = 0, y = y01; j <= n_grid_y; j++, y += dy) {
			double tx, ty;
			double val;

			tx = x;
			ty = -y;
			if (fabs (fstrike) > DBL_EPSILON) coordinates_transform (fstrike, &tx, &ty);

			clear_all_singular_flag ();
			check_singular_point (tx, ty, eps);
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
			val = piezomagnetic_effect (component, u1, u2, u3, tx, ty, z_obs);
			fprintf (stream, "%.4f\t%.4f\t%.8f\n", x, y, val);
		}
	}
	return;
}
