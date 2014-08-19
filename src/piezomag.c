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

/*c*********************************************************************
 * calculates specified component of seismo-magnetic field
 * on obs. point (xobs, yobs, zobs)
 * int	component:		X_COMP(0), Y_COMP(1), Z_COMP(2) or TOTAL_FORCE(3)
 *c*********************************************************************/
static double
piezomagnetic_effect_component (int component, const fault_params *fault, const magnetic_params *mag, double xobs, double yobs, double zobs)
{
	double	val = 0.0;

	if (fabs (fault->u1) > DBL_EPSILON) val += fault->u1 * strike_slip (component, fault, mag, xobs, yobs, zobs);
	if (fabs (fault->u2) > DBL_EPSILON) val += fault->u2 * dip_slip (component, fault, mag, xobs, yobs, zobs);
	if (fabs (fault->u3) > DBL_EPSILON) val += fault->u3 * tensile_opening (component, fault, mag, xobs, yobs, zobs);

	return val;
}

/*c******************************************************************
 * calculates specified component of seismo-magnetic field
 * on obs. point (xobs, yobs, zobs)
 * int	component:	X_COMP(0), Y_COMP(1), Z_COMP(2) or TOTAL_FORCE(3)
 *c******************************************************************/
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
		if (fabs (fault->fstrike) > DBL_EPSILON) rotate (-fault->fstrike, &hx, &hy);
		val = total_force (hx, hy, hz, mag->exf_inc, mag->exf_dec);
	} else val = piezomagnetic_effect_component (component, fault, mag, xobs, yobs, zobs);

	return val;
}

/*c**********************************************************************
 * calculates and outputs seismo-magnetic effect
 * in the range [xobs1 : dx : xobs2], [yobs1 : dy : yobs2] and z = zobs
 *c**********************************************************************/
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
			if (fabs (fault->fstrike) > DBL_EPSILON) rotate (fault->fstrike, &tx, &ty);

			clear_all_singular_flags ();
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

