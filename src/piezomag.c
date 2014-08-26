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

#ifdef MIN
#undef MIN
#endif
#define PIEZOMAG_MIN(x, y) ((x) <= (y)) ? (x) : (y)

/* allowable distance between obs. and singular point.
 * default = 1.e-4 (km) */
double	eps_dist = 1.e-4;

/*c*********************************************************************
 * calculates specified component of seismomagnetic field
 * on obs. point (tx, ty, zobs) of internal fault coordinate system.
 * calculated magnetic component is also in fault coordinate system
  *c*********************************************************************/
static double
seismomagnetic_component_in_fault_coordinate (int component, const fault_params *fault, const magnetic_params *mag, double tx, double ty, double zobs)
{
	double	val = 0.0;

	if (fabs (fault->u1) > DBL_EPSILON) val += fault->u1 * strike_slip (component, fault, mag, tx, ty, zobs);
	if (fabs (fault->u2) > DBL_EPSILON) val += fault->u2 * dip_slip (component, fault, mag, tx, ty, zobs);
	if (fabs (fault->u3) > DBL_EPSILON) val += fault->u3 * tensile_opening (component, fault, mag, tx, ty, zobs);

	return val;
}

/*c*****************************************************************************
 * calculates specified component of seismomagnetic field
 * on obs. point (xobs, yobs, zobs)
 * int	component:	X_COMP(0:NS), Y_COMP(1:EW), Z_COMP(2:DwonUp) or TOTAL_FORCE(3)
 *c*****************************************************************************/
bool
seismomagnetic_field (int component, const fault_params *fault, const magnetic_params *mag, double xobs, double yobs, double zobs, double *val)
{
	bool	status = true;	// if obs. point is singular point, status is set to false
	double	tx, ty;	// obs. point on fault coordinate system

	if (component < 0 || component >= 4) {
		fprintf (stderr, "ERROR: seismomagnetic_field: component must be X_COMP, Y_COMP, Z_COMP or TOTAL_FORCE\n");
		exit (1);
	}

	// rotate coordinate system
	tx = xobs;
	ty = yobs;
	if (fabs (fault->fstrike) > DBL_EPSILON) rotate (fault->fstrike, &tx, &ty);

	clear_all_singular_flags ();
	check_singular_point (fault, tx, ty, eps_dist);
	if (is_singular_point (singular_R)) status = false;
	if (is_singular_point (singular_RE)) status = false;

	if (component == Z_COMP) {
		*val = seismomagnetic_component_in_fault_coordinate (Z_COMP, fault, mag, tx, ty, zobs);
	} else if (component == X_COMP) {
		double	hx = seismomagnetic_component_in_fault_coordinate (X_COMP, fault, mag, tx, ty, zobs);
		if (fabs (fault->fstrike) > DBL_EPSILON) {
			double	hy = seismomagnetic_component_in_fault_coordinate (Y_COMP, fault, mag, tx, ty, zobs);
			rotate (-fault->fstrike, &hx, &hy);
		}
		*val = hx;
	} else if (component == Y_COMP) {
		double	hy = seismomagnetic_component_in_fault_coordinate (Y_COMP, fault, mag, tx, ty, zobs);
		if (fabs (fault->fstrike) > DBL_EPSILON) {
			double	hx = seismomagnetic_component_in_fault_coordinate (X_COMP, fault, mag, tx, ty, zobs);
			rotate (-fault->fstrike, &hx, &hy);
		}
		*val = hy;
	} else if (component == TOTAL_FORCE) {
		double	hx = seismomagnetic_component_in_fault_coordinate (X_COMP, fault, mag, tx, ty, zobs);
		double	hy = seismomagnetic_component_in_fault_coordinate (Y_COMP, fault, mag, tx, ty, zobs);
		double hz = seismomagnetic_component_in_fault_coordinate (Z_COMP, fault, mag, tx, ty, zobs);
		if (fabs (fault->fstrike) > DBL_EPSILON) rotate (-fault->fstrike, &hx, &hy);
		*val = total_force (hx, hy, hz, mag->exf_inc, mag->exf_dec);
	}

	return status;
}

/*c**********************************************************************
 * calculates and outputs seismomagnetic field
 * in the range [xobs1:dx:xobs2](NS), [yobs1:dy:yobs2](EW) and z = zobs
 *c**********************************************************************/
void
fprintf_seismomagnetic_field (FILE *stream, int component, const fault_params *fault, const magnetic_params *mag,
		double xobs1, double xobs2, double dx, double yobs1, double yobs2, double dy, double zobs)
{
	int		i, j;
	double x, y;
	int		n_grid_x = (int) floor ((xobs2 - xobs1) / dx);
	int		n_grid_y = (int) floor ((yobs2 - yobs1) / dy);
	bool	status;

	/* eps_dist is set to 2 * (minimum grid interval) */
	eps_dist = 2. * ((double) PIEZOMAG_MIN (dx, dy));
//	fprintf (stderr, "1: eps = %.4e, %f, %f\n", eps_dist, dx, dy);
	for (i = 0, x = xobs1; i <= n_grid_x; i++, x += dx) {
		for (j = 0, y = yobs1; j <= n_grid_y; j++, y += dy) {
			double val;

			clear_all_singular_flags ();
			status = seismomagnetic_field (component, fault, mag, x, y, zobs, &val);
			if (!status && is_singular_point (singular_R)) {
				if (verbos) {
					fprintf (stderr, "## SINGULAR: R: ");
					fprintf (stderr, "x = %f, y = %f", x, y);
					fprintf (stderr, " evaluation skiped ##\n");
				}
				continue;
			}
			if (!status && is_singular_point (singular_RE)) {
				if (verbos) {
					fprintf (stderr, "## SINGULAR: R + ETA: ");
					fprintf (stderr, "x = %f, y = %f", x, y);
					fprintf (stderr, " evaluation skiped ##\n");
				}
				continue;
			}
			fprintf (stream, "%.4f\t%.4f\t%.8f\n", x, y, val);
		}
	}
	return;
}
