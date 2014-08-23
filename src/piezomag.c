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
 * calculates specified component of seismomagnetic field
 * on obs. point (tx, ty, zobs) of internal fault coordinates
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

/*c******************************************************************
 * calculates specified component of seismomagnetic field
 * on obs. point (tx, ty, zobs) of internal fault coordinates
 *c******************************************************************/
static double
seismomagnetic_field_in_fault_coordinate (int component, const fault_params *fault, const magnetic_params *mag, double tx, double ty, double zobs)
{
	double	val;
	if (component == TOTAL_FORCE) {
		double	hx = seismomagnetic_component_in_fault_coordinate (X_COMP, fault, mag, tx, ty, zobs);
		double	hy = seismomagnetic_component_in_fault_coordinate (Y_COMP, fault, mag, tx, ty, zobs);
		double hz = seismomagnetic_component_in_fault_coordinate (Z_COMP, fault, mag, tx, ty, zobs);
		if (fabs (fault->fstrike) > DBL_EPSILON) rotate (-fault->fstrike, &hx, &hy);
		val = total_force (hx, hy, hz, mag->exf_inc, mag->exf_dec);
	} else val = seismomagnetic_component_in_fault_coordinate (component, fault, mag, tx, ty, zobs);

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
	bool	status = true;	// obs. point is singular point, status is set to false
	double	tx, ty;	// obs. point on fault coordinate system

	if (component < 0 || component >= 4) {
		fprintf (stderr, "ERROR: seismomagnetic_field: component must be X_COMP, Y_COMP, Z_COMP or TOTAL_FORCE\n");
		exit (1);
	}

	tx = yobs;
	ty = - xobs;
	if (fabs (fault->fstrike) > DBL_EPSILON) rotate (fault->fstrike, &tx, &ty);

	clear_all_singular_flags ();
	check_singular_point (fault, tx, ty, 1.e-4);
	if (is_singular_point (singular_R)) status = false;
	if (is_singular_point (singular_RE)) status = false;

	*val = seismomagnetic_field_in_fault_coordinate (component, fault, mag, tx, ty, zobs);

	return status;
}

/*c******************************************************************************
 * calculates and outputs seismomagnetic field
 * in the range [xobs1 : dx : xobs2](NS), [yobs1 : dy : yobs2](EW) and z = zobs
 *c******************************************************************************/
void
fprintf_seismomagnetic_field (FILE *stream, int component, const fault_params *fault, const magnetic_params *mag,
		double xobs1, double xobs2, double dx, double yobs1, double yobs2, double dy, double zobs)
{
	int		i, j;
	double x, y;
	int		n_grid_x = (int) floor ((xobs2 - xobs1) / dx);
	int		n_grid_y = (int) floor ((yobs2 - yobs1) / dy);
	double eps = 2. * MIN (dx, dy);

	for (i = 0, x = xobs1; i <= n_grid_x; i++, x += dx) {
		for (j = 0, y = yobs1; j <= n_grid_y; j++, y += dy) {
			double tx, ty;
			double val;

			tx = x;
			ty = y;
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
			val = seismomagnetic_field_in_fault_coordinate (component, fault, mag, tx, ty, zobs);
			fprintf (stream, "%.4f\t%.4f\t%.8f\n", x, y, val);
		}
	}
	return;
}
