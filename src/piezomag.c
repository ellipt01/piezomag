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

#define PIEZOMAG_MIN(x, y) ((x) <= (y)) ? (x) : (y)

/*c*********************************************************************
 * calculates specified component of seismomagnetic field
 * on obs. point (tx, ty, zobs) of internal fault coordinate system.
 * calculated magnetic component is also in fault coordinate system
  *c*********************************************************************/
static double
seismomagnetic_component_in_fault_coordinate (MagComp component, SeismoMagTerm term, const fault_params *fault, const magnetic_params *mag, double tx, double ty, double zobs)
{
	double	val = 0.0;

	if (fabs (fault->u1) > DBL_EPSILON) val += fault->u1 * strike_slip (component, term, fault, mag, tx, ty, zobs);
	if (fabs (fault->u2) > DBL_EPSILON) val += fault->u2 * dip_slip (component, term, fault, mag, tx, ty, zobs);
	if (fabs (fault->u3) > DBL_EPSILON) val += fault->u3 * tensile_opening (component, term, fault, mag, tx, ty, zobs);

	return val;
}

/*** public functions ***/

/*c***************************************************************************
 * calculates specified component and term of seismomagnetic field
 * on obs. point (xobs, yobs, zobs)
 * MagComp component: output magnetic component
 *         MAG_COMP_X(1:NS)
 *         MAG_COMP_Y(2:EW)
 *         MAG_COMP_Z(3:DwonUp)
 *         MAG_COMP_F(0)
 * SeismoMagTerm term: output seismomagnetic term
 *         SEISMO_MAG_MAIN      (0: main term (0))
 *         SEISMO_MAG_MIRROR    (1: mirror image term (H0))
 *         SEISMO_MAG_SUBMIRROR (2: sub-mirror image term (HI, HIII or HII))
 *         SEISMO_MAG_TOTAL     (3: total seismomagnetic field)
 *c***************************************************************************/
bool
seismomagnetic_field_term (MagComp component, SeismoMagTerm term, const fault_params *fault, const magnetic_params *mag, double xobs, double yobs, double zobs, double *val)
{
	bool	status = true;	// if obs. point is singular point, status is set to false
	double	tx, ty;	// obs. point on fault coordinate system

	// z_obs must be < 0, i.e. outside of medium
	if (zobs >= 0.) {
		fprintf (stderr, "ERROR: seismomagnetic_field_term: z_obs must be < 0.\n");
		fprintf (stderr, "observation point must be located outside the medium.\n");
		return false;
	}
	// component must be X_COMP(0:NS), Y_COMP(1:EW), Z_COMP(2:DownUp) or TOTAL_FORCE(3)
	if (!check_mag_component (component)) {
		fprintf (stderr, "ERROR: seismomagnetic_field_term: component must be MAG_COMP_X, MAG_COMP_Y, MAG_COMP_Z or MAG_COMP_F\n");
		return false;
	}
	// term must be SEISMO_MAG_MAIN, SEISMO_MAG_MIRROR, SEISMO_MAG_SUBMIRROR or SEISMO_MAG_TOTAL
	if (!check_seismo_mag_term (term)) {
		fprintf (stderr, "ERROR: seismomagnetic_field_term: term must be SEISMO_MAG_MAIN, SEISMO_MAG_MIRROR, SEISMO_MAG_SUBMIRROR\n");
		return false;
	}

	// rotate coordinate system
	tx = xobs;
	ty = yobs;
	if (fabs (fault->fstrike) > DBL_EPSILON) rotate (fault->fstrike, &tx, &ty);

	clear_all_singular_flags ();

	if (component == MAG_COMP_Z) {
		*val = seismomagnetic_component_in_fault_coordinate (MAG_COMP_Z, term, fault, mag, tx, ty, zobs);
	} else if (component == MAG_COMP_X) {
		double	hx = seismomagnetic_component_in_fault_coordinate (MAG_COMP_X, term, fault, mag, tx, ty, zobs);
		if (fabs (fault->fstrike) > DBL_EPSILON) {
			double	hy = seismomagnetic_component_in_fault_coordinate (MAG_COMP_Y, term, fault, mag, tx, ty, zobs);
			rotate (-fault->fstrike, &hx, &hy);
		}
		*val = hx;
	} else if (component == MAG_COMP_Y) {
		double	hy = seismomagnetic_component_in_fault_coordinate (MAG_COMP_Y, term, fault, mag, tx, ty, zobs);
		if (fabs (fault->fstrike) > DBL_EPSILON) {
			double	hx = seismomagnetic_component_in_fault_coordinate (MAG_COMP_X, term, fault, mag, tx, ty, zobs);
			rotate (-fault->fstrike, &hx, &hy);
		}
		*val = hy;
	} else if (component == MAG_COMP_F) {
		double	hx = seismomagnetic_component_in_fault_coordinate (MAG_COMP_X, term, fault, mag, tx, ty, zobs);
		double	hy = seismomagnetic_component_in_fault_coordinate (MAG_COMP_Y, term, fault, mag, tx, ty, zobs);
		double hz = seismomagnetic_component_in_fault_coordinate (MAG_COMP_Z, term, fault, mag, tx, ty, zobs);
		if (fabs (fault->fstrike) > DBL_EPSILON) rotate (-fault->fstrike, &hx, &hy);
		*val = total_force (hx, hy, hz, mag->exf_inc, mag->exf_dec);
	}
	if (is_singular_point ()) status = false;

	return status;
}

/*c*******************************************************
 * calculates specified component of seismomagnetic field
 * on obs. point (xobs, yobs, zobs)
 * MagComp component: output magnetic component
 *         MAG_COMP_X(1:NS)
 *         MAG_COMP_Y(2:EW)
 *         MAG_COMP_Z(3:DwonUp)
 *         MAG_COMP_F(0)
 *c*******************************************************/
bool
seismomagnetic_field (MagComp component, const fault_params *fault, const magnetic_params *mag, double xobs, double yobs, double zobs, double *val)
{
	// z_obs must be < 0, i.e. outside of medium
	if (zobs >= 0.) {
		fprintf (stderr, "ERROR: seismomagnetic_field: z_obs must be < 0.\n");
		fprintf (stderr, "observation point must be located outside the medium.\n");
		return false;
	}
	// component must be X_COMP(0:NS), Y_COMP(1:EW), Z_COMP(2:DownUp) or TOTAL_FORCE(3)
	if (!check_mag_component (component)) {
		fprintf (stderr, "ERROR: seismomagnetic_field: component must be MAG_COMP_X, MAG_COMP_Y, MAG_COMP_Z or MAG_COMP_F\n");
		return false;
	}
	return seismomagnetic_field_term (component, SEISMO_MAG_TOTAL, fault, mag, xobs, yobs, zobs, val);
}

/*c****************************************************************************
 * calculates and outputs specified component and term of seismomagnetic field
 * in the range [xobs1:dx:xobs2](NS), [yobs1:dy:yobs2](EW) and z = zobs
 *c****************************************************************************/
void
fprintf_seismomagnetic_field_term (FILE *stream, MagComp component, SeismoMagTerm term,
		const fault_params *fault, const magnetic_params *mag,
		double xobs1, double xobs2, double dx, double yobs1, double yobs2, double dy, double zobs)
{
	int		i, j;
	double x, y;
	int		n_grid_x, n_grid_y;
	bool	status;

	// z_obs must be < 0, i.e. outside of medium
	if (zobs >= 0.) {
		fprintf (stderr, "ERROR: fprintf_seismomagnetic_field_term: z_obs must be < 0.\n");
		fprintf (stderr, "observation point must be located outside the medium.\n");
		return;
	}
	// component must be X_COMP(0:NS), Y_COMP(1:EW), Z_COMP(2:DownUp) or TOTAL_FORCE(3)
	if (!check_mag_component (component)) {
		fprintf (stderr, "ERROR: fprintf_seismomagnetic_field_term: component must be MAG_COMP_X, MAG_COMP_Y, MAG_COMP_Z or MAG_COMP_F\n");
		return;
	}
	// term must be SEISMO_MAG_MAIN, SEISMO_MAG_MIRROR, SEISMO_MAG_SUBMIRROR or SEISMO_MAG_TOTAL
	if (!check_seismo_mag_term (term)) {
		fprintf (stderr, "ERROR: fprintf_seismomagnetic_field_term: term must be SEISMO_MAG_MAIN, SEISMO_MAG_MIRROR, SEISMO_MAG_SUBMIRROR\n");
		return;
	}

	n_grid_x = (int) floor ((xobs2 - xobs1) / dx);
	n_grid_y = (int) floor ((yobs2 - yobs1) / dy);

	for (i = 0, x = xobs1; i <= n_grid_x; i++, x += dx) {
		for (j = 0, y = yobs1; j <= n_grid_y; j++, y += dy) {
			double val;

			clear_all_singular_flags ();
			status = seismomagnetic_field_term (component, term, fault, mag, x, y, zobs, &val);
			if (!status) continue;
			fprintf (stream, "%.4f\t%.4f\t%.8f\n", x, y, val);
		}
	}
	return;
}

/*c**********************************************************************
 * calculates and outputs specified component of seismomagnetic field
 * in the range [xobs1:dx:xobs2](NS), [yobs1:dy:yobs2](EW) and z = zobs
 *c**********************************************************************/
void
fprintf_seismomagnetic_field (FILE *stream, MagComp component,
		const fault_params *fault, const magnetic_params *mag,
		double xobs1, double xobs2, double dx, double yobs1, double yobs2, double dy, double zobs)
{
	// z_obs must be < 0, i.e. outside of medium
	if (zobs >= 0.) {
		fprintf (stderr, "ERROR: fprintf_seismomagnetic_field: z_obs must be < 0.\n");
		fprintf (stderr, "observation point must be located outside the medium.\n");
		return;
	}
	// component must be X_COMP(0:NS), Y_COMP(1:EW), Z_COMP(2:DownUp) or TOTAL_FORCE(3)
	if (!check_mag_component (component)) {
		fprintf (stderr, "ERROR: fprintf_seismomagnetic_field: component must be MAG_COMP_X, MAG_COMP_Y, MAG_COMP_Z or MAG_COMP_F\n");
		return;
	}

	fprintf_seismomagnetic_field_term (stream, component, SEISMO_MAG_TOTAL, fault, mag, xobs1, xobs2, dx, yobs1, yobs2, dy, zobs);
	return;
}
