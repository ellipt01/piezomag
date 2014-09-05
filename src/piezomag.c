/*
 * piezomag.c
 *
 *  Created on: 2014/08/08
 *      Author: utsugi
 */

#include <stdio.h>
#include <math.h>
#include <float.h>

#include "piezomag.h"
#include "private.h"

#define PIEZOMAG_MIN(x, y) ((x) <= (y)) ? (x) : (y)

typedef double (seismo_mag_term_func)
		(MagComp component,
		 const fault_params *fault, const magnetic_params *mag,
		 double xi, double et, double qq, double y, double z);

/*c*****************************************************
 * returns specified component of seismomagnetic term
 * calculated by seismomag_term_func *func.
 * resultant component is on geodetic coordinate system
 *c*****************************************************/
static double
seismomagnetic_component (MagComp component,
		const fault_params *fault, const magnetic_params *mag,
		double xi, double et, double qq, double y, double z,
		seismo_mag_term_func *func)
{
	double	val = 0.0;
	if (component == MAG_COMP_Z) {
		val = func (MAG_COMP_Z, fault, mag, xi, et, qq, y, z);
	} else if (component == MAG_COMP_X) {
		val = func (MAG_COMP_X, fault, mag, xi, et, qq, y, z);
		if (fabs (fault->fstrike) > DBL_EPSILON) {
			double	hy = func (MAG_COMP_Y, fault, mag, xi, et, qq, y, z);
			rotate (-fault->fstrike, &val, &hy);
		}
	} else if (component == MAG_COMP_Y) {
		val = func (MAG_COMP_Y, fault, mag, xi, et, qq, y, z);
		if (fabs (fault->fstrike) > DBL_EPSILON) {
			double	hx = func (MAG_COMP_X, fault, mag, xi, et, qq, y, z);
			rotate (-fault->fstrike, &hx, &val);
		}
	} else if (component == MAG_COMP_F) {
		double	hx = func (MAG_COMP_X, fault, mag, xi, et, qq, y, z);
		double	hy = func (MAG_COMP_Y, fault, mag, xi, et, qq, y, z);
		double	hz = func (MAG_COMP_Z, fault, mag, xi, et, qq, y, z);
		rotate (-fault->fstrike, &hx, &hy);
		val = total_force (hx, hy, hz, mag->exf_inc, mag->exf_dec);
	}
	return val;
}

/* main source */
static double
seismomagnetic0 (MagComp component,
		const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
{
	int		i;
	double res[4] = {0., 0., 0., 0.};
	double p, q;

	p = y * cd - (d[1] - z) * sd;
	q = y * sd + (d[1] - z) * cd;

	for (i = 0; i < 4; i++) {
		double xi, et;
		double sign = 1.0;
		double	val;

		xi = x;
		if (i < 2) xi += fault->flength1;
		else xi -= fault->flength2;

		et = p;
		if (i % 2 != 0) et += fault->fwidth1;
		else et -= fault->fwidth2;

		calc_geometry_variables (sign, xi, et, q);
		if (fabs (fault->u1) > DBL_EPSILON) {
			val = seismomagnetic_component (component, fault, mag, xi, et, q, y, z, &strike0);
			res[i] += fault->u1 * val;
		}
		if (fabs (fault->u2) > DBL_EPSILON) {
			val = seismomagnetic_component (component, fault, mag, xi, et, q, y, z, &dip0);
			res[i] += fault->u2 * val;
		}
		if (fabs (fault->u3) > DBL_EPSILON) {
			val = seismomagnetic_component (component, fault, mag, xi, et, q, y, z, &tensile0);
			res[i] += fault->u3 * val;
		}
	}
	return (res[0] + res[3]) - (res[1] + res[2]);
}

/* mirror image */
static double
seismomagneticH0 (MagComp component,
		const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
{
	int		i;
	double res[4] = {0., 0., 0., 0.};
	double p, q;

	p = y * cd - (d[3] - z) * sd;
	q = y * sd + (d[3] - z) * cd;

	for (i = 0; i < 4; i++) {
		double xi, et;
		double sign = 1.0;
		double	val;

		xi = x;
		if (i < 2) xi += fault->flength1;
		else xi -= fault->flength2;

		et = p;
		if (i % 2 != 0) et += fault->fwidth1;
		else et -= fault->fwidth2;

		calc_geometry_variables (sign, xi, et, q);
		if (fabs (fault->u1) > DBL_EPSILON) {
			val = seismomagnetic_component (component, fault, mag, xi, et, q, y, z, &strikeH0);
			res[i] += fault->u1 * val;
		}
		if (fabs (fault->u2) > DBL_EPSILON) {
			val = seismomagnetic_component (component, fault, mag, xi, et, q, y, z, &dipH0);
			res[i] += fault->u2 * val;
		}
		if (fabs (fault->u3) > DBL_EPSILON) {
			val = seismomagnetic_component (component, fault, mag, xi, et, q, y, z, &tensileH0);
			res[i] += fault->u3 * val;
		}
	}
	return (res[0] + res[3]) - (res[1] + res[2]);
}

/* sub-mirror image: type I */
double
seismomagneticHI (MagComp component,
		const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
{
	int		i;
	double res[4] = {0., 0., 0., 0.};
	double p, q;
	double	val;

	p = y * cd - (d[2] + z) * sd;
	q = y * sd + (d[2] + z) * cd;

	for (i = 0; i < 4; i++) {
		double xi, et;
		double sign = -1.0;

		xi = x;
		if (i < 2) xi += fault->flength1;
		else xi -= fault->flength2;

		et = p;
		if (i % 2 != 0) et += fault->fwidth1;
		else et -= fault->fwidth2;

		calc_geometry_variables (sign, xi, et, q);
		if (fabs (fault->u1) > DBL_EPSILON) {
			val = seismomagnetic_component (component, fault, mag, xi, et, q, y, z, &strikeHI);
			res[i] += fault->u1 * val;
		}
		if (fabs (fault->u2) > DBL_EPSILON) {
			val = seismomagnetic_component (component, fault, mag, xi, et, q, y, z, &dipHI);
			res[i] += fault->u2 * val;
		}
		if (fabs (fault->u3) > DBL_EPSILON) {
			val = seismomagnetic_component (component, fault, mag, xi, et, q, y, z, &tensileHI);
			res[i] += fault->u3 * val;
		}
	}
	return (res[0] + res[3]) - (res[1] + res[2]);
}

/* sub-mirror image: type III */
double
seismomagneticHIII (MagComp component,
		const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
{
	int	i;
	double res[4] = {0., 0., 0., 0.};
	double p, q;
	double	val;

	p = y * cd - (d[1] - z) * sd;
	q = y * sd + (d[1] - z) * cd;

	for (i = 0; i < 4; i++) {
		double xi, et;
		double sign = 1.0;

		xi = x;
		if (i < 2) xi += fault->flength1;
		else xi -= fault->flength2;

		et = p;
		if (i % 2 != 0) et += fault->fwidth1;
		else et -= fault->fwidth2;

		calc_geometry_variables (sign, xi, et, q);
		if (fabs (fault->u1) > DBL_EPSILON) {
			val = seismomagnetic_component (component, fault, mag, xi, et, q, y, z, &strikeHIII);
			res[i] += fault->u1 * val;
		}
		if (fabs (fault->u2) > DBL_EPSILON) {
			val = seismomagnetic_component (component, fault, mag, xi, et, q, y, z, &dipHIII);
			res[i] += fault->u2 * val;
		}
		if (fabs (fault->u3) > DBL_EPSILON) {
			val = seismomagnetic_component (component, fault, mag, xi, et, q, y, z, &tensileHIII);
			res[i] += fault->u3 * val;
		}
	}
	return (res[0] + res[3]) - (res[1] + res[2]);
}

/* sub-mirror image: type II */
double
seismomagneticHII (MagComp component,
		const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
{
	int	i;
	double res[4] = {0., 0., 0., 0.};
	double p, q;
	double w = (mag->dcurier - fault->fdepth) / sd;
	double	val;

	p = y * cd - (d[2] + z) * sd;
	q = y * sd + (d[2] + z) * cd;

	for (i = 0; i < 4; i++) {
		double xi, et;
		double	sign = -1.0;

		xi = x;
		if (i < 2) xi += fault->flength1;
		else xi -= fault->flength2;

		et = p;
		if (i % 2 != 0) et += fault->fwidth1;
		else et -= w;

		calc_geometry_variables (sign, xi, et, q);
		if (fabs (fault->u1) > DBL_EPSILON) {
			val = seismomagnetic_component (component, fault, mag, xi, et, q, y, z, &strikeHI);
			res[i] += fault->u1 * val;
		}
		if (fabs (fault->u2) > DBL_EPSILON) {
			val = seismomagnetic_component (component, fault, mag, xi, et, q, y, z, &dipHI);
			res[i] += fault->u2 * val;
		}
		if (fabs (fault->u3) > DBL_EPSILON) {
			val = seismomagnetic_component (component, fault, mag, xi, et, q, y, z, &tensileHI);
			res[i] += fault->u3 * val;
		}
	}

	p = y * cd - (d[1] - z) * sd;
	q = y * sd + (d[1] - z) * cd;

	for (i = 0; i < 4; i++) {
		double xi, et;
		double	sign = 1.0;

		xi = x;
		if (i < 2) xi += fault->flength1;
		else xi -= fault->flength2;

		et = p;
		if (i % 2 != 0) et -= w;
		else et -= fault->fwidth2;

		calc_geometry_variables (sign, xi, et, q);
		if (fabs (fault->u1) > DBL_EPSILON) {
			val = seismomagnetic_component (component, fault, mag, xi, et, q, y, z, &strikeHIII);
			res[i] += fault->u1 * val;
		}
		if (fabs (fault->u2) > DBL_EPSILON) {
			val = seismomagnetic_component (component, fault, mag, xi, et, q, y, z, &dipHIII);
			res[i] += fault->u2 * val;
		}
		if (fabs (fault->u3) > DBL_EPSILON) {
			val = seismomagnetic_component (component, fault, mag, xi, et, q, y, z, &tensileHIII);
			res[i] += fault->u3 * val;
		}
	}
	return (res[0] + res[3]) - (res[1] + res[2]);
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
seismomagnetic_field_term (MagComp component, SeismoMagTerm term,
		const fault_params *fault, const magnetic_params *mag, double xobs, double yobs, double zobs, double *val)
{
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

	*val = 0.0;
	if (term & SEISMO_MAG_MAIN) *val += seismomagnetic0 (component, fault, mag, tx, ty, zobs);
	if (term & SEISMO_MAG_MIRROR) *val += seismomagneticH0 (component, fault, mag, tx, ty, zobs);
	if (term & SEISMO_MAG_SUBMIRROR) {
		if (fault->fdepth + fault->fwidth2 * sd <= mag->dcurier) *val += seismomagneticHI (component, fault, mag, tx, ty, zobs);
		else if (fault->fdepth - fault->fwidth1 * sd >= mag->dcurier) *val += seismomagneticHIII (component, fault, mag, tx, ty, zobs);
		else *val += seismomagneticHII (component, fault, mag, tx, ty, zobs);
	}

	return (is_singular_point ()) ? false : true;
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
seismomagnetic_field (MagComp component,
		const fault_params *fault, const magnetic_params *mag,
		double xobs, double yobs, double zobs, double *val)
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
