#include <stdio.h>
#include <math.h>
#include <float.h>

#include "piezomag.h"
#include "private.h"

/*** main term ***/
static double
strikex0 (MagComp component, double xi, double et, double qq)
{
	return 2.0 * K1 (component, 1.0, xi, et, qq);
}

static double
strikey0 (MagComp component, double xi, double et, double qq)
{
	return 2.0 * K2 (component, 1.0, xi, et, qq);
}

static double
strikez0 (MagComp component, double xi, double et, double qq)
{
	return 2.0 * K3 (component, 1.0, xi, et, qq);
}

static double
strike0 (MagComp component, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
{
	int		i;
	double res[4];
	double p, q;

	p = y * cd - (d[1] - z) * sd;
	q = y * sd + (d[1] - z) * cd;

	for (i = 0; i < 4; i++) {
		double xi, et;
		double hx = 0., hy = 0., hz = 0.;
		double sign = 1.0;

		xi = x + fault->flength1;
		et = p + fault->fwidth1;

		if (i >= 2)		 xi = x - fault->flength2;
		if (i % 2 == 0) et = p - fault->fwidth2;

		calc_geometry_variables (sign, xi, et, q);
		if (fabs (mag->cx) > DBL_EPSILON) hx = strikex0 (component, xi, et, q);
		if (fabs (mag->cy) > DBL_EPSILON) hy = strikey0 (component, xi, et, q);
		if (fabs (mag->cz) > DBL_EPSILON) hz = strikez0 (component, xi, et, q);

		res[i] = mag->cx * hx + mag->cy * hy + mag->cz * hz;
	}
	return (res[0] + res[3]) - (res[1] + res[2]);
}

/*** contributions from the mirror image H0 ***/
static double
strikexH0 (MagComp component, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K1_val = K1 (component, 1.0, xi, et, qq);
	double atan_xe_qr_val = atan_xe_qr (component, 1.0, xi, et, qq);
	double J1_val = J1 (component, 1.0, xi, et, qq);
	double M1_val = M1 (component, 1.0, xi, et, qq);
	double L1_val = L1 (component, 1.0, xi, et, qq);
	double M1y_val = M1y (component, 1.0, xi, et, qq);
	double M1z_val = M1z (component, 1.0, xi, et, qq);

	val = - (2.0 - alpha4) * K1_val
		- alpha1 * atan_xe_qr_val - alpha3 * J1_val
		- alpha3 * (qd * M1_val + (z - h) * L1_val * sd)
		- 2.0 * alpha4 * h * (M1_val * cd - L1_val * sd)
		- 4.0 * alpha1 * h * L1_val * sd
		+ 2.0 * alpha2 * h * ((qd + h * cd) * M1z_val - (z - 2.0 * h) * M1y_val * sd);

	return val;
}

static double
strikeyH0 (MagComp component, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K2_val = K2 (component, 1.0, xi, et, qq);
	double log_re_val = log_re (component, 1.0, xi, et, qq);
	double J2_val = J2 (component, 1.0, xi, et, qq);
	double L2_val = L2 (component, 1.0, xi, et, qq);
	double M2_val = M2 (component, 1.0, xi, et, qq);
	double M3_val = M3 (component, 1.0, xi, et, qq);
	double M2y_val = M2y (component, 1.0, xi, et, qq);
	double M2z_val = M2z (component, 1.0, xi, et, qq);

	val = - (2.0 - alpha4) * K2_val
		+ alpha4 * log_re_val * sd + alpha3 * J2_val
		- alpha3 * (qd * M2_val + (z - h) * L2_val * sd * cd)
		+ 2.0 * alpha4 * h * (M2_val * cd - L2_val * sd * cd)
		- 4.0 * alpha1 * h * L2_val * sd * cd
		+ 2.0 * alpha2 * h * M3_val * sd
		+ 2.0 * alpha2 * h * ((qd + h * cd) * M2z_val - (z - 2.0 * h) * M2y_val * sd);

	return val;
}

static double
strikezH0 (MagComp component, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K3_val = K3 (component, 1.0, xi, et, qq);
	double log_re_val = log_re (component, 1.0, xi, et, qq);
	double M2_val = M2 (component, 1.0, xi, et, qq);
	double M3_val = M3 (component, 1.0, xi, et, qq);
	double M3y_val = M2z (component, 1.0, xi, et, qq);
	double M3z_val = M2y (component, 1.0, xi, et, qq);

	val = - (2.0 + alpha5) * K3_val
		- fault->alpha * log_re_val * cd
		- alpha3 * (qd * M3_val - (z - h) * M2_val * sd)
		+ 2.0 * fault->alpha * h * (M3_val * cd + M2_val * sd)
		- 4.0 * alpha1 * h * M2_val * sd
		- 2.0 * alpha2 * h * M2_val * sd
		- 2.0 * alpha2 * h * ((qd + h * cd) * M3z_val - (z - 2.0 * h) * M3y_val * sd);

	return val;
}

static double
strikeH0 (MagComp component, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
{
	int		i;
	double res[4];
	double p, q;

	p = y * cd - (d[3] - z) * sd;
	q = y * sd + (d[3] - z) * cd;

	for (i = 0; i < 4; i++) {
		double xi, et;
		double hx = 0., hy = 0., hz = 0.;
		double sign = 1.0;

		xi = x + fault->flength1;
		et = p + fault->fwidth1;

		if (i >= 2)		 xi = x - fault->flength2;
		if (i % 2 == 0) et = p - fault->fwidth2;

		calc_geometry_variables (sign, xi, et, q);
		if (fabs (mag->cx) > DBL_EPSILON) hx = strikexH0 (component, fault, mag, xi, et, q, y, z);
		if (fabs (mag->cy) > DBL_EPSILON) hy = strikeyH0 (component, fault, mag, xi, et, q, y, z);
		if (fabs (mag->cz) > DBL_EPSILON) hz = strikezH0 (component, fault, mag, xi, et, q, y, z);

		res[i] = mag->cx * hx + mag->cy * hy + mag->cz * hz;
	}
	return (res[0] + res[3]) - (res[1] + res[2]);
}

/*** contributions from the mirror image HI ***/
static double
strikexHI (MagComp component, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K1_val = K1 (component, -1.0, xi, et, qq);
	double atan_xe_qr_val = atan_xe_qr (component, -1.0, xi, et, qq);
	double J1_val = J1 (component, -1.0, xi, et, qq);
	double M1_val = M1 (component, -1.0, xi, et, qq);
	double L1_val = L1 (component, -1.0, xi, et, qq);

	val = alpha4 * K1_val
		- alpha1 * atan_xe_qr_val - alpha3 * J1_val
		- alpha2 * (qd * M1_val + (z - h) * L1_val * sd);

	return val;
}

static double
strikeyHI (MagComp component, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K2_val = K2 (component, -1.0, xi, et, qq);
	double log_re_val = log_re (component, -1.0, xi, et, qq);
	double J2_val = J2 (component, -1.0, xi, et, qq);
	double L2_val = L2 (component, -1.0, xi, et, qq);
	double M2_val = M2 (component, -1.0, xi, et, qq);

	val = alpha4 * K2_val + alpha6 * log_re_val * sd - alpha3 * J2_val
		- alpha2 * (qd * M2_val + (z - h) * L2_val * sd * cd);

	return val;
}

static double
strikezHI (MagComp component, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K3_val = K3 (component, -1.0, xi, et, qq);
	double log_re_val = log_re (component, -1.0, xi, et, qq);
	double M2_val = M2 (component, -1.0, xi, et, qq);
	double M3_val = M3 (component, -1.0, xi, et, qq);

	val = alpha5 * K3_val + fault->alpha * log_re_val * cd
		+ alpha2 * (qd * M3_val - (z - h) * M2_val * sd);

	return val;
}

static double
strikeHI (MagComp component, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
{
	int		i;
	double res[4];
	double p, q;

	p = y * cd - (d[2] + z) * sd;
	q = y * sd + (d[2] + z) * cd;

	for (i = 0; i < 4; i++) {
		double xi, et;
		double hx = 0., hy = 0., hz = 0.;
		double sign = -1.0;

		xi = x + fault->flength1;
		et = p + fault->fwidth1;

		if (i >= 2)		 xi = x - fault->flength2;
		if (i % 2 == 0) et = p - fault->fwidth2;

		calc_geometry_variables (sign, xi, et, q);
		if (fabs (mag->cx) > DBL_EPSILON) hx = strikexHI (component, fault, mag, xi, et, q, y, z);
		if (fabs (mag->cy) > DBL_EPSILON) hy = strikeyHI (component, fault, mag, xi, et, q, y, z);
		if (fabs (mag->cz) > DBL_EPSILON) hz = strikezHI (component, fault, mag, xi, et, q, y, z);

		res[i] = mag->cx * hx + mag->cy * hy + mag->cz * hz;
	}
	return (res[0] + res[3]) - (res[1] + res[2]);
}

/*** contributions from the mirror image HIII ***/
static double
strikexHIII (MagComp component, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K1_val = K1 (component, 1.0, xi, et, qq);
	double atan_xe_qr_val = atan_xe_qr (component, 1.0, xi, et, qq);
	double J1_val = J1 (component, 1.0, xi, et, qq);
	double M1_val = M1 (component, 1.0, xi, et, qq);
	double L1_val = L1 (component, 1.0, xi, et, qq);

	val = - alpha4 * K1_val
		+ alpha1 * atan_xe_qr_val + alpha3 * J1_val
		+ alpha3 * (qd * M1_val + (z - h) * L1_val * sd);

	return val;
}

static double
strikeyHIII (MagComp component, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K2_val = K2 (component, 1.0, xi, et, qq);
	double log_re_val = log_re (component, 1.0, xi, et, qq);
	double J2_val = J2 (component, 1.0, xi, et, qq);
	double L2_val = L2 (component, 1.0, xi, et, qq);
	double M2_val = M2 (component, 1.0, xi, et, qq);

	val = - alpha4 * K2_val - alpha4 * log_re_val * sd - alpha3 * J2_val
		+ alpha3 * (qd * M2_val + (z - h) * L2_val * sd * cd);

	return val;
}

static double
strikezHIII (MagComp component, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K3_val = K3 (component, 1.0, xi, et, qq);
	double log_re_val = log_re (component, 1.0, xi, et, qq);
	double M2_val = M2 (component, 1.0, xi, et, qq);
	double M3_val = M3 (component, 1.0, xi, et, qq);

	val = alpha5 * K3_val + fault->alpha * log_re_val * cd
		+ alpha3 * (qd * M3_val - (z - h) * M2_val * sd);

	return val;
}

static double
strikeHIII (MagComp component, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
{
	int		i;
	double res[4];
	double p, q;

	p = y * cd - (d[1] - z) * sd;
	q = y * sd + (d[1] - z) * cd;

	for (i = 0; i < 4; i++) {
		double xi, et;
		double hx = 0., hy = 0., hz = 0.;
		double sign = 1.0;

		xi = x + fault->flength1;
		et = p + fault->fwidth1;

		if (i >= 2)		 xi = x - fault->flength2;
		if (i % 2 == 0) et = p - fault->fwidth2;

		calc_geometry_variables (sign, xi, et, q);
		if (fabs (mag->cx) > DBL_EPSILON) hx = strikexHIII (component, fault, mag, xi, et, q, y, z);
		if (fabs (mag->cy) > DBL_EPSILON) hy = strikeyHIII (component, fault, mag, xi, et, q, y, z);
		if (fabs (mag->cz) > DBL_EPSILON) hz = strikezHIII (component, fault, mag, xi, et, q, y, z);

		res[i] = mag->cx * hx + mag->cy * hy + mag->cz * hz;
	}
	return (res[0] + res[3]) - (res[1] + res[2]);
}


static double
strikeHII (MagComp component, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
{
	int		i;
	double res[4];
	double p, q;
	double w = (mag->dcurier - fault->fdepth) / sd;

	p = y * cd - (d[2] + z) * sd;
	q = y * sd + (d[2] + z) * cd;

	for (i = 0; i < 4; i++) {
		double xi, et;
		double hx = 0., hy = 0., hz = 0.;
		double	sign = -1.0;

		xi = x + fault->flength1;
		et = p + fault->fwidth1;

		if (i >= 2)		 xi = x - fault->flength2;
		if (i % 2 == 0) et = p - w;

		calc_geometry_variables (sign, xi, et, q);
		if (fabs (mag->cx) > DBL_EPSILON) hx = strikexHI (component, fault, mag, xi, et, q, y, z);
		if (fabs (mag->cy) > DBL_EPSILON) hy = strikeyHI (component, fault, mag, xi, et, q, y, z);
		if (fabs (mag->cz) > DBL_EPSILON) hz = strikezHI (component, fault, mag, xi, et, q, y, z);

		res[i] = mag->cx * hx + mag->cy * hy + mag->cz * hz;
	}

	p = y * cd - (d[1] - z) * sd;
	q = y * sd + (d[1] - z) * cd;

	for (i = 0; i < 4; i++) {
		double xi, et;
		double hx = 0., hy = 0., hz = 0.;
		double	sign = 1.0;

		xi = x + fault->flength1;
		et = p - w;

		if (i >= 2)		 xi = x - fault->flength2;
		if (i % 2 == 0) et = p - fault->fwidth2;

		calc_geometry_variables (sign, xi, et, q);
		if (fabs (mag->cx) > DBL_EPSILON) hx = strikexHIII (component, fault, mag, xi, et, q, y, z);
		if (fabs (mag->cy) > DBL_EPSILON) hy = strikeyHIII (component, fault, mag, xi, et, q, y, z);
		if (fabs (mag->cz) > DBL_EPSILON) hz = strikezHIII (component, fault, mag, xi, et, q, y, z);

		res[i] += mag->cx * hx + mag->cy * hy + mag->cz * hz;
	}
	return (res[0] + res[3]) - (res[1] + res[2]);
}

static double
strike_slip_main (MagComp component, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
{
	return strike0 (component, fault, mag, x, y, z);
}

static double
strike_slip_mirror_image (MagComp component, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
{
	return strikeH0 (component, fault, mag, x, y, z);
}

static double
strike_slip_submirror_image (MagComp component, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
{
	double	val;
	if (fault->fdepth + fault->fwidth2 * sd <= mag->dcurier) val = strikeHI (component, fault, mag, x, y, z);
	else if (fault->fdepth - fault->fwidth1 * sd >= mag->dcurier) val = strikeHIII (component, fault, mag, x, y, z);
	else val = strikeHII (component, fault, mag, x, y, z);
	return val;
}

/*** public functions ***/
double
strike_slip (MagComp component, SeismoMagTerm term, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
{
	double res = 0.;
	if (!check_mag_component (component)) {
		fprintf (stderr, "ERROR: strike_slip: component must be MAG_COMP_X, MAG_COMP_Y, MAG_COMP_Z or MAG_COMP_F\n");
		return 0.0;
	}
	if (!check_seismo_mag_term (term)) {
		fprintf (stderr, "ERROR: strike_slip: term must be SEISMO_MAG_MAIN, SEISMO_MAG_MIRROR, SEISMO_MAG_SUBMIRROR or SEISMO_MAG_TOTAL\n");
		return 0.0;
	}
	if (term & SEISMO_MAG_MAIN) res += strike_slip_main (component, fault, mag, x, y, z);
	if (term & SEISMO_MAG_MIRROR) res += strike_slip_mirror_image (component, fault, mag, x, y, z);
	if (term & SEISMO_MAG_SUBMIRROR) res += strike_slip_submirror_image (component, fault, mag, x, y, z);
	return res;
}
