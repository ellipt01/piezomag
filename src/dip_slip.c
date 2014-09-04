#include <stdio.h>
#include <math.h>
#include <float.h>

#include "piezomag.h"
#include "private.h"

/*** main term ***/
static double
dipx0 (MagComp component, double xi, double et, double qq)
{
	return -2.0 * K4 (component, 1.0, xi, et, qq);
}

static double
dipy0 (MagComp component, double xi, double et, double qq)
{
	return 2.0 * K5 (component, 1.0, xi, et, qq);
}

static double
dipz0 (MagComp component, double xi, double et, double qq)
{
	return 2.0 * K6 (component, 1.0, xi, et, qq);
}

static double
dip0 (MagComp component, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
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
		if (fabs (mag->cx) > DBL_EPSILON) hx = dipx0 (component, xi, et, q);
		if (fabs (mag->cy) > DBL_EPSILON) hy = dipy0 (component, xi, et, q);
		if (fabs (mag->cz) > DBL_EPSILON) hz = dipz0 (component, xi, et, q);

		res[i] = mag->cx * hx + mag->cy * hy + mag->cz * hz;
	}
	return (res[0] + res[3]) - (res[1] + res[2]);
}

/*** contributions from the mirror image H0 ***/
static double
dipxH0 (MagComp component, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K4_val = K4 (component, 1.0, xi, et, qq);
	double log_re_val = log_re (component, 1.0, xi, et, qq);
	double J2_val = J2 (component, 1.0, xi, et, qq) * (fault_is_vertical ? 1. : secd);
	double O1_val = O1 (component, 1.0, xi, et, qq);
	double N2_val = N2 (component, 1.0, xi, et, qq);
	double M2_val = M2 (component, 1.0, xi, et, qq);
	double O2z_val = O1z (component, 1.0, xi, et, qq);
	double N2z_val = N2z (component, 1.0, xi, et, qq);

	val =	(2.0 - alpha4) * K4_val
		+ alpha3 * (log_re_val * s2d - J2_val * c2d)
		+ alpha3 * (qd * O1_val + (z - h) * N2_val * sd)
		+ 2.0 * alpha5 * h * (O1_val * cd - N2_val * sd)
		+ 4.0 * alpha1 * h * M2_val
		- 2.0 * alpha2 * h * ((qd + h * cd) * O2z_val + (z - 2.0 * h) * N2z_val * sd);

	return val;
}

static double
dipyH0 (MagComp component, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K5_val = K5 (component, 1.0, xi, et, qq);
	double log_rx_val = log_rx (component, 1.0, xi, et, qq);
	double atan_xe_qr_val = atan_xe_qr (component, 1.0, xi, et, qq);
	double J1_val = J1 (component, 1.0, xi, et, qq) * (fault_is_vertical ? 1. : secd);
	double O2_val = O2 (component, 1.0, xi, et, qq);
	double O3_val = O3 (component, 1.0, xi, et, qq);
	double N1_val = N1 (component, 1.0, xi, et, qq);
	double L1_val = L1 (component, 1.0, xi, et, qq);
	double O2z_val = O2z (component, 1.0, xi, et, qq);
	double O3z_val = O3z (component, 1.0, xi, et, qq);
	double N1z_val = N1z (component, 1.0, xi, et, qq);

	val = - (2.0 - alpha4) * K5_val
		- fault->alpha * log_rx_val * sd + alpha4 * atan_xe_qr_val * cd
		- alpha3 * J1_val * c2d
		+ alpha3 * (qd * O2_val - (z - h) * (N1_val - O3_val) * sd)
		+ 2.0 * alpha4 * h * (O2_val * cd + (N1_val - O3_val) * sd)
		- 4.0 * alpha1 * h * L1_val * sd * cd
		- 2.0 * alpha2 * h * ((qd + h * cd) * O2z_val
				- (z - 2.0 * h) * (N1z_val - O3z_val) * sd
				+ O3_val * sd);

	return val;
}

static double
dipzH0 (MagComp component, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K6_val = K6 (component, 1.0, xi, et, qq);
	double log_rx_val = log_rx (component, 1.0, xi, et, qq);
	double atan_xe_qr_val = atan_xe_qr (component, 1.0, xi, et, qq);
	double O2_val = O2 (component, 1.0, xi, et, qq);
	double O3_val = O3 (component, 1.0, xi, et, qq);
	double M1_val = M1 (component, 1.0, xi, et, qq);
	double O2z_val = O2z (component, 1.0, xi, et, qq);
	double O3z_val = O3z (component, 1.0, xi, et, qq);

	val = - (2.0 + alpha5) * K6_val
		+ fault->alpha * log_rx_val * cd + alpha4 * atan_xe_qr_val * sd
		+ alpha3 * (qd * O3_val - (z - h) * O2_val * sd)
		- 2.0 * fault->alpha * h * (O3_val * cd + O2_val * sd)
		- 4.0 * alpha1 * h * M1_val * sd * cd
		+ 2.0 * alpha2 * h * ((qd + h * cd) * O3z_val
				- (z - 2.0 * h) * O2z_val * sd - O2_val * sd);

	return val;
}

static double
dipH0 (MagComp component, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
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
		if (fabs (mag->cx) > DBL_EPSILON) hx = dipxH0 (component, fault, mag, xi, et, q, y, z);
		if (fabs (mag->cy) > DBL_EPSILON) hy = dipyH0 (component, fault, mag, xi, et, q, y, z);
		if (fabs (mag->cz) > DBL_EPSILON) hz = dipzH0 (component, fault, mag, xi, et, q, y, z);

		res[i] = mag->cx * hx + mag->cy * hy + mag->cz * hz;
	}
	return (res[0] + res[3]) - (res[1] + res[2]);
}

/*** contributions from the mirror image HI ***/
static double
dipxHI (MagComp component, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K4_val = K4 (component, -1.0, xi, et, qq);
	double log_re_val = log_re (component, -1.0, xi, et, qq);
	double J2_val = J2 (component, -1.0, xi, et, qq) * (fault_is_vertical ? 1. : secd);
	double O1_val = O1 (component, -1.0, xi, et, qq);
	double N2_val = N2 (component, -1.0, xi, et, qq);

	val = - alpha4 * K4_val + alpha3 * (log_re_val * s2d + J2_val * c2d)
		- alpha2 * (qd * O1_val + (z - h) * N2_val * sd);

	return val;
}

static double
dipyHI (MagComp component, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K5_val = K5 (component, -1.0, xi, et, qq);
	double log_rx_val = log_rx (component, -1.0, xi, et, qq);
	double atan_xe_qr_val = atan_xe_qr (component, -1.0, xi, et, qq);
	double J1_val = J1 (component, -1.0, xi, et, qq) * (fault_is_vertical ? 1. : secd);
	double O2_val = O2 (component, -1.0, xi, et, qq);
	double O3_val = O3 (component, -1.0, xi, et, qq);
	double N1_val = N1 (component, -1.0, xi, et, qq);

	val = alpha4 * K5_val + fault->alpha * log_rx_val * sd
		+ alpha4 * atan_xe_qr_val * cd
		- alpha3 * J1_val * c2d
		+ alpha2 * (qd * O2_val - (z - h) * (N1_val - O3_val) * sd);

	return val;
}

static double
dipzHI (MagComp component, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K6_val = K6 (component, -1.0, xi, et, qq);
	double log_rx_val = log_rx (component, -1.0, xi, et, qq);
	double atan_xe_qr_val = atan_xe_qr (component, -1.0, xi, et, qq);
	double O2_val = O2 (component, -1.0, xi, et, qq);
	double O3_val = O3 (component, -1.0, xi, et, qq);

	val = - alpha5 * K6_val
		- fault->alpha * log_rx_val * cd + alpha4 * atan_xe_qr_val * sd
		- alpha2 * (qd * O3_val - (z - h) * O2_val * sd);

	return val;
}

static double
dipHI (MagComp component, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
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
		if (fabs (mag->cx) > DBL_EPSILON) hx = dipxHI (component, fault, mag, xi, et, q, y, z);
		if (fabs (mag->cy) > DBL_EPSILON) hy = dipyHI (component, fault, mag, xi, et, q, y, z);
		if (fabs (mag->cz) > DBL_EPSILON) hz = dipzHI (component, fault, mag, xi, et, q, y, z);

		res[i] = mag->cx * hx + mag->cy * hy + mag->cz * hz;
	}
	return (res[0] + res[3]) - (res[1] + res[2]);
}

/*** contributions from the mirror image HIII ***/
static double
dipxHIII (MagComp component, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K4_val = K4 (component, 1.0, xi, et, qq);
	double log_re_val = log_re (component, 1.0, xi, et, qq);
	double J2_val = J2 (component, 1.0, xi, et, qq) * (fault_is_vertical ? 1. : secd);
	double O1_val = O1 (component, 1.0, xi, et, qq);
	double N2_val = N2 (component, 1.0, xi, et, qq);

	val = alpha4 * K4_val - alpha3 * (log_re_val * s2d - J2_val * c2d)
		- alpha3 * (qd * O1_val + (z - h) * N2_val * sd);

	return val;
}

static double
dipyHIII (MagComp component, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K5_val = K5 (component, 1.0, xi, et, qq);
	double log_rx_val = log_rx (component, 1.0, xi, et, qq);
	double atan_xe_qr_val = atan_xe_qr (component, 1.0, xi, et, qq);
	double J1_val = J1 (component, 1.0, xi, et, qq) * (fault_is_vertical ? 1. : secd);
	double O2_val = O2 (component, 1.0, xi, et, qq);
	double O3_val = O3 (component, 1.0, xi, et, qq);
	double N1_val = N1 (component, 1.0, xi, et, qq);

	val = - alpha4 * K5_val + fault->alpha * log_rx_val * sd
		- alpha4 * atan_xe_qr_val * cd
		+ alpha3 * J1_val * c2d
		- alpha3 * (qd * O2_val - (z - h) * (N1_val - O3_val) * sd);

	return val;
}

static double
dipzHIII (MagComp component, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K6_val = K6 (component, 1.0, xi, et, qq);
	double log_rx_val = log_rx (component, 1.0, xi, et, qq);
	double atan_xe_qr_val = atan_xe_qr (component, 1.0, xi, et, qq);
	double O2_val = O2 (component, 1.0, xi, et, qq);
	double O3_val = O3 (component, 1.0, xi, et, qq);

	val = alpha5 * K6_val
		- fault->alpha * log_rx_val * cd - alpha4 * atan_xe_qr_val * sd
		- alpha3 * (qd * O3_val - (z - h) * O2_val * sd);

	return val;
}

static double
dipHIII (MagComp component, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
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
		if (fabs (mag->cx) > DBL_EPSILON) hx = dipxHIII (component, fault, mag, xi, et, q, y, z);
		if (fabs (mag->cy) > DBL_EPSILON) hy = dipyHIII (component, fault, mag, xi, et, q, y, z);
		if (fabs (mag->cz) > DBL_EPSILON) hz = dipzHIII (component, fault, mag, xi, et, q, y, z);

		res[i] = mag->cx * hx + mag->cy * hy + mag->cz * hz;
	}
	return (res[0] + res[3]) - (res[1] + res[2]);
}


static double
dipHII (MagComp component, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
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
		double sign = -1.0;

		xi = x + fault->flength1;
		et = p + fault->fwidth1;

		if (i >= 2)		 xi = x - fault->flength2;
		if (i % 2 == 0) et = p - w;

		calc_geometry_variables (sign, xi, et, q);
		if (fabs (mag->cx) > DBL_EPSILON) hx = dipxHI (component, fault, mag, xi, et, q, y, z);
		if (fabs (mag->cy) > DBL_EPSILON) hy = dipyHI (component, fault, mag, xi, et, q, y, z);
		if (fabs (mag->cz) > DBL_EPSILON) hz = dipzHI (component, fault, mag, xi, et, q, y, z);

		res[i] = mag->cx * hx + mag->cy * hy + mag->cz * hz;
	}

	p = y * cd - (d[1] - z) * sd;
	q = y * sd + (d[1] - z) * cd;

	for (i = 0; i < 4; i++) {
		double xi, et;
		double hx = 0., hy = 0., hz = 0.;
		double sign = 1.0;

		xi = x + fault->flength1;
		et = p - w;

		if (i >= 2)		 xi = x - fault->flength2;
		if (i % 2 == 0) et = p - fault->fwidth2;

		calc_geometry_variables (sign, xi, et, q);
		if (fabs (mag->cx) > DBL_EPSILON) hx = dipxHIII (component, fault, mag, xi, et, q, y, z);
		if (fabs (mag->cy) > DBL_EPSILON) hy = dipyHIII (component, fault, mag, xi, et, q, y, z);
		if (fabs (mag->cz) > DBL_EPSILON) hz = dipzHIII (component, fault, mag, xi, et, q, y, z);

		res[i] += mag->cx * hx + mag->cy * hy + mag->cz * hz;
	}
	return (res[0] + res[3]) - (res[1] + res[2]);
}

static double
dip_slip_main (MagComp component, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
{
	return dip0 (component, fault, mag, x, y, z);
}

static double
dip_slip_mirror_image (MagComp component, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
{
	return dipH0 (component, fault, mag, x, y, z);
}

static double
dip_slip_submirror_image (MagComp component, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
{
	double	val;
	if (fault->fdepth + fault->fwidth2 * sd <= mag->dcurier) val = dipHI (component, fault, mag, x, y, z);
	else if (fault->fdepth - fault->fwidth1 * sd >= mag->dcurier) val = dipHIII (component, fault, mag, x, y, z);
	else val = dipHII (component, fault, mag, x, y, z);
	return val;
}

/*** public functions ***/
double
dip_slip (MagComp component, SeismoMagTerm term, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
{
	double res = 0.0;
	if (!check_mag_component (component)) {
		fprintf (stderr, "ERROR: dip_slip: component must be MAG_COMP_X, MAG_COMP_Y, MAG_COMP_Z or MAG_COMP_F\n");
		return 0.0;
	}
	if (!check_seismo_mag_term (term)) {
		fprintf (stderr, "ERROR: dip_slip: term must be SEISMO_MAG_MAIN, SEISMO_MAG_MIRROR, SEISMO_MAG_SUBMIRROR or SEISMO_MAG_TOTAL\n");
		return 0.0;
	}
	if (term & SEISMO_MAG_MAIN) res += dip_slip_main (component, fault, mag, x, y, z);
	if (term & SEISMO_MAG_MIRROR) res += dip_slip_mirror_image (component, fault, mag, x, y, z);
	if (term & SEISMO_MAG_SUBMIRROR) res += dip_slip_submirror_image (component, fault, mag, x, y, z);
	return res;
}
