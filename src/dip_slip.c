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

/*** seismomagnetic terms on fault coordinate system due to dip-slip fault ***/

/* The following functions dose not calculate geometry variables such as r, rx, re, ..., etc,
 * only refer the global variables. So, before calling the following functions,
 * calc_geometry_variables() function must be called. */

/* main source */
double
dip0 (MagComp component, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double	val = 0.;
	if (fabs (mag->cx) > DBL_EPSILON) val += mag->cx * dipx0 (component, xi, et, qq);
	if (fabs (mag->cy) > DBL_EPSILON) val += mag->cy * dipy0 (component, xi, et, qq);
	if (fabs (mag->cz) > DBL_EPSILON) val += mag->cz * dipz0 (component, xi, et, qq);
	return val;
}

/* mirror image */
double
dipH0 (MagComp component, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double	val = 0.;
	if (fabs (mag->cx) > DBL_EPSILON) val += mag->cx * dipxH0 (component, fault, mag, xi, et, qq, y, z);
	if (fabs (mag->cy) > DBL_EPSILON) val += mag->cy * dipyH0 (component, fault, mag, xi, et, qq, y, z);
	if (fabs (mag->cz) > DBL_EPSILON) val += mag->cz * dipzH0 (component, fault, mag, xi, et, qq, y, z);
	return val;
}

/* sub-mirror image: type I */
double
dipHI (MagComp component, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double	val = 0.;
	if (fabs (mag->cx) > DBL_EPSILON) val += mag->cx * dipxHI (component, fault, mag, xi, et, qq, y, z);
	if (fabs (mag->cy) > DBL_EPSILON) val += mag->cy * dipyHI (component, fault, mag, xi, et, qq, y, z);
	if (fabs (mag->cz) > DBL_EPSILON) val += mag->cz * dipzHI (component, fault, mag, xi, et, qq, y, z);
	return val;
}

/* sub-mirror image type: III */
double
dipHIII (MagComp component, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double	val = 0.;
	if (fabs (mag->cx) > DBL_EPSILON) val += mag->cx * dipxHIII (component, fault, mag, xi, et, qq, y, z);
	if (fabs (mag->cy) > DBL_EPSILON) val += mag->cy * dipyHIII (component, fault, mag, xi, et, qq, y, z);
	if (fabs (mag->cz) > DBL_EPSILON) val += mag->cz * dipzHIII (component, fault, mag, xi, et, qq, y, z);
	return val;
}
