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

/*** seismomagnetic terms on fault coordinate system due to strike-slip fault ***/

/* The following functions dose not calculate geometry variables such as r, rx, re, ..., etc,
 * only refer the global variables. So, before calling the following functions,
 * calc_geometry_variables() function must be called. */

/* main source */
double
strike0 (MagComp component, const magnetic_params *mag, double xi, double et, double qq)
{
	double	val = 0.;
	if (fabs (mag->cx) > DBL_EPSILON) val += mag->cx * strikex0 (component, xi, et, qq);
	if (fabs (mag->cy) > DBL_EPSILON) val += mag->cy * strikey0 (component, xi, et, qq);
	if (fabs (mag->cz) > DBL_EPSILON) val += mag->cz * strikez0 (component, xi, et, qq);
	return val;
}

/* mirror image */
double
strikeH0 (MagComp component, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double	val = 0.;
	if (fabs (mag->cx) > DBL_EPSILON) val += mag->cx * strikexH0 (component, fault, mag, xi, et, qq, y, z);
	if (fabs (mag->cy) > DBL_EPSILON) val += mag->cy * strikeyH0 (component, fault, mag, xi, et, qq, y, z);
	if (fabs (mag->cz) > DBL_EPSILON) val += mag->cz * strikezH0 (component, fault, mag, xi, et, qq, y, z);
	return val;
}

/* sub-mirror image: type I */
double
strikeHI (MagComp component, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double	val = 0.;
	if (fabs (mag->cx) > DBL_EPSILON) val += mag->cx * strikexHI (component, fault, mag, xi, et, qq, y, z);
	if (fabs (mag->cy) > DBL_EPSILON) val += mag->cy * strikeyHI (component, fault, mag, xi, et, qq, y, z);
	if (fabs (mag->cz) > DBL_EPSILON) val += mag->cz * strikezHI (component, fault, mag, xi, et, qq, y, z);
	return val;
}

/* sub-mirror image type: III */
double
strikeHIII (MagComp component, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double	val = 0.;
	if (fabs (mag->cx) > DBL_EPSILON) val += mag->cx * strikexHIII (component, fault, mag, xi, et, qq, y, z);
	if (fabs (mag->cy) > DBL_EPSILON) val += mag->cy * strikeyHIII (component, fault, mag, xi, et, qq, y, z);
	if (fabs (mag->cz) > DBL_EPSILON) val += mag->cz * strikezHIII (component, fault, mag, xi, et, qq, y, z);
	return val;
}
