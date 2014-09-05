#include <stdio.h>
#include <math.h>
#include <float.h>

#include "piezomag.h"
#include "private.h"

/*** main term ***/
double
tensilex0 (MagComp component, double xi, double et, double qq)
{
	return -2.0 * K7 (component, 1.0, xi, et, qq);
}

double
tensiley0 (MagComp component, double xi, double et, double qq)
{
	return 2.0 * K8 (component, 1.0, xi, et, qq);
}

double
tensilez0 (MagComp component, double xi, double et, double qq)
{
	return 2.0 * K9 (component, 1.0, xi, et, qq);
}

/*** contributions from the mirror image H0 ***/
static double
tensilexH0 (MagComp component, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K7_val = K7 (component, 1.0, xi, et, qq);
	double log_re_val = log_re (component, 1.0, xi, et, qq);
	double log_rc_val = log_rc (component, 1.0, xi, et, qq);
	double P1_val = P1 (component, 1.0, xi, et, qq);
	// todo: M2*td and N2*td are not evaluated correctly when detla = 90
	double M2_val = M2 (component, 1.0, xi, et, qq) * (fault_is_vertical ? 1. : td);
	double N2_val = N2 (component, 1.0, xi, et, qq) * (fault_is_vertical ? 1. : td);
	double M3_val = M3 (component, 1.0, xi, et, qq);
	double P1z_val = P1z (component, 1.0, xi, et, qq);
	double P1y_val = P1y (component, 1.0, xi, et, qq);

	val =	(2.0 - alpha4) * K7_val
		+ 2.0 * alpha3 * log_rc_val * sd - fault->alpha * log_re_val
		+ alpha3 * (qd * P1_val - (z - h) * (M2_val + N2_val * sd))
//		+ alpha3 * (qd * P1_val - (z - h) * (M2_val + N2_val * sd) * td)
		+ 2.0 * alpha5 * h * (P1_val * cd + (M2_val + N2_val * sd))
//		+ 2.0 * alpha5 * h * (P1_val * cd + (M2_val + N2_val * sd) * td)
		+ 6.0 * alpha3 * h * M3_val
		- 2.0 * alpha2 * h * ((qd + h * cd) * P1z_val - (z - 2.0 * h) * P1y_val * sd);

	return val;
}

static double
tensileyH0 (MagComp component, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K8_val = K8 (component, 1.0, xi, et, qq);
	double log_rx_val = log_rx (component, 1.0, xi, et, qq);
	double atan_xe_qr_val = atan_xe_qr (component, 1.0, xi, et, qq);
	double J1_val = J1 (component, 1.0, xi, et, qq);
	double O2_val = O2 (component, 1.0, xi, et, qq);
	double O3_val = O3 (component, 1.0, xi, et, qq);
	double M1_val = M1 (component, 1.0, xi, et, qq);
	double L1_val = L1 (component, 1.0, xi, et, qq);
	double O2y_val = O2y (component, 1.0, xi, et, qq);
	double O2z_val = O2z (component, 1.0, xi, et, qq);
	double M1y_val = M1y (component, 1.0, xi, et, qq);

	val = - (2.0 - alpha4) * K8_val
		- alpha4 * atan_xe_qr_val * sd - fault->alpha * log_rx_val * cd
		+ 2.0 * alpha3 * J1_val * sd
		+ alpha3 * (qd * (O3_val + M1_val * sd) - (z - h) * (O2_val - L1_val * sd) * sd)
		- 2.0 * fault->alpha * h * ((O3_val+ M1_val * sd) * cd + (O2_val - L1_val * sd) * sd)
		+ 6.0 * alpha3 * h * L1_val * sd2
		+ 2.0 * alpha2 * h * ((qd + h * cd) * (O2y_val + M1y_val *cd)
				+ (z - 2.0 * h) * (O2z_val + M1y_val * sd)
				+ M1_val * sd2 * cd);

	return val;
}

static double
tensilezH0 (MagComp component, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K9_val = K9 (component, 1.0, xi, et, qq);
	double log_rx_val = log_rx (component, 1.0, xi, et, qq);
	double atan_xe_qr_val = atan_xe_qr (component, 1.0, xi, et, qq);
	double P2_val = P2 (component, 1.0, xi, et, qq);
	double P3_val = P3 (component, 1.0, xi, et, qq);
	double M1_val = M1 (component, 1.0, xi, et, qq);
	double P3y_val = P3y (component, 1.0, xi, et, qq);
	double P3z_val = P3z (component, 1.0, xi, et, qq);

	val = - (2.0 + alpha5) * K9_val
		- fault->alpha * log_rx_val * sd + alpha4 * atan_xe_qr_val * cd
		+ alpha3 * (qd * P3_val - (z - h) * P2_val * sd)
		+ 2.0 * alpha4 * h * (P3_val * cd + P2_val * sd)
		+ 4.0 * alpha1 * h * M1_val * sd2
		+ 2.0 * alpha2 * h * ((qd + h * cd) * P3z_val
		// todo: check here (P3_val * sd -> P3_val * sd * cd is correct ?)
		// 		- (z - 2.0 * h) * P3y_val + 2.0 * P3_val * sd);
				- (z - 2.0 * h) * P3y_val * sd + 2.0 * P2_val * sd * cd);

	return val;
}

/*** contributions from the mirror image HI ***/
static double
tensilexHI (MagComp component, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K7_val = K7 (component, -1.0, xi, et, qq);
	double log_re_val = log_re (component, -1.0, xi, et, qq);
	double log_rc_val = log_rc (component, -1.0, xi, et, qq);
	double P1_val = P1 (component, -1.0, xi, et, qq);
	// todo: M2*td and N2*td are not evaluated correctly when detla = 90
	double M2_val = M2 (component, -1.0, xi, et, qq) * (fault_is_vertical ? 1. : td);
	double N2_val = N2 (component, -1.0, xi, et, qq) * (fault_is_vertical ? 1. : td);

	val =	- alpha4 * K7_val
		- 2.0 * alpha3 * log_rc_val * sd + fault->alpha * log_re_val
		+ alpha2 * (qd * P1_val + (z - h) * (M2_val - N2_val * sd));
	//	+ alpha2 * (qd * P1_val + (z - h) * (M2_val - N2_val * sd) * td);

	return val;
}

static double
tensileyHI (MagComp component, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K8_val = K8 (component, -1.0, xi, et, qq);
	double log_rx_val = log_rx (component, -1.0, xi, et, qq);
	double atan_xe_qr_val = atan_xe_qr (component, -1.0, xi, et, qq);
	double J1_val = J1 (component, -1.0, xi, et, qq);
	double O2_val = O2 (component, -1.0, xi, et, qq);
	double O3_val = O3 (component, -1.0, xi, et, qq);
	double M1_val = M1 (component, -1.0, xi, et, qq);
	double L1_val = L1 (component, -1.0, xi, et, qq);

	val = alpha4 * K8_val
		- alpha4 * atan_xe_qr_val * sd + fault->alpha * log_rx_val * cd
		+ 2.0 * alpha3 * J1_val * sd
		- alpha2 * (qd * (O3_val - M1_val * sd) - (z - h) * (O2_val + L1_val * sd) * sd);

	return val;
}

static double
tensilezHI (MagComp component, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K9_val = K9 (component, -1.0, xi, et, qq);
	double log_rx_val = log_rx (component, -1.0, xi, et, qq);
	double atan_xe_qr_val = atan_xe_qr (component, -1.0, xi, et, qq);
	double P2_val = P2 (component, -1.0, xi, et, qq);
	double P3_val = P3 (component, -1.0, xi, et, qq);

	// todo: check sign of K9_val
	val = (fault_is_vertical ? 1. : -1.) * alpha5 * K9_val
		+ fault->alpha * log_rx_val * sd + alpha4 * atan_xe_qr_val * cd
		- alpha2 * (qd * P3_val - (z - h) * P2_val * sd);

	return val;
}

/*** contributions from the mirror image HIII ***/
static double
tensilexHIII (MagComp component, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K7_val = K7 (component, 1.0, xi, et, qq);
	double log_re_val = log_re (component, 1.0, xi, et, qq);
	double log_rc_val = log_rc (component, 1.0, xi, et, qq);
	double P1_val = P1 (component, 1.0, xi, et, qq);
	// todo: M2*td and N2*td are not evaluated correctly when detla = 90
	double M2_val = M2 (component, 1.0, xi, et, qq) * (fault_is_vertical ? 1. : td);
	double N2_val = N2 (component, 1.0, xi, et, qq) * (fault_is_vertical ? 1. : td);

	val =	alpha4 * K7_val
		- 2.0 * alpha3 * log_rc_val * sd + fault->alpha * log_re_val
		- alpha3 * (qd * P1_val - (z - h) * (M2_val + N2_val * sd));
	//	- alpha3 * (qd * P1_val - (z - h) * (M2_val + N2_val * sd) * td);

	return val;
}

static double
tensileyHIII (MagComp component, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K8_val = K8 (component, 1.0, xi, et, qq);
	double log_rx_val = log_rx (component, 1.0, xi, et, qq);
	double atan_xe_qr_val = atan_xe_qr (component, 1.0, xi, et, qq);
	double J1_val = J1 (component, 1.0, xi, et, qq);
	double O2_val = O2 (component, 1.0, xi, et, qq);
	double O3_val = O3 (component, 1.0, xi, et, qq);
	double M1_val = M1 (component, 1.0, xi, et, qq);
	double L1_val = L1 (component, 1.0, xi, et, qq);

	val = - alpha4 * K8_val
		+ alpha4 * atan_xe_qr_val * sd + fault->alpha * log_rx_val * cd
		- 2.0 * alpha3 * J1_val * sd
		- alpha3 * (qd * (O3_val + M1_val * sd) - (z - h) * (O2_val - L1_val * sd) * sd);

	return val;
}

static double
tensilezHIII (MagComp component, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K9_val = K9 (component, 1.0, xi, et, qq);
	double log_rx_val = log_rx (component, 1.0, xi, et, qq);
	double atan_xe_qr_val = atan_xe_qr (component, 1.0, xi, et, qq);
	double P2_val = P2 (component, 1.0, xi, et, qq);
	double P3_val = P3 (component, 1.0, xi, et, qq);

	val = alpha5 * K9_val
		+ fault->alpha * log_rx_val * sd - alpha4 * atan_xe_qr_val * cd
		- alpha3 * (qd * P3_val - (z - h) * P2_val * sd);

	return val;
}

/*** seismomagnetic terms on fault coordinate system due to tensile-opening fault ***/

/* The following functions dose not calculate geometry variables such as r, rx, re, ..., etc,
 * only refer the global variables. So, before calling the following functions,
 * calc_geometry_variables() function must be called. */

/* main source */
double
tensile0 (MagComp component, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double	val = 0.;
	if (fabs (mag->cx) > DBL_EPSILON) val += mag->cx * tensilex0 (component, xi, et, qq);
	if (fabs (mag->cy) > DBL_EPSILON) val += mag->cy * tensiley0 (component, xi, et, qq);
	if (fabs (mag->cz) > DBL_EPSILON) val += mag->cz * tensilez0 (component, xi, et, qq);
	return val;
}

/* mirror image on fault coordinate system */
double
tensileH0 (MagComp component, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double	val = 0.;
	if (fabs (mag->cx) > DBL_EPSILON) val += mag->cx * tensilexH0 (component, fault, mag, xi, et, qq, y, z);
	if (fabs (mag->cy) > DBL_EPSILON) val += mag->cy * tensileyH0 (component, fault, mag, xi, et, qq, y, z);
	if (fabs (mag->cz) > DBL_EPSILON) val += mag->cz * tensilezH0 (component, fault, mag, xi, et, qq, y, z);
	return val;
}

/* sub-mirror image: type I */
double
tensileHI (MagComp component, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double	val = 0.;
	if (fabs (mag->cx) > DBL_EPSILON) val += mag->cx * tensilexHI (component, fault, mag, xi, et, qq, y, z);
	if (fabs (mag->cy) > DBL_EPSILON) val += mag->cy * tensileyHI (component, fault, mag, xi, et, qq, y, z);
	if (fabs (mag->cz) > DBL_EPSILON) val += mag->cz * tensilezHI (component, fault, mag, xi, et, qq, y, z);
	return val;
}

/* sub-mirror image: type III */
double
tensileHIII (MagComp component, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double	val = 0.;
	if (fabs (mag->cx) > DBL_EPSILON) val += mag->cx * tensilexHIII (component, fault, mag, xi, et, qq, y, z);
	if (fabs (mag->cy) > DBL_EPSILON) val += mag->cy * tensileyHIII (component, fault, mag, xi, et, qq, y, z);
	if (fabs (mag->cz) > DBL_EPSILON) val += mag->cz * tensilezHIII (component, fault, mag, xi, et, qq, y, z);
	return val;
}
