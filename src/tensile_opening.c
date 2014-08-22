#include <stdio.h>
#include <math.h>

#include "piezomag.h"
#include "private.h"

/*** main term ***/
double
tensilex0 (int flag, double xi, double et, double qq)
{
	return -2.0 * K7 (flag, 1.0, xi, et, qq);
}

double
tensiley0 (int flag, double xi, double et, double qq)
{
	return 2.0 * K8 (flag, 1.0, xi, et, qq);
}

double
tensilez0 (int flag, double xi, double et, double qq)
{
	return 2.0 * K9 (flag, 1.0, xi, et, qq);
}

static double
tensile0 (int flag, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
{
	int		i;
	double res[4];
	double p, q;
	double sign = 1.0;
	double hx, hy, hz;

	p = y * cd - d[1] * sd;
	q = y * sd + d[1] * cd;

	for (i = 0; i < 4; i++) {
		double xi, et;

		xi = x + fault->flength1;
		et = p + fault->fwidth1;

		if (i >= 2)		 xi = x - fault->flength2;
		if (i % 2 == 0) et = p - fault->fwidth2;

		set_geometry_variables (sign, xi, et, q);
		hx = tensilex0 (flag, xi, et, q);
		hy = tensiley0 (flag, xi, et, q);
		hz = tensilez0 (flag, xi, et, q);

		res[i] = mag->cx * hx + mag->cy * hy + mag->cz * hz;
	}
	return (res[0] + res[3]) - (res[1] + res[2]);
}

/*** contributions from the mirror image H0 ***/
static double
tensilexH0 (int flag, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K7_val = K7 (flag, 1.0, xi, et, qq);
	double log_re_val = log_re (flag, 1.0, xi, et, qq);
	double log_rc_val = log_rc (flag, 1.0, xi, et, qq);
	double P1_val = P1 (flag, 1.0, xi, et, qq);
	double M2_val = M2 (flag, 1.0, xi, et, qq) * (fault_is_vertical ? 1. : td);
	double N2_val = N2 (flag, 1.0, xi, et, qq) * (fault_is_vertical ? 1. : td);
	double M3_val = M3 (flag, 1.0, xi, et, qq);
	double P1z_val = P1z (flag, 1.0, xi, et, qq);
	double P1y_val = P1y (flag, 1.0, xi, et, qq);

	val =	(2.0 - alpha4) * K7_val
		+ 2.0 * alpha3 * log_rc_val * sd - fault->alpha * log_re_val
		+ alpha3 * (qd * P1_val - (z - h) * (M2_val + N2_val * sd))
		+ 2.0 * alpha5 * h * (P1_val * cd + (M2_val + N2_val * sd))
		//		+ alpha3 * (qd * P1_val - (z - h) * (M2_val + N2_val * sd) * td)
		//		+ 2.0 * alpha5 * h * (P1_val * cd + (M2_val + N2_val * sd) * td)
		+ 6.0 * alpha3 * h * M3_val
		- 2.0 * alpha2 * h * ((qd + h * cd) * P1z_val - (z - 2.0 * h) * P1y_val * sd);

	return val;
}

static double
tensileyH0 (int flag, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K8_val = K8 (flag, 1.0, xi, et, qq);
	double log_rx_val = log_rx (flag, 1.0, xi, et, qq);
	double atan_xe_qr_val = atan_xe_qr (flag, 1.0, xi, et, qq);
	double J1_val = J1 (flag, 1.0, xi, et, qq);
	double O2_val = O2 (flag, 1.0, xi, et, qq);
	double O3_val = O3 (flag, 1.0, xi, et, qq);
	double M1_val = M1 (flag, 1.0, xi, et, qq);
	double L1_val = L1 (flag, 1.0, xi, et, qq);
	double O2y_val = O2y (flag, 1.0, xi, et, qq);
	double O2z_val = O2z (flag, 1.0, xi, et, qq);
	double M1y_val = M1y (flag, 1.0, xi, et, qq);

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
tensilezH0 (int flag, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K9_val = K9 (flag, 1.0, xi, et, qq);
	double log_rx_val = log_rx (flag, 1.0, xi, et, qq);
	double atan_xe_qr_val = atan_xe_qr (flag, 1.0, xi, et, qq);
	double P2_val = P2 (flag, 1.0, xi, et, qq);
	double P3_val = P3 (flag, 1.0, xi, et, qq);
	double M1_val = M1 (flag, 1.0, xi, et, qq);
	double P3y_val = P3y (flag, 1.0, xi, et, qq);
	double P3z_val = P3z (flag, 1.0, xi, et, qq);

	val = - (2.0 + alpha5) * K9_val
		- fault->alpha * log_rx_val * sd + alpha4 * atan_xe_qr_val * cd
		+ alpha3 * (qd * P3_val - (z - h) * P2_val * sd)
		+ 2.0 * alpha4 * h * (P3_val * cd + P2_val * sd)
		+ 4.0 * alpha1 * h * M1_val * sd2
		+ 2.0 * alpha2 * h * ((qd + h * cd) * P3z_val
				- (z - 2.0 * h) * P3y_val + 2.0 * P3_val * sd);

	return val;
}

static double
tensileH0 (int flag, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
{
	int		i;
	double res[4];
	double p, q;
	double sign = 1.0;
	double hx, hy, hz;

	p = y * cd - d[3] * sd;
	q = y * sd + d[3] * cd;

	for (i = 0; i < 4; i++) {
		double xi, et;

		xi = x + fault->flength1;
		et = p + fault->fwidth1;

		if (i >= 2)		 xi = x - fault->flength2;
		if (i % 2 == 0) et = p - fault->fwidth2;

		set_geometry_variables (sign, xi, et, q);
		hx = tensilexH0 (flag, fault, mag, xi, et, q, y, z);
		hy = tensileyH0 (flag, fault, mag, xi, et, q, y, z);
		hz = tensilezH0 (flag, fault, mag, xi, et, q, y, z);

		res[i] = mag->cx * hx + mag->cy * hy + mag->cz * hz;
	}
	return res[0] - res[1] - res[2] + res[3];
}

/*** contributions from the mirror image HI ***/
static double
tensilexHI (int flag, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K7_val = K7 (flag, -1.0, xi, et, qq);
	double log_re_val = log_re (flag, -1.0, xi, et, qq);
	double log_rc_val = log_rc (flag, -1.0, xi, et, qq);
	double P1_val = P1 (flag, -1.0, xi, et, qq);
	double M2_val = M2 (flag, -1.0, xi, et, qq) * (fault_is_vertical ? 1. : td);
	double N2_val = N2 (flag, -1.0, xi, et, qq) * (fault_is_vertical ? 1. : td);

	val =	- alpha4 * K7_val
		- 2.0 * alpha3 * log_rc_val * sd + fault->alpha * log_re_val
		+ alpha2 * (qd * P1_val + (z - h) * (M2_val - N2_val * sd));
	//	+ alpha2 * (qd * P1_val + (z - h) * (M2_val - N2_val * sd) * td);

	return val;
}

static double
tensileyHI (int flag, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K8_val = K8 (flag, -1.0, xi, et, qq);
	double log_rx_val = log_rx (flag, -1.0, xi, et, qq);
	double atan_xe_qr_val = atan_xe_qr (flag, -1.0, xi, et, qq);
	double J1_val = J1 (flag, -1.0, xi, et, qq);
	double O2_val = O2 (flag, -1.0, xi, et, qq);
	double O3_val = O3 (flag, -1.0, xi, et, qq);
	double M1_val = M1 (flag, -1.0, xi, et, qq);
	double L1_val = L1 (flag, -1.0, xi, et, qq);

	val = alpha4 * K8_val
		- alpha4 * atan_xe_qr_val * sd + fault->alpha * log_rx_val * cd
		+ 2.0 * alpha3 * J1_val * sd
		- alpha2 * (qd * (O3_val - M1_val * sd) - (z - h) * (O2_val + L1_val * sd) * sd);

	return val;
}

static double
tensilezHI (int flag, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K9_val = K9 (flag, -1.0, xi, et, qq);
	double log_rx_val = log_rx (flag, -1.0, xi, et, qq);
	double atan_xe_qr_val = atan_xe_qr (flag, -1.0, xi, et, qq);
	double P2_val = P2 (flag, -1.0, xi, et, qq);
	double P3_val = P3 (flag, -1.0, xi, et, qq);

	val = (fault_is_vertical ? 1. : -1.) * alpha5 * K9_val
		+ fault->alpha * log_rx_val * sd + alpha4 * atan_xe_qr_val * cd
		- alpha2 * (qd * P3_val - (z - h) * P2_val * sd);

	return val;
}

static double
tensileHI (int flag, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
{
	int		i;
	double res[4];
	double p, q;
	double sign = -1.0;
	double hx, hy, hz;

	p = y * cd - d[2] * sd;
	q = y * sd + d[2] * cd;

	for (i = 0; i < 4; i++) {
		double xi, et;

		xi = x + fault->flength1;
		et = p + fault->fwidth1;

		if (i >= 2)		 xi = x - fault->flength2;
		if (i % 2 == 0) et = p - fault->fwidth2;

		set_geometry_variables (sign, xi, et, q);
		hx = tensilexHI (flag, fault, mag, xi, et, q, y, z);
		hy = tensileyHI (flag, fault, mag, xi, et, q, y, z);
		hz = tensilezHI (flag, fault, mag, xi, et, q, y, z);

		res[i] = mag->cx * hx + mag->cy * hy + mag->cz * hz;
	}
	return (res[0] + res[3]) - (res[1] + res[2]);
}

/*** contributions from the mirror image HIII ***/
static double
tensilexHIII (int flag, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K7_val = K7 (flag, 1.0, xi, et, qq);
	double log_re_val = log_re (flag, 1.0, xi, et, qq);
	double log_rc_val = log_rc (flag, 1.0, xi, et, qq);
	double P1_val = P1 (flag, 1.0, xi, et, qq);
	double M2_val = M2 (flag, 1.0, xi, et, qq) * (fault_is_vertical ? 1. : td);
	double N2_val = N2 (flag, 1.0, xi, et, qq) * (fault_is_vertical ? 1. : td);

	val =	alpha4 * K7_val
		- 2.0 * alpha3 * log_rc_val * sd + fault->alpha * log_re_val
		- alpha3 * (qd * P1_val - (z - h) * (M2_val + N2_val * sd));
	//	- alpha3 * (qd * P1_val - (z - h) * (M2_val + N2_val * sd) * td);

	return val;
}

static double
tensileyHIII (int flag, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K8_val = K8 (flag, 1.0, xi, et, qq);
	double log_rx_val = log_rx (flag, 1.0, xi, et, qq);
	double atan_xe_qr_val = atan_xe_qr (flag, 1.0, xi, et, qq);
	double J1_val = J1 (flag, 1.0, xi, et, qq);
	double O2_val = O2 (flag, 1.0, xi, et, qq);
	double O3_val = O3 (flag, 1.0, xi, et, qq);
	double M1_val = M1 (flag, 1.0, xi, et, qq);
	double L1_val = L1 (flag, 1.0, xi, et, qq);

	val = - alpha4 * K8_val
		+ alpha4 * atan_xe_qr_val * sd + fault->alpha * log_rx_val * cd
		- 2.0 * alpha3 * J1_val * sd
		- alpha3 * (qd * (O3_val + M1_val * sd) - (z - h) * (O2_val - L1_val * sd) * sd);

	return val;
}

static double
tensilezHIII (int flag, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K9_val = K9 (flag, 1.0, xi, et, qq);
	double log_rx_val = log_rx (flag, 1.0, xi, et, qq);
	double atan_xe_qr_val = atan_xe_qr (flag, 1.0, xi, et, qq);
	double P2_val = P2 (flag, 1.0, xi, et, qq);
	double P3_val = P3 (flag, 1.0, xi, et, qq);

	val = alpha5 * K9_val
		+ fault->alpha * log_rx_val * sd - alpha4 * atan_xe_qr_val * cd
		- alpha3 * (qd * P3_val - (z - h) * P2_val * sd);

	return val;
}

static double
tensileHIII (int flag, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
{
	int		i;
	double res[4];
	double p, q;
	double sign = 1.0;
	double hx, hy, hz;

	p = y * cd - d[1] * sd;
	q = y * sd + d[1] * cd;

	for (i = 0; i < 4; i++) {
		double xi, et;

		xi = x + fault->flength1;
		et = p + fault->fwidth1;

		if (i >= 2)		 xi = x - fault->flength2;
		if (i % 2 == 0) et = p - fault->fwidth2;

		set_geometry_variables (sign, xi, et, q);
		hx = tensilexHIII (flag, fault, mag, xi, et, q, y, z);
		hy = tensileyHIII (flag, fault, mag, xi, et, q, y, z);
		hz = tensilezHIII (flag, fault, mag, xi, et, q, y, z);

		res[i] = mag->cx * hx + mag->cy * hy + mag->cz * hz;
	} 
	return (res[0] + res[3]) - (res[1] + res[2]);
}


static double
tensileHII (int flag, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
{
	int		i;
	double res[4];
	double p, q;
	double sign;
	double w = (mag->dcurier - fault->fdepth) / sd;
	double hx, hy, hz;

	p = y * cd - d[2] * sd;
	q = y * sd + d[2] * cd;

	for (i = 0; i < 4; i++) {
		double xi, et;

		xi = x + fault->flength1;
		et = p + fault->fwidth1;

		if (i >= 2)		 xi = x - fault->flength2;
		if (i % 2 == 0) et = p - w;

		sign = - 1.0;
		set_geometry_variables (sign, xi, et, q);
		hx = tensilexHI (flag, fault, mag, xi, et, q, y, z);
		hy = tensileyHI (flag, fault, mag, xi, et, q, y, z);
		hz = tensilezHI (flag, fault, mag, xi, et, q, y, z);

		res[i] = mag->cx * hx + mag->cy * hy + mag->cz * hz;
	}

	p = y * cd - d[1] * sd;
	q = y * sd + d[1] * cd;

	for (i = 0; i < 4; i++) {
		double xi, et;

		xi = x + fault->flength1;
		et = p - w;

		if (i >= 2)		 xi = x - fault->flength2;
		if (i % 2 == 0) et = p - fault->fwidth2;

		sign = 1.0;
		set_geometry_variables (sign, xi, et, q);
		hx = tensilexHIII (flag, fault, mag, xi, et, q, y, z);
		hy = tensileyHIII (flag, fault, mag, xi, et, q, y, z);
		hz = tensilezHIII (flag, fault, mag, xi, et, q, y, z);

		res[i] += mag->cx * hx + mag->cy * hy + mag->cz * hz;
	}
	return (res[0] + res[3]) - (res[1] + res[2]);
}


double
tensile_opening (int flag, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
{
	double res;

	res = tensile0 (flag, fault, mag, x, y, z);
	res += tensileH0 (flag, fault, mag, x, y, z);
	if (fault->fdepth + fault->fwidth2 * sd < mag->dcurier) res += tensileHI (flag, fault, mag, x, y, z);
	else if (fault->fdepth - fault->fwidth1 * sd > mag->dcurier) res += tensileHIII (flag, fault, mag, x, y, z);
	else res += tensileHII (flag, fault, mag, x, y, z);

	return res;
}
