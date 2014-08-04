#include <stdio.h>
#include <math.h>

#include "constants.h"
#include "piez.h"

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

double
tensile0 (int flag, double x, double y, double z)
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

		xi = x + flength1;
		et = p + fwidth1;

		if (i >= 2)		 xi = x - flength2;
		if (i % 2 == 0) et = p - fwidth2;

		set_geometry_variables (sign, xi, et, q);
		hx = tensilex0 (flag, xi, et, q);
		hy = tensiley0 (flag, xi, et, q);
		hz = tensilez0 (flag, xi, et, q);

		res[i] = cx * hx + cy * hy + cz * hz;
	}
	return (res[0] + res[3]) - (res[1] + res[2]);
}

/*** contributions from the mirror image H0 ***/
double
tensilexH0 (int flag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = dcurier;
	double qd = y * sd + (fdepth - h) * cd;
	double K7_val = K7 (flag, 1.0, xi, et, qq);
	double log_re_val = log_re (flag, 1.0, xi, et, qq);
	double log_rc_val = log_rc (flag, 1.0, xi, et, qq);
	double P1_val = P1 (flag, 1.0, xi, et, qq);
	double M2_val = M2 (flag, 1.0, xi, et, qq);
	double N2_val = N2 (flag, 1.0, xi, et, qq);
	double M3_val = M3 (flag, 1.0, xi, et, qq);
	double P1z_val = P1z (flag, 1.0, xi, et, qq);
	double P1y_val = P1y (flag, 1.0, xi, et, qq);

	val =	(2.0 - alpha4) * K7_val
		+ 2.0 * alpha3 * log_rc_val * sd - alpha * log_re_val
		+ alpha3 * (qd * P1_val - (z - h) * (M2_val + N2_val * sd) * td)
		+ 2.0 * alpha5 * h * (P1_val * cd + (M2_val + N2_val * sd) * td)
		+ 6.0 * alpha3 * h * M3_val
		- 2.0 * alpha2 * h * ((qd + h * cd) * P1z_val - (z - 2.0 * h) * P1y_val * sd);

	return val;
}

double
tensileyH0 (int flag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = dcurier;
	double qd = y * sd + (fdepth - h) * cd;
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
		- alpha4 * atan_xe_qr_val * sd - alpha * log_rx_val * cd
		+ 2.0 * alpha3 * J1_val * sd
		+ alpha3 * (qd * (O3_val + M1_val * sd) - (z - h) * (O2_val - L1_val * td) * sd)
		- 2.0 * alpha * h * ((O3_val+ M1_val * sd) * cd + (O2_val - L1_val * td) * sd)
		+ 6.0 * alpha3 * h * L1_val * sd * td
		+ 2.0 * alpha2 * h * ((qd + h * cd) * (O2y_val + M1y_val *cd)
				+ (z - 2.0 * h) * (O2z_val + M1y_val * sd)
				+ M1_val * sd2 * cd);

	return val;
}

double
tensilezH0 (int flag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = dcurier;
	double qd = y * sd + (fdepth - h) * cd;
	double K9_val = K9 (flag, 1.0, xi, et, qq);
	double log_rx_val = log_rx (flag, 1.0, xi, et, qq);
	double atan_xe_qr_val = atan_xe_qr (flag, 1.0, xi, et, qq);
	double P2_val = P2 (flag, 1.0, xi, et, qq);
	double P3_val = P3 (flag, 1.0, xi, et, qq);
	double M1_val = M1 (flag, 1.0, xi, et, qq);
	double P3y_val = P3y (flag, 1.0, xi, et, qq);
	double P3z_val = P3z (flag, 1.0, xi, et, qq);

	val = - (2.0 + alpha5) * K9_val
		- alpha * log_rx_val * sd + alpha4 * atan_xe_qr_val * cd
		+ alpha3 * (qd * P3_val - (z - h) * P2_val * sd)
		+ 2.0 * alpha4 * h * (P3_val * cd + P2_val * sd)
		+ 4.0 * alpha1 * h * M1_val * sd2
		+ 2.0 * alpha2 * h * ((qd + h * cd) * P3z_val
				- (z - 2.0 * h) * P3y_val + 2.0 * P3_val * sd);

	return val;
}

double
tensileH0 (int flag, double x, double y, double z)
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

		xi = x + flength1;
		et = p + fwidth1;

		if (i >= 2)		 xi = x - flength2;
		if (i % 2 == 0) et = p - fwidth2;

		set_geometry_variables (sign, xi, et, q);
		hx = tensilexH0 (flag, xi, et, q, y, z);
		hy = tensileyH0 (flag, xi, et, q, y, z);
		hz = tensilezH0 (flag, xi, et, q, y, z);

		res[i] = cx * hx + cy * hy + cz * hz;
	}
	return res[0] - res[1] - res[2] + res[3];
}

/*** contributions from the mirror image HI ***/
double
tensilexHI (int flag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = dcurier;
	double qd = y * sd + (fdepth - h) * cd;
	double K7_val = K7 (flag, -1.0, xi, et, qq);
	double log_re_val = log_re (flag, -1.0, xi, et, qq);
	double log_rc_val = log_rc (flag, -1.0, xi, et, qq);
	double P1_val = P1 (flag, -1.0, xi, et, qq);
	double M2_val = M2 (flag, -1.0, xi, et, qq);
	double N2_val = N2 (flag, -1.0, xi, et, qq);

	val =	- alpha4 * K7_val
		- 2.0 * alpha3 * log_rc_val * sd + alpha * log_re_val
		+ alpha2 * (qd * P1_val + (z - h) * (M2_val - N2_val * sd) * td);

	return val;
}

double
tensileyHI (int flag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = dcurier;
	double qd = y * sd + (fdepth - h) * cd;
	double K8_val = K8 (flag, -1.0, xi, et, qq);
	double log_rx_val = log_rx (flag, -1.0, xi, et, qq);
	double atan_xe_qr_val = atan_xe_qr (flag, -1.0, xi, et, qq);
	double J1_val = J1 (flag, -1.0, xi, et, qq);
	double O2_val = O2 (flag, -1.0, xi, et, qq);
	double O3_val = O3 (flag, -1.0, xi, et, qq);
	double M1_val = M1 (flag, -1.0, xi, et, qq);
	double L1_val = L1 (flag, -1.0, xi, et, qq);

	val = alpha4 * K8_val
		- alpha4 * atan_xe_qr_val * sd + alpha * log_rx_val * cd
		+ 2.0 * alpha3 * J1_val * sd
		- alpha2 * (qd * (O3_val - M1_val * sd) - (z - h) * (O2_val + L1_val * td) * sd);

	return val;
}

double
tensilezHI (int flag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = dcurier;
	double qd = y * sd + (fdepth - h) * cd;
	double K9_val = K9 (flag, -1.0, xi, et, qq);
	double log_rx_val = log_rx (flag, -1.0, xi, et, qq);
	double atan_xe_qr_val = atan_xe_qr (flag, -1.0, xi, et, qq);
	double P2_val = P2 (flag, -1.0, xi, et, qq);
	double P3_val = P3 (flag, -1.0, xi, et, qq);

	val = - alpha5 * K9_val
		+ alpha * log_rx_val * sd + alpha4 * atan_xe_qr_val * cd
		- alpha2 * (qd * P3_val - (z - h) * P2_val * sd);

	return val;
}

double
tensileHI (int flag, double x, double y, double z)
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

		xi = x + flength1;
		et = p + fwidth1;

		if (i >= 2)		 xi = x - flength2;
		if (i % 2 == 0) et = p - fwidth2;

		set_geometry_variables (sign, xi, et, q);
		hx = tensilexHI (flag, xi, et, q, y, z);
		hy = tensileyHI (flag, xi, et, q, y, z);
		hz = tensilezHI (flag, xi, et, q, y, z);

		res[i] = cx * hx + cy * hy + cz * hz;
	}
	return (res[0] + res[3]) - (res[1] + res[2]);
}

/*** contributions from the mirror image HIII ***/
double
tensilexHIII (int flag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = dcurier;
	double qd = y * sd + (fdepth - h) * cd;
	double K7_val = K7 (flag, 1.0, xi, et, qq);
	double log_re_val = log_re (flag, 1.0, xi, et, qq);
	double log_rc_val = log_rc (flag, 1.0, xi, et, qq);
	double P1_val = P1 (flag, 1.0, xi, et, qq);
	double M2_val = M2 (flag, 1.0, xi, et, qq);
	double N2_val = N2 (flag, 1.0, xi, et, qq);

	val =	alpha4 * K7_val
		- 2.0 * alpha3 * log_rc_val * sd + alpha * log_re_val
		- alpha3 * (qd * P1_val - (z - h) * (M2_val + N2_val * sd) * td);

	return val;
}

double
tensileyHIII (int flag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = dcurier;
	double qd = y * sd + (fdepth - h) * cd;
	double K8_val = K8 (flag, 1.0, xi, et, qq);
	double log_rx_val = log_rx (flag, 1.0, xi, et, qq);
	double atan_xe_qr_val = atan_xe_qr (flag, 1.0, xi, et, qq);
	double J1_val = J1 (flag, 1.0, xi, et, qq);
	double O2_val = O2 (flag, 1.0, xi, et, qq);
	double O3_val = O3 (flag, 1.0, xi, et, qq);
	double M1_val = M1 (flag, 1.0, xi, et, qq);
	double L1_val = L1 (flag, 1.0, xi, et, qq);

	val = - alpha4 * K8_val
		+ alpha4 * atan_xe_qr_val * sd + alpha * log_rx_val * cd
		- 2.0 * alpha3 * J1_val * sd
		- alpha3 * (qd * (O3_val + M1_val * sd) - (z - h) * (O2_val - L1_val * td) * sd);

	return val;
}

double
tensilezHIII (int flag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = dcurier;
	double qd = y * sd + (fdepth - h) * cd;
	double K9_val = K9 (flag, 1.0, xi, et, qq);
	double log_rx_val = log_rx (flag, 1.0, xi, et, qq);
	double atan_xe_qr_val = atan_xe_qr (flag, 1.0, xi, et, qq);
	double P2_val = P2 (flag, 1.0, xi, et, qq);
	double P3_val = P3 (flag, 1.0, xi, et, qq);

	val = alpha5 * K9_val
		+ alpha * log_rx_val * sd - alpha4 * atan_xe_qr_val * cd
		- alpha3 * (qd * P3_val - (z - h) * P2_val * sd);

	return val;
}

double
tensileHIII (int flag, double x, double y, double z)
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

		xi = x + flength1;
		et = p + fwidth1;

		if (i >= 2)		 xi = x - flength2;
		if (i % 2 == 0) et = p - fwidth2;

		set_geometry_variables (sign, xi, et, q);
		hx = tensilexHIII (flag, xi, et, q, y, z);
		hy = tensileyHIII (flag, xi, et, q, y, z);
		hz = tensilezHIII (flag, xi, et, q, y, z);

		res[i] = cx * hx + cy * hy + cz * hz;
	} 
	return (res[0] + res[3]) - (res[1] + res[2]);
}


double
tensileHII (int flag, double x, double y, double z)
{
	int		i;
	double res[4];
	double p, q;
	double sign;
	double w = (dcurier - fdepth) / sd;
	double hx, hy, hz;

	p = y * cd - d[2] * sd;
	q = y * sd + d[2] * cd;

	for (i = 0; i < 4; i++) {
		double xi, et;

		xi = x + flength1;
		et = p + fwidth1;

		if (i >= 2)		 xi = x - flength2;
		if (i % 2 == 0) et = p - w;

		sign = - 1.0;
		set_geometry_variables (sign, xi, et, q);
		hx = tensilexHI (flag, xi, et, q, y, z);
		hy = tensileyHI (flag, xi, et, q, y, z);
		hz = tensilezHI (flag, xi, et, q, y, z);

		res[i] = cx * hx + cy * hy + cz * hz;
	}

	p = y * cd - d[1] * sd;
	q = y * sd + d[1] * cd;

	for (i = 0; i < 4; i++) {
		double xi, et;

		xi = x + flength1;
		et = p - w;

		if (i >= 2)		 xi = x - flength2;
		if (i % 2 == 0) et = p - fwidth2;

		sign = 1.0;
		set_geometry_variables (sign, xi, et, q);
		hx = tensilexHIII (flag, xi, et, q, y, z);
		hy = tensileyHIII (flag, xi, et, q, y, z);
		hz = tensilezHIII (flag, xi, et, q, y, z);

		res[i] += cx * hx + cy * hy + cz * hz;
	}
	return (res[0] + res[3]) - (res[1] + res[2]);
}


double
tensile_opening (int flag, double x, double y, double z)
{
	double res;

	res = tensile0 (flag, x, y, z);
	res += tensileH0 (flag, x, y, z);
	if (fdepth + fwidth2 * sd < dcurier) res += tensileHI (flag, x, y, z);
	else if (fdepth - fwidth1 * sd > dcurier) res += tensileHIII (flag, x, y, z);
	else res += tensileHII (flag, x, y, z);

	return res;
}
