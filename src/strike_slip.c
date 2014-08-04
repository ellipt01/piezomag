#include <stdio.h>
#include <math.h>

#include "constants.h"
#include "piez.h"

/*** main term ***/
double
strikex0 (int flag, double xi, double et, double qq)
{
	return 2.0 * K1 (flag, 1.0, xi, et, qq);
}

double
strikey0 (int flag, double xi, double et, double qq)
{
	return 2.0 * K2 (flag, 1.0, xi, et, qq);
}

double
strikez0 (int flag, double xi, double et, double qq)
{
	return 2.0 * K3 (flag, 1.0, xi, et, qq);
}

double
strike0 (int flag, double x, double y, double z)
{
	int		i;
	double res[4];
	double p, q;
	double sign = 1.0;
	double hx, hy, hz;

	set_singular_flag (1);
	p = y * cd - d[1] * sd;
	q = y * sd + d[1] * cd;

	for (i = 0; i < 4; i++) {
		double xi, et;

		xi = x + flength1;
		et = p + fwidth1;

		if (i >= 2)		 xi = x - flength2;
		if (i % 2 == 0) et = p - fwidth2;

		set_geometry_variables (sign, xi, et, q);
		hx = strikex0 (flag, xi, et, q);
		hy = strikey0 (flag, xi, et, q);
		hz = strikez0 (flag, xi, et, q);

		res[i] = cx * hx + cy * hy + cz * hz;
	}
	clear_singular_flag (0);
	return (res[0] + res[3]) - (res[1] + res[2]);
}

/*** contributions from the mirror image H0 ***/
double
strikexH0 (int flag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = dcurier;
	double qd = y * sd + (fdepth - h) * cd;
	double K1_val = K1 (flag, 1.0, xi, et, qq);
	double atan_xe_qr_val = atan_xe_qr (flag, 1.0, xi, et, qq);
	double J1_val = J1 (flag, 1.0, xi, et, qq);
	double M1_val = M1 (flag, 1.0, xi, et, qq);
	double L1_val = L1 (flag, 1.0, xi, et, qq);
	double M1y_val = M1y (flag, 1.0, xi, et, qq);
	double M1z_val = M1z (flag, 1.0, xi, et, qq);

	val = - (2.0 - alpha4) * K1_val
		- alpha1 * atan_xe_qr_val - alpha3 * J1_val
		- alpha3 * (qd * M1_val + (z - h) * L1_val * td)
		- 2.0 * alpha4 * h * (M1_val * cd - L1_val * td)
		- 4.0 * alpha1 * h * L1_val * td
		+ 2.0 * alpha2 * h * ((qd + h * cd) * M1z_val - (z - 2.0 * h) * M1y_val * sd);

	return val;
}

double
strikeyH0 (int flag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = dcurier;
	double qd = y * sd + (fdepth - h) * cd;
	double K2_val = K2 (flag, 1.0, xi, et, qq);
	double log_re_val = log_re (flag, 1.0, xi, et, qq);
	double J2_val = J2 (flag, 1.0, xi, et, qq);
	double L2_val = L2 (flag, 1.0, xi, et, qq);
	double M2_val = M2 (flag, 1.0, xi, et, qq);
	double M3_val = M3 (flag, 1.0, xi, et, qq);
	double M2y_val = M2y (flag, 1.0, xi, et, qq);
	double M2z_val = M2z (flag, 1.0, xi, et, qq);

	val = - (2.0 - alpha4) * K2_val
		+ alpha4 * log_re_val * sd + alpha3 * J2_val
		- alpha3 * (qd * M2_val + (z - h) * L2_val * sd)
		+ 2.0 * alpha4 * h * (M2_val * cd - L2_val * sd)
		- 4.0 * alpha1 * h * L2_val * sd
		+ 2.0 * alpha2 * h * M3_val * sd
		+ 2.0 * alpha2 * h * ((qd + h * cd) * M2z_val - (z - 2.0 * h) * M2y_val * sd);

	return val;
}

double
strikezH0 (int flag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = dcurier;
	double qd = y * sd + (fdepth - h) * cd;
	double K3_val = K3 (flag, 1.0, xi, et, qq);
	double log_re_val = log_re (flag, 1.0, xi, et, qq);
	double M2_val = M2 (flag, 1.0, xi, et, qq);
	double M3_val = M3 (flag, 1.0, xi, et, qq);
	double M3y_val = M2z (flag, 1.0, xi, et, qq);
	double M3z_val = M2y (flag, 1.0, xi, et, qq);

	val = - (2.0 + alpha5) * K3_val
		- alpha * log_re_val * cd
		- alpha3 * (qd * M3_val - (z - h) * M2_val * sd)
		+ 2.0 * alpha * h * (M3_val * cd + M2_val * sd)
		- 4.0 * alpha1 * h * M2_val * sd
		- 2.0 * alpha2 * h * M2_val * sd
		- 2.0 * alpha2 * h * ((qd + h * cd) * M3z_val - (z - 2.0 * h) * M3y_val * sd);

	return val;
}

double
strikeH0 (int flag, double x, double y, double z)
{
	int		i;
	double res[4];
	double p, q;
	double sign = 1.0;
	double hx, hy, hz;

	set_singular_flag (3);
	p = y * cd - d[3] * sd;
	q = y * sd + d[3] * cd;

	for (i = 0; i < 4; i++) {
		double xi, et;

		xi = x + flength1;
		et = p + fwidth1;

		if (i >= 2)		 xi = x - flength2;
		if (i % 2 == 0) et = p - fwidth2;

		set_geometry_variables (sign, xi, et, q);
		hx = strikexH0 (flag, xi, et, q, y, z);
		hy = strikeyH0 (flag, xi, et, q, y, z);
		hz = strikezH0 (flag, xi, et, q, y, z);

		res[i] = cx * hx + cy * hy + cz * hz;
	}
	clear_singular_flag (0);
	return (res[0] + res[3]) - (res[1] + res[2]);
}

/*** contributions from the mirror image HI ***/
double
strikexHI (int flag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = dcurier;
	double qd = y * sd + (fdepth - h) * cd;
	double K1_val = K1 (flag, -1.0, xi, et, qq);
	double atan_xe_qr_val = atan_xe_qr (flag, -1.0, xi, et, qq);
	double J1_val = J1 (flag, -1.0, xi, et, qq);
	double M1_val = M1 (flag, -1.0, xi, et, qq);
	double L1_val = L1 (flag, -1.0, xi, et, qq);

	val = alpha4 * K1_val
		- alpha1 * atan_xe_qr_val - alpha3 * J1_val
		- alpha2 * (qd * M1_val + (z - h) * L1_val * td);

	return val;
}

double
strikeyHI (int flag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = dcurier;
	double qd = y * sd + (fdepth - h) * cd;
	double K2_val = K2 (flag, -1.0, xi, et, qq);
	double log_re_val = log_re (flag, -1.0, xi, et, qq);
	double J2_val = J2 (flag, -1.0, xi, et, qq);
	double L2_val = L2 (flag, -1.0, xi, et, qq);
	double M2_val = M2 (flag, -1.0, xi, et, qq);

	val = alpha4 * K2_val + alpha6 * log_re_val * sd - alpha3 * J2_val
		- alpha2 * (qd * M2_val + (z - h) * L2_val * sd);

	return val;
}

double
strikezHI (int flag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = dcurier;
	double qd = y * sd + (fdepth - h) * cd;
	double K3_val = K3 (flag, -1.0, xi, et, qq);
	double log_re_val = log_re (flag, -1.0, xi, et, qq);
	double M2_val = M2 (flag, -1.0, xi, et, qq);
	double M3_val = M3 (flag, -1.0, xi, et, qq);

	val = alpha5 * K3_val + alpha * log_re_val * cd
		+ alpha2 * (qd * M3_val - (z - h) * M2_val * sd);

	return val;
}

double
strikeHI (int flag, double x, double y, double z)
{
	int		i;
	double res[4];
	double p, q;
	double sign = -1.0;
	double hx, hy, hz;

	set_singular_flag (2);
	p = y * cd - d[2] * sd;
	q = y * sd + d[2] * cd;

	for (i = 0; i < 4; i++) {
		double xi, et;

		xi = x + flength1;
		et = p + fwidth1;

		if (i >= 2)		 xi = x - flength2;
		if (i % 2 == 0) et = p - fwidth2;

		set_geometry_variables (sign, xi, et, q);
		hx = strikexHI (flag, xi, et, q, y, z);
		hy = strikeyHI (flag, xi, et, q, y, z);
		hz = strikezHI (flag, xi, et, q, y, z);

		res[i] = cx * hx + cy * hy + cz * hz;
	}
	clear_singular_flag (0);
	return (res[0] + res[3]) - (res[1] + res[2]);
}

/*** contributions from the mirror image HIII ***/
double
strikexHIII (int flag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = dcurier;
	double qd = y * sd + (fdepth - h) * cd;
	double K1_val = K1 (flag, 1.0, xi, et, qq);
	double atan_xe_qr_val = atan_xe_qr (flag, 1.0, xi, et, qq);
	double J1_val = J1 (flag, 1.0, xi, et, qq);
	double M1_val = M1 (flag, 1.0, xi, et, qq);
	double L1_val = L1 (flag, 1.0, xi, et, qq);

	val = - alpha4 * K1_val
		+ alpha1 * atan_xe_qr_val + alpha3 * J1_val
		+ alpha3 * (qd * M1_val + (z - h) * L1_val * td);

	return val;
}

double
strikeyHIII (int flag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = dcurier;
	double qd = y * sd + (fdepth - h) * cd;
	double K2_val = K2 (flag, 1.0, xi, et, qq);
	double log_re_val = log_re (flag, 1.0, xi, et, qq);
	double J2_val = J2 (flag, 1.0, xi, et, qq);
	double L2_val = L2 (flag, 1.0, xi, et, qq);
	double M2_val = M2 (flag, 1.0, xi, et, qq);

	val = - alpha4 * K2_val - alpha4 * log_re_val * sd - alpha3 * J2_val
		+ alpha3 * (qd * M2_val + (z - h) * L2_val * sd);

	return val;
}

double
strikezHIII (int flag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = dcurier;
	double qd = y * sd + (fdepth - h) * cd;
	double K3_val = K3 (flag, 1.0, xi, et, qq);
	double log_re_val = log_re (flag, 1.0, xi, et, qq);
	double M2_val = M2 (flag, 1.0, xi, et, qq);
	double M3_val = M3 (flag, 1.0, xi, et, qq);

	val = alpha5 * K3_val + alpha * log_re_val * cd
		+ alpha3 * (qd * M3_val - (z - h) * M2_val * sd);

	return val;
}

double
strikeHIII (int flag, double x, double y, double z)
{
	int		i;
	double res[4];
	double p, q;
	double sign = 1.0;
	double hx, hy, hz;

	set_singular_flag (1);
	p = y * cd - d[1] * sd;
	q = y * sd + d[1] * cd;

	for (i = 0; i < 4; i++) {
		double xi, et;

		xi = x + flength1;
		et = p + fwidth1;

		if (i >= 2)		 xi = x - flength2;
		if (i % 2 == 0) et = p - fwidth2;

		set_geometry_variables (sign, xi, et, q);
		hx = strikexHIII (flag, xi, et, q, y, z);
		hy = strikeyHIII (flag, xi, et, q, y, z);
		hz = strikezHIII (flag, xi, et, q, y, z);

		res[i] = cx * hx + cy * hy + cz * hz;
	}
	clear_singular_flag (0);
	return (res[0] + res[3]) - (res[1] + res[2]);
}


double
strikeHII (int flag, double x, double y, double z)
{
	int		i;
	double res[4];
	double p, q;
	double sign;
	double w = (dcurier - fdepth) / sd;
	double hx, hy, hz;

	set_singular_flag (2);
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
		hx = strikexHI (flag, xi, et, q, y, z);
		hy = strikeyHI (flag, xi, et, q, y, z);
		hz = strikezHI (flag, xi, et, q, y, z);

		res[i] = cx * hx + cy * hy + cz * hz;
	}

	set_singular_flag (1);
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
		hx = strikexHIII (flag, xi, et, q, y, z);
		hy = strikeyHIII (flag, xi, et, q, y, z);
		hz = strikezHIII (flag, xi, et, q, y, z);

		res[i] += cx * hx + cy * hy + cz * hz;
	}
	set_singular_flag (0);
	return (res[0] + res[3]) - (res[1] + res[2]);
}


double
strike_slip (int flag, double x, double y, double z)
{
	double res;

	res = strike0 (flag, x, y, z);
	res += strikeH0 (flag, x, y, z);
	if (fdepth + fwidth2 * sd < dcurier) res += strikeHI (flag, x, y, z);
	else if (fdepth - fwidth1 * sd > dcurier) res += strikeHIII (flag, x, y, z);
	else res += strikeHII (flag, x, y, z);

	return res;
}
