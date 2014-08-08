#include <stdio.h>
#include <math.h>

#include "piezomag.h"
#include "private.h"

/*** main term ***/
static double
dipx0 (int flag, double xi, double et, double qq)
{
	return -2.0 * K4 (flag, 1.0, xi, et, qq);
}

static double
dipy0 (int flag, double xi, double et, double qq)
{
	return 2.0 * K5 (flag, 1.0, xi, et, qq);
}

static double
dipz0 (int flag, double xi, double et, double qq)
{
	return 2.0 * K6 (flag, 1.0, xi, et, qq);
}

static double
dip0 (int flag, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
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
		hx = dipx0 (flag, xi, et, q);
		hy = dipy0 (flag, xi, et, q);
		hz = dipz0 (flag, xi, et, q);

		res[i] = mag->cx * hx + mag->cy * hy + mag->cz * hz;
	}
	return (res[0] + res[3]) - (res[1] + res[2]);
}

/*** contributions from the mirror image H0 ***/
static double
dipxH0 (int flag, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K4_val = K4 (flag, 1.0, xi, et, qq);
	double log_re_val = log_re (flag, 1.0, xi, et, qq);
	double J2_val = J2 (flag, 1.0, xi, et, qq);
	double O1_val = O1 (flag, 1.0, xi, et, qq);
	double N2_val = N2 (flag, 1.0, xi, et, qq);
	double M2_val = M2 (flag, 1.0, xi, et, qq);
	double O2z_val = O1z (flag, 1.0, xi, et, qq);
	double N2z_val = N2z (flag, 1.0, xi, et, qq);

	val =	(2.0 - alpha4) * K4_val
		+ alpha3 * (log_re_val * s2d - J2_val * c2d * secd)
		+ alpha3 * (qd * O1_val + (z - h) * N2_val * sd)
		+ 2.0 * alpha5 * h * (O1_val * cd - N2_val * sd)
		+ 4.0 * alpha1 * h * M2_val
		- 2.0 * alpha2 * h * ((qd + h * cd) * O2z_val + (z - 2.0 * h) * N2z_val * sd);

	return val;
}

static double
dipyH0 (int flag, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K5_val = K5 (flag, 1.0, xi, et, qq);
	double log_rx_val = log_rx (flag, 1.0, xi, et, qq);
	double atan_xe_qr_val = atan_xe_qr (flag, 1.0, xi, et, qq);
	double J1_val = J1 (flag, 1.0, xi, et, qq);
	double O2_val = O2 (flag, 1.0, xi, et, qq);
	double O3_val = O3 (flag, 1.0, xi, et, qq);
	double N1_val = N1 (flag, 1.0, xi, et, qq);
	double L1_val = L1 (flag, 1.0, xi, et, qq);
	double O2z_val = O2z (flag, 1.0, xi, et, qq);
	double O3z_val = O3z (flag, 1.0, xi, et, qq);
	double N1z_val = N1z (flag, 1.0, xi, et, qq);

	val = - (2.0 - alpha4) * K5_val
		- fault->alpha * log_rx_val * sd + alpha4 * atan_xe_qr_val * cd
		- alpha3 * J1_val * c2d * secd
		+ alpha3 * (qd * O2_val - (z - h) * (N1_val - O3_val) * sd)
		+ 2.0 * alpha4 * h * (O2_val * cd + (N1_val - O3_val) * sd)
		- 4.0 * alpha1 * h * L1_val * sd
		- 2.0 * alpha2 * h * ((qd + h * cd) * O2z_val
				- (z - 2.0 * h) * (N1z_val - O3z_val) * sd
				+ O3_val * sd);

	return val;
}

static double
dipzH0 (int flag, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K6_val = K6 (flag, 1.0, xi, et, qq);
	double log_rx_val = log_rx (flag, 1.0, xi, et, qq);
	double atan_xe_qr_val = atan_xe_qr (flag, 1.0, xi, et, qq);
	double O2_val = O2 (flag, 1.0, xi, et, qq);
	double O3_val = O3 (flag, 1.0, xi, et, qq);
	double M1_val = M1 (flag, 1.0, xi, et, qq);
	double O2z_val = O2z (flag, 1.0, xi, et, qq);
	double O3z_val = O3z (flag, 1.0, xi, et, qq);

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
dipH0 (int flag, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
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
		hx = dipxH0 (flag, fault, mag, xi, et, q, y, z);
		hy = dipyH0 (flag, fault, mag, xi, et, q, y, z);
		hz = dipzH0 (flag, fault, mag, xi, et, q, y, z);

		res[i] = mag->cx * hx + mag->cy * hy + mag->cz * hz;
	}
	return (res[0] + res[3]) - (res[1] + res[2]);
}

/*** contributions from the mirror image HI ***/
static double
dipxHI (int flag, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K4_val = K4 (flag, -1.0, xi, et, qq);
	double log_re_val = log_re (flag, -1.0, xi, et, qq);
	double J2_val = J2 (flag, -1.0, xi, et, qq);
	double O1_val = O1 (flag, -1.0, xi, et, qq);
	double N2_val = N2 (flag, -1.0, xi, et, qq);

	val = - alpha4 * K4_val + alpha3 * (log_re_val * s2d + J2_val * c2d * secd)
		- alpha2 * (qd * O1_val + (z - h) * N2_val * sd);

	return val;
}

static double
dipyHI (int flag, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K5_val = K5 (flag, -1.0, xi, et, qq);
	double log_rx_val = log_rx (flag, -1.0, xi, et, qq);
	double atan_xe_qr_val = atan_xe_qr (flag, -1.0, xi, et, qq);
	double J1_val = J1 (flag, -1.0, xi, et, qq);
	double O2_val = O2 (flag, -1.0, xi, et, qq);
	double O3_val = O3 (flag, -1.0, xi, et, qq);
	double N1_val = N1 (flag, -1.0, xi, et, qq);

	val = alpha4 * K5_val + fault->alpha * log_rx_val * sd
		+ alpha4 * atan_xe_qr_val * cd
		- alpha3 * J1_val * c2d * secd
		+ alpha2 * (qd * O2_val - (z - h) * (N1_val - O3_val) * sd);

	return val;
}

static double
dipzHI (int flag, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K6_val = K6 (flag, -1.0, xi, et, qq);
	double log_rx_val = log_rx (flag, -1.0, xi, et, qq);
	double atan_xe_qr_val = atan_xe_qr (flag, -1.0, xi, et, qq);
	double O2_val = O2 (flag, -1.0, xi, et, qq);
	double O3_val = O3 (flag, -1.0, xi, et, qq);

	val = - alpha5 * K6_val
		- fault->alpha * log_rx_val * cd + alpha4 * atan_xe_qr_val * sd
		- alpha2 * (qd * O3_val - (z - h) * O2_val * sd);

	return val;
}

static double
dipHI (int flag, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
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
		hx = dipxHI (flag, fault, mag, xi, et, q, y, z);
		hy = dipyHI (flag, fault, mag, xi, et, q, y, z);
		hz = dipzHI (flag, fault, mag, xi, et, q, y, z);

		res[i] = mag->cx * hx + mag->cy * hy + mag->cz * hz;
	}
	return (res[0] + res[3]) - (res[1] + res[2]);
}

/*** contributions from the mirror image HIII ***/
static double
dipxHIII (int flag, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K4_val = K4 (flag, 1.0, xi, et, qq);
	double log_re_val = log_re (flag, 1.0, xi, et, qq);
	double J2_val = J2 (flag, 1.0, xi, et, qq);
	double O1_val = O1 (flag, 1.0, xi, et, qq);
	double N2_val = N2 (flag, 1.0, xi, et, qq);

	val = alpha4 * K4_val - alpha3 * (log_re_val * s2d - J2_val * c2d * secd)
		- alpha3 * (qd * O1_val + (z - h) * N2_val * sd);

	return val;
}

static double
dipyHIII (int flag, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K5_val = K5 (flag, 1.0, xi, et, qq);
	double log_rx_val = log_rx (flag, 1.0, xi, et, qq);
	double atan_xe_qr_val = atan_xe_qr (flag, 1.0, xi, et, qq);
	double J1_val = J1 (flag, 1.0, xi, et, qq);
	double O2_val = O2 (flag, 1.0, xi, et, qq);
	double O3_val = O3 (flag, 1.0, xi, et, qq);
	double N1_val = N1 (flag, 1.0, xi, et, qq);

	val = - alpha4 * K5_val + fault->alpha * log_rx_val * sd
		- alpha4 * atan_xe_qr_val * cd
		+ alpha3 * J1_val * c2d * secd
		- alpha3 * (qd * O2_val - (z - h) * (N1_val - O3_val) * sd);

	return val;
}

static double
dipzHIII (int flag, const fault_params *fault, const magnetic_params *mag, double xi, double et, double qq, double y, double z)
{
	double val;
	double h = mag->dcurier;
	double qd = y * sd + (fault->fdepth - h) * cd;
	double K6_val = K6 (flag, 1.0, xi, et, qq);
	double log_rx_val = log_rx (flag, 1.0, xi, et, qq);
	double atan_xe_qr_val = atan_xe_qr (flag, 1.0, xi, et, qq);
	double O2_val = O2 (flag, 1.0, xi, et, qq);
	double O3_val = O3 (flag, 1.0, xi, et, qq);

	val = alpha5 * K6_val
		- fault->alpha * log_rx_val * cd - alpha4 * atan_xe_qr_val * sd
		- alpha3 * (qd * O3_val - (z - h) * O2_val * sd);

	return val;
}

static double
dipHIII (int flag, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
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
		hx = dipxHIII (flag, fault, mag, xi, et, q, y, z);
		hy = dipyHIII (flag, fault, mag, xi, et, q, y, z);
		hz = dipzHIII (flag, fault, mag, xi, et, q, y, z);

		res[i] = mag->cx * hx + mag->cy * hy + mag->cz * hz;
	}
	return (res[0] + res[3]) - (res[1] + res[2]);
}


static double
dipHII (int flag, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
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
		hx = dipxHI (flag, fault, mag, xi, et, q, y, z);
		hy = dipyHI (flag, fault, mag, xi, et, q, y, z);
		hz = dipzHI (flag, fault, mag, xi, et, q, y, z);

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
		hx = dipxHIII (flag, fault, mag, xi, et, q, y, z);
		hy = dipyHIII (flag, fault, mag, xi, et, q, y, z);
		hz = dipzHIII (flag, fault, mag, xi, et, q, y, z);

		res[i] += mag->cx * hx + mag->cy * hy + mag->cz * hz;
	}
	return (res[0] + res[3]) - (res[1] + res[2]);
}


double
dip_slip (int flag, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
{
	double res;

	res = dip0 (flag, fault, mag, x, y, z);
	res += dipH0 (flag, fault, mag, x, y, z);
	if (fault->fdepth + fault->fwidth2 * sd < mag->dcurier) res += dipHI (flag, fault, mag, x, y, z);
	else if (fault->fdepth - fault->fwidth1 * sd > mag->dcurier) res += dipHIII (flag, fault, mag, x, y, z);
	else res += dipHII (flag, fault, mag, x, y, z);

	return res;
}
