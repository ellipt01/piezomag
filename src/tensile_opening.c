#include <stdio.h>
#include <math.h>

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

static double
tensile0 (MagComp component, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
{
	int		i;
	double res[4];
	double p, q;
	double sign = 1.0;
	double hx, hy, hz;

	p = y * cd - (d[1] - z) * sd;
	q = y * sd + (d[1] - z) * cd;

	for (i = 0; i < 4; i++) {
		double xi, et;

		xi = x + fault->flength1;
		et = p + fault->fwidth1;

		if (i >= 2)		 xi = x - fault->flength2;
		if (i % 2 == 0) et = p - fault->fwidth2;

		set_geometry_variables (sign, xi, et, q);
		hx = tensilex0 (component, xi, et, q);
		hy = tensiley0 (component, xi, et, q);
		hz = tensilez0 (component, xi, et, q);

		res[i] = mag->cx * hx + mag->cy * hy + mag->cz * hz;
	}
	return (res[0] + res[3]) - (res[1] + res[2]);
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

static double
tensileH0 (MagComp component, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
{
	int		i;
	double res[4];
	double p, q;
	double sign = 1.0;
	double hx, hy, hz;

	p = y * cd - (d[3] - z) * sd;
	q = y * sd + (d[3] - z) * cd;

	for (i = 0; i < 4; i++) {
		double xi, et;

		xi = x + fault->flength1;
		et = p + fault->fwidth1;

		if (i >= 2)		 xi = x - fault->flength2;
		if (i % 2 == 0) et = p - fault->fwidth2;

		set_geometry_variables (sign, xi, et, q);
		hx = tensilexH0 (component, fault, mag, xi, et, q, y, z);
		hy = tensileyH0 (component, fault, mag, xi, et, q, y, z);
		hz = tensilezH0 (component, fault, mag, xi, et, q, y, z);

		res[i] = mag->cx * hx + mag->cy * hy + mag->cz * hz;
	}
	return res[0] - res[1] - res[2] + res[3];
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

static double
tensileHI (MagComp component, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
{
	int		i;
	double res[4];
	double p, q;
	double sign = -1.0;
	double hx, hy, hz;

	p = y * cd - (d[2] + z) * sd;
	q = y * sd + (d[2] + z) * cd;

	for (i = 0; i < 4; i++) {
		double xi, et;

		xi = x + fault->flength1;
		et = p + fault->fwidth1;

		if (i >= 2)		 xi = x - fault->flength2;
		if (i % 2 == 0) et = p - fault->fwidth2;

		set_geometry_variables (sign, xi, et, q);
		hx = tensilexHI (component, fault, mag, xi, et, q, y, z);
		hy = tensileyHI (component, fault, mag, xi, et, q, y, z);
		hz = tensilezHI (component, fault, mag, xi, et, q, y, z);

		res[i] = mag->cx * hx + mag->cy * hy + mag->cz * hz;
	}
	return (res[0] + res[3]) - (res[1] + res[2]);
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

static double
tensileHIII (MagComp component, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
{
	int		i;
	double res[4];
	double p, q;
	double sign = 1.0;
	double hx, hy, hz;

	p = y * cd - (d[1] - z) * sd;
	q = y * sd + (d[1] - z) * cd;

	for (i = 0; i < 4; i++) {
		double xi, et;

		xi = x + fault->flength1;
		et = p + fault->fwidth1;

		if (i >= 2)		 xi = x - fault->flength2;
		if (i % 2 == 0) et = p - fault->fwidth2;

		set_geometry_variables (sign, xi, et, q);
		hx = tensilexHIII (component, fault, mag, xi, et, q, y, z);
		hy = tensileyHIII (component, fault, mag, xi, et, q, y, z);
		hz = tensilezHIII (component, fault, mag, xi, et, q, y, z);

		res[i] = mag->cx * hx + mag->cy * hy + mag->cz * hz;
	} 
	return (res[0] + res[3]) - (res[1] + res[2]);
}


static double
tensileHII (MagComp component, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
{
	int		i;
	double res[4];
	double p, q;
	double sign;
	double w = (mag->dcurier - fault->fdepth) / sd;
	double hx, hy, hz;

	p = y * cd - (d[2] + z) * sd;
	q = y * sd + (d[2] + z) * cd;

	for (i = 0; i < 4; i++) {
		double xi, et;

		xi = x + fault->flength1;
		et = p + fault->fwidth1;

		if (i >= 2)		 xi = x - fault->flength2;
		if (i % 2 == 0) et = p - w;

		sign = - 1.0;
		set_geometry_variables (sign, xi, et, q);
		hx = tensilexHI (component, fault, mag, xi, et, q, y, z);
		hy = tensileyHI (component, fault, mag, xi, et, q, y, z);
		hz = tensilezHI (component, fault, mag, xi, et, q, y, z);

		res[i] = mag->cx * hx + mag->cy * hy + mag->cz * hz;
	}

	p = y * cd - (d[1] - z) * sd;
	q = y * sd + (d[1] - z) * cd;

	for (i = 0; i < 4; i++) {
		double xi, et;

		xi = x + fault->flength1;
		et = p - w;

		if (i >= 2)		 xi = x - fault->flength2;
		if (i % 2 == 0) et = p - fault->fwidth2;

		sign = 1.0;
		set_geometry_variables (sign, xi, et, q);
		hx = tensilexHIII (component, fault, mag, xi, et, q, y, z);
		hy = tensileyHIII (component, fault, mag, xi, et, q, y, z);
		hz = tensilezHIII (component, fault, mag, xi, et, q, y, z);

		res[i] += mag->cx * hx + mag->cy * hy + mag->cz * hz;
	}
	return (res[0] + res[3]) - (res[1] + res[2]);
}

static double
tensile_opening_main (MagComp component, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
{
	return tensile0 (component, fault, mag, x, y, z);
}

static double
tensile_opening_mirror_image (MagComp component, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
{
	return tensileH0 (component, fault, mag, x, y, z);
}

static double
tensile_opening_submirror_image (MagComp component, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
{
	double	val;
	if (fault->fdepth + fault->fwidth2 * sd < mag->dcurier) val = tensileHI (component, fault, mag, x, y, z);
	else if (fault->fdepth - fault->fwidth1 * sd > mag->dcurier) val = tensileHIII (component, fault, mag, x, y, z);
	else val = tensileHII (component, fault, mag, x, y, z);
	return val;
}

/*** public functions ***/
double
tensile_opening (MagComp component, SeismoMagTerm term, const fault_params *fault, const magnetic_params *mag, double x, double y, double z)
{
	double res = 0.0;
	if (!check_mag_component (component)) {
		fprintf (stderr, "ERROR: tensile_opening: component must be MAG_COMP_X, MAG_COMP_Y, MAG_COMP_Z or MAG_COMP_F\n");
		return 0.0;
	}
	if (!check_seismo_mag_term (term)) {
		fprintf (stderr, "ERROR: tensile_opening: term must be SEISMO_MAG_MAIN, SEISMO_MAG_MIRROR, SEISMO_MAG_SUBMIRROR\n");
		return 0.0;
	}
	if (term & SEISMO_MAG_MAIN) res += tensile_opening_main (component, fault, mag, x, y, z);
	if (term & SEISMO_MAG_MIRROR) res += tensile_opening_mirror_image (component, fault, mag, x, y, z);
	if (term & SEISMO_MAG_SUBMIRROR) res += tensile_opening_submirror_image (component, fault, mag, x, y, z);
	return res;
}
