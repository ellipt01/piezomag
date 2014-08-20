/*
 * dipole.c
 *
 *  Created on: 2014/08/20
 *      Author: utsugi
 */

#include <stdio.h>
#include <math.h>

#include "piezomag.h"
#include "private.h"

#define DERIV_X 1
#define DERIV_Y 2
#define DERIV_Z 3

/*c**************************
 *c******** dipoles *********
 *c**************************/

/* log (R_i + xi) */
double
log_rx_x (double sign, double xi, double et, double qq)
{
	return ir;
}

double
log_rx_y (double sign, double xi, double et, double qq)
{
	return yy * (ir * irx);
}

double
log_rx_z (double sign, double xi, double et, double qq)
{
	return - cc * (ir * irx);
}

/* log (R_i + eta_i) */
double
log_re_x (double sign, double xi, double et, double qq)
{
	return xi * ir * ire;
}

double
log_re_y (double sign, double xi, double et, double qq)
{
	return cd * ire + yy * ir * ire;
}

double
log_re_z (double sign, double xi, double et, double qq)
{
	return sign * sd * ire - cc * ir * ire;
}

/* log (R_i + c_i) */
double
log_rc_x (double sign, double xi, double et, double qq)
{
	return xi * ir * irc;
}

double
log_rc_y (double sign, double xi, double et, double qq)
{
	return yy * ir * irc;
}

double
log_rc_z (double sign, double xi, double et, double qq)
{
	return - ir;
}

/* atan( (xi eta_i)/(q_i R_i) ) */
double
atan_xe_qr_x (double sign, double xi, double et, double qq)
{
	return - qq * (ir * ire);
}

double
atan_xe_qr_y (double sign, double xi, double et, double qq)
{
	return xi * sd * (ir * ire) - sign * cc * (ir * irx);
}

double
atan_xe_qr_z (double sign, double xi, double et, double qq)
{
	return - sign * xi * cd * (ir * ire) - sign * yy * (ir * irx);
}

/* J_1 */
double
J1_x (double sign, double xi, double et, double qq)
{
	//	double res = qq * (ir * ire) + sign * yy * (ir * irc);
	return ir * (qq * ire + sign * yy * irc);
}

double
J1_y (double sign, double xi, double et, double qq)
{
	//	double res = - xi * sd * (ir * ire) - sign * xi * (ir * irc);
	return - xi * ir * (sd * ire + sign * irc);
}

double
J1_z (double sign, double xi, double et, double qq)
{
	return sign * xi * cd * (ir * ire);
}

/* J_2 */
double
J2_x (double sign, double xi, double et, double qq)
{
	//	double res = xi * (ir * irc) + sign * xi * sd * (ir * ire);
	return xi * ir * (irc + sign * sd * ire);
}

double
J2_y (double sign, double xi, double et, double qq)
{
	/*
	double res = yy * (ir * irc)
		+ sign * sd * cd * ire + sign * yy * sd * (ir * ire);
	*/
	return yy * (ir * irc) + sign * sd * ire * (cd + yy * ir);
}

double
J2_z (double sign, double xi, double et, double qq)
{
	//	double res = - ir + sd * sd * ire - sign * cc * sd * (ir * ire);
	return - ir + sd * ire * (sd - sign * cc * ir);
}

double
K1_x (double sign, double xi, double et, double qq)
{
	double r_J1_x = J1_x (sign, xi, et, qq);
	return td * (irc - xi * xi * (ir * irc2) + r_J1_x * td);
}

double
K1_y (double sign, double xi, double et, double qq)
{
	double r_J1_y = J1_y (sign, xi, et, qq);
	return td * (- xi * yy * (ir * irc2) + r_J1_y * td);
}

double
K1_z (double sign, double xi, double et, double qq)
{
	double r_J1_z = J1_z (sign, xi, et, qq);
	return td * (xi * irc2 + xi * cc * (ir * irc2) + r_J1_z * td);
}

/* vertical fault */
static double
K2_x_0 (double sign, double xi, double et, double qq)
{
	return - 0.5 * xi * et * (ir * irc2)
			+ sign * xi * pow (qq, 2.) * (ir * irc3)
			+ 0.5 * xi * (ir * ire);
}

static double
K2_y_0 (double sign, double xi, double et, double qq)
{
	return - 0.5 * et * qq * (ir * irc2)
			- sign * qq * irc2
			+ sign * pow (qq, 3.) * (ir * irc3)
			+ 0.5 * qq * (ir * ire);
}

static double
K2_z_0 (double sign, double xi, double et, double qq)
{
	return irc - pow (qq, 2.) * (ir * irc2);
}

/* inclined fault */
static double
K2_x_1 (double sign, double xi, double et, double qq)
{
	double r_J2_x = J2_x (sign, xi, et, qq);
	return td * (- xi * yy * (ir * irc2) - sign * r_J2_x * td);
}

static double
K2_y_1 (double sign, double xi, double et, double qq)
{
	double r_J2_y = J2_y (sign, xi, et, qq);
	return td * (irc - yy * yy * (ir * irc2) - sign * r_J2_y * td);
}

static double
K2_z_1 (double sign, double xi, double et, double qq)
{
	double r_J2_z = J2_z (sign, xi, et, qq);
	return td * (yy * irc2 + yy * cc * (ir * irc2) - sign * r_J2_z * td);
}


double
K3_x (double sign, double xi, double et, double qq)
{
	double r_J2_x = J2_x (sign, xi, et, qq);
	return r_J2_x * td;
}

double
K3_y (double sign, double xi, double et, double qq)
{
	double r_J2_y = J2_y (sign, xi, et, qq);
	return r_J2_y * td;
}

double
K3_z (double sign, double xi, double et, double qq)
{
	double r_J2_z = J2_z (sign, xi, et, qq);
	return r_J2_z * td;
}

double
K4_x_0 (double sign, double xi, double et, double qq)
{
	double r_K2_x = K2_x_0 (sign, xi, et, qq);
	return cd * (r_K2_x - xi * sd * (ir * ire));
}

double
K4_y_0 (double sign, double xi, double et, double qq)
{
	double r_K2_y = K2_y_0 (sign, xi, et, qq);
	return cd * (r_K2_y - sd * (cd * ire + yy * (ir * ire)));
}

double
K4_z_0 (double sign, double xi, double et, double qq)
{
	double r_K2_z = K2_z_0 (sign, xi, et, qq);
	return cd * (r_K2_z - sd * (sign * sd * ire - cc * (ir * ire)));
}

double
K4_x_1 (double sign, double xi, double et, double qq)
{
	double r_K2_x = K2_x_1 (sign, xi, et, qq);
	return cd * (r_K2_x - xi * sd * (ir * ire));
}

double
K4_y_1 (double sign, double xi, double et, double qq)
{
	double r_K2_y = K2_y_1 (sign, xi, et, qq);
	return cd * (r_K2_y - sd * (cd * ire + yy * (ir * ire)));
}

double
K4_z_1 (double sign, double xi, double et, double qq)
{
	double r_K2_z = K2_z_1 (sign, xi, et, qq);
	return cd * (r_K2_z - sd * (sign * sd * ire - cc * (ir * ire)));
}

double
K5_x (double sign, double xi, double et, double qq)
{
	double r_K1_x = K1_x (sign, xi, et, qq);
	return r_K1_x * cd;
}

double
K5_y (double sign, double xi, double et, double qq)
{
	double r_K1_y = K1_y (sign, xi, et, qq);
	return r_K1_y * cd;
}

double
K5_z (double sign, double xi, double et, double qq)
{
	double r_K1_z = K1_z (sign, xi, et, qq);
	return r_K1_z * cd;
}

double
K6_x (double sign, double xi, double et, double qq)
{
	double r_J1_x = J1_x (sign, xi, et, qq);
	return - r_J1_x * sd;
}

double
K6_y (double sign, double xi, double et, double qq)
{
	double r_J1_y = J1_y (sign, xi, et, qq);
	return - r_J1_y * sd;
}

double
K6_z (double sign, double xi, double et, double qq)
{
	double r_J1_z = J1_z (sign, xi, et, qq);
	return - r_J1_z * sd;
}

double
K7_x_0 (double sign, double xi, double et, double qq)
{
	double r_K4_x = K4_x_0 (sign, xi, et, qq);
	return - r_K4_x * td;
}

double
K7_y_0 (double sign, double xi, double et, double qq)
{
	double r_K4_y = K4_y_0 (sign, xi, et, qq);
	return - r_K4_y * td;
}

double
K7_z_0 (double sign, double xi, double et, double qq)
{
	double r_K4_z = K4_z_0 (sign, xi, et, qq);
	return - r_K4_z * td;
}

double
K7_x_1 (double sign, double xi, double et, double qq)
{
	double r_K4_x = K4_x_1 (sign, xi, et, qq);
	return - r_K4_x * td;
}

double
K7_y_1 (double sign, double xi, double et, double qq)
{
	double r_K4_y = K4_y_1 (sign, xi, et, qq);
	return - r_K4_y * td;
}

double
K7_z_1 (double sign, double xi, double et, double qq)
{
	double r_K4_z = K4_z_1 (sign, xi, et, qq);
	return - r_K4_z * td;
}

double
K8_x (double sign, double xi, double et, double qq)
{
	double r_K5_x = K5_x (sign, xi, et, qq);
	return - r_K5_x * td;
}

double
K8_y (double sign, double xi, double et, double qq)
{
	double r_K5_y = K5_y (sign, xi, et, qq);
	return - r_K5_y * td;
}

double
K8_z (double sign, double xi, double et, double qq)
{
	double r_K5_z = K5_z (sign, xi, et, qq);
	return - r_K5_z * td;
}

double
K9_x (double sign, double xi, double et, double qq)
{
	double r_K6_x = K6_x (sign, xi, et, qq);
	return - r_K6_x * td;
}

double
K9_y (double sign, double xi, double et, double qq)
{
	double r_K6_y = K6_y (sign, xi, et, qq);
	return - r_K6_y * td;
}

double
K9_z (double sign, double xi, double et, double qq)
{
	double r_K6_z = K6_z (sign, xi, et, qq);
	return - r_K6_z * td;
}

/*****************************************************/
double
log_rx (int flag, double sign, double xi, double et, double qq)
{
	if (flag == DERIV_X) return log_rx_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return log_rx_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return log_rx_z (sign, xi, et, qq);
	return 0.;
}

double
log_re (int flag, double sign, double xi, double et, double qq)
{
	if (flag == DERIV_X) return log_re_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return log_re_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return log_re_z (sign, xi, et, qq);
	return 0.;
}

double
log_rc (int flag, double sign, double xi, double et, double qq)
{
	if (flag == DERIV_X) return log_rc_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return log_rc_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return log_rc_z (sign, xi, et, qq);
	return 0.;
}

double
atan_xe_qr (int flag, double sign, double xi, double et, double qq)
{
	if (flag == DERIV_X) return atan_xe_qr_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return atan_xe_qr_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return atan_xe_qr_z (sign, xi, et, qq);
	return 0.;
}

double
J1 (int flag, double sign, double xi, double et, double qq)
{
	if (flag == DERIV_X) return J1_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return J1_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return J1_z (sign, xi, et, qq);
	return 0.;
}

double
J2 (int flag, double sign, double xi, double et, double qq)
{
	if (flag == DERIV_X) return J2_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return J2_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return J2_z (sign, xi, et, qq);
	return 0.;
}

double
K1 (int flag, double sign, double xi, double et, double qq)
{
	return 0.; /*****/
	if (flag == DERIV_X) return K1_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return K1_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return K1_z (sign, xi, et, qq);
	return 0.;
}

double
K2 (int flag, double sign, double xi, double et, double qq)
{
	if (fault_is_vertical) {
		if (flag == DERIV_X) return K2_x_0 (sign, xi, et, qq);
		else if (flag == DERIV_Y) return K2_y_0 (sign, xi, et, qq);
		else if (flag == DERIV_Z) return K2_z_0 (sign, xi, et, qq);
	} else {
		if (flag == DERIV_X) return K2_x_1 (sign, xi, et, qq);
		else if (flag == DERIV_Y) return K2_y_1 (sign, xi, et, qq);
		else if (flag == DERIV_Z) return K2_z_1 (sign, xi, et, qq);
	}
	return 0.;
}

double
K3 (int flag, double sign, double xi, double et, double qq)
{
	return 0.; /*****/
	if (flag == DERIV_X) return K3_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return K3_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return K3_z (sign, xi, et, qq);
	return 0.;
}

double
K4 (int flag, double sign, double xi, double et, double qq)
{
	if (fault_is_vertical) {
		if (flag == DERIV_X) return K4_x_0 (sign, xi, et, qq);
		else if (flag == DERIV_Y) return K4_y_0 (sign, xi, et, qq);
		else if (flag == DERIV_Z) return K4_z_0 (sign, xi, et, qq);
	} else {
		if (flag == DERIV_X) return K4_x_1 (sign, xi, et, qq);
		else if (flag == DERIV_Y) return K4_y_1 (sign, xi, et, qq);
		else if (flag == DERIV_Z) return K4_z_1 (sign, xi, et, qq);
	}
	return 0.;
}

double
K5 (int flag, double sign, double xi, double et, double qq)
{
	if (flag == DERIV_X) return K5_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return K5_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return K5_z (sign, xi, et, qq);
	return 0.;
}

double
K6 (int flag, double sign, double xi, double et, double qq)
{
	if (flag == DERIV_X) return K6_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return K6_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return K6_z (sign, xi, et, qq);
	return 0.;
}

double
K7 (int flag, double sign, double xi, double et, double qq)
{
	if (fault_is_vertical) {
		if (flag == DERIV_X) return K7_x_0 (sign, xi, et, qq);
		else if (flag == DERIV_Y) return K7_y_0 (sign, xi, et, qq);
		else if (flag == DERIV_Z) return K7_z_0 (sign, xi, et, qq);
	} else {
		if (flag == DERIV_X) return K7_x_1 (sign, xi, et, qq);
		else if (flag == DERIV_Y) return K7_y_1 (sign, xi, et, qq);
		else if (flag == DERIV_Z) return K7_z_1 (sign, xi, et, qq);
	}
	return 0.;
}

double
K8 (int flag, double sign, double xi, double et, double qq)
{
	return 0.; /*****/
	if (flag == DERIV_X) return K8_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return K8_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return K8_z (sign, xi, et, qq);
	return 0.;
}

double
K9 (int flag, double sign, double xi, double et, double qq)
{
	return 0.; /*****/
	if (flag == DERIV_X) return K9_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return K9_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return K9_z (sign, xi, et, qq);
	return 0.;
}
