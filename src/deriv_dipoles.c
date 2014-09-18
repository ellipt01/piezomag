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

/*c**************************
 *c******** dipoles *********
 *c**************************/

/* log (R_i + xi) */
static double
log_rx_x (double sign, double xi, double et, double qq)
{
	return ir;
}

static double
log_rx_y (double sign, double xi, double et, double qq)
{
	return yy * (ir * irx);
}

static double
log_rx_z (double sign, double xi, double et, double qq)
{
	return - cc * (ir * irx);
}

/* log (R_i + eta_i) */
static double
log_re_x (double sign, double xi, double et, double qq)
{
	return xi * ir * ire;
}

static double
log_re_y (double sign, double xi, double et, double qq)
{
	return cd * ire + yy * ir * ire;
}

static double
log_re_z (double sign, double xi, double et, double qq)
{
	return sign * sd * ire - cc * ir * ire;
}

/* log (R_i + c_i) */
static double
log_rc_x (double sign, double xi, double et, double qq)
{
	return xi * ir * irc;
}

static double
log_rc_y (double sign, double xi, double et, double qq)
{
	return yy * ir * irc;
}

static double
log_rc_z (double sign, double xi, double et, double qq)
{
	return - ir;
}

/* atan( (xi eta_i)/(q_i R_i) ) */
static double
atan_xe_qr_x (double sign, double xi, double et, double qq)
{
	return - qq * (ir * ire);
}

static double
atan_xe_qr_y (double sign, double xi, double et, double qq)
{
	return xi * sd * (ir * ire) - sign * cc * (ir * irx);
}

static double
atan_xe_qr_z (double sign, double xi, double et, double qq)
{
	return - sign * xi * cd * (ir * ire) - sign * yy * (ir * irx);
}

/* J_1 */
static double
J1_x (double sign, double xi, double et, double qq)
{
	return (fault_is_vertical) ? 0. : ir * (qq * ire + sign * yy * irc);
}

static double
J1_y (double sign, double xi, double et, double qq)
{
	return (fault_is_vertical) ? 0. : - xi * ir * (sd * ire + sign * irc);
}

static double
J1_z (double sign, double xi, double et, double qq)
{
	return (fault_is_vertical) ? 0. : sign * xi * cd * (ir * ire);
}

/* J_2 */
static double
J2_x (double sign, double xi, double et, double qq)
{
	return (fault_is_vertical) ? 0. : xi * ir * (irc + sign * sd * ire);
}

static double
J2_y (double sign, double xi, double et, double qq)
{
	return (fault_is_vertical) ? 0. : yy * (ir * irc) + sign * sd * ire * (cd + yy * ir);
}

static double
J2_z (double sign, double xi, double et, double qq)
{
	return (fault_is_vertical) ? 0. : - ir + sd * ire * (sd - sign * cc * ir);
}

static double
K1_x (double sign, double xi, double et, double qq)
{
	double	val;
	if (fault_is_vertical) val = - 0.5 * sign * qq * irc2 + sign * pow (xi, 2.) * qq * (ir * irc3);
	else {
		double r_J1_x = J1_x (sign, xi, et, qq);
		val = (irc - xi * xi * (ir * irc2) + r_J1_x * td) * td;
	}
	return val;
}

static double
K1_y (double sign, double xi, double et, double qq)
{
	double	val;
	if (fault_is_vertical) val = - 0.5 * sign * xi * irc2 + sign * xi * pow (qq, 2.) * (ir * irc3);
	else {
		double r_J1_y = J1_y (sign, xi, et, qq);
		val = (- xi * yy * (ir * irc2) + r_J1_y * td) * td;
	}
	return val;
}

static double
K1_z (double sign, double xi, double et, double qq)
{
	double	val;
	if (fault_is_vertical) val = - sign * xi * qq * (ir * irc2);
	else {
		double r_J1_z = J1_z (sign, xi, et, qq);
		val = (xi * irc2 + xi * cc * (ir * irc2) + r_J1_z * td) * td;
	}
	return val;
}

static double
K2_x (double sign, double xi, double et, double qq)
{
	double	val;
	if (fault_is_vertical) {
		val = - 0.5 * xi * et * (ir * irc2)
				+ sign * xi * pow (qq, 2.) * (ir * irc3)
				+ 0.5 * xi * (ir * ire);
	} else {
		double r_J2_x = J2_x (sign, xi, et, qq);
		val = td * (- xi * yy * (ir * irc2) - sign * r_J2_x * td);
	}
	return val;
}

static double
K2_y (double sign, double xi, double et, double qq)
{
	double	val;
	if (fault_is_vertical) {
		val = - 0.5 * et * qq * (ir * irc2)
					- sign * qq * irc2
					+ sign * pow (qq, 3.) * (ir * irc3)
					+ 0.5 * qq * (ir * ire);
	} else {
		double r_J2_y = J2_y (sign, xi, et, qq);
		val = td * (irc - yy * yy * (ir * irc2) - sign * r_J2_y * td);
	}
	return val;
}

static double
K2_z (double sign, double xi, double et, double qq)
{
	double	val;
	if (fault_is_vertical) val = sign * (irc - pow (qq, 2.) * (ir * irc2));
	else {
		double r_J2_z = J2_z (sign, xi, et, qq);
		val = td * (yy * irc2 + yy * cc * (ir * irc2) - sign * r_J2_z * td);
	}
	return val;
}

static double
K3_x (double sign, double xi, double et, double qq)
{
	double	val;
	if (fault_is_vertical) val = - sign * xi * qq * (ir * irc2);
	else {
		double r_J2_x = J2_x (sign, xi, et, qq);
		val = r_J2_x * td;
	}
	return val;
}

static double
K3_y (double sign, double xi, double et, double qq)
{
	double	val;
	if (fault_is_vertical) val = sign * irc - sign * pow (qq, 2.) * (ir * irc2);
	else {
		double r_J2_y = J2_y (sign, xi, et, qq);
		val = r_J2_y * td;
	}
	return val;
}

static double
K3_z (double sign, double xi, double et, double qq)
{
	double	val;
	if (fault_is_vertical) val = sign * qq * (ir * irc);
	else {
		double r_J2_z = J2_z (sign, xi, et, qq);
		val = r_J2_z * td;
	}
	return val;
}

static double
K4_x (double sign, double xi, double et, double qq)
{
	double	val;
	if (fault_is_vertical) val = 0.;
	else {
		double r_K2_x = K2_x (sign, xi, et, qq);
		val = cd * (r_K2_x - xi * sd * (ir * ire));
	}
	return val;
}

static double
K4_y (double sign, double xi, double et, double qq)
{
	double	val;
	if (fault_is_vertical) val = 0.;
	else {
		double r_K2_y = K2_y (sign, xi, et, qq);
		val = cd * (r_K2_y - sd * (cd * ire + yy * (ir * ire)));
	}
	return val;
}

static double
K4_z (double sign, double xi, double et, double qq)
{
	double	val;
	if (fault_is_vertical) val = 0.;
	else {
		double r_K2_z = K2_z (sign, xi, et, qq);
		val = cd * (r_K2_z - sd * (sign * sd * ire - cc * (ir * ire)));
	}
	return val;
}

static double
K5_x (double sign, double xi, double et, double qq)
{
	double	val;
	if (fault_is_vertical) val = 0.;
	else {
		double r_K1_x = K1_x (sign, xi, et, qq);
		val = r_K1_x * cd;
	}
	return val;
}

static double
K5_y (double sign, double xi, double et, double qq)
{
	double	val;
	if (fault_is_vertical) val = 0.;
	else {
		double r_K1_y = K1_y (sign, xi, et, qq);
		val = r_K1_y * cd;
	}
	return val;
}

static double
K5_z (double sign, double xi, double et, double qq)
{
	double	val;
	if (fault_is_vertical) val = 0.;
	else {
		double r_K1_z = K1_z (sign, xi, et, qq);
		val = r_K1_z * cd;
	}
	return val;
}

static double
K6_x (double sign, double xi, double et, double qq)
{
	double	val;
	if (fault_is_vertical) val = 0.;
	else {
		double r_J1_x = J1_x (sign, xi, et, qq);
		val = - r_J1_x * sd;
	}
	return val;
}

static double
K6_y (double sign, double xi, double et, double qq)
{
	double	val;
	if (fault_is_vertical) val = 0.;
	else {
		double r_J1_y = J1_y (sign, xi, et, qq);
		val = - r_J1_y * sd;
	}
	return val;
}

static double
K6_z (double sign, double xi, double et, double qq)
{
	double	val;
	if (fault_is_vertical) val = 0.;
	else {
		double r_J1_z = J1_z (sign, xi, et, qq);
		val = - r_J1_z * sd;
	}
	return val;
}

/* inclined fault */
static double
K7_x (double sign, double xi, double et, double qq)
{
	double	val;
	if (fault_is_vertical) {
		double r_K2_x = K2_x (sign, xi, et, qq);
		val = - r_K2_x + xi * ir * ire;
	} else {
		double r_K4_x = K4_x (sign, xi, et, qq);
		val = - r_K4_x * td;
	}
	return val;
}

static double
K7_y (double sign, double xi, double et, double qq)
{
	double	val;
	if (fault_is_vertical) {
		double r_K2_y = K2_y (sign, xi, et, qq);
		val = - r_K2_y + qq * ir * ire;
	} else {
		double r_K4_y = K4_y (sign, xi, et, qq);
		val = - r_K4_y * td;
	}
	return val;
}

static double
K7_z (double sign, double xi, double et, double qq)
{
	double	val;
	if (fault_is_vertical) {
		double r_K2_z = K2_z (sign, xi, et, qq);
		val = - r_K2_z + sign * ir;
	} else {
		double r_K4_z = K4_z (sign, xi, et, qq);
		val = - r_K4_z * td;
	}
	return val;
}

static double
K8_x (double sign, double xi, double et, double qq)
{
	double	val;
	if (fault_is_vertical) {
		double r_K1_x = K1_x (sign, xi, et, qq);
		val = - r_K1_x * sd;
	} else {
		double r_K5_x = K5_x (sign, xi, et, qq);
		val = - r_K5_x * td;
	}
	return val;
}

static double
K8_y (double sign, double xi, double et, double qq)
{
	double	val;
	if (fault_is_vertical) {
		double r_K1_y = K1_y (sign, xi, et, qq);
		val = - r_K1_y * sd;
	} else {
		double r_K5_y = K5_y (sign, xi, et, qq);
		val = - r_K5_y * td;
	}
	return val;
}

static double
K8_z (double sign, double xi, double et, double qq)
{
	double	val;
	if (fault_is_vertical) {
		double r_K1_z = K1_z (sign, xi, et, qq);
		val = - r_K1_z * sd;
	} else {
		double r_K5_z = K5_z (sign, xi, et, qq);
		val = - r_K5_z * td;
	}
	return val;
}

static double
K9_x (double sign, double xi, double et, double qq)
{
	double	val;
	if (fault_is_vertical) val = - (sign * irc - sign * pow (xi, 2.) * (ir * irc2));
	else {
		double r_K6_x = K6_x (sign, xi, et, qq);
		val = - r_K6_x * td;
	}
	return val;
}

static double
K9_y (double sign, double xi, double et, double qq)
{
	double	val;
	if (fault_is_vertical) val = sign * xi * qq * (ir * irc2);
	else {
		double r_K6_y = K6_y (sign, xi, et, qq);
		val = - r_K6_y * td;
	}
	return val;
}

static double
K9_z (double sign, double xi, double et, double qq)
{
	double	val;
	if (fault_is_vertical) val = - sign * xi * (ir * irc);
	else {
		double r_K6_z = K6_z (sign, xi, et, qq);
		val = - r_K6_z * td;
	}
	return val;
}

/*****************************************************/
double
log_rx (MagComp component, double sign, double xi, double et, double qq)
{
	if (component == MAG_COMP_X) return log_rx_x (sign, xi, et, qq);
	if (component == MAG_COMP_Y) return log_rx_y (sign, xi, et, qq);
	if (component == MAG_COMP_Z) return log_rx_z (sign, xi, et, qq);
	return _PIEZOMAG_DUMMY_;
}

double
log_re (MagComp component, double sign, double xi, double et, double qq)
{
	if (component == MAG_COMP_X) return log_re_x (sign, xi, et, qq);
	if (component == MAG_COMP_Y) return log_re_y (sign, xi, et, qq);
	if (component == MAG_COMP_Z) return log_re_z (sign, xi, et, qq);
	return _PIEZOMAG_DUMMY_;
}

double
log_rc (MagComp component, double sign, double xi, double et, double qq)
{
	if (component == MAG_COMP_X) return log_rc_x (sign, xi, et, qq);
	if (component == MAG_COMP_Y) return log_rc_y (sign, xi, et, qq);
	if (component == MAG_COMP_Z) return log_rc_z (sign, xi, et, qq);
	return _PIEZOMAG_DUMMY_;
}

double
atan_xe_qr (MagComp component, double sign, double xi, double et, double qq)
{
	if (component == MAG_COMP_X) return atan_xe_qr_x (sign, xi, et, qq);
	if (component == MAG_COMP_Y) return atan_xe_qr_y (sign, xi, et, qq);
	if (component == MAG_COMP_Z) return atan_xe_qr_z (sign, xi, et, qq);
	return _PIEZOMAG_DUMMY_;
}

double
J1 (MagComp component, double sign, double xi, double et, double qq)
{
	if (component == MAG_COMP_X) return J1_x (sign, xi, et, qq);
	if (component == MAG_COMP_Y) return J1_y (sign, xi, et, qq);
	if (component == MAG_COMP_Z) return J1_z (sign, xi, et, qq);
	return _PIEZOMAG_DUMMY_;
}

double
J2 (MagComp component, double sign, double xi, double et, double qq)
{
	if (component == MAG_COMP_X) return J2_x (sign, xi, et, qq);
	if (component == MAG_COMP_Y) return J2_y (sign, xi, et, qq);
	if (component == MAG_COMP_Z) return J2_z (sign, xi, et, qq);
	return _PIEZOMAG_DUMMY_;
}

double
K1 (MagComp component, double sign, double xi, double et, double qq)
{
	if (component == MAG_COMP_X) return K1_x (sign, xi, et, qq);
	if (component == MAG_COMP_Y) return K1_y (sign, xi, et, qq);
	if (component == MAG_COMP_Z) return K1_z (sign, xi, et, qq);
	return _PIEZOMAG_DUMMY_;
}

double
K2 (MagComp component, double sign, double xi, double et, double qq)
{
	if (component == MAG_COMP_X) return K2_x (sign, xi, et, qq);
	if (component == MAG_COMP_Y) return K2_y (sign, xi, et, qq);
	if (component == MAG_COMP_Z) return K2_z (sign, xi, et, qq);
	return _PIEZOMAG_DUMMY_;
}

double
K3 (MagComp component, double sign, double xi, double et, double qq)
{
	if (component == MAG_COMP_X) return K3_x (sign, xi, et, qq);
	if (component == MAG_COMP_Y) return K3_y (sign, xi, et, qq);
	if (component == MAG_COMP_Z) return K3_z (sign, xi, et, qq);
	return _PIEZOMAG_DUMMY_;
}

double
K4 (MagComp component, double sign, double xi, double et, double qq)
{
	if (component == MAG_COMP_X) return K4_x (sign, xi, et, qq);
	if (component == MAG_COMP_Y) return K4_y (sign, xi, et, qq);
	if (component == MAG_COMP_Z) return K4_z (sign, xi, et, qq);
	return _PIEZOMAG_DUMMY_;
}

double
K5 (MagComp component, double sign, double xi, double et, double qq)
{
	if (component == MAG_COMP_X) return K5_x (sign, xi, et, qq);
	if (component == MAG_COMP_Y) return K5_y (sign, xi, et, qq);
	if (component == MAG_COMP_Z) return K5_z (sign, xi, et, qq);
	return _PIEZOMAG_DUMMY_;
}

double
K6 (MagComp component, double sign, double xi, double et, double qq)
{
	if (component == MAG_COMP_X) return K6_x (sign, xi, et, qq);
	if (component == MAG_COMP_Y) return K6_y (sign, xi, et, qq);
	if (component == MAG_COMP_Z) return K6_z (sign, xi, et, qq);
	return _PIEZOMAG_DUMMY_;
}

double
K7 (MagComp component, double sign, double xi, double et, double qq)
{
	if (component == MAG_COMP_X) return K7_x (sign, xi, et, qq);
	if (component == MAG_COMP_Y) return K7_y (sign, xi, et, qq);
	if (component == MAG_COMP_Z) return K7_z (sign, xi, et, qq);
	return _PIEZOMAG_DUMMY_;
}

double
K8 (MagComp component, double sign, double xi, double et, double qq)
{
	if (component == MAG_COMP_X) return K8_x (sign, xi, et, qq);
	if (component == MAG_COMP_Y) return K8_y (sign, xi, et, qq);
	if (component == MAG_COMP_Z) return K8_z (sign, xi, et, qq);
	return _PIEZOMAG_DUMMY_;
}

double
K9 (MagComp component, double sign, double xi, double et, double qq)
{
	if (component == MAG_COMP_X) return K9_x (sign, xi, et, qq);
	if (component == MAG_COMP_Y) return K9_y (sign, xi, et, qq);
	if (component == MAG_COMP_Z) return K9_z (sign, xi, et, qq);
	return _PIEZOMAG_DUMMY_;
}
