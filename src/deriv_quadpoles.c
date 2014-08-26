/*
 * quadpole.c
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

/*c****************************
 *c******** quadpoles *********
 *c****************************/

static double
L0 (double sign, double xi, double et, double qq)
{
	double	val;
	if (fault_is_vertical) val = - sign * qq * ir * irc2;
	else val = (ir * irc + sign * sd * (ir * ire)) * secd;
	return val;
}

static double
L0_x (double sign, double xi, double et, double qq)
{
	double	val;
	if (fault_is_vertical) val = sign * xi * qq * (2. * r + rc) * ir3c2 * irc;
	else val = (- xi * (r + rc) * ir3c2 - sign * xi * (r + re) * sd * ir3e2) * secd;
	return val;
}

static double
L0_y (double sign, double xi, double et, double qq)
{
	double	val;
	if (fault_is_vertical) val = - sign * ir * irc2 + sign * pow (qq, 2.) * (2. * r + rc) * ir3c2 * irc;
	else val = (- yy * (r + rc) * ir3c2
		- sign * yy * sd * (r + re) * ir3e2 - sign * sd * cd * (ir * ire2)) * secd;
	return val;
}

static double
L0_z (double sign, double xi, double et, double qq)
{
	double	val;
	if (fault_is_vertical) val = - sign * qq * (r + rc) * ir3c2;
	else val = (ir3 + sign * cc * sd * (r + re) * ir3e2 - sd2 * (ir * ire2)) * secd;
	return val;
}

static double
L1_x (double sign, double xi, double et, double qq)
{
	double l0 = L0 (sign, xi, et, qq);
	double l0_x = L0_x (sign, xi, et, qq);
	return l0 + xi * l0_x;
}

static double
L1_y (double sign, double xi, double et, double qq)
{
	double l0_y = L0_y (sign, xi, et, qq);
	return xi * l0_y;
}

static double
L1_z (double sign, double xi, double et, double qq)
{
	double l0_z = L0_z (sign, xi, et, qq);
	return xi * l0_z;
}


static double
L2_x (double sign, double xi, double et, double qq)
{
	double	r_l0_x = L0_x (sign, xi, et, qq);
	return - sign * sd * xi * ir * ire2 + yy * r_l0_x;
}

static double
L2_y (double sign, double xi, double et, double qq)
{
	double	val;
	double	r_l0 = L0 (sign, xi, et, qq);
	double	r_l0_y = L0_y (sign, xi, et, qq);
	val = - sign * sd * yy * ir * ire2 + r_l0 + yy * r_l0_y;
	if (!fault_is_vertical) val += - sign * sd * cd * ire2 + yy;
	return val;
}

static double
L2_z (double sign, double xi, double et, double qq)
{
	double r_l0_z = L0_z (sign, xi, et, qq);
	return - sd * (sd * ire2 - sign * cc * ir * ire2) + yy * r_l0_z;
}

static double
M1_x (double sign, double xi, double et, double qq)
{
	return (ir * ire) - xi * xi * (r + re) * ir3e2;
}

static double
M1_y (double sign, double xi, double et, double qq)
{
	double	val = - xi * yy * (r + re) * ir3e2;
	if (!fault_is_vertical) val -= xi * cd * (ir * ire2);
	return val;
}

static double
M1_z (double sign, double xi, double et, double qq)
{
	double	val;
	if (fault_is_vertical) val = - sign * xi * ir3;
	else val = - sign * xi * sd * (ir * ire2) + xi * cc * (r + re) * ir3e2;
	return val;
}

static double
M2_x (double sign, double xi, double et, double qq)
{
	return M1_y (sign, xi, et, qq);
}

static double
M2_y (double sign, double xi, double et, double qq)
{
	double	val;
	if (fault_is_vertical) val = ir * ire - pow (qq, 2.) * (r + re) * ir3e2;
	else val = sd2 * ire2 - (yy * cd + sign * cc * sd) * (ir * ire2)
		- yy * yy * (r + re) * ir3e2;
	return val;
}

static double
M2_z (double sign, double xi, double et, double qq)
{
	double	val;
	if (fault_is_vertical) val = - sign * qq * ir3;
	else val = - sign * sd * cd * ire2 + (cc * cd - sign * yy * sd) * (ir * ire2)
		+ yy * cc * (r + re) * ir3e2;
	return val;
}

static double
M3_x (double sign, double xi, double et, double qq)
{
	return M1_z (sign, xi, et, qq);
}

static double
M3_y (double sign, double xi, double et, double qq)
{
	return M2_z (sign, xi, et, qq);
}

static double
M3_z (double sign, double xi, double et, double qq)
{
	double	val;
	if (fault_is_vertical) val = - et * ir3;
	else val = cd2 * ire2 + (yy * cd + sign * cc * sd) * (ir * ire2)
		- cc * cc * (r + re) * ir3e2;
	return val;
}

static double
N1_x (double sign, double xi, double et, double qq)
{
	return (ir * irc) - xi * xi * (r + rc) * ir3c2;
}

static double
N1_y (double sign, double xi, double et, double qq)
{
	return - xi * yy * (r + rc) * ir3c2;
}

static double
N1_z (double sign, double xi, double et, double qq)
{
	return xi * ir3;
}

static double
N2_x (double sign, double xi, double et, double qq)
{
	return N1_y (sign, xi, et, qq);
}

static double
N2_y (double sign, double xi, double et, double qq)
{
	return (ir * irc) - yy * yy * (r + rc) * ir3c2;
}

static double
N2_z (double sign, double xi, double et, double qq)
{
	return yy * ir3;
}

static double
O1_x (double sign, double xi, double et, double qq)
{
	return - xi * ir3;
}

static double
O1_y (double sign, double xi, double et, double qq)
{
	return - yy * ir3;
}

static double
O1_z (double sign, double xi, double et, double qq)
{
	return cc * ir3;
}

static double
O2_x (double sign, double xi, double et, double qq)
{
	return O1_y (sign, xi, et, qq);
}

static double
O2_y (double sign, double xi, double et, double qq)
{
	return (ir * irx) - yy * yy * (r + rx) * ir3x2;
}

static double
O2_z (double sign, double xi, double et, double qq)
{
	return yy * cc * (r + rx) * ir3x2;
}

static double
O3_x (double sign, double xi, double et, double qq)
{
	return O1_z (sign, xi, et, qq);
}

static double
O3_y (double sign, double xi, double et, double qq)
{
	return O2_z (sign, xi, et, qq);
}

static double
O3_z (double sign, double xi, double et, double qq)
{
	return 1.0 * (ir * irx) - cc * cc * (r + rx) * ir3x2;
}

static double
P1_x (double sign, double xi, double et, double qq)
{
	return xi * qq * (r + re) * ir3e2;
}

static double
P1_y (double sign, double xi, double et, double qq)
{
	return sign * cc * ir3 + sd * (ir * ire)
		- pow (xi, 2.) * (r + re) * sd * ir3e2;
}

static double
P1_z (double sign, double xi, double et, double qq)
{
	double	val = sign * yy * ir3;
	if (!fault_is_vertical)
		val += sign * (- (ir * ire) + pow (xi, 2.) * (r + re) * ir3e2) * cd;
	return val;
}

static double
P2_x (double sign, double xi, double et, double qq)
{
	return P1_y (sign, xi, et, qq);
}

static double
P2_y (double sign, double xi, double et, double qq)
{
	double	val = sign * yy * cc * (r + rx) * ir3x2
		- xi * yy * (r + re) * sd * ir3e2;
	if (!fault_is_vertical) val -= xi * sd * cd * (ir * ire2);
	return val;
}

static double
P2_z (double sign, double xi, double et, double qq)
{
	double	val = sign * ((ir * irx) - pow (cc, 2.) * (r + rx) * ir3x2);
	if (fault_is_vertical) val += - sign * xi * ir3;
	else val += - sign * xi * sd2 * (ir * ire2) + xi * cc * (r + re) * sd * ir3e2;
	return val;
}

static double
P3_x (double sign, double xi, double et, double qq)
{
	return P1_z (sign, xi, et, qq);
}

static double
P3_y (double sign, double xi, double et, double qq)
{
	return P2_z (sign, xi, et, qq);
}

static double
P3_z (double sign, double xi, double et, double qq)
{
	double	val = - sign * yy * cc * (r + rx) * ir3x2;
	if (!fault_is_vertical)
		val += xi * (sd * (ir * ire2) - sign * cc * (r + re) * ir3e2) * cd;
	return val;
}

/*****************************************************/
double
L1 (int flag, double sign, double xi, double et, double qq)
{
	if (flag == DERIV_X) return L1_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return L1_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return L1_z (sign, xi, et, qq);
	return 0.;
}

double
L2 (int flag, double sign, double xi, double et, double qq)
{
	if (flag == DERIV_X) return L2_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return L2_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return L2_z (sign, xi, et, qq);
	return 0.;
}

double
M1 (int flag, double sign, double xi, double et, double qq)
{
	if (flag == DERIV_X) return M1_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return M1_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return M1_z (sign, xi, et, qq);
	return 0.;
}

double
M2 (int flag, double sign, double xi, double et, double qq)
{
	if (flag == DERIV_X) return M2_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return M2_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return M2_z (sign, xi, et, qq);
	return 0.;
}

double
M3 (int flag, double sign, double xi, double et, double qq)
{
	if (flag == DERIV_X) return M3_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return M3_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return M3_z (sign, xi, et, qq);
	return 0.;
}

double
N1 (int flag, double sign, double xi, double et, double qq)
{
	if (flag == DERIV_X) return N1_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return N1_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return N1_z (sign, xi, et, qq);
	return 0.;
}

double
N2 (int flag, double sign, double xi, double et, double qq)
{
	if (flag == DERIV_X) return N2_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return N2_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return N2_z (sign, xi, et, qq);
	return 0.;
}

double
O1 (int flag, double sign, double xi, double et, double qq)
{
	if (flag == DERIV_X) return O1_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return O1_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return O1_z (sign, xi, et, qq);
	return 0.;
}

double
O2 (int flag, double sign, double xi, double et, double qq)
{
	if (flag == DERIV_X) return O2_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return O2_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return O2_z (sign, xi, et, qq);
	return 0.;
}

double
O3 (int flag, double sign, double xi, double et, double qq)
{
	if (flag == DERIV_X) return O3_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return O3_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return O3_z (sign, xi, et, qq);
	return 0.;
}

double
P1 (int flag, double sign, double xi, double et, double qq)
{
	if (flag == DERIV_X) return P1_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return P1_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return P1_z (sign, xi, et, qq);
	return 0.;
}

double
P2 (int flag, double sign, double xi, double et, double qq)
{
	if (flag == DERIV_X) return P2_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return P2_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return P2_z (sign, xi, et, qq);
	return 0.;
}

double
P3 (int flag, double sign, double xi, double et, double qq)
{
	if (flag == DERIV_X) return P3_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return P3_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return P3_z (sign, xi, et, qq);
	return 0.;
}
