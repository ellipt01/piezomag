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

double
L0 (double sign, double xi, double et, double qq)
{
	return 1.0 * (ir * irc) + sign * sd * (ir * ire);
}

double
L0_x (double sign, double xi, double et, double qq)
{
	return - xi * (r + rc) * ir3c2 - sign * xi * (r + re) * sd * ir3e2;
}

double
L0_y (double sign, double xi, double et, double qq)
{
	return - yy * (r + rc) * ir3c2
		- sign * yy * sd * (r + re) * ir3e2 - sign * sd * cd * (ir * ire2);
}

double
L0_z (double sign, double xi, double et, double qq)
{
	return ir3 + sign * cc * sd * (r + re) * ir3e2 - sd2 * (ir * ire2);
}

double
L1_x (double sign, double xi, double et, double qq)
{
	double l0 = L0 (sign, xi, et, qq);
	double l0_x = L0_x (sign, xi, et, qq);
	return l0 + xi * l0_x;
}

double
L1_y (double sign, double xi, double et, double qq)
{
	double l0_y = L0_y (sign, xi, et, qq);
	return xi * l0_y;
}

double
L1_z (double sign, double xi, double et, double qq)
{
	double l0_z = L0_z (sign, xi, et, qq);
	return xi * l0_z;
}

double
L2_x (double sign, double xi, double et, double qq)
{
	double res = L1_y (sign, xi, et, qq);
	res *= secd;
	return res;
}

double
L2_y (double sign, double xi, double et, double qq)
{
	double l0_y = L0_y (sign, xi, et, qq);
	return sign * (sd * ire2 - sign * cc * (ir * ire2)) * sd * td
		+ ((ir * irc) + yy * l0_y) * secd;
}

double
L2_z (double sign, double xi, double et, double qq)
{
	double l0_z = L0_z (sign, xi, et, qq);
	return - (sd * ire2 - sign * cc * (ir * ire2)) * sd
		+ yy * l0_z * secd;
}

double
M1_x (double sign, double xi, double et, double qq)
{
	return (ir * ire) - xi * xi * (r + re) * ir3e2;
}

double
M1_y (double sign, double xi, double et, double qq)
{
	return - xi * cd * (ir * ire2) - xi * yy * (r + re) * ir3e2;
}

double
M1_z (double sign, double xi, double et, double qq)
{
	return - sign * xi * sd * (ir * ire2) + xi * cc * (r + re) * ir3e2;
}

double
M2_x (double sign, double xi, double et, double qq)
{
	return M1_y (sign, xi, et, qq);
}

double
M2_y (double sign, double xi, double et, double qq)
{
	return sd2 * ire2 - (yy * cd + sign * cc * sd) * (ir * ire2)
		- yy * yy * (r + re) * ir3e2;
}

double
M2_z (double sign, double xi, double et, double qq)
{
	return - sign * sd * cd * ire2 + (cc * cd - sign * yy * sd) * (ir * ire2)
		+ yy * cc * (r + re) * ir3e2;
}

double
M3_x (double sign, double xi, double et, double qq)
{
	return M1_z (sign, xi, et, qq);
}

double
M3_y (double sign, double xi, double et, double qq)
{
	return M2_z (sign, xi, et, qq);
}

double
M3_z (double sign, double xi, double et, double qq)
{
	return cd2 * ire2 + (yy * cd + sign * cc * sd) * (ir * ire2)
		- cc * cc * (r + re) * ir3e2;
}

double
N1_x (double sign, double xi, double et, double qq)
{
	return (ir * irc) - xi * xi * (r + rc) * ir3c2;
}

double
N1_y (double sign, double xi, double et, double qq)
{
	return - xi * yy * (r + rc) * ir3c2;
}

double
N1_z (double sign, double xi, double et, double qq)
{
	return xi * ir3;
}

double
N2_x (double sign, double xi, double et, double qq)
{
	return N1_y (sign, xi, et, qq);
}

double
N2_y (double sign, double xi, double et, double qq)
{
	return (ir * irc) - yy * yy * (r + rc) * ir3c2;
}

double
N2_z (double sign, double xi, double et, double qq)
{
	return yy * ir3;
}

double
O1_x (double sign, double xi, double et, double qq)
{
	return - xi * ir3;
}

double
O1_y (double sign, double xi, double et, double qq)
{
	return - yy * ir3;
}

double
O1_z (double sign, double xi, double et, double qq)
{
	return cc * ir3;
}

double
O2_x (double sign, double xi, double et, double qq)
{
	return O1_y (sign, xi, et, qq);
}

double
O2_y (double sign, double xi, double et, double qq)
{
	return (ir * irx) - yy * yy * (r + rx) * ir3x2;
}

double
O2_z (double sign, double xi, double et, double qq)
{
	return yy * cc * (r + rx) * ir3x2;
}

double
O3_x (double sign, double xi, double et, double qq)
{
	return O1_z (sign, xi, et, qq);
}

double
O3_y (double sign, double xi, double et, double qq)
{
	return O2_z (sign, xi, et, qq);
}

double
O3_z (double sign, double xi, double et, double qq)
{
	return 1.0 * (ir * irx) - cc * cc * (r + rx) * ir3x2;
}

double
P1_x (double sign, double xi, double et, double qq)
{
	return xi * qq * (r + re) * ir3e2;
}

double
P1_y (double sign, double xi, double et, double qq)
{
	return sign * cc * ir3 + sd * (ir * ire)
		- xi * xi * (r + re) * sd * ir3e2;
}

double
P1_z (double sign, double xi, double et, double qq)
{
	return sign * (yy * ir3 - cd * (ir * ire)
					 + xi * xi * (r + re) * cd * ir3e2);
}

double
P2_x (double sign, double xi, double et, double qq)
{
	return P1_y (sign, xi, et, qq);
}

double
P2_y (double sign, double xi, double et, double qq)
{
	return sign * yy * cc * (r + rx) * ir3x2
		- xi * sd * cd * (ir * ire2)
		- xi * yy * (r + re) * sd * ir3e2;
}

double
P2_z (double sign, double xi, double et, double qq)
{
	return sign * ((ir * irx) - cc * cc * (r + rx) * ir3x2)
		- sign * xi * sd2 * (ir * ire2) + xi * cc * (r + re) * sd * ir3e2;
}

double
P3_x (double sign, double xi, double et, double qq)
{
	return P1_z (sign, xi, et, qq);
}

double
P3_y (double sign, double xi, double et, double qq)
{
	return P2_z (sign, xi, et, qq);
}

double
P3_z (double sign, double xi, double et, double qq)
{
	return - sign * yy * cc * (r + rx) * ir3x2
		+ xi * sd * cd * (ir * ire2)
		- sign * xi * cc * (r + re) * cd * ir3e2;
}

/*****************************************************/
double
L1 (int flag, double sign, double xi, double et, double qq)
{
	return 0.; /*****/
	if (flag == DERIV_X) return L1_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return L1_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return L1_z (sign, xi, et, qq);
	return 0.;
}

double
L2 (int flag, double sign, double xi, double et, double qq)
{
	return 0.; /*****/
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
	return 0.; /*****/
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
	return 0.; /*****/
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
