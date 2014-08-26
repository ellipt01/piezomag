/*
 * octpoles.c
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

/*c***************************
 *c******** octpoles *********
 *c***************************/
static double
M1y_x (double sign, double xi, double et, double qq)
{
	return xi * xi * cd * (3.0 * r + et) * (ir3 * ire3)
		- yy * (r + re) * ir3e2
		+ xi * xi * yy * (8.0 * r2 + 9.0 * r * et + 3 * et * et) * ir5e3
		- cd * (ir * ire2);
}

static double
M1y_y (double sign, double xi, double et, double qq)
{
	return - 2.0 * xi * sd2 * (ir * ire3)
		+ xi * (yy * cd + sign * cc * sd) * (3.0 * r + et) * (ir3 * ire3)
		+ xi * yy * yy * (8.0 * r2 + 9.0 * r * et + 3 * et * et) * ir5e3;
}

static double
M1y_z (double sign, double xi, double et, double qq)
{
	double res =
		xi * ire3 * (sign * 2.0 * sd * cd * ir
		 - (cc * cd - sign * yy * sd) * (3.0 * r + et) * ir3
		 - yy * cc * (8.0 * r2 + 9.0 * r * et + 3 * et * et) * ir5);
	res = 0.0;
	return res;
}

static double
M1z_x (double sign, double xi, double et, double qq)
{
	double res =
		- (sign * sd - cc * (r + re) * ir * ir) * ir * ire2
		+ xi * xi * (sign * sd * (3.0 * r + et)
		 - cc * 8.0 - cc * 9.0 * et * ir
		 - cc * 3 * et * et * ir * ir) * ir3 * ire3;
	return res;
}

static double
M1z_y (double sign, double xi, double et, double qq)
{
	return M1y_z (sign, xi, et, qq);
}

static double
M1z_z (double sign, double xi, double et, double qq)
{
	double res =
		xi * ire3 * (- 2.0 * cd2 * ir
		 - (yy * cd + sign * cc * sd) * (3.0 * r + et) * ir3
		 + cc * cc * (8.0 * r2 + 9.0 * r * et + 3 * et * et) * ir5);
	res = 0.0;
	return res;
}

static double
M2y_x (double sign, double xi, double et, double qq)
{
	return - 2.0 * xi * sd2 * (ir * ire3)
		+ xi * (yy * cd + sign * cc * sd) * (3.0 * r + et) * (ir3 * ire3)
		+ xi * yy * yy * (8.0 * r2 + 9.0 * r * et + 3 * et * et) * ir5e3;
}

static double
M2y_y (double sign, double xi, double et, double qq)
{
	double res =
		- ire2 * (cd + 2.0 * yy * (r + re) * ir * ir) * ir
		+ ire3 * (- 2.0 * sd2 * cd
				+ 2.0 * ir * (yy * c2d + sign * cc * sd * cd)
				+ yy * ir3 * ( sign * cc * sd * (3.0 * r + et)
					+ 2.0 * yy * cd * (3.0 * r + et)
					+ 8.0 * yy * yy)
				+ yy * yy * yy * et * (9.0 * r + 3.0 * et) * ir5);
	return res;
}

static double
M2y_z (double sign, double xi, double et, double qq)
{
	double res = 2.0 * cc * sd2 * (ir * ire3) + sign * sd * (ir * ire2)
		- cc * (yy * cd + sign * cc * sd) * (3.0 * r + et) * (ir3 * ire3)
		- cc * yy * yy * (8.0 * r2 + 9.0 * et * r + 3.0 * et * et) * ir5e3
		- sign * 2.0 * sd * sd2 * ire3
		+ sign * 2.0 * (yy * cd + sign * cc * sd) * sd * (ir * ire3)
		+ sign * yy * yy * (3.0 * r + et) * sd * (ir3 * ire3);
	return res;
}

static double
M2z_x (double sign, double xi, double et, double qq)
{
	return M1z_y (sign, xi, et, qq);
}

static double
M2z_y (double sign, double xi, double et, double qq)
{
	return M2y_z (sign, xi, et, qq);
}

static double
M2z_z (double sign, double xi, double et, double qq)
{
	double res = - 2.0 * yy * cd2 * (ir * ire3) + cd * (ir * ire2)
		- yy * (yy * cd + sign * cc * sd) * (3.0 * r + et) * (ir3 * ire3)
		+ yy * cc * cc * (8.0 * r2 + 9.0 * et * r + 3.0 * et * et) * ir5e3
		- 2.0 * cd2 * cd * ire3
		- 2.0 * (yy * cd + sign * cc * sd) * cd * (ir * ire3)
		+ cc * cc * (3.0 * r + et) * cd * (ir3 * ire3);
	return res;
}

static double
M3y_x (double sign, double xi, double et, double qq)
{
	return M2z_x (sign, xi, et, qq);
}

static double
M3y_y (double sign, double xi, double et, double qq)
{
	return M2z_y (sign, xi, et, qq);
}

static double
M3y_z (double sign, double xi, double et, double qq)
{
	return M2z_z (sign, xi, et, qq);
}

static double
M3z_x (double sign, double xi, double et, double qq)
{
	return M1z_z (sign, xi, et, qq);
}

static double
M3z_y (double sign, double xi, double et, double qq)
{
	return M2z_z (sign, xi, et, qq);
}

static double
M3z_z (double sign, double xi, double et, double qq)
{
	double res = 2.0 * cc * cd2 * (ir * ire3)
		- sign * sd * (ir * ire2)
		+ cc * (yy * cd + sign * cc * sd) * (3.0 * r + et) * (ir3 * ire3)
		+ 2.0 * cc * (r + re) * ir3e2
		- cc * cc * cc * (8.0 * r2 + 9.0 * r * et + 3 * et * et) * ir5e3
		- sign * 2.0 * sd * cd2 * ire3
		- sign * 2.0 * (yy * cd + sign * cc * sd) * sd * (ir * ire3)
		+ sign * cc * cc * (3.0 * r + et) * sd * (ir3 * ire3);
	return res;
}

static double
N1z_x (double sign, double xi, double et, double qq)
{
	return 1.0 * ir3 - 3.0 * xi * xi * ir5;
}

static double
N1z_y (double sign, double xi, double et, double qq)
{
	return - 3.0 * xi * yy * ir5;
}

static double
N1z_z (double sign, double xi, double et, double qq)
{
	return 3.0 * xi * cc * ir5;
}

static double
N2z_x (double sign, double xi, double et, double qq)
{
	return N1z_y (sign, xi, et, qq);
}

static double
N2z_y (double sign, double xi, double et, double qq)
{
	return ir3 - 3.0 * yy * yy * ir5;
}

static double
N2z_z (double sign, double xi, double et, double qq)
{
	return 3.0 * yy * cc * ir5;
}

static double
O1z_x (double sign, double xi, double et, double qq)
{
	return -3.0 * xi * cc * ir5;
}

static double
O1z_y (double sign, double xi, double et, double qq)
{
	return -3.0 * yy * cc * ir5;
}

static double
O1z_z (double sign, double xi, double et, double qq)
{
	return - ir3 + 3.0 * cc * cc * ir5;
}

static double
O2y_x (double sign, double xi, double et, double qq)
{
	return - ir3 + 3.0 * yy * yy * ir5;
}

static double
O2y_y (double sign, double xi, double et, double qq)
{
	return - 3.0 * yy * (r + rx) * ir3x2
		+ yy * yy * yy * (8.0 * r2 + 9.0 * xi * r + 3.0 * xi * xi) * ir5x3;
}

static double
O2y_z (double sign, double xi, double et, double qq)
{
	return cc * (r + rx) * ir3x2
		- yy * yy * cc * (8.0 * r2 + 9.0 * xi * r + 3.0 * xi * xi) * ir5x3;
}

static double
O2z_x (double sign, double xi, double et, double qq)
{
	return - 3.0 * yy * cc * ir5;
}

static double
O2z_y (double sign, double xi, double et, double qq)
{
	return O2y_z (sign, xi, et, qq);
}

static double
O2z_z (double sign, double xi, double et, double qq)
{
	return - yy * (r + rx) * ir3x2
		+ yy * cc * cc * (8.0 * r2 + 9.0 * xi * r + 3.0 * xi * xi) * ir5x3;
}

static double
O3z_x (double sign, double xi, double et, double qq)
{
	return - ir3 + 3.0 * cc * cc * ir5;
}

static double
O3z_y (double sign, double xi, double et, double qq)
{
	return O2z_z (sign, xi, et, qq);
}

static double
O3z_z (double sign, double xi, double et, double qq)
{
	return 3.0 * cc * (r + rx) * ir3x2
		- cc * cc * cc * (8.0 * r2 + 9.0 * xi * r + 3.0 * xi * xi) * ir5x3;
}

static double
P1y_x (double sign, double xi, double et, double qq)
{
	return - sign * 3.0 * xi * cc * ir5
		- 3.0 * xi * (r + re) * sd * ir3e2
		+ xi * xi * xi * (8.0 * r2 + 9.0 * et * r + 3.0 * et * et) * sd * ir5e3;
}

static double
P1y_y (double sign, double xi, double et, double qq)
{
	double res = - sign * 3.0 * yy * cc * ir5
		- yy * (r + re) * sd * ir3e2
		+ xi * xi * yy * (8.0 * r2 + 9.0 * et * r + 3.0 * et * et) * sd * ir5e3
		- sd * cd * (ir * ire2)
		+ xi * xi * (3.0 * r + et) * sd * cd * (ir3 * ire3);
	return res;
}

static double
P1y_z (double sign, double xi, double et, double qq)
{
	double res = - sign * ir3 + sign * 3.0 * cc * cc * ir5
		+ cc * (r + re) * sd * ir3e2
		- xi * xi * cc * (8.0 * r2 + 9.0 * et * r + 3.0 * et * et) * sd * ir5e3
		- sign * sd2 * (ir * ire2)
		+ sign *	xi * xi * (3.0 * r + et) * sd2 * (ir3 * ire3);
	return res;
}

static double
P1z_x (double sign, double xi, double et, double qq)
{
	double res =
		sign * (- 3.0 * xi * yy * ir5
			+ 3.0 * xi * (r + re) * cd * ir3e2
			- xi * xi * xi *
			(8.0 * r2 + 9.0 * et * r + 3.0 * et * et) * cd * ir5e3);
	return res;
}

static double
P1z_y (double sign, double xi, double et, double qq)
{
	return P1y_z (sign, xi, et, qq);
}

static double
P1z_z (double sign, double xi, double et, double qq)
{
	double res = sign * 3.0 * yy * cc * ir5
		- sign * cc * (r + re) * cd * ir3 * ire2
		+ sign * xi * xi * cc *
		(8.0 * r2 + 9.0 * et * r + 3.0 * et * et) * cd * ir5e3
		+ sd * cd * ir * ire2
		- xi * xi * (3.0 * r + et) * sd * cd * ir3 * ire3;
	return res;
}

static double
P3y_x (double sign, double xi, double et, double qq)
{
	return P1y_z (sign, xi, et, qq);
}

static double
P3y_y (double sign, double xi, double et, double qq)
{
	double res = - sign * yy * (r + rx) * ir3x2
		+ sign * yy * cc * cc *
		(8.0 * r2 + 9.0 * xi * r + 3.0 * xi * xi) * ir5x3
		+ sign * xi * yy * (3.0 * r + et) * sd2 * (ir3 * ire3)
		- xi * yy * cc *
		(8.0 * r2 + 9.0 * et * r + 3.0 * et * et) * sd * ir5e3
		+ sign * 2.0 * xi * sd2 * cd * (ir * ire3)
		- xi * cc * (3.0 * r + et) * sd * cd * (ir3 * ire3);
	return res;
}

static double
P3y_z (double sign, double xi, double et, double qq)
{
	double res = sign * 3.0 * cc * (r + rx) * ir3x2
		- sign * cc * cc * cc *
		(8.0 * r2 + 9.0 * xi * r + 3.0 * xi * xi) * ir5x3
		- sign * xi * cc * (3.0 * r + et) * sd2 * (ir3 * ire3)
		- xi * (r + re) * sd * ir3e2
		+ xi * cc * cc *
		(8.0 * r2 + 9.0 * et * r + 3.0 * et * et) * sd * ir5e3
		+ 2.0 * xi * sd2 * sd * (ir * ire3)
		- sign * xi * cc * (3.0 * r + et) * sd2 * (ir3 * ire3);
	return res;
}

static double
P3z_x (double sign, double xi, double et, double qq)
{
	return P1z_z (sign, xi, et, qq);
}

static double
P3z_y (double sign, double xi, double et, double qq)
{
	return P3y_z (sign, xi, et, qq);
}

static double
P3z_z (double sign, double xi, double et, double qq)
{
	double res = sign * yy * (r + rx) * ir3x2
		- sign * yy * cc * cc *
		(8.0 * r2 + 9.0 * xi * r + 3.0 * xi * xi) * ir5x3
		+ xi * cc * (3.0 * r + et) * sd * cd * (ir3 * ire3)
		+ sign * xi * (r + re) * cd * ir3e2
		- sign * xi * cc * cc *
		(8.0 * r2 + 9.0 * et * r + 3.0 * et * et) * cd * ir5e3
		- sign * 2.0 * xi * sd2 * cd * (ir * ire3)
		+ xi * cc * (3.0 * r + et) * sd * cd * (ir3 * ire3);
	return res;
}

/*****************************************************/
double
M1y (int flag, double sign, double xi, double et, double qq)
{
	if (singular_RE[0]) return 0.;
	if (flag == DERIV_X) return M1y_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return M1y_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return M1y_z (sign, xi, et, qq);
	return 0.;
}

double
M1z (int flag, double sign, double xi, double et, double qq)
{
	if (singular_RE[0]) return 0.;
	if (flag == DERIV_X) return M1z_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return M1z_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return M1z_z (sign, xi, et, qq);
	return 0.;
}

double
M2y (int flag, double sign, double xi, double et, double qq)
{
	if (singular_RE[0]) return 0.;
	if (flag == DERIV_X) return M2y_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return M2y_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return M2y_z (sign, xi, et, qq);
	return 0.;
}

double
M2z (int flag, double sign, double xi, double et, double qq)
{
	if (singular_RE[0]) return 0.;
	if (flag == DERIV_X) return M2z_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return M2z_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return M2z_z (sign, xi, et, qq);
	return 0.;
}

double
M3y (int flag, double sign, double xi, double et, double qq)
{
	if (singular_RE[0]) return 0.;
	if (flag == DERIV_X) return M3y_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return M3y_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return M3y_z (sign, xi, et, qq);
	return 0.;
}

double
M3z (int flag, double sign, double xi, double et, double qq)
{
	if (singular_RE[0]) return 0.;
	if (flag == DERIV_X) return M3z_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return M3z_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return M3z_z (sign, xi, et, qq);
	return 0.;
}

double
N1z (int flag, double sign, double xi, double et, double qq)
{
	if (flag == DERIV_X) return N1z_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return N1z_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return N1z_z (sign, xi, et, qq);
	return 0.;
}

double
N2z (int flag, double sign, double xi, double et, double qq)
{
	if (flag == DERIV_X) return N2z_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return N2z_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return N2z_z (sign, xi, et, qq);
	return 0.;
}

double
O1z (int flag, double sign, double xi, double et, double qq)
{
	if (flag == DERIV_X) return O1z_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return O1z_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return O1z_z (sign, xi, et, qq);
	return 0.;
}

double
O2y (int flag, double sign, double xi, double et, double qq)
{
	if (flag == DERIV_X) return O2y_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return O2y_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return O2y_z (sign, xi, et, qq);
	return 0.;
}

double
O2z (int flag, double sign, double xi, double et, double qq)
{
	if (flag == DERIV_X) return O2z_x (sign, xi, et, qq);
	if (flag == DERIV_Y) return O2z_y (sign, xi, et, qq);
	if (flag == DERIV_Z) return O2z_z (sign, xi, et, qq);
	return 0.;
}

double
O3z (int flag, double sign, double xi, double et, double qq)
{
	if (flag == DERIV_X) return O3z_x (sign, xi, et, qq);
	if (flag == DERIV_Y) return O3z_y (sign, xi, et, qq);
	if (flag == DERIV_Z) return O3z_z (sign, xi, et, qq);
	return 0.;
}

double
P1y (int flag, double sign, double xi, double et, double qq)
{
	if (flag == DERIV_X) return P1y_x (sign, xi, et, qq);
	if (flag == DERIV_Y) return P1y_y (sign, xi, et, qq);
	if (flag == DERIV_Z) return P1y_z (sign, xi, et, qq);
	return 0.;
}

double
P1z (int flag, double sign, double xi, double et, double qq)
{
	if (flag == DERIV_X) return P1z_x (sign, xi, et, qq);
	if (flag == DERIV_Y) return P1z_y (sign, xi, et, qq);
	if (flag == DERIV_Z) return P1z_z (sign, xi, et, qq);
	return 0.;
}

double
P3y (int flag, double sign, double xi, double et, double qq)
{
	if (flag == DERIV_X) return P3y_x (sign, xi, et, qq);
	if (flag == DERIV_Y) return P3y_y (sign, xi, et, qq);
	if (flag == DERIV_Z) return P3y_z (sign, xi, et, qq);
	return 0.;
}

double
P3z (int flag, double sign, double xi, double et, double qq)
{
	if (flag == DERIV_X) return P3z_x (sign, xi, et, qq);
	if (flag == DERIV_Y) return P3z_y (sign, xi, et, qq);
	if (flag == DERIV_Z) return P3z_z (sign, xi, et, qq);
	return 0.;
}

