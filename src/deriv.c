#include <stdio.h>
#include <math.h>

#include "constants.h"
#include "piez.h"

#define DERIV_X 0
#define DERIV_Y 1
#define DERIV_Z 2

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

double
K2_x (double sign, double xi, double et, double qq)
{
	double r_J2_x = J2_x (sign, xi, et, qq);
	return td * (- xi * yy * (ir * irc2) - sign * r_J2_x * td);
}

double
K2_y (double sign, double xi, double et, double qq)
{
	double r_J2_y = J2_y (sign, xi, et, qq);
	return td * (irc - yy * yy * (ir * irc2) - sign * r_J2_y * td);
}

double
K2_z (double sign, double xi, double et, double qq)
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
K4_x (double sign, double xi, double et, double qq)
{
	double r_K2_x = K2_x (sign, xi, et, qq);
	return cd * (r_K2_x - xi * sd * (ir * ire));
}

double
K4_y (double sign, double xi, double et, double qq)
{
	double r_K2_y = K2_y (sign, xi, et, qq);
	return cd * (r_K2_y - sd * (cd * ire + yy * (ir * ire)));
}

double
K4_z (double sign, double xi, double et, double qq)
{
	double r_K2_z = K2_z (sign, xi, et, qq);
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
K7_x (double sign, double xi, double et, double qq)
{
	double r_K4_x = K4_x (sign, xi, et, qq);
	return - r_K4_x * td;
}

double
K7_y (double sign, double xi, double et, double qq)
{
	double r_K4_y = K4_y (sign, xi, et, qq);
	return - r_K4_y * td;
}

double
K7_z (double sign, double xi, double et, double qq)
{
	double r_K4_z = K4_z (sign, xi, et, qq);
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

/********************************************************
 ****************** quad-poles ******************
 ********************************************************/
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

double
M1y_x (double sign, double xi, double et, double qq)
{
	return xi * xi * cd * (3.0 * r + et) * (ir3 * ire3)
		- yy * (r + re) * ir3e2
		+ xi * xi * yy * (8.0 * r2 + 9.0 * r * et + 3 * et * et) * ir5e3
		- cd * (ir * ire2);
}

double
M1y_y (double sign, double xi, double et, double qq)
{
	return - 2.0 * xi * sd2 * (ir * ire3)
		+ xi * (yy * cd + sign * cc * sd) * (3.0 * r + et) * (ir3 * ire3)
		+ xi * yy * yy * (8.0 * r2 + 9.0 * r * et + 3 * et * et) * ir5e3;
}

double
M1y_z (double sign, double xi, double et, double qq)
{
	/*
	double res = sign * 2.0 * xi * sd * cd * (ir * ire3)
		- xi * (cc * cd - sign * yy * sd) * (3.0 * r + et) * (ir3 * ire3)
		- xi * yy * cc * (8.0 * r2 + 9.0 * r * et + 3 * et * et) * ir5e3;
	*/
	double res =
		xi * ire3 * (sign * 2.0 * sd * cd * ir
		 - (cc * cd - sign * yy * sd) * (3.0 * r + et) * ir3
		 - yy * cc * (8.0 * r2 + 9.0 * r * et + 3 * et * et) * ir5);
	res = 0.0;
	return res;
}

double
M1z_x (double sign, double xi, double et, double qq)
{
	/*
	double res = sign * xi * xi * sd * (3.0 * r + et) * (ir3 * ire3)
		+ cc * (r + re) * ir3e2
		- xi * xi * cc * (8.0 * r2 + 9.0 * r * et + 3 * et * et) * ir5e3
		- sign * sd * (ir * ire2);
	*/
	double res =
		- (sign * sd - cc * (r + re) * ir * ir) * ir * ire2
		+ xi * xi * (sign * sd * (3.0 * r + et)
		 - cc * 8.0 - cc * 9.0 * et * ir
		 - cc * 3 * et * et * ir * ir) * ir3 * ire3;
	return res;
}

double
M1z_y (double sign, double xi, double et, double qq)
{
	return M1y_z (sign, xi, et, qq);
}

double
M1z_z (double sign, double xi, double et, double qq)
{
	/*
	double res = - 2.0 * xi * cd2 * (ir * ire3)
		- xi * (yy * cd + sign * cc * sd) * (3.0 * r + et) * (ir3 * ire3)
		+ xi * cc * cc * (8.0 * r2 + 9.0 * r * et + 3 * et * et) * ir5e3;
	*/
	double res =
		xi * ire3 * (- 2.0 * cd2 * ir
		 - (yy * cd + sign * cc * sd) * (3.0 * r + et) * ir3
		 + cc * cc * (8.0 * r2 + 9.0 * r * et + 3 * et * et) * ir5);
	res = 0.0;
	return res;
}

double
M2y_x (double sign, double xi, double et, double qq)
{
	return - 2.0 * xi * sd2 * (ir * ire3)
		+ xi * (yy * cd + sign * cc * sd) * (3.0 * r + et) * (ir3 * ire3)
		+ xi * yy * yy * (8.0 * r2 + 9.0 * r * et + 3 * et * et) * ir5e3;
}

double
M2y_y (double sign, double xi, double et, double qq)
{
	/*
	double res = - 2.0 * yy * sd2 * (ir * ire3) - cd * (ir * ire2)
		+ yy * (yy * cd + sign * cc * sd) * (3.0 * r + et) * (ir3 * ire3)
		- 2.0 * yy * (r + re) * (ir3 * ire2)
		+ yy * yy * yy * (8.0 * r2 + 9.0 * et * r + 3.0 * et * et) * ir5e3
		- 2.0 * sd2 * cd * ire3
		+ 2.0 * (yy * cd + sign * cc * sd) * cd * (ir * ire3)
		+ yy * yy * (3.0 * r + et) * cd * (ir3 * ire3);
	*/
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

double
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

double
M2z_x (double sign, double xi, double et, double qq)
{
	return M1z_y (sign, xi, et, qq);
}

double
M2z_y (double sign, double xi, double et, double qq)
{
	return M2y_z (sign, xi, et, qq);
}

double
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

double
M3y_x (double sign, double xi, double et, double qq)
{
	return M2z_x (sign, xi, et, qq);
}

double
M3y_y (double sign, double xi, double et, double qq)
{
	return M2z_y (sign, xi, et, qq);
}

double
M3y_z (double sign, double xi, double et, double qq)
{
	return M2z_z (sign, xi, et, qq);
}

double
M3z_x (double sign, double xi, double et, double qq)
{
	return M1z_z (sign, xi, et, qq);
}

double
M3z_y (double sign, double xi, double et, double qq)
{
	return M2z_z (sign, xi, et, qq);
}

double
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

double
N1z_x (double sign, double xi, double et, double qq)
{
	return 1.0 * ir3 - 3.0 * xi * xi * ir5;
}

double
N1z_y (double sign, double xi, double et, double qq)
{
	return - 3.0 * xi * yy * ir5;
}

double
N1z_z (double sign, double xi, double et, double qq)
{
	return 3.0 * xi * cc * ir5;
}

double
N2z_x (double sign, double xi, double et, double qq)
{
	return N1z_y (sign, xi, et, qq);
}

double
N2z_y (double sign, double xi, double et, double qq)
{
	return ir3 - 3.0 * yy * yy * ir5;
}

double
N2z_z (double sign, double xi, double et, double qq)
{
	return 3.0 * yy * cc * ir5;
}

double
O1z_x (double sign, double xi, double et, double qq)
{
	return -3.0 * xi * cc * ir5;
}

double
O1z_y (double sign, double xi, double et, double qq)
{
	return -3.0 * yy * cc * ir5;
}

double
O1z_z (double sign, double xi, double et, double qq)
{
	return - ir3 + 3.0 * cc * cc * ir5;
}

double
O2y_x (double sign, double xi, double et, double qq)
{
	return - ir3 + 3.0 * yy * yy * ir5;
}

double
O2y_y (double sign, double xi, double et, double qq)
{
	return - 3.0 * yy * (r + rx) * ir3x2
		+ yy * yy * yy * (8.0 * r2 + 9.0 * xi * r + 3.0 * xi * xi) * ir5x3;
}

double
O2y_z (double sign, double xi, double et, double qq)
{
	return cc * (r + rx) * ir3x2
		- yy * yy * cc * (8.0 * r2 + 9.0 * xi * r + 3.0 * xi * xi) * ir5x3;
}

double
O2z_x (double sign, double xi, double et, double qq)
{
	return - 3.0 * yy * cc * ir5;
}

double
O2z_y (double sign, double xi, double et, double qq)
{
	return O2y_z (sign, xi, et, qq);
}

double
O2z_z (double sign, double xi, double et, double qq)
{
	return - yy * (r + rx) * ir3x2
		+ yy * cc * cc * (8.0 * r2 + 9.0 * xi * r + 3.0 * xi * xi) * ir5x3;
}

double
O3z_x (double sign, double xi, double et, double qq)
{
	return - ir3 + 3.0 * cc * cc * ir5;
}

double
O3z_y (double sign, double xi, double et, double qq)
{
	return O2z_z (sign, xi, et, qq);
}

double
O3z_z (double sign, double xi, double et, double qq)
{
	return 3.0 * cc * (r + rx) * ir3x2
		- cc * cc * cc * (8.0 * r2 + 9.0 * xi * r + 3.0 * xi * xi) * ir5x3;
}

double
P1y_x (double sign, double xi, double et, double qq)
{
	return - sign * 3.0 * xi * cc * ir5
		- 3.0 * xi * (r + re) * sd * ir3e2
		+ xi * xi * xi * (8.0 * r2 + 9.0 * et * r + 3.0 * et * et) * sd * ir5e3;
}

double
P1y_y (double sign, double xi, double et, double qq)
{
	double res = - sign * 3.0 * yy * cc * ir5
		- yy * (r + re) * sd * ir3e2
		+ xi * xi * yy * (8.0 * r2 + 9.0 * et * r + 3.0 * et * et) * sd * ir5e3
		- sd * cd * (ir * ire2)
		+ xi * xi * (3.0 * r + et) * sd * cd * (ir3 * ire3);
	return res;
}

double
P1y_z (double sign, double xi, double et, double qq)
{
	double res = - sign * ir3 + sign * 3.0 * cc * cc * ir5
		+ cc * (r + re) * sd * ir3e2
		- xi * xi * cc * (8.0 * r2 + 9.0 * et * r + 3.0 * et * et) * sd * ir5e3
		- sign * sd2 * (ir * ire2)
		+ sign *	xi * xi * (3.0 * r + et) * sd2 * (ir3 * ire3);
	return res;
}

double
P1z_x (double sign, double xi, double et, double qq)
{
	double res =
		sign * (- 3.0 * xi * yy * ir5
			+ 3.0 * xi * (r + re) * cd * ir3e2
			- xi * xi * xi *
			(8.0 * r2 + 9.0 * et * r + 3.0 * et * et) * cd * ir5e3);
	return res;
}

double
P1z_y (double sign, double xi, double et, double qq)
{
	return P1y_z (sign, xi, et, qq);
}

double
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

double
P3y_x (double sign, double xi, double et, double qq)
{
	return P1y_z (sign, xi, et, qq);
}

double
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

double
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

double
P3z_x (double sign, double xi, double et, double qq)
{
	return P1z_z (sign, xi, et, qq);
}

double
P3z_y (double sign, double xi, double et, double qq)
{
	return P3y_z (sign, xi, et, qq);
}

double
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
	if (flag == DERIV_X) return K1_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return K1_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return K1_z (sign, xi, et, qq);
	return 0.;
}

double
K2 (int flag, double sign, double xi, double et, double qq)
{
	if (flag == DERIV_X) return K2_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return K2_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return K2_z (sign, xi, et, qq);
	return 0.;
}

double
K3 (int flag, double sign, double xi, double et, double qq)
{
	if (flag == DERIV_X) return K3_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return K3_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return K3_z (sign, xi, et, qq);
	return 0.;
}

double
K4 (int flag, double sign, double xi, double et, double qq)
{
	if (flag == DERIV_X) return K4_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return K4_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return K4_z (sign, xi, et, qq);
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
	if (flag == DERIV_X) return K7_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return K7_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return K7_z (sign, xi, et, qq);
	return 0.;
}

double
K8 (int flag, double sign, double xi, double et, double qq)
{
	if (flag == DERIV_X) return K8_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return K8_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return K8_z (sign, xi, et, qq);
	return 0.;
}

double
K9 (int flag, double sign, double xi, double et, double qq)
{
	if (flag == DERIV_X) return K9_x (sign, xi, et, qq);
	else if (flag == DERIV_Y) return K9_y (sign, xi, et, qq);
	else if (flag == DERIV_Z) return K9_z (sign, xi, et, qq);
	return 0.;
}

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
