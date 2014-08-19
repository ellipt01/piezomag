/*
 * private.h
 *
 *  Created on: 2014/08/08
 *      Author: utsugi
 */

#ifndef _PRIVATE_H_
#define _PRIVATE_H_

/*** constants and functions which are used internally ***/

double	sd;		// sin (delta)
double	cd;		// cos (delta)
double	td;		// tan (delta)
double	secd;	// sec (delta)
double	sd2;	// sin^2 (delta)
double	cd2;	// cos^2 (delta)
double	s2d;	// sin (2 * delta)
double	c2d;	// cos (2 * delta)

// d1, d2, d3
// d0 = dummy (not used)
// d1 = fdepth - z_obs
// d2 = fdepth - 2 * dcurier + z_obs
// d3 = fdepth + 2 * dcurier - z_obs
double	d[4];

// alpha;
double	alpha0;	// 4 * alpha - 1
double	alpha1;	// 3 * alpha / alpha0
double	alpha2;	// 6 * alpha^2 / alpha0
double	alpha3;	// 2 * alpha * (1 - alpha) / alpha0
double	alpha4;	// alpha * (2 * alpha + 1) / alpha0
double	alpha5;	// alpha * (2 * alpha - 5) / alpha0
double	alpha6;	// 3 * alpha * (1 - 2 * alpha) / alpha0

double	yy;	// yy = eta * cos(delta) + qq * sin(delta)
double	cc;	// cc = sign * (qq * sin(delta) - eta * sin(delta))

// r, r^2: r = sqrt (xi^2 + eta^2 + q^2)
double	r;
double	r2;

// r + xi, r + eta and r + c
double	rx;
double	re;
double	rc;

// 1 / r, 1 / r^3 and 1 / r^5
double	ir;
double	ir3;
double	ir5;

// 1 / (r + xi), 1 / (r + eta) and 1 / (r + c)
double	irx;
double	ire;
double	irc;

// 1 / (r + xi)^2, 1 / (r + eta)^2 and 1 / (r + c)^2
double	irx2;
double	ire2;
double	irc2;

// 1 / (r + xi)^3, 1 / (r + eta)^3 and 1 / (r + c)^3
double	irx3;
double	ire3;
double	irc3;

// 1 / r^3 * (r + xi)^2, 1 / r^3 * (r + eta)^2 and 1 / r^3 * (r + c)^2
double	ir3x2;
double	ir3e2;
double	ir3c2;

// 1 / r^5 * (r + xi)^3, 1 / r^5 * (r + eta)^3 and 1 / r^5 * (r + c)^3
double	ir5x3;
double	ir5e3;
double	ir5c3;

/********* deriv.c *********/
/***** dipole terms *****/
double log_rx (int flag, double sign, double xi, double et, double qq);
double log_re (int flag, double sign, double xi, double et, double qq);
double log_rc (int flag, double sign, double xi, double et, double qq);
double atan_xe_qr (int flag, double sign, double xi, double et, double qq);
double J1 (int flag, double sign, double xi, double et, double qq);
double J2 (int flag, double sign, double xi, double et, double qq);
double K1 (int flag, double sign, double xi, double et, double qq);
double K2 (int flag, double sign, double xi, double et, double qq);
double K3 (int flag, double sign, double xi, double et, double qq);
double K4 (int flag, double sign, double xi, double et, double qq);
double K5 (int flag, double sign, double xi, double et, double qq);
double K6 (int flag, double sign, double xi, double et, double qq);
double K7 (int flag, double sign, double xi, double et, double qq);
double K8 (int flag, double sign, double xi, double et, double qq);
double K9 (int flag, double sign, double xi, double et, double qq);
/***** quad-pole terms *****/
double L1 (int flag, double sign, double xi, double et, double qq);
double L2 (int flag, double sign, double xi, double et, double qq);
double M1 (int flag, double sign, double xi, double et, double qq);
double M2 (int flag, double sign, double xi, double et, double qq);
double M3 (int flag, double sign, double xi, double et, double qq);
double N1 (int flag, double sign, double xi, double et, double qq);
double N2 (int flag, double sign, double xi, double et, double qq);
double O1 (int flag, double sign, double xi, double et, double qq);
double O2 (int flag, double sign, double xi, double et, double qq);
double O3 (int flag, double sign, double xi, double et, double qq);
double P1 (int flag, double sign, double xi, double et, double qq);
double P2 (int flag, double sign, double xi, double et, double qq);
double P3 (int flag, double sign, double xi, double et, double qq);
/***** oct-pole terms *****/
double M1y (int flag, double sign, double xi, double et, double qq);
double M1z (int flag, double sign, double xi, double et, double qq);
double M2y (int flag, double sign, double xi, double et, double qq);
double M2z (int flag, double sign, double xi, double et, double qq);
double M3y (int flag, double sign, double xi, double et, double qq);
double M3z (int flag, double sign, double xi, double et, double qq);
double N1z (int flag, double sign, double xi, double et, double qq);
double N2z (int flag, double sign, double xi, double et, double qq);
double O1z (int flag, double sign, double xi, double et, double qq);
double O2y (int flag, double sign, double xi, double et, double qq);
double O2z (int flag, double sign, double xi, double et, double qq);
double O3z (int flag, double sign, double xi, double et, double qq);
double P1y (int flag, double sign, double xi, double et, double qq);
double P1z (int flag, double sign, double xi, double et, double qq);
double P3y (int flag, double sign, double xi, double et, double qq);
double P3z (int flag, double sign, double xi, double et, double qq);

void	set_geometry_variables (double sign, double xi, double et, double qq);

/********* utils.c *********/
double	total_force (double hx, double hy, double hz, double exf_inc, double exf_dec);
void	rotate (double theta, double *x, double *y);
//void	clear_singular_flag (int i);
void	clear_all_singular_flags (void);
void	set_singular_flag (int i);
bool	is_singular_point (bool *flag);
void	check_singular_point (const fault_params *fault, double x, double y, double eps);

/********* strike.c *********/
double strike_slip (int flag, const fault_params *fault, const magnetic_params *mag, double x, double y, double z);
/********* dip.c *********/
double dip_slip (int flag, const fault_params *fault, const magnetic_params *mag, double x, double y, double z);
/********* tensile.c *********/
double tensile_opening (int flag, const fault_params *fault, const magnetic_params *mag, double x, double y, double z);

#endif /* _PRIVATE_H_ */
