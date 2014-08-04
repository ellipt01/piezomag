#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

// z-coordinate of obs. point
double	z_obs;

double	lambda;
double	mu;
double	beta;

/*** fault parameters ***/
// fault dimension
double	flength1;
double	flength2;
double	fwidth1;
double	fwidth2;

// angle
double	fdepth;
double	fstrike;
double	fdip;

/*** magnetic field parameter ***/
// intensity, inclination and declination
// magnetization
double	mgz_int;
double	mgz_inc;
double	mgz_dec;
// external field
double	exf_inc;
double	exf_dec;

// dislocation vector
double	u1;
double	u2;
double	u3;

// depth of Curier point isotherm
double	dcurier;

double	x01;
double	x02;
double	dx;
double	y01;
double	y02;
double	dy;
int		output_comp;

double	cx;
double	cy;
double	cz;

// sin(delta), cos(delta) and tan(delta)
double	sd;
double	cd;
double	td;
// sec(delta), sin^2(delta), cos^2(delta), sin(2*delta) and cos(2*delta)
double	secd;
double	sd2;
double	cd2;
double	s2d;
double	c2d;

double	d[4];

double	alpha;
double	alpha0;
double	alpha1;
double	alpha2;
double	alpha3;
double	alpha4;
double	alpha5;
double	alpha6;

double	yy;
double	cc;

// r, r^2
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

#endif // _CONSTANTS_H_

