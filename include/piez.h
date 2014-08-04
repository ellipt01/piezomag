#ifndef _PIEZ_H_
#define _PIEZ_H_

#include <stdbool.h>

#ifdef  X_COMP
#undef  X_COMP
#endif
#define X_COMP 0

#ifdef  Y_COMP
#undef  Y_COMP
#endif
#define Y_COMP 1

#ifdef  Z_COMP
#undef  Z_COMP
#endif
#define Z_COMP 2

#ifdef  TOTAL_FORCE
#undef  TOTAL_FORCE
#endif
#define TOTAL_FORCE 3

#ifdef deg2rad
#undef deg2rad
#endif
#define deg2rad(a) ((a) * M_PI / 180.)

bool	singular_R[4];
bool	singular_RE[4];

bool	verbos;

/********* util.c *********/
void	usage (char *toolname);
bool	initialize (int argc, char **argv);
bool	fread_params (FILE *fp);
void	fwrite_params (FILE *stream);
void	coordinates_transform (double theta, double *x, double *y);
void	clear_singular_flag (int i);
void	clear_all_singular_flag (void);
void	set_singular_flag (int i);
bool	is_singular_point (bool *flag);
void	check_singular_point (double x, double y, double eps);


/********* common.c *********/
bool	set_constants (void);
void	set_geometry_variables (double sign, double xi, double et, double qq);
double	piezomagnetic_effect_component (int component,
				       double u1, double u2, double u3,
				       double x, double y, double z);
double	piezomagnetic_effect (int component,
			     double u1, double u2, double u3,
			     double x, double y, double z);
void	fprintf_piezomagnetic_effect (FILE *stream, int component,
				     double x1, double x2, double dx,
				     double y1, double y2, double dy, double z);

/********* dipoles *********/
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
/********* quad-poles *********/
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
/********* oct-poles *********/
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

double strike_slip (int flag, double x, double y, double z);
double dip_slip (int flag, double x, double y, double z);
double tensile_opening (int flag, double x, double y, double z);

#endif	// _PIEZ_H_
