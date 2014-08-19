#ifndef _PIEZOMAG_H_
#define _PIEZOMAG_H_

#include <stdbool.h>

#ifdef MIN
#undef MIN
#endif
#define MIN(x, y) ((x) <= (y)) ? (x) : (y)

#ifdef deg2rad
#undef deg2rad
#endif
#define deg2rad(a) ((a) * M_PI / 180.)

/*** magnetic component ***/
#ifdef  TOTAL_FORCE
#undef  TOTAL_FORCE
#endif
#define TOTAL_FORCE 0

#ifdef  X_COMP
#undef  X_COMP
#endif
#define X_COMP 1

#ifdef  Y_COMP
#undef  Y_COMP
#endif
#define Y_COMP 2

#ifdef  Z_COMP
#undef  Z_COMP
#endif
#define Z_COMP 3

/*** z-coordinate of observation point ***/
double	z_obs;

// specify component of output : X_COMP (1), Y_COMP (2), Z_COMP (3) or TOTAL_FORCE (0)
int		output_comp;

/*** flags ***/
bool	singular_R[4];	// R is singular?
bool	singular_RE[4];	// R + eta is singular?
bool	verbos;			// verbos mode

/*** fault parameters ***/
typedef struct s_fault_params	fault_params;

struct s_fault_params {

	// crustal parameters
	double	lambda;	// lame constants
	double	mu;
	double	alpha;		// = (lambda + mu) / (lambda + 2 * mu)

	// fault dimension
	double	flength1;
	double	flength2;
	double	fwidth1;
	double	fwidth2;

	double	fdepth;	// depth of the origin of fault coordinate
	double	fstrike;	// strike angle
	double	fdip;		// dip angle

	// dislocation vector
	double	u1;	// strike slip
	double	u2;	// dip slip
	double	u3;	// tensile opening

};

/*** crustal and magnetic parameters ***/
typedef struct s_magnetic_params	magnetic_params;

struct s_magnetic_params {

	// magnetic parameters
	double	beta;		// seismo-magnetic sensitivity

	double	mgz_int;	// intensity of magnetization of the uniformly magnetized crust
	double	mgz_inc;	// inclination
	double	mgz_dec;	// declination
	// external field
	double	exf_inc;	// inclination of external field
	double	exf_dec;	// declination

	// depth of Curier point isotherm
	double	dcurier;

	// seismo-magnetic moment vector
	double	cx;
	double	cy;
	double	cz;
	double	c0; // = sqrt (cx^2 + cy^2 + cz^2)

};

/********* util.c *********/
bool	fread_params (FILE *fp, fault_params *fault, magnetic_params *mag);

/********* common.c *********/
double	piezomagnetic_effect (int component, const fault_params *fault, const magnetic_params *mag, double xobs, double yobs, double zobs);
void	fprintf_piezomagnetic_effect (FILE *stream, int component, const fault_params *fault, const magnetic_params *mag,
			double xobs1, double xobs2, double dx, double yobs1, double yobs2, double dy, double zobs);

#endif	// _PIEZOMAG_H_
