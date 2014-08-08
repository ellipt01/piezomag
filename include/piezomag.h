#ifndef _PIEZOMAG_H_
#define _PIEZOMAG_H_

#include <stdbool.h>

#ifdef MIN
#undef MIN
#endif
#define MIN(x, y) ((x) <= (y)) ? (x) : (y)

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

/*** range for calculating magnetic field ***/
double	z_obs;		// z-coordinate of obs. point
double	x_west;	// left end of x coordinate (EW) in the range
double	x_east;	// right edge in range
double	dx;			// x interval
double	y_south;	// lower end of y coordinate (NS) in the range
double	y_north;	// top of the range
double	dy;			// y interval

// specify component of output : X_COMP (0), Y_COMP (1), Z_COMP (2) or TOTAL_FORCE (3)
int		output_comp;

/*** flags ***/
bool	singular_R[4];	// is R is singular
bool	singular_RE[4];	// is R + eta is singular
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
	double	c0;
	double	cx;
	double	cy;
	double	cz;
};

/********* util.c *********/
bool	initialize (int argc, char **argv, fault_params **fault, magnetic_params **mag);
void	coordinates_transform (double theta, double *x, double *y);

void	clear_singular_flag (int i);
void	clear_all_singular_flag (void);
void	set_singular_flag (int i);
bool	is_singular_point (bool *flag);
void	check_singular_point (const fault_params *fault, double x, double y, double eps);

/********* common.c *********/
//bool	set_constants (fault_params *fault, magnetic_params *mag);
void	set_geometry_variables (double sign, double xi, double et, double qq);
double	piezomagnetic_effect_component (int component, const fault_params *fault, const magnetic_params *mag, double x, double y, double z);
double	piezomagnetic_effect (int component, const fault_params *fault, const magnetic_params *mag, double x, double y, double z);
void	fprintf_piezomagnetic_effect (FILE *stream, int component, const fault_params *fault, const magnetic_params *mag,
			double x1, double x2, double dx, double y1, double y2, double dy, double z);

double strike_slip (int flag, const fault_params *fault, const magnetic_params *mag, double x, double y, double z);
double dip_slip (int flag, const fault_params *fault, const magnetic_params *mag, double x, double y, double z);
double tensile_opening (int flag, const fault_params *fault, const magnetic_params *mag, double x, double y, double z);

#endif	// _PIEZOMAG_H_
