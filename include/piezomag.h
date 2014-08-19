#ifndef _PIEZOMAG_H_
#define _PIEZOMAG_H_

#include <stdbool.h>

/**** utilities ***/
#ifdef MIN
#undef MIN
#endif
#define MIN(x, y) ((x) <= (y)) ? (x) : (y)

/* degree -> radian */
#ifdef deg2rad
#undef deg2rad
#endif
#define deg2rad(a) ((a) * M_PI / 180.)

/*** definition of magnetic components ***/
#ifdef  TOTAL_FORCE
#undef  TOTAL_FORCE
#endif
#define TOTAL_FORCE 0	// total force

#ifdef  X_COMP
#undef  X_COMP
#endif
#define X_COMP 1	// x component

#ifdef  Y_COMP
#undef  Y_COMP
#endif
#define Y_COMP 2	// y component

#ifdef  Z_COMP
#undef  Z_COMP
#endif
#define Z_COMP 3	// z component


/*c********************
 *c  global variables
 *c********************/

/* z-coordinate of observation point ***/
double	z_obs;

/* specify component of output : X_COMP (1), Y_COMP (2), Z_COMP (3) or TOTAL_FORCE (0) */
int		output_comp;

/* verbos mide */
bool	verbos;


/*c********************
 *c    structures
 *c********************/

/* fault parameters */
typedef struct s_fault_params	fault_params;

struct s_fault_params {

	// crustal parameters
	double	lambda;	// lame constants
	double	mu;
	double	alpha;		// = (lambda + mu) / (lambda + 2 * mu)

	// fault geometry
	double	flength1;	// length
	double	flength2;
	double	fwidth1;	// width
	double	fwidth2;

	double	fdepth;	// depth of the origin of fault coordinate
	double	fstrike;	// strike angle
	double	fdip;		// dip angle

	// dislocation vector
	double	u1;	// strike slip
	double	u2;	// dip slip
	double	u3;	// tensile opening

};

/* magnetic parameters */
typedef struct s_magnetic_params	magnetic_params;

struct s_magnetic_params {

	double	beta;		// seismo-magnetic sensitivity

	double	mgz_int;	// intensity of magnetization of the uniformly magnetized crust
	double	mgz_inc;	// inclination
	double	mgz_dec;	// declination
	// external field
	double	exf_inc;	// inclination of external field
	double	exf_dec;	// declination

	double	dcurier;	// depth of Curier point isotherm

	// seismomagnetic moment vector
	double	cx;
	double	cy;
	double	cz;
	double	c0; // = sqrt (cx^2 + cy^2 + cz^2)

};


/*c********************
 *c  public functions
 *c********************/

/** util.c **/

/* read parameters from file and store them to global variables
   and pre-allocated structures: fault_params *fault, magnetic_params *mag. */
bool	fread_params (FILE *fp, fault_params *fault, magnetic_params *mag);


/** piezomag.c **/

/* calculate seismomagnetic field on observation point (xobs, yobs, zobs) */
double	seismomagnetic_effect (int component, const fault_params *fault, const magnetic_params *mag, double xobs, double yobs, double zobs);

/* calculate seismomagnetic field on grid x=[xobs1:dx:xobs2], y=[yobs1:dy:yobs2], z=zobs */
void	fprintf_seismomagnetic_effect (FILE *stream, int component, const fault_params *fault, const magnetic_params *mag,
			double xobs1, double xobs2, double dx, double yobs1, double yobs2, double dy, double zobs);

#endif	// _PIEZOMAG_H_
