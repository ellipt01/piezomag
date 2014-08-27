#ifndef _PIEZOMAG_H_
#define _PIEZOMAG_H_

#include <stdbool.h>

/*c**********
 * utilities
 *c**********/

/* degree -> radian */
#define _deg2rad_(a) ((a) * M_PI / 180.)


/*c*******
 *  enums
 *c*******/

/* definition of magnetic component */
typedef enum {
	MAG_COMP_NONE = -1,
	MAG_COMP_F    =  0,	// total force
	MAG_COMP_X    =  1,	// x component
	MAG_COMP_Y    =  2,	// y component
	MAG_COMP_Z    =  3,	// z component
	MAG_COMP_NUM_ITEMS = 4
} MagComponent;

/* definition of seismomagnetic term */
typedef enum {
	SEISMO_MAG_NONE,
	SEISMO_MAG_MAIN       =  1 << 0,	// main term (0)
	SEISMO_MAG_MIRROR     =  1 << 1,	// mirror image (H0)
	SEISMO_MAG_SUBMIRROR  =  1 << 2,	// sub-mirror image (HI, HIII or HII)
	// total seismomagnetic field (0 + H0 + (HI, HIII or HII))
	SEISMO_MAG_TOTAL =  SEISMO_MAG_MAIN | SEISMO_MAG_MIRROR | SEISMO_MAG_SUBMIRROR,
	SEISMO_MAG_NUM_ITEMS
} SeismoMagTerm;


/*c********************
 *c  global variables
 *c********************/

/* z-coordinate of observation point ***/
double			z_obs;

/* specify component of output : X_COMP (1), Y_COMP (2), Z_COMP (3) or TOTAL_FORCE (0) */
MagComponent	output_comp;

/* verbos mide */
bool			verbos;

/* if fdip = 90., i.e. fault is vertical, set true */
bool			fault_is_vertical;

/* allowable distance between obs. and singular point.
 * if |singular_point - obs. point| < eps_dist, evaluation of
 * seismomagnetic field is avoided for this obs. point. */
extern double	eps_dist;


/*c********************
 *c    structures
 *c********************/

/* fault parameters */
typedef struct s_fault_params	fault_params;

struct s_fault_params {

	// crustal parameters
	double	lambda;	// lame constants
	double	mu;			// rigidity
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

	// initial crustal magnetization
	double	mgz_int;	// intensity of magnetization of the uniformly magnetized crust
	double	mgz_inc;	// inclination
	double	mgz_dec;	// declination

	// external field
	double	exf_inc;	// inclination of external field
	double	exf_dec;	// declination

	double	dcurier;	// depth of Curier point isotherm

	// seismomagnetic moment
	// moment vector in x(NS)-y(EW)-z(DownUp) coordinate system
	double	cx;	// x(NS) component
	double	cy;	// y(EW) component
	double	cz;	// z(DownUp) component
	double	c0; // = sqrt (cx^2 + cy^2 + cz^2)

};


/*c********************
 *c  public functions
 *c********************/

/** util.c **/
/* allocate structures */
fault_params		*fault_params_alloc (void);
magnetic_params	*magnetic_params_alloc (void);

/* read parameters from file FILE *fp and store them to global variables
   and pre-allocated structures: fault_params *fault, magnetic_params *mag. */
bool	fread_params (FILE *fp, fault_params *fault, magnetic_params *mag);

/** piezomag.c **/
/* calculate seismomagnetic field on observation point (xobs, yobs, zobs)
   output specified seismomagnetic term: main (0), mirror image (H0) or sub-mirror image term (HI, HIII or HII) */
bool	seismomagnetic_field_term (MagComponent component, SeismoMagTerm term, const fault_params *fault, const magnetic_params *mag,
			double xobs, double yobs, double zobs, double *val);
bool	seismomagnetic_field (MagComponent component, const fault_params *fault, const magnetic_params *mag,
			double xobs, double yobs, double zobs, double *val);

/* calculate seismomagnetic field on grid x=[xobs1:dx:xobs2], y=[yobs1:dy:yobs2], z=zobs */
void	fprintf_seismomagnetic_field_term (FILE *stream, MagComponent component, SeismoMagTerm term, const fault_params *fault, const magnetic_params *mag,
			double xobs1, double xobs2, double dx, double yobs1, double yobs2, double dy, double zobs);
void	fprintf_seismomagnetic_field (FILE *stream, MagComponent component, const fault_params *fault, const magnetic_params *mag,
			double xobs1, double xobs2, double dx, double yobs1, double yobs2, double dy, double zobs);

/** strike.c **/
double	strike_slip (MagComponent component, SeismoMagTerm term, const fault_params *fault, const magnetic_params *mag, double x, double y, double z);

/** dip.c **/
double	dip_slip (MagComponent component, SeismoMagTerm term, const fault_params *fault, const magnetic_params *mag, double x, double y, double z);

/** tensile.c **/
double	tensile_opening (MagComponent component, SeismoMagTerm term, const fault_params *fault, const magnetic_params *mag, double x, double y, double z);

#endif	// _PIEZOMAG_H_
