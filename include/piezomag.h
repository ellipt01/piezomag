#ifndef _PIEZOMAG_H_
#define _PIEZOMAG_H_

#include <stdbool.h>

/*c***********
 *c   enums
 *c***********/

/* definition of magnetic components */
typedef enum {
	MAG_COMP_NONE = -1,	// dummy
	MAG_COMP_F    =  0,	// total force
	MAG_COMP_X    =  1,	// x component
	MAG_COMP_Y    =  2,	// y component
	MAG_COMP_Z    =  3,	// z component
} MagComp;

/* flags for seismomagnetic terms */
typedef enum {
	SEISMO_MAG_MAIN       =  1 << 0,	// main term (0)
	SEISMO_MAG_MIRROR     =  1 << 1,	// mirror image (H0)
	SEISMO_MAG_SUBMIRROR  =  1 << 2,	// sub-mirror image (HI, HIII or HII)
	// total seismomagnetic field (0 + H0 + (HI, HIII or HII))
	SEISMO_MAG_TOTAL =  SEISMO_MAG_MAIN | SEISMO_MAG_MIRROR | SEISMO_MAG_SUBMIRROR
} SeismoMagTerm;


/*c****************
 *c   structures
 *c****************/

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

	// external geomagnetic field
	double	exf_inc;	// inclination of external field
	double	exf_dec;	// declination

	double	dcurier;	// depth of Curier point isotherm

	// seismomagnetic moment on fault coordinate system
	double	cx;	// x component
	double	cy;	// y component
	double	cz;	// z component
	double	c0; // = sqrt (cx^2 + cy^2 + cz^2)

};


/*c**********************
 *c   public functions
 *c**********************/

/*** utils.c ***/

/* allocate structures */
fault_params		*fault_params_alloc (void);
magnetic_params	*magnetic_params_alloc (void);

/*** params_io.c ***/

/* read parameters from file FILE *fp and store them to
 * pre-allocated structures: fault_params *fault, magnetic_params *mag. */
bool	fread_params (FILE *fp, fault_params *fault, magnetic_params *mag);

/*** piezomag.c ***/

/* calculate seismomagnetic term (main(0), mirror image(H0) or sub-mirror image(HI, HIII or HII))
 * on observation point (xobs, yobs, zobs) */
bool	seismomagnetic_field_term (MagComp component, SeismoMagTerm term,
			const fault_params *fault, const magnetic_params *mag,
			double xobs, double yobs, double zobs, double *val);

/* calculate seismomagnetic field on observation point (xobs, yobs, zobs) */
bool	seismomagnetic_field (MagComp component,
			const fault_params *fault, const magnetic_params *mag,
			double xobs, double yobs, double zobs, double *val);

/* calculate and fprintf seismomagnetic term on grid x=[xobs1:dx:xobs2], y=[yobs1:dy:yobs2], z=zobs */
void	fprintf_seismomagnetic_field_term (FILE *stream, MagComp component, SeismoMagTerm term,
			const fault_params *fault, const magnetic_params *mag,
			double xobs1, double xobs2, double dx, double yobs1, double yobs2, double dy, double zobs);

/* calculate and fprintf seismomagnetic field on grid x=[xobs1:dx:xobs2], y=[yobs1:dy:yobs2], z=zobs */
void	fprintf_seismomagnetic_field (FILE *stream, MagComp component,
			const fault_params *fault, const magnetic_params *mag,
			double xobs1, double xobs2, double dx, double yobs1, double yobs2, double dy, double zobs);

#endif	// _PIEZOMAG_H_
