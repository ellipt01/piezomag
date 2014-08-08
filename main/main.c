#include <stdio.h>
#include <stdlib.h>

#include "piezomag.h"

int
main (int argc, char **argv)
{
	fault_params		*fault;
	magnetic_params	*mag;

	if (!initialize (argc, argv, &fault, &mag)) exit (1);

  fprintf_piezomagnetic_effect (stdout, output_comp, fault, mag, x_west, x_east, dx, y_south, y_north, dy, z_obs);

  free (fault);
  free (mag);
  return EXIT_SUCCESS;
}
