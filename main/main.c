#include <stdio.h>
#include <stdlib.h>

#include "constants.h"
#include "piez.h"

int
main (int argc, char **argv)
{
	if (!initialize (argc, argv)) exit (1);

  fprintf_piezomagnetic_effect (stdout, output_comp, x01, x02, dx, y01, y02, dy, z_obs);

  return EXIT_SUCCESS;
}
