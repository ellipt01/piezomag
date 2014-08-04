#include <stdio.h>
#include <math.h>


int
main (void)
{
	double sign = 1.0;
	double xi = 0.;
	double et = -10.0;
	double qq = 0.01;

	double r2 = pow (xi, 2.0) + pow (et, 2.0) + pow (qq, 2.0);
	double r = sqrt (r2);

	double re = r + et;
	double ir		= 1.0 / r;
	double ir3	 = 1.0 / pow (r, 3.0);
	double ire2	= 1.0 / pow (re, 2.0);
	double ire3	= 1.0 / pow (re, 3.0);
	double ir3e2 = 1.0 / (pow (r, 3.0) * pow (re, 2.0));
	double ir5e3 = 1.0 / (pow (r, 5.0) * pow (re, 3.0));

	double cd = 1.0 / sqrt (2.0);
	double sd = 1.0 / sqrt (2.0);
	double sd2 = sd * sd;
	double yy = et * cd + qq * sd;
	double cc = qq * cd - et * sd;

	double val, val1, val2, val3, val4, val5, val6, val7;

	val1 = - 2.0 * yy * sd2 * (ir * ire3) - cd * (ir * ire2);
	val6 = 2.0 * (yy * cd + sign * cc * sd) * cd * (ir * ire3);

	val2 = yy * (yy * cd + sign * cc * sd) * (3.0 * r + et) * (ir3 * ire3);
	val3 = - 2.0 * yy * (r + re) * (ir3 * ire2);
	val4 = yy * yy * yy * (8.0 * r2 + 9.0 * et * r + 3.0 * et * et) * ir5e3;
	val5 = - 2.0 * sd2 * cd * ire3;
	val7 = yy * yy * (3.0 * r + et) * cd * (ir3 * ire3);

	val =	val1 + val2 + val3 + val4 + val5 + val6 + val7;

	fprintf (stdout, "%.3e %.3e %.3e %.3e\n", val1, val2, val3, val4);
	fprintf (stdout, "%.3e %.3e %.3e %.3e\n", val5, val6, val7, val);

	return 0;
}
