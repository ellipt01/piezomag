# piezomag

A program and library for calculating seismomagnetic field caused by inclined rectangular fault.

## Description
This package contains library for calculating seismomagnetic field due to fault motion, strike-slip, dip-slip and tensile opening, on inclined rectangular fault plane which located inside the perfectly elastic half space.
Any dip angle of fault plane are available including 90 degrees (vertical fault plane).

## Installation

Run the following commands on the top source directry:
```
$ ./autogen.sh
$ ./configure
$ make
$ make install
```

By above, shared library (libpiezomag.so) is installed in your system. By default, this library and relevant files, such as include files, are installed under /usr/local.
To change this, please specify prefix as follows:
```
$ ./configure --prefix=hoge
```

## Public functions
This library provides the following public functions:
```
bool
fread_params (FILE *fp, fault_params *fault, magnetic_params *mag);
```
This function reads user-defined parameters from input file and store them to global variables and pre-allocated structures:

```fault_params *fault: ``` structure which stores fault parameters.

```magnetic_params *mag:``` structure which stores crustal magnetic properties.

These structures are allocated by the following functions, respectivly:

```fault_params *fault_params_alloc (void)```

```magnetic_params *magnetic_params_alloc (void)```

```
double
seismomagnetic_effect (int component, const fault_params *fault, const magnetic_params *mag,
    double xobs, double yobs, double zobs);
```
This function calculates the seismomagnetic field on obervation point ```(xobs, yobs, zobs)```.
```int component``` specifies the output magnetic component and it takes ```X_COMP (=0)```, ```Y_COMP (=1)```, ```Z_COMP (=3)``` or ```TOTAL_FORCE (=0)```.

```
void
fprintf_seismomagnetic_effect (FILE *stream, int component, const fault_params *fault, const magnetic_params *mag,
    double xobs1, double xobs2, double dx, double yobs1, double yobs2, double dy, double zobs);
```
This function calculates the seismomagnetic field on the grid ```x=[xobs1:dx:xobs2]```, ```y=[yobs1:dy:yobs2]``` and ```z=zobs```.


## Sample program

An example program "main/piez" is also contained in this package.
To run this program, do
```
$ ./main/piez -f <parameter file name>
```

For more details, please see work/example_params.dat (example of parameter file) and work/calcomp.sh (example script to run main/piez).

## Licence
LGPL
