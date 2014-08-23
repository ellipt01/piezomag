# piezomag

A program and library for calculating seismomagnetic field, i.e. geomagnetic field change caused by earthquake
due to piezomagnetic effect.

## Description
This package contains library for calculating seismomagnetic field due to a vertical or inclined
rectangular fault located inside the perfectly elastic half space.
Using this library, you can calculate seismomagnetic field easily,
and any dip angle and all types of fault motion, strike-slip, dip-slip and tensile opening is available
for this library.

## Installation

Run the following commands on the top source directry:
```
$ ./autogen.sh
$ ./configure
$ make
$ make install
```

By above, shared library (libpiezomag.so) is installed in your system. By default, this library and relevant files, such as include files, are installed under /usr/local.
To change this, specify prefix as follows:
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

---
```
bool
seismomagnetic_field (int component, const fault_params *fault, const magnetic_params *mag,
    double xobs, double yobs, double zobs, double *val);
```
This function calculates the seismomagnetic field on obervation point ```(xobs, yobs, zobs)```.
```int component``` specifies the output magnetic component and
it takes ```X_COMP (=0:N+S-)```, ```Y_COMP (=1:E+W-)```, ```Z_COMP (=3:Down+Up-)``` or ```TOTAL_FORCE (=0)```.
If obervation point is on the singular point, return ```false```.

---
```
void
fprintf_seismomagnetic_field (FILE *stream, int component, const fault_params *fault, const magnetic_params *mag,
    double xobs1, double xobs2, double dx, double yobs1, double yobs2, double dy, double zobs);
```
This function calculates the seismomagnetic field on the grid ```x=[xobs1:dx:xobs2]```, ```y=[yobs1:dy:yobs2]``` and ```z=zobs```.
Output format is ``` X(NS)  Y(EW)  VAL```.

## Sample program

An example program "main/piez" is also contained in this package.
To run this program, do
```
$ ./main/piez -f <parameter file name>
```

For more details, see work/example_params.dat (example of parameter file) and work/calcomp.sh (example script to run main/piez).

## Licence
LGPL
