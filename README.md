# piezomag

A program and library for calculating seismomagnetic field caused by inclined rectangular fault.

## Description
This package contains library for calculating seismomagnetic field due to fault motion, strike-slip, dip-slip and tensile opening, on inclined rectangular fault plane which located inside the perfectly elastic half space.
Any dip angle of fault plane are available except for 90 degrees (vertical fault plane).

## Installation

Run the following commands on the top source directry:
```
$ ./autogen.sh
$ ./configure
$ make
$ make install
```

By above, shared library (libpiezomag.so) will be installed in your system. By default, this library and relevant files, such as include files, are installed under /usr/local.
To change this, please specify prefix as follows:
```
$ ./configure --prefix=hoge
```

## Sample program

An example program "main/piez" is also contained in this package.
To run this program, do
```
$ ./main/piez -f <parameter file name>
```

For more details, please see work/example_params.dat (example of parameter file) and work/calcomp.sh (example script to run main/piez).

## Licence
LGPL



