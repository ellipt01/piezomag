## Process this file with automake to produce Makefile.in

INCLUDES = -I$(top_srcdir)/include

AM_CFLAGS =\
	 -Wall\
	 -g

lib_LTLIBRARIES = libpiezomag.la

libpiezomag_la_SOURCES =\
	utils.c\
	deriv_dipoles.c\
	deriv_quadpoles.c\
	deriv_octpoles.c\
	dip_slip.c\
	strike_slip.c\
	tensile_opening.c\
	piezomag.c

libpiezomag_la_CPPFLAGS = -I$(top_srcdir)/include

libpiezomag_la_LDFLAGS = -L$(top_srcdir)/src

libpiezomag_la_LIBADD =

