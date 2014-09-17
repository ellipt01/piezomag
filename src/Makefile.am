## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir)/include

AM_CFLAGS =\
	 -Wall\
	 -g

lib_LTLIBRARIES = libpiezomag.la

libpiezomag_la_SOURCES =\
	private.c\
	params_io.c\
	utils.c\
	deriv_dipoles.c\
	deriv_quadpoles.c\
	deriv_octpoles.c\
	dip_slip.c\
	strike_slip.c\
	tensile_opening.c\
	piezomag.c

libpiezomag_la_CPPFLAGS = -I$(top_srcdir)/include

libpiezomag_la_LDFLAGS = -L$(top_srcdir)/src -version-info $(PIEZOMAG_LT_VERSION)

libpiezomag_la_LIBADD =

