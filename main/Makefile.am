## Process this file with automake to produce Makefile.in

# sample program
bin_PROGRAMS = piez
piez_SOURCES = main.c
piez_CPPFLAGS =\
	-I$(top_srcdir)/include\
	-I$(top_srcdir)/main

piez_LDFLAGS = -L$(top_srcdir)/src

piez_LDADD = -lpiezomag -lm
