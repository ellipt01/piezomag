dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(piezomag, 0.4.1)
AC_PROG_CC
AC_CONFIG_HEADERS([config.h])

AC_CANONICAL_SYSTEM

dnl Initialize automake
AM_INIT_AUTOMAKE

dnl this allows us specify individual liking flags for each target
AM_PROG_CC_C_O 

dnl Initialize Libtool
LT_INIT

dnl Check if Libtool is present
dnl Libtool is used for building share libraries 
AC_PROG_LIBTOOL

# pkgconfig
PIEZOMAG_DESCRIPTION="Library for calculating seismo-magnetic effect due to an inclined rectangular fault"
AC_SUBST([PIEZOMAG_DESCRIPTION])
AC_CONFIG_FILES([piezomag.pc:pkgconfig/piezomag.pc.in])

AC_CONFIG_FILES(Makefile
                include/Makefile
                src/Makefile
                main/Makefile)
AC_OUTPUT
