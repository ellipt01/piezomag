dnl Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(piezomag, 0.5.2)

PIEZOMAG_CURRENT=1
PIEZOMAG_REVISION=0
PIEZOMAG_AGE=0
PIEZOMAG_LT_VERSION="${PIEZOMAG_CURRENT}:${PIEZOMAG_REVISION}:${PIEZOMAG_AGE}"
AC_SUBST(PIEZOMAG_LT_VERSION)

AC_PROG_CC
AC_CONFIG_HEADERS([config.h])

AC_CANONICAL_SYSTEM

dnl Initialize automake
AM_INIT_AUTOMAKE([foreign subdir-objects -Wall])

AM_PROG_AR([ACT-IF-FAIL])

dnl this allows us specify individual liking flags for each target
AM_PROG_CC_C_O 

dnl Initialize Libtool
LT_INIT

dnl Check if Libtool is present
dnl Libtool is used for building share libraries 
AC_PROG_LIBTOOL

# pkgconfig
PIEZOMAG_DESCRIPTION="A library for calculating seismomagnetic effect due to a rectangular fault"
AC_SUBST([PIEZOMAG_DESCRIPTION])
AC_CONFIG_FILES([piezomag.pc:pkgconfig/piezomag.pc.in])

AC_CONFIG_FILES(Makefile
                include/Makefile
                src/Makefile
                main/Makefile)
AC_OUTPUT
