ACLOCAL_AMFLAGS = -I m4

#
# DejaGNU configuration
#

AUTOMAKE_OPTIONS=dejagnu
RUNTESTFLAGS=--tool ale ALE=./ale --srcdir=./testsuite

#
# Gnulib support functions (generating lib/libgnu.la).
# 

SUBDIRS = lib

#
# Flags
# 

AM_CPPFLAGS = @OPTION_CFLAGS@ @OPTION_CPPFLAGS@ -I$(top_srcdir)/lib
AM_LDFLAGS = @OPTION_LDFLAGS@

#
# Executable target.
#

bin_PROGRAMS = ale
ale_SOURCES = ale.cc ui/ui.cc ui/input.cc d2.cc d3.cc
ale_LDADD = lib/libgnu.la @OPTION_LIBS@
