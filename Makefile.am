ACLOCAL_AMFLAGS = -I m4

#
# Gnulib support functions (generating lib/libgnu.la).
# 

SUBDIRS = lib
AM_CPPFLAGS = -I$(top_srcdir)/lib

#
# Executable target.
#

bin_PROGRAMS = ale
ale_SOURCES = ale.cc ui/ui.cc ui/input.cc d2.cc d3.cc
ale_LDADD = lib/libgnu.la
