ACLOCAL_AMFLAGS = -I m4

#
# Local targets.
#

clean-local: dejagnu-clean

#
# DejaGNU configuration
#

AUTOMAKE_OPTIONS=dejagnu
RUNTESTFLAGS=--tool ale ALE=./ale --srcdir=./testsuite

dejagnu-clean:
	rm -f ale.sum ale.log
	find testsuite -name "*temp.*" -print0 | xargs -0 rm -rf
	find testsuite -name "*.output.*" -print0 | xargs -0 rm -rf

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
