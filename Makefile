# Copyright 2002, 2003 David Hilvert <dhilvert@auricle.dyndns.org>,
#                                    <dhilvert@ugcs.caltech.edu>

#  This file is part of the Anti-Lamenessing Engine.
#
#  The Anti-Lamenessing Engine is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  The Anti-Lamenessing Engine is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with Anti-Lamenessing Engine; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#
# Compile options (defaults)
#

IMAGEMAGICK=0
FFTW=0
DEBUG=0
COLORS=SINGLE
COORDINATES=SINGLE

#
# Variable that is false when ImageMagick is not being used
#

use_imagemagick:=$(subst 0,,$(IMAGEMAGICK))

#
# Compiler flags
#

DEBUG_CFLAGS:=$(if $(subst 0,,$(DEBUG)),-DDEBUG,-DNDEBUG)
IMAGEMAGICK_CFLAGS:=$(if $(use_imagemagick),-DUSE_MAGICK $(shell Magick-config --cflags --cppflags),)
IMAGEMAGICK_LDFLAGS:=$(if $(use_imagemagick),$(shell Magick-config --ldflags --libs),)
FFTW_CFLAGS:=$(if $(subst 0,,$(FFTW)),-DUSE_FFTW,)
PRECISION_CFLAGS:=$(if $(subst SINGLE,,$(COLORS)),,-DALE_COLORS=SINGLE)\
                  $(if $(subst DOUBLE,,$(COLORS)),,-DALE_COLORS=DOUBLE)\
                  $(if $(subst HALF,,$(COLORS)),,-DALE_COLORS=HALF)\
                  $(if $(subst SINGLE,,$(COORDINATES)),,-DALE_COORDINATES=SINGLE)\
                  $(if $(subst DOUBLE,,$(COORDINATES)),,-DALE_COORDINATES=DOUBLE)
FFTW_LDFLAGS:=$(if $(subst 0,,$(FFTW)),-lfftw3,)

CFLAGS:= $(DEBUG_CFLAGS) $(FFTW_CFLAGS) $(PRECISION_CFLAGS) \
         $(if $(use_imagemagick),$(IMAGEMAGICK_CFLAGS),-Wall -O2) 
LDFLAGS:=$(if $(use_imagemagick),$(IMAGEMAGICK_LDFLAGS)) $(FFTW_LDFLAGS) -lm

#
# We're using 'ale-phony' because we're using make for configuration,
# and we want to make sure that the output is updated if the configuration
# changes.  
# 
# XXX: is there a better way to do this?
#

all: ale-phony

clean:
	rm -f ale

ale-phony: ale.cc d2.cc *.h d2/*.h d2/render/*.h d2/render/psf/*.h
	g++ -o ale $(CFLAGS) ale.cc d2.cc $(LDFLAGS)

# The following approach to building a Windows binary is probably very
# dependent on the host platform configuration.  The above target may be a
# better place to start for Windows users.

ale.exe: ale.cc d2.cc *.h d2/*.h d2/render/*.h d2/render/psf/*.h
	i586-mingw32msvc-g++ -Wall $(PRECISION_CFLAGS) $(DEBUG_CFLAGS) -o ale.exe -O2 ale.cc d2.cc -lm 
