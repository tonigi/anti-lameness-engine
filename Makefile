# Copyright 2002 David Hilvert <dhilvert@ugcs.caltech.edu>

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

BITS=8
IMAGEMAGICK=0
PLAIN=0

#
# Compilation flags passed to GCC
#

BITSFLAGS:=$(if $(subst 8,,$(BITS)), -DBITS_16,)
USE_IMAGEMAGICK:=$(subst 0,,$(IMAGEMAGICK))
PPMFLAGS:=$(if $(subst 0,,$(PLAIN)), -DPPM_PLAIN)
CFLAGS:=-DNDEBUG $(BITSFLAGS) $(if $(USE_IMAGEMAGICK), -DUSE_MAGICK $(shell Magick-config --cflags --cppflags), $(PPMFLAGS) -Wall -Os) 
LDFLAGS:=-lm $(if $(USE_IMAGEMAGICK), $(shell Magick-config --ldflags --libs))

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

ale-phony: ale.cc align.cc image_rw.cc *.h render/*.h render/ipc/*.h
	g++ -o ale $(CFLAGS) ale.cc align.cc image_rw.cc $(LDFLAGS)

# XXX: This approach to building a Windows binary is probably very dependent on
# the build configuration.  The above target may be a better place to start for
# Windows users.

ale.exe: ale.cc align.cc image_rw.cc *.h render/*.h render/ipc/*.h
	i586-mingw32msvc-g++ -Wall $(BITSFLAGS) $(PPMFLAGS) -DNDEBUG -o ale.exe -Os ale.cc align.cc image_rw.cc -lm 
