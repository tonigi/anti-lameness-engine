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

all: ale

clean:
	rm -f ale

ale: ale.c *.h
	gcc -o $@ -Wall -DNDEBUG -O3 ale.c -lm

# XXX: this isn't even remotely portable.  It's how I make the Windows
# binary, though, so I might as well put it here.

ale.exe: ale.c *.h
	i586-mingw32msvc-gcc -Wall -DNDEBUG -o ale.exe -O3 ale.c -lm 
