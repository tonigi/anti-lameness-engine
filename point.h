// Copyright 2002 David Hilvert <dhilvert@ugcs.caltech.edu>

/*  This file is part of the Anti-Lamenessing Engine.

    The Anti-Lamenessing Engine is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    The Anti-Lamenessing Engine is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the Anti-Lamenessing Engine; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef __point_h__
#define __point_h__

#include "my_real.h"

/*
 * Structure to describe a point
 */
class point {
private:	
	my_real x[2];

public:
	point() {
	}

	point(my_real x0, my_real x1) {
		x[0] = x0;
		x[1] = x1;
	}

	my_real &operator[](int i) {
		assert (i >= 0);
		assert (i < 2);

		return x[i];
	}
};

#endif
