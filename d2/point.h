// Copyright 2002 David Hilvert <dhilvert@auricle.dyndns.org>,
//                              <dhilvert@ugcs.caltech.edu>

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

#ifndef __d2point_h__
#define __d2point_h__

/*
 * Structure to describe a point
 */

class point {
private:	
	ale_pos x[2];

public:
	point() {
	}

	point(ale_pos x0, ale_pos x1) {
		x[0] = x0;
		x[1] = x1;
	}

	const ale_pos &operator[](int i) const {
		assert (i >= 0);
		assert (i < 2);

		return x[i];
	}

	ale_pos &operator[](int i) {
		assert (i >= 0);
		assert (i < 2);

		return x[i];
	}

	point operator+(point p) const {
		return point(p[0] + x[0], p[1] + x[1]);
	}

	point operator-(point p) const {
		return point(x[0] - p[0], x[1] - p[1]);
	}
};

#endif
