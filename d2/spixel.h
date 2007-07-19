// Copyright 2002 David Hilvert <dhilvert@auricle.dyndns.org>,
//                              <dhilvert@ugcs.caltech.edu>

/*  This file is part of the Anti-Lamenessing Engine.

    The Anti-Lamenessing Engine is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    The Anti-Lamenessing Engine is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the Anti-Lamenessing Engine; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifndef __spixel_h__
#define __spixel_h__

/*
 * Structure to describe a pixel to be used in storage.
 */

#define HALF 3

#if ALE_COLORS != HALF

typedef pixel spixel;

#else 

class spixel {
private:	
	ale_sreal x[3];

public:
	spixel() {
		x[0] = 0;
		x[1] = 0;
		x[2] = 0;
	}

	spixel(pixel p) {
		x[0] = p[0];
		x[1] = p[1];
		x[2] = p[2];
	}

	spixel(ale_sreal x0, ale_sreal x1, ale_sreal x2) {
		x[0] = x0;
		x[1] = x1;
		x[2] = x2;
	}

	operator pixel() const {
		pixel result;

		result[0] = x[0];
		result[1] = x[1];
		result[2] = x[2];

		return result;
	}

	const ale_sreal &operator[](int i) const {
		assert (i >= 0);
		assert (i < 3);

		return x[i];
	}

	ale_sreal &operator[](int i) {
		assert (i >= 0);
		assert (i < 3);

		return x[i];
	}

	spixel operator+=(pixel p) {
		return spixel(x[0] += p[0], x[1] += p[1], x[2] += p[2]);
	}

	spixel operator*=(pixel p) {
		return spixel(x[0] *= p[0], x[1] *= p[1], x[2] *= p[2]);
	}

	spixel operator*=(ale_real d) {
		return spixel(x[0] *= d, x[1] *= d, x[2] *= d);
	}

	spixel operator/=(pixel p) {
		return spixel(x[0] /= p[0], x[1] /= p[1], x[2] /= p[2]);
	}

	spixel operator/=(ale_real d) {
		return spixel(x[0] /= d, x[1] /= d, x[2] /= d);
	}
};

#endif

#undef HALF

#endif
