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

#ifndef __pixel_accum_h__
#define __pixel_accum_h__

#include "pixel.h"

/*
 * Structure to accumulate values over many pixels.
 */

class pixel_accum {
private:	
	ale_accum x[3];

public:
	pixel_accum() {
		x[0] = 0;
		x[1] = 0;
		x[2] = 0;
	}

	pixel_accum(ale_accum x0, ale_accum x1, ale_accum x2) {
		x[0] = x0;
		x[1] = x1;
		x[2] = x2;
	}

	pixel_accum(pixel p) {
		x[0] = p[0];
		x[1] = p[1];
		x[2] = p[2];
	}

	operator pixel() {
		pixel result;
		
		result[0] = x[0];
		result[1] = x[1];
		result[2] = x[2];

		return result;
	}

//	Due to automatic typecasts and automatic int <==> ale_accum *
//	conversions, this can cause some really weird bugs.
//
//	pixel_accum(ale_accum *_x) {
//		x[0] = _x[0];
//		x[1] = _x[1];
//		x[2] = _x[2];
//	}

	const ale_accum &operator[](int i) const {
		assert (i >= 0);
		assert (i < 3);

		return x[i];
	}

	ale_accum &operator[](int i) {
		assert (i >= 0);
		assert (i < 3);

		return x[i];
	}

	pixel_accum operator+(pixel_accum p) const {
		return pixel_accum(p[0] + x[0], p[1] + x[1], p[2] + x[2]);
	}

	pixel_accum operator-(pixel_accum p) const {
		return pixel_accum(x[0] - p[0], x[1] - p[1], x[2] - p[2]);
	}

	pixel_accum operator/(pixel_accum p) const {
		return pixel_accum(x[0] / p[0], x[1] / p[1], x[2] / p[2]);
	}

	pixel_accum operator/(ale_accum d) const {
		return pixel_accum(x[0] / d, x[1] / d, x[2] / d);
	}

	pixel_accum operator*(pixel_accum p) const {
		return pixel_accum(x[0] * p[0], x[1] * p[1], x[2] * p[2]);
	}

	pixel_accum operator*(ale_accum d) const {
		return pixel_accum(x[0] * d, x[1] * d, x[2] * d);
	}

	pixel_accum operator+=(pixel_accum p) {
		return pixel_accum(x[0] += p[0], x[1] += p[1], x[2] += p[2]);
	}

	pixel_accum operator*=(pixel_accum p) {
		return pixel_accum(x[0] *= p[0], x[1] *= p[1], x[2] *= p[2]);
	}

	pixel_accum operator*=(ale_accum d) {
		return pixel_accum(x[0] *= d, x[1] *= d, x[2] *= d);
	}

	pixel_accum operator/=(pixel_accum p) {
		return pixel_accum(x[0] /= p[0], x[1] /= p[1], x[2] /= p[2]);
	}

	pixel_accum operator/=(ale_accum d) {
		return pixel_accum(x[0] /= d, x[1] /= d, x[2] /= d);
	}
};

inline pixel_accum operator*(float d, const pixel_accum &p) {
	return p * d;
}

inline pixel_accum operator*(double d, const pixel_accum &p) {
	return p * d;
}

inline std::ostream &operator<<(std::ostream &o, const pixel_accum &p) {
	o << "[" << p[0] << " " << p[1] << " " << p[2] << "]";
	return o;
}

inline pixel_accum ppow(pixel_accum p, float d) {
	return pixel_accum(
			pow((ale_accum) p[0], d),
			pow((ale_accum) p[1], d),
			pow((ale_accum) p[2], d));
}

inline pixel_accum ppow(pixel_accum p, double d) {
	return pixel_accum(
			pow((ale_accum) p[0], d),
			pow((ale_accum) p[1], d),
			pow((ale_accum) p[2], d));
}

#endif
