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

#ifndef __pixel_h__
#define __pixel_h__

/*
 * Structure to describe a pixel
 */

class pixel {
private:	
	ale_real x[3];

public:
	pixel() {
		x[0] = 0;
		x[1] = 0;
		x[2] = 0;
	}

	pixel(ale_real x0, ale_real x1, ale_real x2) {
		x[0] = x0;
		x[1] = x1;
		x[2] = x2;
	}

//	Due to automatic typecasts and automatic int <==> ale_real *
//	conversions, this can cause some really weird bugs.
//
//	pixel(ale_real *_x) {
//		x[0] = _x[0];
//		x[1] = _x[1];
//		x[2] = _x[2];
//	}

	const ale_real &operator[](int i) const {
		assert (i >= 0);
		assert (i < 3);

		return x[i];
	}

	ale_real &operator[](int i) {
		assert (i >= 0);
		assert (i < 3);

		return x[i];
	}

	pixel operator+(pixel p) const {
		return pixel(p[0] + x[0], p[1] + x[1], p[2] + x[2]);
	}

	pixel operator-(pixel p) const {
		return pixel(x[0] - p[0], x[1] - p[1], x[2] - p[2]);
	}

	pixel operator/(pixel p) const {
		return pixel(x[0] / p[0], x[1] / p[1], x[2] / p[2]);
	}

	pixel operator/(ale_real d) const {
		return pixel(x[0] / d, x[1] / d, x[2] / d);
	}

	pixel operator*(pixel p) const {
		return pixel(x[0] * p[0], x[1] * p[1], x[2] * p[2]);
	}

	pixel mult(ale_real d) const {
		return pixel(x[0] * d, x[1] * d, x[2] * d);
	}

	pixel operator+=(pixel p) {
		return pixel(x[0] += p[0], x[1] += p[1], x[2] += p[2]);
	}

	pixel operator*=(pixel p) {
		return pixel(x[0] *= p[0], x[1] *= p[1], x[2] *= p[2]);
	}

	pixel operator*=(ale_real d) {
		return pixel(x[0] *= d, x[1] *= d, x[2] *= d);
	}

	pixel operator/=(pixel p) {
		return pixel(x[0] /= p[0], x[1] /= p[1], x[2] /= p[2]);
	}

	pixel operator/=(ale_real d) {
		return pixel(x[0] /= d, x[1] /= d, x[2] /= d);
	}

	pixel clamp() const {
		pixel result;

		for (int i = 0; i < 3; i++)
			if (x[i] > 1.0)
				result[i] = 1.0;
			else if (x[i] < 0.0)
				result[i] = 0.0;
			else
				result[i] = x[i];

		return result;
	}
		

	pixel abs() {
		return pixel(fabs(x[0]), fabs(x[1]), fabs(x[2]));
	}

	ale_real normsq() {
		return x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
	}

	ale_real norm() {
		return sqrt(normsq());
	}

	ale_real lnorm() {
		return x[0] + x[1] + x[2];
	}

	ale_real maxabs_norm() {
		ale_real m = fabs(x[0]);
		if (fabs(x[1]) > m)
			m = fabs(x[1]);
		if (fabs(x[2]) > m)
			m = fabs(x[2]);

		return m;
	}

	ale_real minabs_norm() {
		ale_real m = fabs(x[0]);
		if (fabs(x[1]) < m)
			m = fabs(x[1]);
		if (fabs(x[2]) < m)
			m = fabs(x[2]);

		return m;
	}

	ale_real min_norm() {
		ale_real m = x[0];
		if (x[1] < m)
			m = x[1];
		if (x[2] < m)
			m = x[2];

		return m;
	}

	ale_real max_norm() {
		ale_real m = x[0];
		if (x[1] > m)
			m = x[1];
		if (x[2] > m)
			m = x[2];

		return m;
	}

	static pixel zero() {
		return pixel(0, 0, 0);
	}

	int operator==(const pixel &p) {
		return x[0] == p[0]
		    && x[1] == p[1]
		    && x[2] == p[2];
	}

	int operator!=(const pixel &p) {
		return !operator==(p);
	}
};

inline pixel operator*(const pixel &p, float d) {
	return p.mult(d);
}
inline pixel operator*(const pixel &p, double d) {
	return p.mult(d);
}
inline pixel operator*(double d, const pixel &p) {
	return p.mult(d);
}
inline pixel operator*(float d, const pixel &p) {
	return p.mult(d);
}

inline std::ostream &operator<<(std::ostream &o, const pixel &p) {
	o << "[" << p[0] << " " << p[1] << " " << p[2] << "]";
	return o;
}

inline pixel ppow(pixel p, double d) {
	return pixel(
			pow((double) p[0], d),
			pow((double) p[1], d),
			pow((double) p[2], d));
}

inline pixel ppow(pixel p, float d) {
	return pixel(
			pow((double) p[0], d),
			pow((double) p[1], d),
			pow((double) p[2], d));
}

#endif
