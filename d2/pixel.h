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

	pixel(const pixel &p) {
		x[0] = p[0];
		x[1] = p[1];
		x[2] = p[2];
	}

	pixel &operator=(const pixel &p) {
		x[0] = p[0];
		x[1] = p[1];
		x[2] = p[2];

		return (*this);
	}


//	Due to automatic typecasts and automatic int <==> ale_real *
//	conversions, this can cause some really weird bugs.
//
//	pixel(ale_real *_x) {
//		x[0] = _x[0];
//		x[1] = _x[1];
//		x[2] = _x[2];
//	}

	const ale_real &operator[](unsigned int i) const {
#if 0
		/*
		 * This may be expensive.
		 */
		assert (i < 3);
#endif
		return x[i];
	}

	ale_real &operator[](unsigned int i) {
#if 0
		/*
		 * This may be expensive.
		 */
		assert (i < 3);
#endif
		return x[i];
	}

	pixel operator+(pixel p) const {
		return pixel(p[0] + x[0], p[1] + x[1], p[2] + x[2]);
	}

	pixel operator-(pixel p) const {
		return pixel(x[0] - p[0], x[1] - p[1], x[2] - p[2]);
	}

	pixel operator-() const {
		return pixel(-x[0], -x[1], -x[2]);
	}

	pixel operator/(pixel p) const {
		return pixel(x[0] / p[0], x[1] / p[1], x[2] / p[2]);
	}

	pixel operator/(ale_real d) const {
		return pixel(x[0] / d, x[1] / d, x[2] / d);
	}

	pixel mult(pixel p) const {
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

	ale_real min_norm() const {
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

	static pixel one() {
		return pixel(1, 1, 1);
	}

	int operator==(const pixel &p) {
		return x[0] == p[0]
		    && x[1] == p[1]
		    && x[2] == p[2];
	}

	int operator!=(const pixel &p) {
		return !operator==(p);
	}

	int finite() {
		return ::finite(x[0]) && ::finite(x[1]) && ::finite(x[2]);
	}

	static pixel undefined() {
		ale_real zero = 0;
		return pixel(zero / zero, zero / zero, zero / zero);

	}
};

inline pixel operator*(const pixel &p, const pixel &q) {
	return p.mult(q);
}

template<typename T>
inline pixel operator*(T d, const pixel &p) {
	return p.mult(d);
}
template<typename T>
inline pixel operator*(const pixel &p, T d) {
	return p.mult(d);
}

inline std::ostream &operator<<(std::ostream &o, const pixel &p) {
	o << "[" << (double) p[0] << " " << (double) p[1] << " " << (double) p[2] << "]";
	return o;
}

template<typename T>
inline pixel ppow(pixel p, T d) {
	return pixel(
			pow(p[0], d),
			pow(p[1], d),
			pow(p[2], d));
}

inline pixel pexp(pixel p) {
	return pixel(
			exp((double) p[0]),
			exp((double) p[1]),
			exp((double) p[2]));
}

inline pixel psqrt(pixel p) {
	return pixel(
			sqrt(p[0]),
			sqrt(p[1]),
			sqrt(p[2]));
}

#endif
