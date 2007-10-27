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

	const ale_pos &operator[](unsigned int i) const {
#if 0
		/*
		 * This may be expensive.
		 */
		assert (i < 2);
#endif

		return x[i];
	}

	ale_pos &operator[](unsigned int i) {
#if 0
		/*
		 * This may be expensive.
		 */
		assert (i < 2);
#endif

		return x[i];
	}

	point operator+(point p) const {
		return point(p[0] + x[0], p[1] + x[1]);
	}

	point operator-(point p) const {
		return point(x[0] - p[0], x[1] - p[1]);
	}

	point operator-() const {
		return point(-x[0], -x[1]);
	}

	point operator+=(point p) { 
		(*this) = (*this) + p;

		return *this;
	}

	point operator-=(point p) { 
		(*this) = (*this) - p;

		return *this;
	}

	point mult(ale_pos d) const {
		return point(x[0] * d, x[1] * d);
	}

	point operator*(point p) const {
		/*
		 * element-wise multiplication
		 */
		return point(x[0] * p[0], x[1] * p[1]);
	}
	
	point operator *=(ale_pos d) {
		(*this) = mult(d);
		return *this;
	}

	point operator/(ale_pos d) const {
		return point(x[0] / d, x[1] / d);
	}

	ale_pos normsq() const {
		return x[0] * x[0] + x[1] * x[1];
	}

	ale_pos norm() const {
		return sqrt(normsq());
	}

	ale_pos absmaxnorm() const {
		ale_pos a = fabs(x[0]);
		ale_pos b = fabs(x[1]);

		return (a > b) ? a : b;
	}

	ale_pos lengthtosq(point p) const {
		point diff = operator-(p);

		return diff[0] * diff[0] + diff[1] * diff[1];
	}
	ale_pos lengthto(point p) const {
		return sqrt(lengthtosq(p));
	}

	ale_pos dproduct(point p) const {
		return (x[0] * p[0] + x[1] * p[1]);
	}

	ale_pos anglebetw(point p, point q) {
		/*
		 * by the law of cosines, the cosine is equal to:
		 *
		 *  	(lengthtosq(p) + lengthtosq(q) - p.lengthtosq(q))
		 *    / (2 * lengthto(p) * lengthto(q))
		 */

		ale_pos to_p = lengthtosq(p);
		ale_pos to_q = lengthtosq(q);

		ale_pos cos_of = (double) (to_p + to_q - p.lengthtosq(q))
			       / (2 * sqrt(to_p) * sqrt(to_q));

		return acos(cos_of);
	}

	static point posinf() {
		ale_pos a = +1;
		ale_pos z = +0;

		a = a / z;

		assert (isinf(a) && a > 0);

		return point(a, a);
	}

	static point neginf() {
		point n = -posinf();

		assert (isinf(n[0]) && n[0] < 0);

		return n;
	}

	void accumulate_max(point p) {
		for (int d = 0; d < 2; d++)
			if (p[d] > x[d])
				x[d] = p[d];
	}

	void accumulate_min(point p) {
		for (int d = 0; d < 2; d++)
			if (p[d] < x[d])
				x[d] = p[d];
	}

	static point undefined() {
		double a = 0;

		point p(0, 0);

		return p / a;
	}

	int defined() const {
		return (!isnan(x[0])
		     && !isnan(x[1]));
	}

	int finite() const {
		return (::finite(x[0])
		     && ::finite(x[1]));
	}

	static int defined(const point &p) {
		return p.defined();
	}
};

inline point operator*(const point &p, double d) {
	return p.mult(d);
}
inline point operator*(double d, const point &p) {
	return p.mult(d);
}
inline point operator*(float d, const point &p) {
	return p.mult(d);
}
#endif
