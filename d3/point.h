// Copyright 2002, 2004 David Hilvert <dhilvert@auricle.dyndns.org>,
//                                    <dhilvert@ugcs.caltech.edu>

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

#ifndef __d3point_h__
#define __d3point_h__

/*
 * Structure to describe a point in three dimensions.
 */
class point {
private:	
	ale_pos x[3];

public:
	point() {
		x[0] = 0;
		x[1] = 0;
		x[2] = 0;
	}

	point(ale_pos x0, ale_pos x1, ale_pos x2) {
		x[0] = x0;
		x[1] = x1;
		x[2] = x2;
	}

	/*
	 * Z-values of zero are almost never the right thing to do, but
	 * for cases when they are ...
	 */
	point(d2::point p) {
		x[0] = p[0];
		x[1] = p[1];
		x[2] = 0;
	}

	const ale_pos &operator[](int i) const {
		assert (i >= 0);
		assert (i < 3);

		return x[i];
	}

	ale_pos &operator[](int i) {
		assert (i >= 0);
		assert (i < 3);

		return x[i];
	}

	d2::point xy() {
		d2::point result;

		result[0] = x[0];
		result[1] = x[1];

		return result;
	}

	point operator-(point p) const {
		return point(x[0] - p[0], x[1] - p[1], x[2] - p[2]);
	}

	point operator /=(ale_real r) {
		x[0] /= r;
		x[1] /= r;
		x[2] /= r;

		return *this;
	}

	point operator *=(ale_real r) {
		x[0] *= r;
		x[1] *= r;
		x[2] *= r;

		return *this;
	}

	int operator !=(point p) {
		return (x[0] != p[0]
		     || x[1] != p[1]
		     || x[2] != p[2]);
	}

	ale_pos lengthtosq(point p) const {
		point diff = operator-(p);

		return diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2];
	}
	ale_pos lengthto(point p) const {
		return sqrt(lengthtosq(p));
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

		ale_pos cos_of = (to_p + to_q - p.lengthtosq(q))
			       / (2 * sqrt(to_p) * sqrt(to_q));

		return acos(cos_of);
	}
};

#endif