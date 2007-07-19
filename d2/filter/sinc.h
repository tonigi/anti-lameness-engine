// Copyright 2004 David Hilvert <dhilvert@auricle.dyndns.org>,
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

#ifndef __sinc_h__
#define __sinc_h__

/*
 * A filtering class for the sinc function.
 *
 * Notes on Lanczos and Sinc by Dave Martindale et alia are available here:
 *
 * 	http://www.binbooks.com/books/photo/i/l/57186AF8DE
 * 	http://www.binbooks.com/books/photo/i/l/57186AE7E6&orig=1
 */

class sinc : public filter {
public:

	/*
	 * Sinc filter.
	 *
	 * Width is infinite.
	 */
	static double _sinc(ale_pos p) {
		if (fabs(p) < 0.001)
			return 1;
		return sin(M_PI * p) / (M_PI * p);
	}

	static double _sinc(point p) {
		return _sinc(p[0]) * _sinc(p[1]);
	}

	virtual int equals(const filter *f) const {
		if (typeid(*f) == typeid(*this))
			return 1;
		return 0;
	}

	/*
	 * Size of filter support, in number of half-cycles to each side of the
	 * filter center.
	 */
	virtual ale_real support() const {
		double zero = 0;
		double infinity = 1 / zero;
		assert (!isnan(infinity));
		assert (1 < infinity);
		return (ale_real) infinity;
	}

	/*
	 * Response of filter at point p
	 */
	virtual ale_real response(point p) const {
		return _sinc(p);
	}
};
#endif
