// Copyright 2004 David Hilvert <dhilvert@auricle.dyndns.org>,
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

#ifndef __lanczos_h__
#define __lanczos_h__

/*
 * A lanczos filter class.
 *
 * Notes on Lanczos and Sinc by Dave Martindale et alia are available here:
 *
 * 	http://www.binbooks.com/books/photo/i/l/57186AF8DE
 * 	http://www.binbooks.com/books/photo/i/l/57186AE7E6&orig=1
 */

class lanczos : public filter {
private:
	double half_width;

	/*
	 * Lanczos function
	 *
	 * The lanczos function is the central lobe of the sinc function.
	 */
	ale_real _lanczos(ale_pos p) const {
		if (fabs(p) >= half_width)
			return 0;

		return sinc::_sinc(p / half_width);
	}

	ale_real _lanczos(point p) const {
		return _lanczos(p[0]) * _lanczos(p[1]);
	}

public:

	/*
	 * Size of filter support, in number of half-cycles to each side of the
	 * filter center.
	 */
	ale_real support() const {
		return half_width;
	}

	virtual int equals(const filter *f) const {
		if (typeid(*f) == typeid(*this))
			return ((lanczos *)f)->half_width == half_width;
		return 0;
	}

	/*
	 * Response of filter at point p
	 */
	virtual ale_real response(point p) const {
		return _lanczos(p);
	}

	lanczos(double half_width) {
		this->half_width = half_width;
	}

};
#endif
