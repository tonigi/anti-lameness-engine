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

#ifndef __mult_h__
#define __mult_h__

/*
 * A class for pointwise multiplication of filters.
 */

class mult : public filter {
private:
	ale_real _support;
	filter *f1, *f2;
public:

	/*
	 * Size of filter support, in number of half-cycles to each side of the
	 * filter center.
	 */
	virtual ale_real support() const {
		return _support;
	}

	virtual int equals(const filter *f) const {
		if (typeid(*f) != typeid(*this))
			return 0;

		const mult *m = (const mult *)f;

		/*
		 * XXX: if we wished, we could recognize commutativity.
		 */
		return (m->f1->equals(f1)
		     && m->f2->equals(f2));
	}

	/*
	 * Response of filter at point p
	 */
	virtual ale_real response(point p) const {
		return f1->response(p) * f2->response(p);
	}

	mult(filter *f1, filter *f2) {
		this->f1 = f1;
		this->f2 = f2;
		assert (f1 != NULL);
		assert (f2 != NULL);

		_support = f1->support();

		if (f2->support() < _support)
			_support = f2->support();
	}

};
#endif
